"""
PyMOL plugin for loading UniProt-linked PDB structures and AlphaFold models.

Install this directory as a PyMOL plugin, or run:
    run H:/My Drive/VSCode_Github/Pymol/UNIPROT_plugin/__init__.py

Then use:
    uniprot_structures P06730
    uniprot_structures IF4E_HUMAN
"""

from __future__ import annotations

import json
import csv
import queue
import os
import re
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

try:
    from pymol import cmd
except Exception:  # pragma: no cover - lets editors import the file outside PyMOL
    cmd = None


UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
USER_AGENT = "PyMOL-UniProt-Structures-Plugin/1.5"
_OPEN_DIALOGS = []
_RUN_LOCK = threading.Lock()


@dataclass
class PDBReference:
    pdb_id: str
    chains: List[str]
    method: str = ""
    resolution: str = ""


@dataclass
class UniProtEntry:
    accession: str
    entry_name: str
    pdb_refs: List[PDBReference]


def _http_json(url: str, timeout: int = 30) -> object:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


class Cancelled(Exception):
    """Cooperative cancellation; completed files and objects are retained."""


def _check_cancel(cancel):
    if cancel is not None and cancel.is_set():
        raise Cancelled("Cancelled; completed files and objects have been kept.")


def _report(progress, message, completed=0, total=0):
    print(f"[UniProt plugin] {message}", flush=True)
    if progress is not None:
        progress({"message": message, "completed": completed, "total": total})


def _cached(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _download(url, path, timeout=60, progress=None, cancel=None):
    # Bounded buffers and atomic replacement prevent partial files becoming cache hits.
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    fd, partial = tempfile.mkstemp(prefix=".uniprot-", suffix=".part", dir=os.path.dirname(os.path.abspath(path)))
    try:
        with os.fdopen(fd, "wb") as handle:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                size = int(response.headers.get("Content-Length") or 0)
                received = 0
                last_update = 0.0
                while True:
                    _check_cancel(cancel)
                    chunk = response.read(256 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    received += len(chunk)
                    now = time.monotonic()
                    if now - last_update >= 1:
                        fraction = f" ({received / size:.0%})" if size else ""
                        _report(progress, f"Downloading {os.path.basename(path)}: {received / 1048576:.1f} MB{fraction}")
                        last_update = now
                if received == 0 or (size and received != size):
                    raise IOError("Empty or incomplete download")
        _check_cancel(cancel)
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.remove(partial)


def _download_pdb_or_cif(pdb_id, cache_dir, progress=None, cancel=None):
    attempts = [(RCSB_PDB_URL.format(pdb_id=pdb_id), os.path.join(cache_dir, f"{pdb_id}.pdb")),
                (RCSB_CIF_URL.format(pdb_id=pdb_id), os.path.join(cache_dir, f"{pdb_id}.cif"))]
    for _, path in attempts:
        if _cached(path):
            return path
    for url, path in attempts:
        _check_cancel(cancel)
        try:
            _download(url, path, progress=progress, cancel=cancel)
            return path
        except urllib.error.HTTPError as exc:
            if exc.code != 404:
                raise
            _report(progress, f"{pdb_id}: {os.path.basename(path)} unavailable; trying next format.")
    raise IOError(f"{pdb_id}: neither PDB nor mmCIF file is available")


def _safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip())
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "UniProt"


def _parse_chain_property(value: str) -> List[str]:
    """Parse UniProt PDB cross-reference chain strings like 'A/B=1-217, C=5-80'."""
    chains: List[str] = []
    for chunk in re.split(r"\s*,\s*", value or ""):
        left = chunk.split("=", 1)[0].strip()
        if not left:
            continue
        for chain in re.split(r"[/\s]+", left):
            chain = chain.strip()
            if chain and chain not in chains:
                chains.append(chain)
    return chains


def _xref_property(properties: Sequence[dict], key: str) -> str:
    for item in properties or []:
        if item.get("key") == key:
            return str(item.get("value", ""))
    return ""


def _entry_from_uniprot_json(data: dict) -> UniProtEntry:
    accession = data.get("primaryAccession") or ""
    entry_name = data.get("uniProtkbId") or accession
    pdb_refs: List[PDBReference] = []

    for xref in data.get("uniProtKBCrossReferences", []) or []:
        if xref.get("database") != "PDB":
            continue
        props = xref.get("properties", []) or []
        chains = _parse_chain_property(_xref_property(props, "Chains"))
        pdb_refs.append(
            PDBReference(
                pdb_id=str(xref.get("id", "")).upper(),
                chains=chains,
                method=_xref_property(props, "Method"),
                resolution=_xref_property(props, "Resolution"),
            )
        )

    return UniProtEntry(accession=accession, entry_name=entry_name, pdb_refs=pdb_refs)


def parse_uniprot_queries(text):
    """Split pasted entry names/accessions on whitespace, commas or semicolons."""
    return list(dict.fromkeys(value for value in re.split(r"[\s,;]+", text.strip()) if value))


def resolve_uniprot(query: str) -> UniProtEntry:
    query = query.strip()
    if not query:
        raise ValueError("Please provide a UniProt accession or entry name.")
    if len(parse_uniprot_queries(query)) != 1 or re.search(r"[\s,;]", query):
        raise ValueError("Expected one UniProt entry. For several entries, paste them in the GUI or use uniprot_structures_file.")

    quoted = urllib.parse.quote(query)
    direct_url = f"{UNIPROT_BASE}/{quoted}.json"
    try:
        return _entry_from_uniprot_json(_http_json(direct_url))
    except urllib.error.HTTPError as exc:
        if exc.code != 404:
            raise

    search = f'(accession:"{query}") OR (id:"{query}")'
    params = urllib.parse.urlencode({"query": search, "format": "json", "size": "1"})
    data = _http_json(f"{UNIPROT_BASE}/search?{params}")
    results = data.get("results", []) if isinstance(data, dict) else []
    if not results:
        raise ValueError(f"No UniProt entry found for '{query}'.")
    return _entry_from_uniprot_json(results[0])


def get_alphafold_pdb_url(accession: str) -> Optional[str]:
    try:
        data = _http_json(ALPHAFOLD_API.format(accession=urllib.parse.quote(accession)))
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise

    if not isinstance(data, list) or not data:
        return None
    first = data[0]
    return first.get("pdbUrl") or first.get("cifUrl")


def _chain_selection(object_name: str, chains: Sequence[str]) -> str:
    if not chains:
        return object_name
    chain_terms = " or ".join(f"chain {chain}" for chain in chains)
    return f"({object_name}) and ({chain_terms})"


def _unique_object_name(base: str) -> str:
    if cmd is None:
        return base
    name = base
    index = 2
    existing = set(cmd.get_names("all"))
    while name in existing:
        name = f"{base}_{index}"
        index += 1
    return name

def _method_subgroup(method: str) -> str:
    normalized = (method or "").lower()
    if "x-ray" in normalized or "xray" in normalized or "diffraction" in normalized:
        label = "X-ray"
    elif "electron" in normalized or normalized == "em" or "microscopy" in normalized:
        label = "EM"
    elif "nmr" in normalized or "magnetic resonance" in normalized:
        label = "NMR"
    else:
        label = "Other_methods"
    return label


def _parse_pdb_methods(value):
    """None/all selects all experimental methods; an empty iterable selects none."""
    all_methods = {"NMR", "X-ray", "EM", "Other_methods"}
    if value is None or (isinstance(value, str) and value.strip().lower() == "all"):
        return all_methods
    if isinstance(value, str):
        value = re.split(r"[+,;|\s]+", value.strip()) if value.strip() else []
    aliases = {"nmr": "NMR", "xray": "X-ray", "x-ray": "X-ray", "em": "EM",
               "other": "Other_methods", "other_methods": "Other_methods"}
    result = set()
    for method in value:
        key = str(method).strip().lower()
        if key == "none":
            continue
        if key not in aliases:
            raise ValueError(f"Unknown PDB category {method!r}; use NMR, X-ray, EM, Other, all, or none.")
        result.add(aliases[key])
    return result


def _update_method_group(entry_group, prefix, label, current, count):
    """Create/rename only nonempty PyMOL groups; counts describe loaded objects."""
    name = _unique_object_name(f"{prefix}_{label}_{count}")
    if current:
        cmd.set_name(current, name)
    else:
        cmd.group(name)
    cmd.group(entry_group, name)
    return name


def _protein_chain(object_name, mapped_chains=()):
    available = cmd.get_chains(f"({object_name}) and polymer.protein")
    candidates = [chain for chain in available if not mapped_chains or chain in mapped_chains]
    if not candidates:
        raise ValueError(f"No protein chains match UniProt mapping {list(mapped_chains)}")
    return "A" if "A" in candidates else candidates[0]


def _protein_selection(object_name, chain):
    # Quoted chain identifiers support blank and multi-character mmCIF identifiers.
    return f"({object_name}) and polymer.protein and chain {json.dumps(chain)}"


def _load_related_chains(path, pdb_id, chains, subgroup_name, progress=None):
    temp_name = _unique_object_name(f"_uniprot_{pdb_id}_full")
    final_name = _unique_object_name(pdb_id)
    try:
        cmd.load(path, temp_name, zoom=0)
        chain = _protein_chain(temp_name, chains)
        selection = _protein_selection(temp_name, chain)
        # Keep every non-protein component, including ligands on other/blank chains.
        # Only other protein chains are pruned; the complete object is transformed
        # during fitting so retained ligands keep their binding-site coordinates.
        cmd.remove(f"({temp_name}) and polymer.protein and not ({selection})")
        cmd.set_name(temp_name, final_name)
        cmd.group(subgroup_name, final_name)
        cmd.hide("everything", final_name)
        _report(progress, f"{pdb_id}: keeping mapped protein chain {chain!r}, all ligands, non-protein components and waters (all model states).")
        return final_name
    except BaseException:
        cmd.delete(temp_name)
        cmd.delete(final_name)
        raise


def _show_structure(object_name):
    cmd.show("cartoon", f"({object_name}) and polymer.protein")
    # Sticks display bonded ligands; nonbonded displays ions and water oxygens.
    nonprotein = f"({object_name}) and not polymer.protein"
    cmd.show("sticks", f"({nonprotein}) and not solvent")
    cmd.show("nonbonded", nonprotein)


def _align_object(object_name, reference, method, progress, reference_b_min=None):
    mobile = f"({object_name}) and polymer.protein"
    target_chain = _protein_chain(reference)
    target = _protein_selection(reference, target_chain)
    if reference_b_min is not None:
        target = f"({target}) and b > {float(reference_b_min):g}"
        if cmd.count_atoms(f"({target}) and name CA") < 3:
            raise ValueError(f"AlphaFold has fewer than 3 CA atoms with pLDDT > {reference_b_min:g}; restricted alignment skipped.")
        _report(progress, f"Alignment reference restricted to AlphaFold pLDDT > {reference_b_min:g}; experimental B-factors are not filtered.")
    _report(progress, f"Aligning {object_name} to {reference}, reference chain {target_chain!r}, using {method}; target: {target}")
    if method == "super":
        result = cmd.super(mobile, target, transform=1)
    elif method == "cealign":
        result = cmd.cealign(target, mobile)
    else:
        result = cmd.align(mobile, target, transform=1)
    # API error returns and empty alignments must not be treated as successful fits.
    if isinstance(result, (int, float)) or result is None:
        raise ValueError(f"Alignment returned no fit: {result}")
    count = result.get("alignment_length", 0) if isinstance(result, dict) else result[1]
    if count <= 0:
        raise ValueError("Alignment contains no matched atoms")
    _report(progress, f"Aligned {object_name}: {result}")


def read_uniprot_file(filename):
    """Read TXT lists or a UniProt column in CSV/TSV; preserve order, remove duplicates."""
    with open(os.path.expanduser(filename), encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("The input file contains no UniProt entries.")
    headers = {"uniprot", "uniprot_id", "uniprot_accession", "accession", "entry", "entry_name", "uniprot_entry_name"}
    delimiter = next((d for d in ("\t", ",", ";") if d in lines[0]), None)
    rows = list(csv.reader(lines, delimiter=delimiter)) if delimiter else [line.split() for line in lines]
    normalized = [re.sub(r"[ -]+", "_", field.strip().lower()) for field in rows[0]]
    column = next((i for i, field in enumerate(normalized) if field in headers), None)
    if column is not None:
        rows = rows[1:]
        values = [row[column].strip() for row in rows if len(row) > column]
    elif delimiter:
        values = [row[0].strip() for row in rows if row]
    else:
        values = [value for row in rows for value in row]
    queries = list(dict.fromkeys(value for value in values if value))
    if not queries:
        raise ValueError("The input file contains no UniProt entries.")
    return queries


def _load_uniprot_structures(query, include_pdb=True, include_alphafold=True, align=True,
                             alignment_method="align", cache_dir=None, limit=0,
                             download_only=False, progress=None, cancel=None,
                             pdb_methods=None, af_confidence_only=False, af_plddt_min=50.0):
    if not download_only and cmd is None:
        raise RuntimeError("Loading structures requires PyMOL.")
    if limit < 0:
        raise ValueError("PDB limit must be zero or positive.")
    if alignment_method not in ("align", "super", "cealign"):
        raise ValueError("Alignment method must be align, super, or cealign.")
    methods = _parse_pdb_methods(pdb_methods)
    af_plddt_min = float(af_plddt_min)
    if not 0 <= af_plddt_min <= 100:
        raise ValueError("AlphaFold pLDDT threshold must be between 0 and 100.")
    if af_confidence_only and align and not download_only and not include_alphafold:
        raise ValueError("Confidence-filtered alignment requires AlphaFold. Enable AlphaFold or turn off the confidence filter.")
    _check_cancel(cancel)
    _report(progress, f"Resolving UniProt {query}...")
    entry = resolve_uniprot(query)
    _check_cancel(cancel)
    group_name = _safe_name(entry.entry_name) if download_only else _unique_object_name(_safe_name(entry.entry_name))
    cache_dir = cache_dir or os.path.join(tempfile.gettempdir(), "pymol_uniprot_structures", _safe_name(entry.entry_name))
    os.makedirs(cache_dir, exist_ok=True)
    pdb_refs = [ref for ref in entry.pdb_refs if ref.pdb_id and _method_subgroup(ref.method) in methods] if include_pdb else []
    pdb_refs = list({ref.pdb_id: ref for ref in pdb_refs}.values())[:limit or None]
    _report(progress, f"{entry.entry_name} ({entry.accession}): {len(entry.pdb_refs)} PDB references found; {len(pdb_refs)} selected. Cache: {cache_dir}")
    af_url = None
    errors = []
    if include_alphafold:
        _report(progress, "Looking up AlphaFold model...")
        try:
            af_url = get_alphafold_pdb_url(entry.accession)
            if not af_url:
                _report(progress, "No AlphaFold model available.")
        except MemoryError:
            raise
        except Exception as exc:
            errors.append(f"AlphaFold lookup: {exc}")
            _report(progress, errors[-1])
    jobs = ([None] if af_url else []) + pdb_refs
    total = len(jobs)
    completed = downloaded = cached = failed = 0
    loaded, visible, files = [], [], []
    reference = None
    reference_is_af = False
    method_groups, method_counts = {}, {}
    prefix = _safe_name(entry.entry_name.split("_", 1)[0])
    if not download_only and jobs:
        cmd.group(group_name)

    def status(message):
        _report(progress, f"{entry.entry_name}: {completed}/{total} processed ({completed / total if total else 1:.0%}); "
                f"{total-completed} remaining; {downloaded} downloaded, {cached} cached, {failed} failed; "
                f"{len(loaded)} loaded. {message}", completed, total)

    def activity(event):
        if progress is not None:
            progress(dict(event, completed=completed, total=total))

    status("Starting." if total else "No coordinate files to process.")
    for ref in jobs:
        _check_cancel(cancel)
        is_af = ref is None
        label = "AlphaFold" if is_af else ref.pdb_id
        status(f"Fetching {label}...")
        try:
            if is_af:
                extension = ".cif" if urllib.parse.urlparse(af_url).path.lower().endswith(".cif") else ".pdb"
                path = os.path.join(cache_dir, f"AlphaFold_{entry.accession}{extension}")
                was_cached = _cached(path)
                if not was_cached:
                    _download(af_url, path, progress=activity, cancel=cancel)
            else:
                was_cached = any(_cached(os.path.join(cache_dir, ref.pdb_id + ext)) for ext in (".pdb", ".cif"))
                path = _download_pdb_or_cif(ref.pdb_id, cache_dir, progress=activity, cancel=cancel)
            cached += int(was_cached)
            downloaded += int(not was_cached)
            files.append(path)
        except (Cancelled, MemoryError):
            raise
        except Exception as exc:
            failed += 1
            errors.append(f"{label} download: {exc}")
            completed += 1
            status(errors[-1])
            continue
        _check_cancel(cancel)
        if not download_only:
            obj = None
            try:
                status(f"Loading {label} into PyMOL...")
                if is_af:
                    obj = _load_related_chains(path, f"AF_{_safe_name(entry.entry_name)}", [], group_name, activity)
                else:
                    subgroup_label = _method_subgroup(ref.method)
                    current = method_groups.get(subgroup_label)
                    # Loading into the root first avoids leaving empty method groups on failure.
                    obj = _load_related_chains(path, ref.pdb_id, ref.chains, group_name, activity)
                    count = method_counts.get(subgroup_label, 0) + 1
                    subgroup = _update_method_group(group_name, prefix, subgroup_label, current, count)
                    cmd.group(subgroup, obj)
                    method_groups[subgroup_label] = subgroup
                    method_counts[subgroup_label] = count
                loaded.append(obj)
                _check_cancel(cancel)
                if align and af_confidence_only and not is_af and not reference_is_af:
                    raise ValueError("No AlphaFold reference is loaded; confidence-filtered alignment skipped. Structure retained hidden.")
                if reference is None:
                    reference = obj
                    reference_is_af = is_af
                    if not is_af and align:
                        _report(activity, f"AlphaFold unavailable/disabled; using first PDB {obj} as reference.")
                elif align:
                    _align_object(obj, reference, alignment_method, activity,
                                  reference_b_min=af_plddt_min if af_confidence_only and reference_is_af else None)
                _show_structure(obj)
                visible.append(obj)
            except (Cancelled, MemoryError):
                raise
            except Exception as exc:
                errors.append(f"{label} load/alignment: {exc}")
                _report(activity, errors[-1])
                if obj:
                    cmd.hide("everything", obj)
                    cmd.disable(obj)
        completed += 1
        status(f"Finished {label}.")
    if visible:
        cmd.zoom(" or ".join(f"({obj})" for obj in visible))
    status(f"Done. {len(visible)} visible; {len(errors)} error(s). Files: {cache_dir}")
    return {"accession": entry.accession, "entry_name": entry.entry_name, "group": group_name,
            "loaded": loaded, "visible": visible, "files": files, "cache_dir": cache_dir,
            "pdb_count": len(entry.pdb_refs), "downloaded": downloaded, "cached": cached,
            "failed": failed, "errors": errors, "method_groups": method_groups,
            "method_counts": method_counts, "pdb_methods": sorted(methods)}


def load_uniprot_structures(query, include_pdb=True, include_alphafold=True, align=True,
                            alignment_method="align", cache_dir=None, limit=0,
                            download_only=False, progress=None, cancel=None,
                            pdb_methods=None, af_confidence_only=False, af_plddt_min=50.0):
    if not _RUN_LOCK.acquire(blocking=False):
        raise RuntimeError("A UniProt job is already running. Cancel it or wait for it to finish.")
    try:
        return _load_uniprot_structures(query, include_pdb, include_alphafold, align,
                                        alignment_method, cache_dir, limit, download_only, progress, cancel,
                                        pdb_methods, af_confidence_only, af_plddt_min)
    finally:
        _RUN_LOCK.release()


def load_uniprot_batch(queries, cache_root=None, progress=None, cancel=None, **options):
    if not _RUN_LOCK.acquire(blocking=False):
        raise RuntimeError("A UniProt job is already running. Cancel it or wait for it to finish.")
    summaries, errors = [], []
    try:
        queries = list(dict.fromkeys(q.strip() for q in queries if q.strip()))
        if not queries:
            raise ValueError("Enter a UniProt accession/name or choose an input file.")
        for index, query in enumerate(queries, 1):
            _check_cancel(cancel)
            prefix = f"Protein {index}/{len(queries)} ({query})"
            def report(event):
                if progress is not None:
                    progress(dict(event, message=f"{prefix}: {event['message']}"))
            _report(report, f"Starting {prefix}.")
            try:
                folder = os.path.join(cache_root, _safe_name(query)) if cache_root else None
                summaries.append(_load_uniprot_structures(query, cache_dir=folder, progress=report, cancel=cancel, **options))
            except Cancelled:
                raise
            except MemoryError:
                raise  # Continuing after memory exhaustion can destabilize PyMOL.
            except Exception as exc:
                errors.append(f"{query}: {exc}")
                _report(report, f"Failed: {exc}; continuing with next protein.")
        detail_errors = len(errors) + sum(len(s["errors"]) for s in summaries)
        _report(progress, f"Batch done: {len(summaries)}/{len(queries)} proteins processed; "
                f"{sum(s['downloaded'] for s in summaries)} downloaded, {sum(s['cached'] for s in summaries)} cached; "
                f"{detail_errors} error(s).", len(queries), len(queries))
        return {"summaries": summaries, "errors": errors}
    finally:
        _RUN_LOCK.release()


def uniprot_structures(query, include_pdb="1", include_alphafold="1", align="1",
                       alignment_method="align", limit="0", download_only="0", cache_dir="",
                       pdb_methods="all", af_confidence_only="0", af_plddt_min="50"):
    """Load one UniProt entry; trailing options preserve existing command syntax."""
    return load_uniprot_structures(query, bool(int(include_pdb)), bool(int(include_alphafold)),
                                   bool(int(align)), alignment_method, cache_dir or None, int(limit),
                                   bool(int(download_only)), pdb_methods=pdb_methods,
                                   af_confidence_only=bool(int(af_confidence_only)), af_plddt_min=float(af_plddt_min))


def uniprot_structures_file(filename, limit="0", download_only="0", cache_dir="",
                            pdb_methods="all", include_alphafold="1", af_confidence_only="0", af_plddt_min="50"):
    """Process a TXT/CSV/TSV UniProt list sequentially."""
    return load_uniprot_batch(read_uniprot_file(filename), cache_root=cache_dir or None,
                              limit=int(limit), download_only=bool(int(download_only)), pdb_methods=pdb_methods,
                              include_alphafold=bool(int(include_alphafold)),
                              af_confidence_only=bool(int(af_confidence_only)), af_plddt_min=float(af_plddt_min))


def list_uniprot_groups():
    """Find entry-name groups, including those restored from a saved session."""
    if cmd is None:
        return []
    return sorted(name for name in cmd.get_names("objects")
                  if cmd.get_type(name) == "object:group"
                  and re.fullmatch(r"[A-Z0-9]+_[A-Z][A-Z0-9]*(?:_\d+)?", name)
                  and not re.search(r"_(?:NMR|EM|X-ray|Other_methods)_\d+(?:_\d+)?$", name))


def uniprot_show(group, only="1"):
    """Show a protein and populated subgroups; only=1 hides all other scene objects."""
    if cmd is None:
        raise RuntimeError("This command requires PyMOL.")
    group = group.strip()
    exclusive = bool(int(only))
    if not _RUN_LOCK.acquire(blocking=False):
        raise RuntimeError("Wait for the running UniProt job to finish before changing protein visibility.")
    try:
        if group not in cmd.get_names("objects") or cmd.get_type(group) != "object:group":
            raise ValueError(f"No loaded protein group named {group!r}.")
        objects = cmd.get_object_list(f"({group})")
        if not objects:
            raise ValueError(f"Group {group!r} contains no molecular structures.")
        if exclusive:
            cmd.disable("all")
        # Enabling a parent alone does not enable its children. Enable each molecular
        # member with its ancestors, without serializing a potentially huge session.
        for name in objects:
            cmd.enable(name, parents=1)
        cmd.zoom(group)
        message = f"Showing {'only ' if exclusive else ''}{group}: {len(objects)} structure(s), with their parent folders enabled."
        _report(None, message)
        return objects
    finally:
        _RUN_LOCK.release()


def uniprot_hide(group):
    """Hide a whole protein through its parent, preserving its children's switches."""
    if cmd is None:
        raise RuntimeError("This command requires PyMOL.")
    if not _RUN_LOCK.acquire(blocking=False):
        raise RuntimeError("Wait for the running UniProt job to finish before changing protein visibility.")
    try:
        group = group.strip()
        if group not in cmd.get_names("objects") or cmd.get_type(group) != "object:group":
            raise ValueError(f"No loaded protein group named {group!r}.")
        cmd.disable(group)
    finally:
        _RUN_LOCK.release()


def _start_job(queries, options):
    events = queue.Queue()
    cancel = threading.Event()
    def worker():
        try:
            result = load_uniprot_batch(queries, progress=lambda event: events.put(("progress", event)),
                                        cancel=cancel, **options)
            events.put(("done", result))
        except Cancelled as exc:
            _report(None, str(exc))
            events.put(("cancelled", str(exc)))
        except MemoryError:
            message = "Out of memory. Job stopped; use Download only or a smaller PDB limit. Completed files and objects were kept."
            _report(None, message)
            events.put(("error", message))
        except Exception as exc:
            _report(None, f"Job stopped: {exc}")
            events.put(("error", str(exc)))
    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    return thread, cancel, events


def _open_qt_gui():
    try:
        from pymol import plugins
        from pymol.Qt import QtWidgets, QtCore
    except Exception as exc:
        print(f"[UniProt plugin] Qt GUI is not available: {exc}")
        return False
    parent = None
    try:
        parent = plugins.get_qtwindow()
    except Exception:
        pass
    job = {"thread": None, "cancel": None, "events": None, "started": 0, "message": "Ready"}
    def busy():
        return job["events"] is not None
    class Dialog(QtWidgets.QDialog):
        def reject(self):
            if busy():
                job["cancel"].set()
                job["message"] = "Cancellation requested; waiting for current operation..."
                return
            super().reject()
        def closeEvent(self, event):
            if busy():
                self.reject()
                event.ignore()
            else:
                event.accept()
    dialog = Dialog(parent)
    dialog.setWindowTitle("UniProt Structures")
    dialog.setMinimumWidth(680)
    query_edit = QtWidgets.QLineEdit(dialog)
    query_edit.setPlaceholderText("One or more entries, e.g. BCL2_HUMAN CAB39_HUMAN (spaces, commas or semicolons)")
    file_edit = QtWidgets.QLineEdit(dialog)
    file_edit.setPlaceholderText("Optional TXT, CSV or TSV list")
    file_button = QtWidgets.QPushButton("Browse...", dialog)
    def browse_file():
        path, _ = QtWidgets.QFileDialog.getOpenFileName(dialog, "UniProt list", "", "UniProt lists (*.txt *.csv *.tsv);;All files (*)")
        if path:
            file_edit.setText(path)
    file_button.clicked.connect(browse_file)
    cache_edit = QtWidgets.QLineEdit(dialog)
    cache_edit.setPlaceholderText("Default: system temporary folder")
    cache_button = QtWidgets.QPushButton("Browse...", dialog)
    def browse_cache():
        path = QtWidgets.QFileDialog.getExistingDirectory(dialog, "Download folder")
        if path:
            cache_edit.setText(path)
    cache_button.clicked.connect(browse_cache)
    pdb_check = QtWidgets.QCheckBox("Download PDB structures", dialog)
    pdb_check.setChecked(True)
    category_checks = {}
    for label in ("NMR", "X-ray", "EM", "Other_methods"):
        checkbox = QtWidgets.QCheckBox(label.replace("_methods", ""), dialog)
        checkbox.setObjectName("category_" + label)
        checkbox.setChecked(True)
        category_checks[label] = checkbox
    af_check = QtWidgets.QCheckBox("Download AlphaFold model", dialog)
    af_check.setChecked(True)
    align_check = QtWidgets.QCheckBox("Align structures (AlphaFold reference when available)", dialog)
    align_check.setChecked(True)
    confidence_check = QtWidgets.QCheckBox("Fit only to AlphaFold regions with pLDDT (B-factor) >", dialog)
    confidence_check.setObjectName("af_confidence_only")
    confidence_spin = QtWidgets.QDoubleSpinBox(dialog)
    confidence_spin.setObjectName("af_plddt_min")
    confidence_spin.setRange(0, 100)
    confidence_spin.setValue(50)
    confidence_spin.setEnabled(False)
    confidence_check.toggled.connect(confidence_spin.setEnabled)
    def update_source_controls():
        for check in category_checks.values():
            check.setEnabled(pdb_check.isChecked())
        eligible = af_check.isChecked() and align_check.isChecked()
        if not eligible:
            confidence_check.setChecked(False)
        confidence_check.setEnabled(eligible)
        confidence_spin.setEnabled(eligible and confidence_check.isChecked())
    pdb_check.toggled.connect(update_source_controls)
    af_check.toggled.connect(update_source_controls)
    align_check.toggled.connect(update_source_controls)
    download_check = QtWidgets.QCheckBox("Download only: save files without loading into PyMOL", dialog)
    method_combo = QtWidgets.QComboBox(dialog)
    method_combo.addItems(["align", "super", "cealign"])
    limit_spin = QtWidgets.QSpinBox(dialog)
    limit_spin.setRange(0, 100000)
    limit_spin.setSpecialValueText("All (no limit)")
    status = QtWidgets.QLabel("Ready", dialog)
    status.setWordWrap(True)
    bar = QtWidgets.QProgressBar(dialog)
    bar.setRange(0, 100)
    bar.setValue(0)
    log = QtWidgets.QPlainTextEdit(dialog)
    log.setReadOnly(True)
    log.setMaximumBlockCount(1500)
    form = QtWidgets.QFormLayout()
    form.addRow("UniProt entries", query_edit)
    for title, edit, button in (("List file (takes priority)", file_edit, file_button), ("Download folder", cache_edit, cache_button)):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(edit)
        row.addWidget(button)
        form.addRow(title, row)
    form.addRow("", pdb_check)
    category_row = QtWidgets.QHBoxLayout()
    for check in category_checks.values():
        category_row.addWidget(check)
    form.addRow("PDB categories", category_row)
    for widget in (af_check, align_check, download_check):
        form.addRow("", widget)
    confidence_row = QtWidgets.QHBoxLayout()
    confidence_row.addWidget(confidence_check)
    confidence_row.addWidget(confidence_spin)
    form.addRow("", confidence_row)
    form.addRow("Alignment", method_combo)
    form.addRow("PDB limit per protein", limit_spin)
    form.addRow("", QtWidgets.QLabel("Keeps one mapped protein chain: A if available, otherwise the first mapped chain.", dialog))
    load_button = QtWidgets.QPushButton("Start", dialog)
    cancel_button = QtWidgets.QPushButton("Cancel", dialog)
    cancel_button.setEnabled(False)
    close_button = QtWidgets.QPushButton("Close", dialog)
    controls = [query_edit, file_edit, file_button, cache_edit, cache_button, pdb_check,
                af_check, align_check, download_check, method_combo, limit_spin, load_button, close_button,
                confidence_check, confidence_spin] + list(category_checks.values())
    viewer = QtWidgets.QGroupBox("View loaded proteins (no download needed)", dialog)
    viewer_layout = QtWidgets.QHBoxLayout(viewer)
    protein_combo = QtWidgets.QComboBox(viewer)
    protein_combo.setObjectName("uniprot_protein_selector")
    refresh_button = QtWidgets.QPushButton("Refresh", viewer)
    only_button = QtWidgets.QPushButton("Show only this protein", viewer)
    show_button = QtWidgets.QPushButton("Show", viewer)
    hide_button = QtWidgets.QPushButton("Hide", viewer)
    for widget in (protein_combo, refresh_button, only_button, show_button, hide_button):
        viewer_layout.addWidget(widget)
        controls.append(widget)
    def refresh_proteins():
        selected = protein_combo.currentText()
        protein_combo.clear()
        protein_combo.addItems(list_uniprot_groups())
        index = protein_combo.findText(selected)
        if index >= 0:
            protein_combo.setCurrentIndex(index)
    def view_protein(action):
        try:
            name = protein_combo.currentText()
            if not name:
                raise ValueError("No loaded proteins found. Load structures first, or click Refresh after opening a saved session.")
            if action == "hide":
                uniprot_hide(name)
            else:
                uniprot_show(name, only="1" if action == "only" else "0")
            status.setText(f"{'Hidden' if action == 'hide' else 'Showing only' if action == 'only' else 'Showing'} {name}")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(dialog, "Protein visibility", str(exc))
    refresh_button.clicked.connect(refresh_proteins)
    only_button.clicked.connect(lambda: view_protein("only"))
    show_button.clicked.connect(lambda: view_protein("show"))
    hide_button.clicked.connect(lambda: view_protein("hide"))
    refresh_proteins()
    def set_running(value):
        for control in controls:
            control.setEnabled(not value)
        cancel_button.setEnabled(value)
        if not value:
            update_source_controls()
    def cancel_job():
        if busy():
            job["cancel"].set()
            job["message"] = "Cancellation requested; waiting for current operation..."
    cancel_button.clicked.connect(cancel_job)
    timer = QtCore.QTimer(dialog)
    def poll():
        if not busy():
            return
        finished = False
        while True:
            try:
                kind, data = job["events"].get_nowait()
            except queue.Empty:
                break
            if kind == "progress":
                job["message"] = data["message"]
                log.appendPlainText(data["message"])
                total = data["total"]
                bar.setRange(0, total if total else 0)
                if total:
                    bar.setValue(data["completed"])
            else:
                finished = True
                if kind != "done":
                    job["message"] = data
                    log.appendPlainText(data)
                if kind == "error":
                    QtWidgets.QMessageBox.critical(dialog, "UniProt Structures", data)
        status.setText(f"Elapsed {int(time.monotonic()-job['started'])} s | {job['message']}")
        if finished:
            timer.stop()
            job["events"] = None
            set_running(False)
            refresh_proteins()
            if bar.maximum() == 0:
                bar.setRange(0, 100)
                bar.setValue(0)
    timer.timeout.connect(poll)
    def on_load():
        try:
            queries = read_uniprot_file(file_edit.text().strip()) if file_edit.text().strip() else parse_uniprot_queries(query_edit.text())
            if not any(queries):
                raise ValueError("Enter a UniProt accession/name or browse to a list file.")
            selected_methods = [label for label, check in category_checks.items() if check.isChecked()]
            if not (pdb_check.isChecked() and selected_methods) and not af_check.isChecked():
                raise ValueError("Select at least one PDB category or AlphaFold.")
            options = dict(include_pdb=pdb_check.isChecked(), include_alphafold=af_check.isChecked(),
                           align=align_check.isChecked(), alignment_method=method_combo.currentText(),
                           limit=limit_spin.value(), download_only=download_check.isChecked(),
                           cache_root=cache_edit.text().strip() or None, pdb_methods=selected_methods,
                           af_confidence_only=confidence_check.isChecked(), af_plddt_min=confidence_spin.value())
            set_running(True)
            log.clear()
            bar.setRange(0, 0)
            job["started"] = time.monotonic()
            job["message"] = "Starting..."
            job["thread"], job["cancel"], job["events"] = _start_job(queries, options)
            timer.start(200)
        except Exception as exc:
            set_running(False)
            QtWidgets.QMessageBox.warning(dialog, "UniProt Structures", str(exc))
    load_button.clicked.connect(on_load)
    close_button.clicked.connect(dialog.close)
    buttons = QtWidgets.QHBoxLayout()
    for button in (load_button, cancel_button, close_button):
        buttons.addWidget(button)
    layout = QtWidgets.QVBoxLayout(dialog)
    layout.addLayout(form)
    layout.addWidget(viewer)
    layout.addWidget(status)
    layout.addWidget(bar)
    layout.addWidget(log)
    layout.addLayout(buttons)
    dialog.finished.connect(lambda _code: _OPEN_DIALOGS.remove(dialog) if dialog in _OPEN_DIALOGS else None)
    _OPEN_DIALOGS.append(dialog)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return True


def _open_gui(*_args, **_kwargs):
    if _open_qt_gui():
        return
    # Tk fallback: worker touches only PyMOL and a queue, never Tk widgets.
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox, filedialog
    except Exception as exc:
        print(f"[UniProt plugin] Tkinter is not available: {exc}")
        return
    window = tk.Toplevel()
    window.title("UniProt Structures")
    frame = ttk.Frame(window, padding=12)
    frame.grid(sticky="nsew")
    query = tk.StringVar()
    filename = tk.StringVar()
    cache = tk.StringVar()
    pdb = tk.BooleanVar(value=True)
    af = tk.BooleanVar(value=True)
    align = tk.BooleanVar(value=True)
    download_only = tk.BooleanVar(value=False)
    method = tk.StringVar(value="align")
    limit = tk.StringVar(value="0")
    status = tk.StringVar(value="Ready")
    categories = {label: tk.BooleanVar(value=True) for label in ("NMR", "X-ray", "EM", "Other_methods")}
    confidence = tk.BooleanVar(value=False)
    threshold = tk.StringVar(value="50")
    controls = []
    for row, (label, value) in enumerate((("UniProt entries", query), ("List file (takes priority)", filename), ("Download folder", cache))):
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w")
        edit = ttk.Entry(frame, textvariable=value, width=55)
        edit.grid(row=row, column=1, sticky="ew")
        controls.append(edit)
    def browse_file():
        path = filedialog.askopenfilename(parent=window, filetypes=[("UniProt lists", "*.txt *.csv *.tsv"), ("All files", "*")])
        if path:
            filename.set(path)
    def browse_cache():
        path = filedialog.askdirectory(parent=window)
        if path:
            cache.set(path)
    for row, action in ((1, browse_file), (2, browse_cache)):
        button = ttk.Button(frame, text="Browse...", command=action)
        button.grid(row=row, column=2)
        controls.append(button)
    for row, (label, value) in enumerate((("Download PDB structures", pdb), ("Download AlphaFold", af),
            ("Align structures", align), ("Download only (no PyMOL loading)", download_only)), 3):
        widget = ttk.Checkbutton(frame, text=label, variable=value)
        widget.grid(row=row, column=0, columnspan=3, sticky="w")
        controls.append(widget)
    category_frame = ttk.Frame(frame)
    category_frame.grid(row=13, column=0, columnspan=3, sticky="w")
    ttk.Label(category_frame, text="PDB categories:").pack(side="left")
    for label, variable in categories.items():
        checkbox = ttk.Checkbutton(category_frame, text=label.replace("_methods", ""), variable=variable)
        checkbox.pack(side="left")
        controls.append(checkbox)
    confidence_frame = ttk.Frame(frame)
    confidence_frame.grid(row=14, column=0, columnspan=3, sticky="w")
    checkbox = ttk.Checkbutton(confidence_frame, text="Fit only to AlphaFold pLDDT (B-factor) >", variable=confidence)
    checkbox.pack(side="left")
    controls.append(checkbox)
    threshold_edit = ttk.Entry(confidence_frame, textvariable=threshold, width=6)
    threshold_edit.pack(side="left")
    controls.append(threshold_edit)
    ttk.Label(frame, text="PDB limit per protein (0 = all)").grid(row=7, column=0)
    widget = ttk.Entry(frame, textvariable=limit, width=10)
    widget.grid(row=7, column=1, sticky="w")
    controls.append(widget)
    method_widget = ttk.Combobox(frame, textvariable=method, values=("align", "super", "cealign"), state="readonly")
    method_widget.grid(row=8, column=1, sticky="w")
    ttk.Label(frame, text="Alignment").grid(row=8, column=0)
    ttk.Label(frame, textvariable=status, wraplength=650).grid(row=9, column=0, columnspan=3, sticky="w")
    bar = ttk.Progressbar(frame, maximum=100, length=650)
    bar.grid(row=10, column=0, columnspan=3)
    log = tk.Text(frame, width=90, height=12, state="disabled")
    log.grid(row=11, column=0, columnspan=3)
    job = {"events": None, "cancel": None}
    def append(message):
        log.configure(state="normal")
        log.insert("end", message + "\n")
        if int(log.index("end-1c").split(".")[0]) > 1500:
            log.delete("1.0", "100.0")
        log.see("end")
        log.configure(state="disabled")
    def set_running(value):
        for control in controls:
            control.configure(state="disabled" if value else "normal")
        method_widget.configure(state="disabled" if value else "readonly")
        cancel_button.configure(state="normal" if value else "disabled")
    def poll():
        finished = False
        while True:
            try:
                kind, data = job["events"].get_nowait()
            except queue.Empty:
                break
            if kind == "progress":
                job["message"] = data["message"]
                append(data["message"])
                if data["total"]:
                    bar.configure(value=100 * data["completed"] / data["total"])
            else:
                finished = True
                if kind != "done":
                    job["message"] = data
                    append(data)
                if kind == "error":
                    messagebox.showerror("UniProt Structures", data, parent=window)
        status.set(f"Elapsed {int(time.monotonic()-job['started'])} s | {job['message']}")
        if finished:
            job["events"] = None
            set_running(False)
        else:
            window.after(200, poll)
    def on_load():
        try:
            queries = read_uniprot_file(filename.get().strip()) if filename.get().strip() else parse_uniprot_queries(query.get())
            n = int(limit.get() or "0")
            if n < 0 or not any(queries) or not (pdb.get() or af.get()):
                raise ValueError("Choose UniProt input, at least one structure source, and a nonnegative limit.")
            selected_methods = [label for label, variable in categories.items() if variable.get()]
            if not (pdb.get() and selected_methods) and not af.get():
                raise ValueError("Select at least one experimental category or AlphaFold.")
            if confidence.get() and align.get() and not download_only.get() and not af.get():
                raise ValueError("Confidence-filtered alignment requires AlphaFold.")
            options = dict(include_pdb=pdb.get(), include_alphafold=af.get(), align=align.get(),
                           download_only=download_only.get(), alignment_method=method.get(), limit=n,
                           cache_root=cache.get().strip() or None, pdb_methods=selected_methods,
                           af_confidence_only=confidence.get(), af_plddt_min=float(threshold.get()))
            set_running(True)
            job["started"] = time.monotonic()
            job["message"] = "Starting..."
            job["thread"], job["cancel"], job["events"] = _start_job(queries, options)
            window.after(200, poll)
        except Exception as exc:
            set_running(False)
            messagebox.showwarning("UniProt Structures", str(exc), parent=window)
    def cancel_job():
        if job["events"] is not None:
            job["cancel"].set()
            job["message"] = "Cancellation requested; waiting for current operation..."
    def close():
        if job["events"] is not None:
            cancel_job()
        else:
            window.destroy()
    start_button = ttk.Button(frame, text="Start", command=on_load)
    start_button.grid(row=12, column=0)
    controls.append(start_button)
    cancel_button = ttk.Button(frame, text="Cancel", command=cancel_job, state="disabled")
    cancel_button.grid(row=12, column=1)
    close_button = ttk.Button(frame, text="Close", command=close)
    close_button.grid(row=12, column=2)
    controls.append(close_button)
    window.protocol("WM_DELETE_WINDOW", close)


def uniprot_structures_gui():
    """Open the UniProt Structures GUI."""
    _open_gui()


def _register_commands():
    if cmd is not None:
        cmd.extend("uniprot_structures", uniprot_structures)
        cmd.extend("uniprot_structures_file", uniprot_structures_file)
        cmd.extend("uniprot_structures_gui", uniprot_structures_gui)
        cmd.extend("uniprot_show", uniprot_show)
        cmd.extend("uniprot_hide", uniprot_hide)


def __init_plugin__(app=None):
    _register_commands()
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt("UniProt Structures", _open_gui)
    except Exception:
        try:
            app.rootmenu.addmenuitem("Plugin", "command", "UniProt Structures", label="UniProt Structures", command=_open_gui)
        except Exception:
            pass


_register_commands()
