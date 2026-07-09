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
import os
import re
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

try:
    from pymol import cmd
except Exception:  # pragma: no cover - lets editors import the file outside PyMOL
    cmd = None


UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
USER_AGENT = "PyMOL-UniProt-Structures-Plugin/1.0"
_OPEN_DIALOGS = []


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


def _download(url: str, path: str, timeout: int = 60) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response, open(path, "wb") as handle:
        handle.write(response.read())


def _download_pdb_or_cif(pdb_id: str, cache_dir: str) -> Optional[str]:
    pdb_path = os.path.join(cache_dir, f"{pdb_id}.pdb")
    if os.path.exists(pdb_path):
        return pdb_path

    cif_path = os.path.join(cache_dir, f"{pdb_id}.cif")
    if os.path.exists(cif_path):
        return cif_path

    attempts = (
        (RCSB_PDB_URL.format(pdb_id=pdb_id), pdb_path),
        (RCSB_CIF_URL.format(pdb_id=pdb_id), cif_path),
    )
    for url, path in attempts:
        try:
            print(f"[UniProt plugin] downloading PDB {pdb_id} from {url}...")
            _download(url, path)
            return path
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                print(f"[UniProt plugin] {pdb_id}: {os.path.basename(path)} not found; trying next format.")
                continue
            print(f"[UniProt plugin] {pdb_id}: download failed ({exc}); skipped.")
            return None
        except Exception as exc:
            print(f"[UniProt plugin] {pdb_id}: download failed ({exc}); skipped.")
            return None

    print(f"[UniProt plugin] {pdb_id}: no PDB or mmCIF coordinate file found; skipped.")
    return None


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


def resolve_uniprot(query: str) -> UniProtEntry:
    query = query.strip()
    if not query:
        raise ValueError("Please provide a UniProt accession or entry name.")

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
    existing = set(cmd.get_names("objects"))
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


def _ensure_method_groups(entry_group: str, include_other: bool = False) -> Dict[str, str]:
    labels = ["NMR", "X-ray", "EM"]
    if include_other:
        labels.append("Other_methods")

    groups: Dict[str, str] = {}
    for label in labels:
        subgroup = _unique_object_name(label)
        cmd.group(subgroup)
        cmd.group(entry_group, subgroup)
        groups[label] = subgroup
    return groups


def _load_related_chains(path: str, pdb_id: str, chains: Sequence[str], subgroup_name: str) -> Optional[str]:
    temp_name = _unique_object_name(f"_{_safe_name(subgroup_name)}_{pdb_id}_full")
    final_name = _unique_object_name(pdb_id)
    cmd.load(path, temp_name)

    if chains:
        selection = _chain_selection(temp_name, chains)
        if cmd.count_atoms(selection) == 0:
            cmd.delete(temp_name)
            print(f"[UniProt plugin] {pdb_id}: no atoms found for chains {', '.join(chains)}; skipped.")
            return None
        cmd.create(final_name, selection)
        cmd.delete(temp_name)
    else:
        cmd.set_name(temp_name, final_name)

    cmd.group(subgroup_name, final_name)
    return final_name


def _mobile_alignment_selection(object_name: str, preferred_chain: str = "A") -> str:
    preferred = f"({object_name}) and chain {preferred_chain}"
    if cmd.count_atoms(preferred):
        return preferred
    print(f"[UniProt plugin] {object_name}: chain {preferred_chain} not found; aligning whole object.")
    return object_name

def _align_to_reference(objects: Iterable[str], reference: str, method: str = "align") -> None:
    for object_name in objects:
        if object_name == reference:
            continue
        mobile_selection = _mobile_alignment_selection(object_name, "A")
        try:
            print(f"[UniProt plugin] aligning {mobile_selection} -> {reference} with {method}")
            if method == "super":
                result = cmd.super(mobile_selection, reference, transform=1)
            elif method == "cealign":
                result = cmd.cealign(reference, mobile_selection)
            else:
                result = cmd.align(mobile_selection, reference, transform=1)
            print(f"[UniProt plugin] aligned {object_name} to {reference}: {result}")
        except Exception as exc:
            print(f"[UniProt plugin] alignment failed for {object_name}: {exc}")

def load_uniprot_structures(
    query: str,
    include_pdb: bool = True,
    include_alphafold: bool = True,
    align: bool = True,
    alignment_method: str = "align",
    cache_dir: Optional[str] = None,
    limit: int = 0,
) -> Dict[str, object]:
    if cmd is None:
        raise RuntimeError("This plugin must be run inside PyMOL.")

    entry = resolve_uniprot(query)
    group_name = _safe_name(entry.entry_name)
    cache_dir = cache_dir or os.path.join(tempfile.gettempdir(), "pymol_uniprot_structures", group_name)
    os.makedirs(cache_dir, exist_ok=True)

    cmd.group(group_name)
    method_groups = _ensure_method_groups(group_name)
    loaded: List[str] = []
    reference: Optional[str] = None

    print(f"[UniProt plugin] UniProt: {entry.entry_name} ({entry.accession})")
    print(f"[UniProt plugin] PDB cross-references in UniProt: {len(entry.pdb_refs)}")

    if include_alphafold:
        af_url = get_alphafold_pdb_url(entry.accession)
        if af_url:
            extension = ".cif" if af_url.lower().endswith(".cif") else ".pdb"
            af_path = os.path.join(cache_dir, f"AlphaFold_{entry.accession}{extension}")
            if not os.path.exists(af_path):
                print(f"[UniProt plugin] downloading AlphaFold model...")
                _download(af_url, af_path)
            reference = _unique_object_name(f"AF_{group_name}")
            cmd.load(af_path, reference)
            cmd.group(group_name, reference)
            loaded.append(reference)
        else:
            print(f"[UniProt plugin] no AlphaFold model found for {entry.accession}.")

    pdb_refs = entry.pdb_refs[: limit or None]
    if include_pdb:
        for ref in pdb_refs:
            if not ref.pdb_id:
                continue
            pdb_path = _download_pdb_or_cif(ref.pdb_id, cache_dir)
            if not pdb_path:
                continue
            time.sleep(0.1)
            subgroup_label = _method_subgroup(ref.method)
            if subgroup_label not in method_groups:
                subgroup = _unique_object_name(subgroup_label)
                cmd.group(subgroup)
                cmd.group(group_name, subgroup)
                method_groups[subgroup_label] = subgroup
            subgroup_name = method_groups[subgroup_label]
            object_name = _load_related_chains(pdb_path, ref.pdb_id, ref.chains, subgroup_name)
            if object_name:
                loaded.append(object_name)

    if align and reference:
        _align_to_reference(loaded, reference, alignment_method)
    elif align:
        print("[UniProt plugin] no AlphaFold reference loaded; alignment skipped.")

    cmd.zoom(group_name)
    print(f"[UniProt plugin] loaded {len(loaded)} object(s) in group '{group_name}'.")
    print(f"[UniProt plugin] downloaded files are in: {cache_dir}")
    return {
        "accession": entry.accession,
        "entry_name": entry.entry_name,
        "group": group_name,
        "loaded": loaded,
        "cache_dir": cache_dir,
        "pdb_count": len(entry.pdb_refs),
    }


def uniprot_structures(
    query: str,
    include_pdb: str = "1",
    include_alphafold: str = "1",
    align: str = "1",
    alignment_method: str = "align",
    limit: str = "0",
) -> None:
    """PyMOL command: load UniProt-linked PDB structures and AlphaFold model."""
    load_uniprot_structures(
        query=query,
        include_pdb=bool(int(include_pdb)),
        include_alphafold=bool(int(include_alphafold)),
        align=bool(int(align)),
        alignment_method=alignment_method,
        limit=int(limit),
    )


def _run_in_thread(function, *args, **kwargs) -> None:
    thread = threading.Thread(target=function, args=args, kwargs=kwargs, daemon=True)
    thread.start()


def _open_qt_gui() -> bool:
    try:
        from pymol import plugins
        from pymol.Qt import QtWidgets
    except Exception as exc:
        print(f"[UniProt plugin] Qt GUI is not available: {exc}")
        return False

    parent = None
    try:
        parent = plugins.get_qtwindow()
    except Exception:
        pass

    dialog = QtWidgets.QDialog(parent)
    dialog.setWindowTitle("UniProt Structures")
    dialog.setMinimumWidth(360)

    query_edit = QtWidgets.QLineEdit(dialog)
    query_edit.setPlaceholderText("P06730 or IF4E_HUMAN")

    pdb_check = QtWidgets.QCheckBox("Load PDB structures", dialog)
    pdb_check.setChecked(True)
    af_check = QtWidgets.QCheckBox("Load AlphaFold model", dialog)
    af_check.setChecked(True)
    align_check = QtWidgets.QCheckBox("Align to AlphaFold", dialog)
    align_check.setChecked(True)

    method_combo = QtWidgets.QComboBox(dialog)
    method_combo.addItems(["align", "super", "cealign"])

    limit_spin = QtWidgets.QSpinBox(dialog)
    limit_spin.setRange(0, 10000)
    limit_spin.setValue(0)
    limit_spin.setSpecialValueText("0 = all")

    status_label = QtWidgets.QLabel("Ready", dialog)

    form = QtWidgets.QFormLayout()
    form.addRow("UniProt", query_edit)
    form.addRow("", pdb_check)
    form.addRow("", af_check)
    form.addRow("", align_check)
    form.addRow("Alignment", method_combo)
    form.addRow("PDB limit", limit_spin)
    form.addRow("", status_label)

    load_button = QtWidgets.QPushButton("Load", dialog)
    close_button = QtWidgets.QPushButton("Close", dialog)
    buttons = QtWidgets.QHBoxLayout()
    buttons.addWidget(load_button)
    buttons.addWidget(close_button)

    layout = QtWidgets.QVBoxLayout(dialog)
    layout.addLayout(form)
    layout.addLayout(buttons)

    def on_load() -> None:
        query = query_edit.text().strip()
        if not query:
            QtWidgets.QMessageBox.warning(dialog, "UniProt Structures", "Enter a UniProt accession or entry name.")
            return

        try:
            status_label.setText("Loading...")
            QtWidgets.QApplication.processEvents()
            summary = load_uniprot_structures(
                query=query,
                include_pdb=pdb_check.isChecked(),
                include_alphafold=af_check.isChecked(),
                align=align_check.isChecked(),
                alignment_method=method_combo.currentText(),
                limit=limit_spin.value(),
            )
            status_label.setText(f"Loaded {len(summary['loaded'])} object(s) in {summary['group']}")
        except Exception as exc:
            status_label.setText("Failed")
            QtWidgets.QMessageBox.critical(dialog, "UniProt Structures", str(exc))

    load_button.clicked.connect(on_load)
    close_button.clicked.connect(dialog.close)
    dialog.finished.connect(lambda _code: _OPEN_DIALOGS.remove(dialog) if dialog in _OPEN_DIALOGS else None)
    _OPEN_DIALOGS.append(dialog)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return True

def _open_gui(*_args, **_kwargs) -> None:
    if _open_qt_gui():
        return

    try:
        import tkinter as tk
        from tkinter import messagebox, ttk
    except Exception as exc:
        print(f"[UniProt plugin] Tkinter is not available: {exc}")
        return

    window = tk.Toplevel()
    window.title("UniProt Structures")
    window.resizable(False, False)

    query_var = tk.StringVar(value="")
    pdb_var = tk.BooleanVar(value=True)
    af_var = tk.BooleanVar(value=True)
    align_var = tk.BooleanVar(value=True)
    method_var = tk.StringVar(value="align")
    limit_var = tk.StringVar(value="0")
    status_var = tk.StringVar(value="Ready")

    frame = ttk.Frame(window, padding=12)
    frame.grid(row=0, column=0, sticky="nsew")

    ttk.Label(frame, text="UniProt accession or entry name").grid(row=0, column=0, columnspan=2, sticky="w")
    entry = ttk.Entry(frame, textvariable=query_var, width=34)
    entry.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(2, 8))
    entry.focus_set()

    ttk.Checkbutton(frame, text="Load PDB structures", variable=pdb_var).grid(row=2, column=0, columnspan=2, sticky="w")
    ttk.Checkbutton(frame, text="Load AlphaFold model", variable=af_var).grid(row=3, column=0, columnspan=2, sticky="w")
    ttk.Checkbutton(frame, text="Align to AlphaFold", variable=align_var).grid(row=4, column=0, columnspan=2, sticky="w")

    ttk.Label(frame, text="Alignment").grid(row=5, column=0, sticky="w", pady=(8, 0))
    ttk.Combobox(frame, textvariable=method_var, values=("align", "super", "cealign"), width=12, state="readonly").grid(
        row=5, column=1, sticky="e", pady=(8, 0)
    )

    ttk.Label(frame, text="PDB limit (0 = all)").grid(row=6, column=0, sticky="w", pady=(8, 0))
    ttk.Entry(frame, textvariable=limit_var, width=8).grid(row=6, column=1, sticky="e", pady=(8, 0))

    ttk.Label(frame, textvariable=status_var).grid(row=7, column=0, columnspan=2, sticky="w", pady=(10, 4))

    def on_load() -> None:
        query = query_var.get().strip()
        if not query:
            messagebox.showwarning("UniProt Structures", "Enter a UniProt accession or entry name.")
            return
        try:
            limit = int(limit_var.get().strip() or "0")
        except ValueError:
            messagebox.showwarning("UniProt Structures", "PDB limit must be a number.")
            return

        def worker() -> None:
            try:
                status_var.set("Loading...")
                summary = load_uniprot_structures(
                    query=query,
                    include_pdb=pdb_var.get(),
                    include_alphafold=af_var.get(),
                    align=align_var.get(),
                    alignment_method=method_var.get(),
                    limit=limit,
                )
                status_var.set(f"Loaded {len(summary['loaded'])} object(s) in {summary['group']}")
            except Exception as exc:
                status_var.set("Failed")
                messagebox.showerror("UniProt Structures", str(exc))

        _run_in_thread(worker)

    ttk.Button(frame, text="Load", command=on_load).grid(row=8, column=0, sticky="ew", pady=(4, 0))
    ttk.Button(frame, text="Close", command=window.destroy).grid(row=8, column=1, sticky="ew", pady=(4, 0), padx=(8, 0))
    window.bind("<Return>", lambda _event: on_load())


def uniprot_structures_gui() -> None:
    """PyMOL command: open the UniProt Structures GUI."""
    _open_gui()


def __init_plugin__(app=None) -> None:
    if cmd is not None:
        cmd.extend("uniprot_structures", uniprot_structures)
        cmd.extend("uniprot_structures_gui", uniprot_structures_gui)
    try:
        from pymol.plugins import addmenuitemqt

        addmenuitemqt("UniProt Structures", _open_gui)
    except Exception:
        try:
            app.rootmenu.addmenuitem("Plugin", "command", "UniProt Structures", label="UniProt Structures", command=_open_gui)
        except Exception:
            pass


if cmd is not None:
    cmd.extend("uniprot_structures", uniprot_structures)
    cmd.extend("uniprot_structures_gui", uniprot_structures_gui)
