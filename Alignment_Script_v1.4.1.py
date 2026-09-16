"""Fast selected-model PDB alignment with complete-source streaming output.

Only the requested model records and selected atoms are parsed.  The fitted
source model is rewritten while the rest of the source PDB is copied in large
binary chunks.  Omit --output for stdout; --output may equal --source.

Dependencies: Biopython and its required NumPy dependency.
"""

import argparse
import os
import shutil
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer


BACKBONE = {"N", "CA", "C", "O", "OXT"}
COORD_RECORDS = {b"ATOM", b"HETATM"}


def residue_spec(text, flag, parser):
    values = set()
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item[1:]:
            at = item.index("-", 1)
            try:
                first, last = int(item[:at]), int(item[at + 1 :])
            except ValueError:
                parser.error(f"Invalid range '{item}' in {flag}.")
            if first > last:
                parser.error(f"Invalid descending range '{item}' in {flag}.")
            values.update(range(first, last + 1))
        else:
            try:
                values.add(int(item))
            except ValueError:
                parser.error(f"Invalid residue '{item}' in {flag}.")
    if not values:
        parser.error(f"{flag} is empty.")
    return frozenset(values)


def arguments():
    p = argparse.ArgumentParser(description="Rapidly align one PDB model and stream the complete source PDB.")
    p.add_argument("-s", "--source", type=Path, required=True)
    p.add_argument("-t", "--target", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, help="Output path; omit or use '-' for stdout; may equal source.")
    p.add_argument("--source-frame", type=int, help="1-based source model; defaults to target frame.")
    p.add_argument("--target-frame", type=int, default=1, help="1-based target model (default: 1).")
    p.add_argument("--atom-scope", choices=("ca", "backbone", "full"), default="ca")
    p.add_argument("--chain")
    p.add_argument("--source-chain")
    p.add_argument("--residue", type=int)
    p.add_argument("--source-residue", type=int)
    p.add_argument("--residue-start", type=int)
    p.add_argument("--residue-end", type=int)
    p.add_argument("--source-residue-start", type=int)
    p.add_argument("--source-residue-end", type=int)
    p.add_argument("--residues")
    p.add_argument("--source-residues")
    a = p.parse_args()
    if a.target_frame < 1 or (a.source_frame is not None and a.source_frame < 1):
        p.error("Frame numbers are 1-based and must be positive.")
    for prefix in ("", "source_"):
        single = getattr(a, prefix + "residue")
        start = getattr(a, prefix + "residue_start")
        end = getattr(a, prefix + "residue_end")
        listed = getattr(a, prefix + "residues")
        if single is not None and (start is not None or end is not None or listed is not None):
            p.error(f"Do not combine --{prefix.replace('_','-')}residue with range/list options.")
        if listed is not None and (start is not None or end is not None):
            p.error(f"Do not combine --{prefix.replace('_','-')}residues with range options.")
        if single is not None:
            setattr(a, prefix + "residue_start", single)
            setattr(a, prefix + "residue_end", single)
        setattr(a, prefix + "residue_numbers", residue_spec(listed, "--" + prefix.replace("_", "-") + "residues", p) if listed else None)
    return a


def same_path(a, b):
    return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def atom_from_line(line):
    """Return (residue key, atom name, alternate location, occupancy, xyz)."""
    try:
        chain = line[21:22].decode("ascii")
        number = int(line[22:26])
        insertion = line[26:27].decode("ascii")
        name = line[12:16].decode("ascii").strip()
        altloc = line[16:17].decode("ascii")
        occupancy_text = line[54:60].strip()
        occupancy = float(occupancy_text) if occupancy_text else 0.0
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except (UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"Malformed coordinate record: {line.rstrip()!r}") from exc
    return (chain, number, insertion), name, altloc, occupancy, xyz


def scan_models(path, wanted):
    """Scan once, retaining atom records and byte spans only for wanted model ordinals."""
    wanted = set(wanted)
    atoms = {number: [] for number in wanted}
    spans = {}
    completed = set()
    saw_model = False
    current = None
    count = 0
    start = None
    with open(path, "rb", buffering=1024 * 1024) as handle:
        while True:
            offset = handle.tell()
            line = handle.readline()
            if not line:
                break
            record = line[:6].strip()
            if record == b"MODEL":
                saw_model = True
                count += 1
                current, start = count, offset
            elif not saw_model and record in COORD_RECORDS:
                count, current, start = 1, 1, 0
            if current in wanted and record == b"ATOM":
                atoms[current].append(atom_from_line(line))
            if record == b"ENDMDL":
                if current in wanted:
                    spans[current] = (start, handle.tell())
                    completed.add(current)
                current = None
                if completed == wanted:
                    break
        if not saw_model and 1 in wanted and atoms[1]:
            spans[1] = (0, handle.seek(0, os.SEEK_END))
            completed.add(1)
    missing = wanted - completed
    if missing:
        raise ValueError(f"Frame(s) {sorted(missing)} not found in {path}; detected {count} model(s).")
    return atoms, spans


def in_scope(name, scope):
    return name == "CA" if scope == "ca" else name in BACKBONE if scope == "backbone" else True


def selected_residues(records, chain, start, end, numbers, scope):
    """Ordered residues with one preferred coordinate per atom name."""
    residues = OrderedDict()
    for residue, name, altloc, occupancy, xyz in records:
        if chain is not None and residue[0] != chain:
            continue
        number = residue[1]
        if numbers is not None:
            if number not in numbers:
                continue
        elif (start is not None and number < start) or (end is not None and number > end):
            continue
        if not in_scope(name, scope):
            continue
        atom_map = residues.setdefault(residue, {})
        old = atom_map.get(name)
        rank = (altloc in (" ", "A"), occupancy, altloc == " ")
        if old is None or rank > old[0]:
            atom_map[name] = (rank, xyz)
    return residues


def fit_coordinates(target_records, source_records, a):
    source_chain = a.source_chain if a.source_chain is not None else a.chain
    source_start = a.source_residue_start if a.source_residue_start is not None else a.residue_start
    source_end = a.source_residue_end if a.source_residue_end is not None else a.residue_end
    source_numbers = a.source_residue_numbers if a.source_residue_numbers is not None else a.residue_numbers
    explicit = any(value is not None for value in (
        a.chain, a.source_chain, a.residue_start, a.residue_end, a.residue_numbers,
        a.source_residue_start, a.source_residue_end, a.source_residue_numbers,
    ))
    target = selected_residues(target_records, a.chain, a.residue_start, a.residue_end, a.residue_numbers, a.atom_scope)
    source = selected_residues(source_records, source_chain, source_start, source_end, source_numbers, a.atom_scope)
    fixed, moving = [], []
    if explicit:
        if len(target) != len(source):
            raise ValueError(f"Target/source selections contain {len(target)} and {len(source)} residues.")
        residue_pairs = zip(target.values(), source.values())
        for target_atoms, source_atoms in residue_pairs:
            for name in sorted(set(target_atoms) & set(source_atoms)):
                fixed.append(target_atoms[name][1]); moving.append(source_atoms[name][1])
    else:
        target_atoms = {(residue, name): value[1] for residue, names in target.items() for name, value in names.items()}
        source_atoms = {(residue, name): value[1] for residue, names in source.items() for name, value in names.items()}
        for key in sorted(set(target_atoms) & set(source_atoms)):
            fixed.append(target_atoms[key]); moving.append(source_atoms[key])
    if len(fixed) < 3:
        raise ValueError("Fewer than three matching atoms were selected.")
    fitter = SVDSuperimposer()
    fitter.set(np.asarray(fixed, dtype=np.float64), np.asarray(moving, dtype=np.float64))
    fitter.run()
    rotation, translation = fitter.get_rotran()
    return rotation, translation, fitter.get_rms(), len(fixed)


def transformed_line(line, rotation, translation):
    try:
        xyz = np.asarray((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    except ValueError as exc:
        raise ValueError(f"Malformed coordinate record: {line.rstrip()!r}") from exc
    x, y, z = np.dot(xyz, rotation) + translation
    fields = f"{x:8.3f}{y:8.3f}{z:8.3f}".encode("ascii")
    if len(fields) != 24:
        raise ValueError("Transformed coordinates exceed PDB fixed-width fields.")
    return line[:30] + fields + line[54:]


def copy_bytes(source, destination, count, block=8 * 1024 * 1024):
    while count:
        data = source.read(min(block, count))
        if not data:
            raise OSError("Unexpected end of source PDB.")
        destination.write(data)
        count -= len(data)


def rewrite(source_path, destination, span, rotation, translation):
    start, end = span
    changed = 0
    with open(source_path, "rb", buffering=8 * 1024 * 1024) as source:
        copy_bytes(source, destination, start)
        while source.tell() < end:
            line = source.readline()
            if line[:6].strip() in COORD_RECORDS:
                line = transformed_line(line, rotation, translation)
                changed += 1
            destination.write(line)
        shutil.copyfileobj(source, destination, length=8 * 1024 * 1024)
    return changed


def output_pdb(source, output, span, rotation, translation):
    if output is None or str(output) == "-":
        return rewrite(source, sys.stdout.buffer, span, rotation, translation)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    inplace = same_path(source, output)
    temporary = None
    try:
        if inplace:
            fd, name = tempfile.mkstemp(prefix=f".{source.name}.", suffix=".tmp", dir=source.parent)
            os.close(fd); temporary = Path(name); destination = temporary
        else:
            destination = output
        with open(destination, "wb", buffering=8 * 1024 * 1024) as handle:
            changed = rewrite(source, handle, span, rotation, translation)
        if inplace:
            shutil.copystat(source, temporary)
            os.replace(temporary, source)
            temporary = None
        return changed
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    a = arguments()
    source, target = a.source.resolve(), a.target.resolve()
    if not source.is_file() or not target.is_file():
        sys.exit("Source and target must both be existing PDB files.")
    source_frame = a.source_frame if a.source_frame is not None else a.target_frame
    try:
        if same_path(source, target):
            records, spans = scan_models(source, {source_frame, a.target_frame})
            source_records, target_records = records[source_frame], records[a.target_frame]
        else:
            source_data, spans = scan_models(source, {source_frame})
            target_data, _ = scan_models(target, {a.target_frame})
            source_records, target_records = source_data[source_frame], target_data[a.target_frame]
        rotation, translation, rmsd, matched = fit_coordinates(target_records, source_records, a)
        changed = output_pdb(source, a.output, spans[source_frame], rotation, translation)
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: {exc}")
    destination = "stdout" if a.output is None or str(a.output) == "-" else a.output
    print(f"Matched atoms: {matched}; fit RMSD: {rmsd:.6f} A", file=sys.stderr)
    print(f"Transformed source frame {source_frame}: {changed} coordinate records", file=sys.stderr)
    print(f"Complete source PDB written to: {destination}", file=sys.stderr)


if __name__ == "__main__":
    main()
