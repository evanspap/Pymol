"""Align every model in one source PDB to one fixed target model.

The target may be a model in the source PDB, a model in another multi-model
PDB, or an external single-model PDB.  The target atom selection is cached
once. Source models are then selected, fitted, transformed, and written one at
a time, so memory use is limited to one model rather than the whole trajectory.

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


def parse_residues(text, flag, parser):
    result = set()
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token[1:]:
            split = token.index("-", 1)
            try:
                first, last = int(token[:split]), int(token[split + 1 :])
            except ValueError:
                parser.error(f"Invalid range '{token}' in {flag}.")
            if first > last:
                parser.error(f"Descending range '{token}' is invalid in {flag}.")
            result.update(range(first, last + 1))
        else:
            try:
                result.add(int(token))
            except ValueError:
                parser.error(f"Invalid residue '{token}' in {flag}.")
    if not result:
        parser.error(f"{flag} is empty.")
    return frozenset(result)


def arguments():
    parser = argparse.ArgumentParser(
        description="Align every source PDB model to one fixed target model and write one complete PDB."
    )
    parser.add_argument("-s", "--source", type=Path, required=True, help="Multi-model mobile/source PDB.")
    parser.add_argument(
        "-t", "--target", type=Path, help="Target PDB; defaults to --source when omitted."
    )
    parser.add_argument("--target-frame", type=int, default=1, help="1-based target model (default: 1).")
    parser.add_argument(
        "-o", "--output", type=Path, help="Complete aligned PDB; omit or use '-' for stdout; may equal source."
    )
    parser.add_argument("--atom-scope", choices=("ca", "backbone", "full"), default="ca")
    parser.add_argument("--chain", help="Target chain; default is all chains.")
    parser.add_argument("--source-chain", help="Source chain; defaults to --chain.")
    parser.add_argument("--residue", type=int)
    parser.add_argument("--source-residue", type=int)
    parser.add_argument("--residue-start", type=int)
    parser.add_argument("--residue-end", type=int)
    parser.add_argument("--source-residue-start", type=int)
    parser.add_argument("--source-residue-end", type=int)
    parser.add_argument("--residues", help="Target selection such as 37,50-60,100.")
    parser.add_argument("--source-residues", help="Source selection; defaults to --residues.")
    parser.add_argument(
        "--progress-every", type=int, default=100, help="Report every N processed models; 0 disables progress."
    )
    args = parser.parse_args()
    if args.target_frame < 1:
        parser.error("--target-frame must be at least 1.")
    if args.progress_every < 0:
        parser.error("--progress-every cannot be negative.")
    for prefix in ("", "source_"):
        single = getattr(args, prefix + "residue")
        start = getattr(args, prefix + "residue_start")
        end = getattr(args, prefix + "residue_end")
        listed = getattr(args, prefix + "residues")
        label = "--" + prefix.replace("_", "-")
        if single is not None and (start is not None or end is not None or listed is not None):
            parser.error(f"Do not combine {label}residue with its range/list options.")
        if listed is not None and (start is not None or end is not None):
            parser.error(f"Do not combine {label}residues with its range options.")
        if single is not None:
            setattr(args, prefix + "residue_start", single)
            setattr(args, prefix + "residue_end", single)
        setattr(
            args,
            prefix + "residue_numbers",
            parse_residues(listed, label + "residues", parser) if listed else None,
        )
    return args


def same_path(first, second):
    return os.path.normcase(os.path.abspath(first)) == os.path.normcase(os.path.abspath(second))


def parse_atom(line):
    try:
        residue = (
            line[21:22].decode("ascii"),
            int(line[22:26]),
            line[26:27].decode("ascii"),
        )
        name = line[12:16].decode("ascii").strip()
        altloc = line[16:17].decode("ascii")
        occupancy_text = line[54:60].strip()
        occupancy = float(occupancy_text) if occupancy_text else 0.0
        coordinates = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except (UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"Malformed PDB coordinate record: {line.rstrip()!r}") from exc
    return residue, name, altloc, occupancy, coordinates


def target_model_atoms(path, target_frame):
    """Read only the requested target model's standard ATOM records."""
    atoms = []
    saw_models = False
    current = None
    count = 0
    with open(path, "rb", buffering=1024 * 1024) as handle:
        for line in handle:
            record = line[:6].strip()
            if record == b"MODEL":
                saw_models = True
                count += 1
                current = count
                continue
            if not saw_models and record == b"ATOM":
                count = current = 1
            if current == target_frame and record == b"ATOM":
                atoms.append(parse_atom(line))
            if record == b"ENDMDL":
                if current == target_frame:
                    break
                current = None
    if not atoms:
        raise ValueError(f"Target frame {target_frame} was not found or contains no ATOM records: {path}")
    return atoms


def in_scope(atom_name, scope):
    if scope == "ca":
        return atom_name == "CA"
    if scope == "backbone":
        return atom_name in BACKBONE
    return True


def select_atoms(records, chain, start, end, numbers, scope):
    residues = OrderedDict()
    for residue, name, altloc, occupancy, coordinates in records:
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
        atoms = residues.setdefault(residue, {})
        previous = atoms.get(name)
        rank = (altloc in (" ", "A"), occupancy, altloc == " ")
        if previous is None or rank > previous[0]:
            atoms[name] = (rank, coordinates)
    return residues


def selection_settings(args):
    source_chain = args.source_chain if args.source_chain is not None else args.chain
    source_start = args.source_residue_start if args.source_residue_start is not None else args.residue_start
    source_end = args.source_residue_end if args.source_residue_end is not None else args.residue_end
    source_numbers = (
        args.source_residue_numbers if args.source_residue_numbers is not None else args.residue_numbers
    )
    explicit = any(
        value is not None
        for value in (
            args.chain,
            args.source_chain,
            args.residue_start,
            args.residue_end,
            args.residue_numbers,
            args.source_residue_start,
            args.source_residue_end,
            args.source_residue_numbers,
        )
    )
    return source_chain, source_start, source_end, source_numbers, explicit


def fit_model(target_selection, source_records, args, settings):
    source_chain, source_start, source_end, source_numbers, explicit = settings
    source_selection = select_atoms(
        source_records, source_chain, source_start, source_end, source_numbers, args.atom_scope
    )
    fixed, moving = [], []
    if explicit:
        if len(target_selection) != len(source_selection):
            raise ValueError(
                f"Target/source selections contain {len(target_selection)} and "
                f"{len(source_selection)} residues."
            )
        for fixed_residue, moving_residue in zip(target_selection.values(), source_selection.values()):
            for name in sorted(set(fixed_residue) & set(moving_residue)):
                fixed.append(fixed_residue[name][1])
                moving.append(moving_residue[name][1])
    else:
        fixed_map = {
            (residue, name): value[1]
            for residue, atoms in target_selection.items()
            for name, value in atoms.items()
        }
        moving_map = {
            (residue, name): value[1]
            for residue, atoms in source_selection.items()
            for name, value in atoms.items()
        }
        for key in sorted(set(fixed_map) & set(moving_map)):
            fixed.append(fixed_map[key])
            moving.append(moving_map[key])
    if len(fixed) < 3:
        raise ValueError("Fewer than three matching atoms were selected.")
    fitter = SVDSuperimposer()
    fitter.set(np.asarray(fixed, dtype=np.float64), np.asarray(moving, dtype=np.float64))
    fitter.run()
    rotation, translation = fitter.get_rotran()
    return rotation, translation, fitter.get_rms(), len(fixed)


def transform_line(line, rotation, translation):
    xyz = np.asarray((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    x, y, z = np.dot(xyz, rotation) + translation
    fields = f"{x:8.3f}{y:8.3f}{z:8.3f}".encode("ascii")
    if len(fields) != 24:
        raise ValueError("Transformed coordinates exceed PDB fixed-width fields.")
    return line[:30] + fields + line[54:]


def process_model(lines, target_selection, args, settings):
    source_records = [parse_atom(line) for line in lines if line[:6].strip() == b"ATOM"]
    rotation, translation, rmsd, matched = fit_model(target_selection, source_records, args, settings)
    transformed = [
        transform_line(line, rotation, translation) if line[:6].strip() in COORD_RECORDS else line
        for line in lines
    ]
    return transformed, rmsd, matched


def align_stream(source_path, output, target_selection, args, settings):
    model_lines = None
    model_count = 0
    rmsd_sum = 0.0
    rmsd_min = None
    rmsd_max = None
    matched = None
    saw_models = False
    with open(source_path, "rb", buffering=8 * 1024 * 1024) as source:
        for line in source:
            record = line[:6].strip()
            if record == b"MODEL":
                saw_models = True
                model_lines = [line]
                continue
            if model_lines is not None:
                model_lines.append(line)
                if record == b"ENDMDL":
                    model_count += 1
                    transformed, rmsd, matched = process_model(model_lines, target_selection, args, settings)
                    output.writelines(transformed)
                    rmsd_sum += rmsd
                    rmsd_min = rmsd if rmsd_min is None else min(rmsd_min, rmsd)
                    rmsd_max = rmsd if rmsd_max is None else max(rmsd_max, rmsd)
                    model_lines = None
                    if args.progress_every and model_count % args.progress_every == 0:
                        print(f"Aligned {model_count} models...", file=sys.stderr)
                continue
            if not saw_models and record in COORD_RECORDS:
                model_lines = [line]
            else:
                output.write(line)
        if model_lines is not None:
            model_count = 1
            transformed, rmsd, matched = process_model(model_lines, target_selection, args, settings)
            output.writelines(transformed)
            rmsd_sum = rmsd_min = rmsd_max = rmsd
    if model_count == 0:
        raise ValueError("Source PDB contains no complete models with coordinate records.")
    return model_count, matched, rmsd_min, rmsd_max, rmsd_sum / model_count


def write_output(source, output_path, target_selection, args, settings):
    if output_path is None or str(output_path) == "-":
        return align_stream(source, sys.stdout.buffer, target_selection, args, settings)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    in_place = same_path(source, output_path)
    temporary = None
    try:
        if in_place:
            descriptor, name = tempfile.mkstemp(prefix=f".{source.name}.", suffix=".tmp", dir=source.parent)
            os.close(descriptor)
            temporary = Path(name)
            destination = temporary
        else:
            destination = output_path
        with open(destination, "wb", buffering=8 * 1024 * 1024) as output:
            result = align_stream(source, output, target_selection, args, settings)
        if in_place:
            shutil.copystat(source, temporary)
            os.replace(temporary, source)
            temporary = None
        return result
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main():
    args = arguments()
    source = args.source.resolve()
    target = (args.target or args.source).resolve()
    if not source.is_file() or not target.is_file():
        sys.exit("Source and target must be existing PDB files.")
    try:
        settings = selection_settings(args)
        target_records = target_model_atoms(target, args.target_frame)
        target_selection = select_atoms(
            target_records,
            args.chain,
            args.residue_start,
            args.residue_end,
            args.residue_numbers,
            args.atom_scope,
        )
        if not target_selection:
            raise ValueError("The target atom selection is empty.")
        count, matched, rmsd_min, rmsd_max, rmsd_mean = write_output(
            source, args.output, target_selection, args, settings
        )
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: {exc}")
    destination = "stdout" if args.output is None or str(args.output) == "-" else args.output
    print(f"Aligned models: {count}; matched atoms/model: {matched}", file=sys.stderr)
    print(
        f"RMSD A: min={rmsd_min:.6f}, mean={rmsd_mean:.6f}, max={rmsd_max:.6f}",
        file=sys.stderr,
    )
    print(f"Complete aligned source PDB written to: {destination}", file=sys.stderr)


if __name__ == "__main__":
    main()
