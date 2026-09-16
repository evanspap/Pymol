"""
PDB Alignment Script - v1.4.0
Date: 9-16-2026

Align one selected model from a source PDB onto one selected model from a
target PDB. Only the requested model(s) are parsed for fitting. The complete
source PDB is then streamed to the output, with the selected source model's
coordinates replaced by their aligned coordinates; all other records and
models are copied unchanged.

If source and target are the same file, both requested models are extracted in
one scan. Omitting --output writes the resulting PDB to standard output.
--output may name a new file or the source itself; in-place output is performed
through a temporary file followed by an atomic replacement.

External dependency: Biopython (which installs its required NumPy dependency).
"""

import argparse
import os
import shutil
import sys
import tempfile
from io import StringIO
from pathlib import Path

from Bio.PDB import PDBParser, Superimposer


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
COORDINATE_RECORDS = {"ATOM", "HETATM"}


def parse_residue_list(spec, flag_name, parser):
    """Parse a specification such as 5,10-15,42 into sorted unique numbers."""
    residue_numbers = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token[1:]:
            dash_index = token.index("-", 1)
            try:
                start = int(token[:dash_index])
                end = int(token[dash_index + 1 :])
            except ValueError:
                parser.error(f"Invalid residue range '{token}' in {flag_name}.")
            if start > end:
                parser.error(f"Invalid residue range '{token}' in {flag_name}: start must be <= end.")
            residue_numbers.update(range(start, end + 1))
        else:
            try:
                residue_numbers.add(int(token))
            except ValueError:
                parser.error(f"Invalid residue number '{token}' in {flag_name}.")
    if not residue_numbers:
        parser.error(f"{flag_name} did not contain any residue numbers.")
    return sorted(residue_numbers)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Align one source PDB model onto one target PDB model, then write the complete "
            "source PDB with that source model transformed."
        )
    )
    parser.add_argument("--source", "-s", type=Path, required=True, help="PDB containing the model to move.")
    parser.add_argument("--target", "-t", type=Path, required=True, help="PDB containing the fixed model.")
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output PDB. Omit (or use '-') for stdout. May be the same path as --source for safe in-place output.",
    )
    parser.add_argument("--atom-scope", choices=["full", "backbone", "ca"], default="ca")
    parser.add_argument("--chain", default=None, help="Target chain ID; default is all chains.")
    parser.add_argument("--source-chain", default=None, help="Source chain ID; defaults to --chain.")
    parser.add_argument("--residue-start", type=int, default=None)
    parser.add_argument("--residue-end", type=int, default=None)
    parser.add_argument("--source-residue-start", type=int, default=None)
    parser.add_argument("--source-residue-end", type=int, default=None)
    parser.add_argument("--residue", type=int, default=None)
    parser.add_argument("--source-residue", type=int, default=None)
    parser.add_argument("--residues", default=None, help="Target residues, e.g. 3,10-15,42.")
    parser.add_argument("--source-residues", default=None, help="Source residues; defaults to --residues.")
    parser.add_argument("--target-frame", type=int, default=1, help="1-based target model ordinal (default: 1).")
    parser.add_argument(
        "--source-frame",
        type=int,
        default=None,
        help="1-based source model ordinal (defaults to --target-frame).",
    )
    args = parser.parse_args()

    if args.target_frame < 1:
        parser.error("--target-frame must be at least 1.")
    if args.source_frame is not None and args.source_frame < 1:
        parser.error("--source-frame must be at least 1.")

    if args.residue is not None:
        if args.residue_start is not None or args.residue_end is not None or args.residues is not None:
            parser.error("Use only one of --residue, --residues, or --residue-start/--residue-end.")
        args.residue_start = args.residue_end = args.residue
    if args.source_residue is not None:
        if (
            args.source_residue_start is not None
            or args.source_residue_end is not None
            or args.source_residues is not None
        ):
            parser.error(
                "Use only one of --source-residue, --source-residues, or "
                "--source-residue-start/--source-residue-end."
            )
        args.source_residue_start = args.source_residue_end = args.source_residue

    args.residue_numbers = None
    if args.residues is not None:
        if args.residue_start is not None or args.residue_end is not None:
            parser.error("Use either --residues or the target residue/range flags, not both.")
        args.residue_numbers = parse_residue_list(args.residues, "--residues", parser)

    args.source_residue_numbers = None
    if args.source_residues is not None:
        if args.source_residue_start is not None or args.source_residue_end is not None:
            parser.error("Use either --source-residues or the source residue/range flags, not both.")
        args.source_residue_numbers = parse_residue_list(args.source_residues, "--source-residues", parser)
    return args


def same_path(path_a, path_b):
    """Compare paths without requiring that an output path already exists."""
    return os.path.normcase(os.path.abspath(path_a)) == os.path.normcase(os.path.abspath(path_b))


def extract_model_texts(pdb_path, requested_frames):
    """
    Extract only requested 1-based model ordinals in one sequential scan.

    Files without MODEL/ENDMDL records are treated as a one-model PDB. MODEL
    serial values are not used; frame numbers refer to model order in the file.
    """
    requested = set(requested_frames)
    captured = {frame: [] for frame in requested}
    found = set()
    saw_model = False
    current_frame = None
    model_count = 0

    with open(pdb_path, "r", encoding="ascii", errors="replace", newline="") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                saw_model = True
                model_count += 1
                current_frame = model_count
                if current_frame in requested:
                    captured[current_frame].append(line)
                    found.add(current_frame)
                continue

            if not saw_model and record in COORDINATE_RECORDS:
                model_count = 1
                current_frame = 1
                if 1 in requested:
                    found.add(1)

            if current_frame in requested and record in COORDINATE_RECORDS:
                captured[current_frame].append(line)
            elif current_frame in requested and record in {"TER", "ENDMDL"}:
                captured[current_frame].append(line)

            if record == "ENDMDL":
                current_frame = None

    missing = requested - found
    if missing:
        available = model_count if model_count else 0
        missing_text = ", ".join(str(value) for value in sorted(missing))
        raise ValueError(
            f"Model frame(s) {missing_text} not found in {pdb_path}; file contains {available} model(s)."
        )
    return captured


def parse_model(model_text, structure_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, StringIO("".join(model_text)))
    try:
        return next(structure.get_models())
    except StopIteration as exc:
        raise ValueError(f"Selected {structure_id} model contains no parseable atoms.") from exc


def ordered_residues(model, chain_id, residue_start, residue_end, residue_numbers):
    wanted = set(residue_numbers) if residue_numbers is not None else None
    result = []
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for residue in chain:
            if residue.id[0] != " ":
                continue
            number = residue.id[1]
            if wanted is not None:
                if number not in wanted:
                    continue
            elif (residue_start is not None and number < residue_start) or (
                residue_end is not None and number > residue_end
            ):
                continue
            result.append(residue)
    return result


def atom_is_in_scope(atom_name, atom_scope):
    if atom_scope == "ca":
        return atom_name == "CA"
    if atom_scope == "backbone":
        return atom_name in BACKBONE_ATOM_NAMES
    return True


def select_atom_map(model, atom_scope):
    selected = {}
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            for atom in residue:
                name = atom.get_name()
                if atom_is_in_scope(name, atom_scope):
                    selected[(chain.id, residue.id, name)] = atom
    return selected


def build_exact_pairs(target_model, source_model, atom_scope):
    target_atoms = select_atom_map(target_model, atom_scope)
    source_atoms = select_atom_map(source_model, atom_scope)
    common = sorted(set(target_atoms) & set(source_atoms))
    return [target_atoms[key] for key in common], [source_atoms[key] for key in common]


def build_positional_pairs(target_residues, source_residues, atom_scope):
    if len(target_residues) != len(source_residues):
        raise ValueError(
            "Target and source selections contain different residue counts "
            f"({len(target_residues)} vs {len(source_residues)})."
        )
    fixed = []
    moving = []
    for target_residue, source_residue in zip(target_residues, source_residues):
        target_by_name = {atom.get_name(): atom for atom in target_residue}
        source_by_name = {atom.get_name(): atom for atom in source_residue}
        for name in sorted(set(target_by_name) & set(source_by_name)):
            if atom_is_in_scope(name, atom_scope):
                fixed.append(target_by_name[name])
                moving.append(source_by_name[name])
    return fixed, moving


def transform_coordinate_line(line, rotation, translation):
    """Transform ATOM/HETATM coordinates while preserving the rest of the record."""
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"Cannot parse PDB coordinates from record: {line.rstrip()}") from exc

    new_x = x * rotation[0][0] + y * rotation[1][0] + z * rotation[2][0] + translation[0]
    new_y = x * rotation[0][1] + y * rotation[1][1] + z * rotation[2][1] + translation[1]
    new_z = x * rotation[0][2] + y * rotation[1][2] + z * rotation[2][2] + translation[2]
    coordinates = f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"
    if len(coordinates) != 24:
        raise ValueError("Transformed coordinates exceed the fixed-width PDB coordinate fields.")
    return line[:30] + coordinates + line[54:]


def stream_transformed_source(source_path, output_handle, source_frame, rotation, translation):
    """Copy the entire source, transforming coordinates only in source_frame."""
    saw_model = False
    current_frame = None
    model_count = 0
    transformed_records = 0

    with open(source_path, "r", encoding="ascii", errors="replace", newline="") as source_handle:
        for line in source_handle:
            record = line[:6].strip()
            if record == "MODEL":
                saw_model = True
                model_count += 1
                current_frame = model_count
                output_handle.write(line)
                continue
            if not saw_model and record in COORDINATE_RECORDS:
                model_count = 1
                current_frame = 1
            if current_frame == source_frame and record in COORDINATE_RECORDS:
                output_handle.write(transform_coordinate_line(line, rotation, translation))
                transformed_records += 1
            else:
                output_handle.write(line)
            if record == "ENDMDL":
                current_frame = None
    if transformed_records == 0:
        raise ValueError(f"No coordinate records were transformed for source frame {source_frame}.")
    return transformed_records


def write_output(source_path, output_path, source_frame, rotation, translation):
    if output_path is None or str(output_path) == "-":
        return stream_transformed_source(source_path, sys.stdout, source_frame, rotation, translation)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    in_place = same_path(source_path, output_path)
    temporary_path = None
    try:
        if in_place:
            fd, temporary_name = tempfile.mkstemp(
                prefix=f".{source_path.name}.", suffix=".v1.4.tmp", dir=source_path.parent
            )
            os.close(fd)
            temporary_path = Path(temporary_name)
            destination = temporary_path
        else:
            destination = output_path

        with open(destination, "w", encoding="ascii", errors="strict", newline="") as output_handle:
            transformed = stream_transformed_source(
                source_path, output_handle, source_frame, rotation, translation
            )

        if in_place:
            shutil.copystat(source_path, temporary_path)
            os.replace(temporary_path, source_path)
            temporary_path = None
        return transformed
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass


def main():
    args = parse_args()
    source_path = args.source.resolve()
    target_path = args.target.resolve()
    if not source_path.is_file():
        sys.exit(f"Source PDB not found: {source_path}")
    if not target_path.is_file():
        sys.exit(f"Target PDB not found: {target_path}")

    source_frame = args.source_frame if args.source_frame is not None else args.target_frame
    try:
        if same_path(source_path, target_path):
            texts = extract_model_texts(source_path, {source_frame, args.target_frame})
            source_model = parse_model(texts[source_frame], "source")
            target_model = parse_model(texts[args.target_frame], "target")
        else:
            source_text = extract_model_texts(source_path, {source_frame})[source_frame]
            target_text = extract_model_texts(target_path, {args.target_frame})[args.target_frame]
            source_model = parse_model(source_text, "source")
            target_model = parse_model(target_text, "target")

        source_chain = args.source_chain if args.source_chain is not None else args.chain
        source_start = args.source_residue_start if args.source_residue_start is not None else args.residue_start
        source_end = args.source_residue_end if args.source_residue_end is not None else args.residue_end
        source_numbers = (
            args.source_residue_numbers if args.source_residue_numbers is not None else args.residue_numbers
        )
        selection_requested = any(
            value is not None
            for value in (
                args.chain,
                args.residue_start,
                args.residue_end,
                args.residue_numbers,
                args.source_chain,
                args.source_residue_start,
                args.source_residue_end,
                args.source_residue_numbers,
            )
        )

        if selection_requested:
            target_residues = ordered_residues(
                target_model, args.chain, args.residue_start, args.residue_end, args.residue_numbers
            )
            source_residues = ordered_residues(
                source_model, source_chain, source_start, source_end, source_numbers
            )
            if not target_residues or not source_residues:
                raise ValueError("No residues found in the requested target/source selection.")
            fixed_atoms, moving_atoms = build_positional_pairs(
                target_residues, source_residues, args.atom_scope
            )
        else:
            fixed_atoms, moving_atoms = build_exact_pairs(target_model, source_model, args.atom_scope)

        if len(fixed_atoms) < 3:
            raise ValueError(
                "Fewer than three matching atoms were found; a 3D superposition requires at least three."
            )

        superimposer = Superimposer()
        superimposer.set_atoms(fixed_atoms, moving_atoms)
        rotation, translation = superimposer.rotran
        transformed_records = write_output(
            source_path, args.output, source_frame, rotation, translation
        )
    except (OSError, ValueError) as exc:
        sys.exit(f"Error: {exc}")

    destination = "standard output" if args.output is None or str(args.output) == "-" else str(args.output)
    print(f"Target frame: {args.target_frame}", file=sys.stderr)
    print(f"Source frame transformed: {source_frame}", file=sys.stderr)
    print(f"Matched atoms used for fit: {len(fixed_atoms)}", file=sys.stderr)
    print(f"Fit RMSD: {superimposer.rms:.6f} Angstrom", file=sys.stderr)
    print(f"Coordinate records transformed: {transformed_records}", file=sys.stderr)
    print(f"Complete source PDB written to: {destination}", file=sys.stderr)


if __name__ == "__main__":
    main()
