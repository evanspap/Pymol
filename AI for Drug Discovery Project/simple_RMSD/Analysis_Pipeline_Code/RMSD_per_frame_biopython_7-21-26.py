"""Biopython RMSD calculator for target-vs-reference frame comparisons and matrix-style RMSD tables.
### Version 1.1 (Subject to Updates and Change)
## Author(s): Yagna Devarakonda, Research Volunteer/Assistant, Stony Brook University, NY, USA
## Dr. Evangelos Papadopoulos, Professor, Dana-Farber Cancer Institute 
## Date: 7-21-26

This file is the full-featured version of the RMSD workflow.
It builds on the smaller, barebones helper file named
RMSD_calculation_function_7-19-26.py, which is meant for simple direct RMSD
calculations in the terminal.

**DISCLAIMER**
=============
**Generative AI was used in writing this code, and the authors take
full responsibility for the final content. Please modify code for personal use
if necessary. The authors are not responsible for any errors or issues that
may arise from using this code which originate from team specific
modifications, dependencies, or other factors outside the authors' control.**




PowerShell usage (copy one scenario at a time):

Command-order convention used below (for clarity):
    1) reference first
    2) targets second

In other words, when both are used, write:
    --reference-frame ... --targets ...
or
    --reference-frames ... --targets ...

Very important output rule:
    --reference-frame + --targets
        = simple target-vs-reference list CSV
        = columns: Frame, RMSD_Angstrom

    --reference-frames + --targets
        = reference-target matrix CSV
        = reference frames on rows (vertical axis)
        = target frames on columns (horizontal axis)

If you want Excel-style row/column orientation, always use:
    --reference-frames ... --targets ...
even if you only have one reference frame.

Option A: Run from the folder that contains this script
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>"

Option B: Run from anywhere by using the full script path
    python "<path-to-this-script>\RMSD_per_frame_biopython_7-21-26.py" --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>"

Sample commands:
    # Default: compare every frame to frame 1 using all atoms
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>"

    # Backbone-only comparison:
    # Uses only protein backbone atoms (N, CA, C, O, and OXT when present)
    # and ignores side-chain atoms. This is useful when you want overall
    # fold/motion trends with less side-chain noise.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --atom-scope backbone

    # Use a custom reference and explicit target frames.
    # This creates a simple list CSV, not a matrix.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 5 --targets 12,48

    # 1vN example (reference 1 vs targets 2,4,5,9).
    # This is still simple list output, not matrix output.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 2,4,5,9

    # Use explicit targets plus specific atoms (reference first, targets second)
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 10-25 --atoms CA,CB,N

    # Use only specific residue numbers (with optional atom filter)
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 10-25 --residues 5,10-20 --atoms CA

Atom-name legend (what these labels mean):
    N   = backbone amide nitrogen
    CA  = alpha carbon (backbone)
    C   = backbone carbonyl carbon
    O   = backbone carbonyl oxygen
    OXT = terminal carboxylate oxygen (often only at C-terminus)
    CB  = beta carbon (first side-chain carbon for most amino acids)

Notes:
    - "Backbone-only" mode uses N, CA, C, O, and OXT.
    - Example custom list --atoms CA,CB,N means:
        alpha carbon + beta carbon + backbone nitrogen.
    - Atom names must match PDB atom-name fields used in your structure.

Scenario 1: Compare all frames to frame 1 (easy default)
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>"

Scenario 1b: Same thing, but say it explicitly
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all

Scenario 2: Backbone-only RMSD relative to the first frame
    # "Backbone-only" means only N, CA, C, O, OXT are used for RMSD.
    # Side-chain atoms are excluded.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --atom-scope backbone

Scenario 3: Compare a reference frame against a target range
    # Output type: simple list CSV
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 10-25

Scenario 4: Compare a reference frame against specific target frame numbers
    # Output type: simple list CSV
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 5,12,48,77

Scenario 4b: Same idea, but targets include ranges too
    # Output type: simple list CSV
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 1-10,20,25-30

Scenario 5: Use a different reference frame and custom output folder
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --reference-frame 10 --output-dir "<path-to-results-folder>"

Scenario 6: Compare specific atoms with explicit targets
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 1-10,20 --atoms CA,CB,N

Scenario 7: Build an all-vs-all RMSD matrix for all frames you want to include
    # Use this when you want every chosen frame compared against every other chosen frame.
    # If your file has 1,000 frames, this means all 1,000 can be included.
    # The --compare-all flag means "use every frame in the PDB file."
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all --all-vs-all

Scenario 8: Explicit 1vN using targets list
    # Output type: simple list CSV
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 2,4,5,9

Scenario 9: Nvs1 matrix with references on rows and target on columns
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 2,4,5,9,12 --targets 45

Scenario 10: 1vN matrix with one reference on rows and many targets on columns
    # Use this when you want frame 3 on the vertical axis and many targets on the horizontal axis.
    # Even though it is conceptually "1 vs N", use --reference-frames for matrix output.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 3 --targets 10,25,50,100,250

Scenario 11: Compare only a residue subset
    # Use this when you want RMSD from specific residue numbers only.
    # You can combine --residues with --atom-scope or --atoms.
    python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 1-10 --targets 5 --residues 1-36 --atoms CA

Important path note for teams cloning from GitHub:
    Replace each placeholder with a path from your own machine.
    If you want to avoid typing the full script path, open PowerShell in the script folder first.

What this script does:
- This is the full-featured RMSD script.
- The smaller helper file RMSD_calculation_function_7-19-26.py is the
    barebones version if you only want one simple RMSD value in the terminal.
- Reads a multi-frame PDB trajectory and extracts each model/frame.
- Uses Biopython to parse atoms in the selected reference frame and selected target frames.
- Supports atom filtering by scope (`full` or `backbone`) or by explicit atom names (`--atoms`).
- Supports residue-number filtering with `--residues` (example: `1-36,45,50-60`).
- Compares each selected target frame against one selected reference frame.
- Supports 1v1, 1vN, 2v1, 5v1, and larger Nv1 comparisons as simple list CSV output.
- Supports NvsM reference-target matrix output with reference frames as CSV rows and target frames as CSV columns (`--reference-frames` + `--targets`).
- Optionally supports all-vs-all matrix generation for selected frames (`--all-vs-all`).
- Requires at least two frames in the input file overall.
- Writes results to CSV in one of two formats:
    - Nv1 mode (default): header `Frame,RMSD_Angstrom` and one row per selected target frame.
    - All-vs-all mode: square matrix CSV where rows/columns are selected frame IDs.
"""

import argparse
import csv
import sys
from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
DEFAULT_OUTPUT_NAME = "rmsd_per_frame.csv"


def parse_args():
    """Parse command-line arguments for a portable per-frame RMSD workflow.

    The goal is to let another team run this script without editing the code.
    They only need to supply their own input PDB and choose a reference frame
    plus the target frames to compare against that reference.
    """
    parser = argparse.ArgumentParser(
        description="Compute target-vs-reference RMSD values from a multi-frame PDB file using Biopython.",
        epilog=(
            "Examples:\n"
            "  python RMSD_per_frame_biopython_7-21-26.py --input <path-to-your-pdb-file>\n"
            "  python RMSD_per_frame_biopython_7-21-26.py --input <path-to-your-pdb-file> --atom-scope backbone\n"
            "  python RMSD_per_frame_biopython_7-21-26.py --input <path-to-your-pdb-file> --reference-frame 5 --targets 2,4,5,9 --output-dir <path-to-results-folder>"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Path to the multi-frame PDB file.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=None,
        help="Directory for the output CSV (defaults to the input file folder).",
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help="Output CSV filename.",
    )
    parser.add_argument(
        "--atom-scope",
        choices=["full", "backbone"],
        default="full",
        help="Choose whether to use all atoms or backbone atoms only. Ignored if --atoms is used.",
    )
    parser.add_argument(
        "--atoms",
        type=str,
        default=None,
        help="Optional comma-separated atom names to use (example: CA,CB,N). If given, this overrides --atom-scope.",
    )
    parser.add_argument(
        "--residues",
        type=str,
        default=None,
        help=(
            "Optional residue-number list using commas/ranges "
            "(example: 1-36,45,50-60). "
            "Only atoms from these residue sequence numbers are kept."
        ),
    )
    parser.add_argument(
        "--reference-frame",
        type=int,
        default=1,
        help="Reference frame number used for RMSD comparison (1-based index).",
    )
    parser.add_argument(
        "--reference-frames",
        type=str,
        default=None,
        help=(
            "Optional reference-frame list for NvsM matrix mode "
            "(example: 2,4,5,9 or 2-10,15). "
            "Use together with --targets to produce a single CSV where references are rows and targets are columns."
        ),
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Explicitly compare every frame in the file (same as the default behavior).",
    )
    parser.add_argument(
        "--frames",
        type=str,
        default=None,
        help="Optional frame list using commas and ranges (example: 5,12,48 or 1-10,20,25-30).",
    )
    parser.add_argument(
        "--targets",
        type=str,
        default=None,
        help=(
            "Optional explicit target-frame list for Nv1 mode only "
            "(example: 2,4,5,9 or 2-10,15). "
            "When provided, this overrides --frames/--frame-start/--frame-end/--compare-all in Nv1 mode."
        ),
    )
    parser.add_argument(
        "--frame-start",
        type=int,
        default=None,
        help="Optional first frame number to include (1-based index).",
    )
    parser.add_argument(
        "--frame-end",
        type=int,
        default=None,
        help="Optional last frame number to include (1-based index, inclusive).",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional: only inspect the first N frames from the file before any other selection.",
    )
    parser.add_argument(
        "--all-vs-all",
        action="store_true",
        help="Compute pairwise RMSD matrix across selected frames instead of Nv1 target-vs-reference output.",
    )
    return parser.parse_args()


def split_frames_from_pdb(pdb_file):
    """Split a multi-frame PDB trajectory into separate frame text blocks.

    The file is expected to contain multiple frames separated by END records.
    Only useful PDB records are kept so each frame can be parsed cleanly by Biopython.
    """
    frames = []
    current = []

    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            # END marks the end of one frame in many trajectory-style PDB files.
            if line.startswith("END"):
                if current:
                    frames.append("".join(current + ["END\n"]))
                    current = []
                continue

            # Keep only records that are useful for structural parsing.
            if line.startswith(("ATOM", "HETATM", "TER", "CRYST1", "HEADER", "REMARK", "TITLE", "MODEL", "ENDMDL", "SEQRES", "COMPND")):
                current.append(line)

    # Capture the final frame if the file did not end with END.
    if current:
        frames.append("".join(current + ["END\n"]))

    if len(frames) < 1:
        raise ValueError("No frames were found in the input PDB file.")

    return frames


def parse_frame_selection(args, total_frames):
    """Convert user selection options into a sorted list of frame numbers.

    Selection rules:
    - If the user gives --frames, those exact frame numbers are used.
    - If the user gives --frame-start / --frame-end, that contiguous range is used.
    - If neither is given, all frames are used.

    Frame numbers are 1-based because that is easier for humans to read and
    matches the CSV output labels.
    """
    selected = set()

    if args.compare_all:
        # Explicitly include every frame. This keeps the default behavior obvious.
        selected.update(range(1, total_frames + 1))

    if args.frames:
        # Allow a simple comma-separated list such as "5,12,48" or a mix of
        # numbers and ranges such as "1-10,20,25-30".
        for part in args.frames.split(","):
            item = part.strip()
            if not item:
                continue
            if "-" in item:
                start_text, end_text = item.split("-", 1)
                try:
                    start = int(start_text.strip())
                    end = int(end_text.strip())
                except ValueError as exc:
                    raise ValueError(f"Invalid frame range in --frames: {item}") from exc
                if start < 1 or end < 1:
                    raise ValueError("Frame numbers must be 1 or larger.")
                if start > end:
                    raise ValueError(f"Invalid frame range in --frames: {item}")
                selected.update(range(start, end + 1))
                continue
            try:
                selected.add(int(item))
            except ValueError as exc:
                raise ValueError(f"Invalid frame number in --frames: {item}") from exc

    if args.frame_start is not None or args.frame_end is not None:
        start = args.frame_start if args.frame_start is not None else 1
        end = args.frame_end if args.frame_end is not None else total_frames
        if start < 1 or end < 1:
            raise ValueError("Frame numbers must be 1 or larger.")
        if start > end:
            raise ValueError("--frame-start cannot be greater than --frame-end.")
        selected.update(range(start, end + 1))

    if not selected:
        # No user-specific selection, so use every frame.
        selected.update(range(1, total_frames + 1))

    # Respect the maximum frame cap after the initial file read.
    if args.max_frames is not None:
        selected = {frame_number for frame_number in selected if frame_number <= args.max_frames}

    return sorted(selected)


def parse_frame_list_text(frame_text, total_frames):
    """Parse a user frame list text like '2,4,5,9' or '2-10,15'.

    This helper is used for --targets so users can explicitly define the exact
    frames they want compared against a single reference frame.
    """
    selected = set()

    if frame_text is None or not str(frame_text).strip():
        return []

    for part in str(frame_text).split(","):
        item = part.strip()
        if not item:
            continue

        if "-" in item:
            start_text, end_text = item.split("-", 1)
            try:
                start = int(start_text.strip())
                end = int(end_text.strip())
            except ValueError as exc:
                raise ValueError(f"Invalid range in frame list: {item}") from exc

            if start < 1 or end < 1 or start > end:
                raise ValueError(f"Invalid range in frame list: {item}")

            selected.update(range(start, end + 1))
            continue

        try:
            selected.add(int(item))
        except ValueError as exc:
            raise ValueError(f"Invalid frame number in frame list: {item}") from exc

    out_of_range = [n for n in selected if n > total_frames]
    if out_of_range:
        raise ValueError(
            f"Frame numbers out of range (max frame is {total_frames}): {sorted(out_of_range)}"
        )

    return sorted(selected)


def parse_atom_selection(atom_text):
    """Convert a comma-separated atom list into a clean set of atom names.

    Example input:
    - "CA,CB,N"

    Output:
    - {"CA", "CB", "N"}

    Names are normalized to uppercase so users do not have to worry about case.
    """
    if atom_text is None or not str(atom_text).strip():
        return None

    atom_names = set()
    for part in str(atom_text).split(","):
        name = part.strip().upper()
        if not name:
            continue
        atom_names.add(name)

    if not atom_names:
        return None

    return atom_names


def parse_residue_selection(residue_text):
    """Convert a residue selection text into a set of residue sequence numbers.

    Example inputs:
    - "10,15,20"
    - "1-36,45,50-60"
    """
    if residue_text is None or not str(residue_text).strip():
        return None

    selected = set()
    for part in str(residue_text).split(","):
        item = part.strip()
        if not item:
            continue

        if "-" in item:
            start_text, end_text = item.split("-", 1)
            try:
                start = int(start_text.strip())
                end = int(end_text.strip())
            except ValueError as exc:
                raise ValueError(f"Invalid residue range in --residues: {item}") from exc

            if start > end:
                raise ValueError(f"Invalid residue range in --residues: {item}")

            selected.update(range(start, end + 1))
            continue

        try:
            selected.add(int(item))
        except ValueError as exc:
            raise ValueError(f"Invalid residue number in --residues: {item}") from exc

    if not selected:
        return None

    return selected


def atom_key(atom):
    """Build a stable sort key so atoms are matched consistently across frames.

    Why this matters:
    If atom order is different between frames, the coordinates must still be lined
    up in the same logical order before RMSD can be calculated safely.
    """
    residue = atom.get_parent()
    chain = residue.get_parent() if residue is not None else None
    residue_id = residue.get_id() if residue is not None else ("", 0, "")
    return (
        chain.get_id() if chain is not None else "",
        residue.get_resname() if residue is not None else "",
        residue_id[0],
        residue_id[1],
        residue_id[2],
        atom.get_name().strip().upper(),
        atom.get_altloc() or "",
    )


def parse_frame_with_biopython(
    frame_text,
    frame_index,
    parser,
    atom_scope,
    selected_atoms=None,
    residue_numbers=None,
):
    """Parse one frame with Biopython and return sorted atom keys and coordinates.

    Sorting by atom identity ensures the coordinates line up across frames even if
    the input atom order is inconsistent.

    If selected_atoms is provided, only those atom names are kept.
    If not, the script falls back to atom_scope (full or backbone).
    """
    # Biopython reads the frame text as if it were a small PDB file.
    structure = parser.get_structure(f"frame_{frame_index}", StringIO(frame_text))
    model = next(structure.get_models())

    atoms = []
    for atom in model.get_atoms():
        atom_name = atom.get_name().strip().upper()
        residue = atom.get_parent()
        residue_seq = residue.get_id()[1] if residue is not None else None

        # If the user supplied residue numbers, keep only those residues.
        if residue_numbers is not None and residue_seq not in residue_numbers:
            continue

        # If the user supplied explicit atom names, those win.
        if selected_atoms is not None and atom_name not in selected_atoms:
            continue

        # Otherwise, backbone mode keeps only the standard protein backbone atoms.
        if selected_atoms is None and atom_scope == "backbone" and atom_name not in BACKBONE_ATOM_NAMES:
            continue
        atoms.append(atom)

    if not atoms:
        raise ValueError(f"Frame {frame_index} has no atoms after filtering.")

    atoms_sorted = sorted(atoms, key=atom_key)
    atom_keys = [atom_key(atom) for atom in atoms_sorted]
    coords = np.array([atom.get_coord() for atom in atoms_sorted], dtype=np.float64)
    return atom_keys, coords


def align_to_reference(frame_keys, frame_coords, reference_keys):
    """Reorder a frame's coordinates so they match the reference frame atom order.

    This step is important because RMSD only makes sense if atom A is being
    compared to atom A, atom B to atom B, and so on.
    """
    if frame_keys == reference_keys:
        return frame_coords

    coord_map = {key: coord for key, coord in zip(frame_keys, frame_coords)}
    try:
        aligned = np.array([coord_map[key] for key in reference_keys], dtype=np.float64)
    except KeyError as exc:
        raise ValueError(f"Missing atom name {exc} while aligning a frame to the reference frame.") from exc
    return aligned


def calculate_rmsd(coords1, coords2):
    """Compute direct-coordinate RMSD between two same-shaped coordinate arrays.

    This uses the standard RMSD formula:

        RMSD = sqrt(mean(sum((x - y)^2)))

    where x and y are the aligned atom coordinates.
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape.")
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    return float(np.sqrt(msd))


def compute_rmsd_series(frames_parsed, reference_frame_index):
    """Compute a single RMSD value for each frame against one reference frame.

    This is the key difference from an all-vs-all matrix:
    - Matrix workflow: every frame is compared against every other frame.
    - This workflow: every frame is compared only against the chosen reference.
    """
    if len(frames_parsed) < 2:
        raise ValueError(
            "RMSD analysis requires at least two parsed frames. "
            "Select two or more frames before running RMSD."
        )

    if reference_frame_index < 1 or reference_frame_index > len(frames_parsed):
        raise ValueError(
            f"Reference frame {reference_frame_index} is out of range for {len(frames_parsed)} parsed frames."
        )

    # Pull the chosen reference frame out of the parsed data.
    reference_entry = None
    for frame_number, frame_keys, frame_coords in frames_parsed:
        if frame_number == reference_frame_index:
            reference_entry = (frame_keys, frame_coords)
            break

    if reference_entry is None:
        raise ValueError(
            f"Reference frame {reference_frame_index} was not included in the selected frame set."
        )

    reference_keys, reference_coords = reference_entry
    rmsd_rows = []

    for frame_number, frame_keys, frame_coords in frames_parsed:
        # Make sure the current frame and the reference frame have matching atom order.
        aligned_coords = align_to_reference(frame_keys, frame_coords, reference_keys)
        rmsd_value = calculate_rmsd(aligned_coords, reference_coords)
        rmsd_rows.append((frame_number, rmsd_value))

    return rmsd_rows


def write_rmsd_series_csv(rmsd_rows, output_csv, reference_frame, atom_scope):
    """Write one RMSD value per frame to a CSV file.

    The output is intentionally simple so beginners can open it in Excel,
    LibreOffice, or any text editor.
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Frame", "RMSD_Angstrom"])
        for frame_number, rmsd_value in rmsd_rows:
            writer.writerow([frame_number, f"{rmsd_value:.6f}"])

    # The CSV itself is intentionally simple. The chosen analysis settings are
    # printed in main() so the workflow remains easy to follow.


def compute_rmsd_matrix(frames_parsed):
    """Compute a full pairwise RMSD matrix for parsed frames.

    Returns:
    - frame_numbers: list of 1-based frame IDs in matrix order
    - matrix: NxN numpy array of RMSD values
    """
    if len(frames_parsed) < 2:
        raise ValueError("All-vs-all mode requires at least two selected frames.")

    frame_numbers = [frame_number for frame_number, _, _ in frames_parsed]
    n = len(frames_parsed)
    matrix = np.zeros((n, n), dtype=np.float64)

    for i, (_, ref_keys, ref_coords) in enumerate(frames_parsed):
        for j, (_, cmp_keys, cmp_coords) in enumerate(frames_parsed):
            aligned_cmp = align_to_reference(cmp_keys, cmp_coords, ref_keys)
            matrix[i, j] = calculate_rmsd(ref_coords, aligned_cmp)

    return frame_numbers, matrix


def write_rmsd_matrix_csv(frame_numbers, matrix, output_csv):
    """Write an all-vs-all RMSD matrix CSV with frame IDs as row/column labels."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Frame"] + frame_numbers)
        for frame_number, row in zip(frame_numbers, matrix):
            writer.writerow([frame_number] + [f"{value:.6f}" for value in row])


def write_reference_target_matrix_csv(reference_frames, target_frames, matrix, output_csv):
    """Write a reference-target RMSD matrix CSV.

    CSV orientation:
    - Vertical axis (rows): reference frames
    - Horizontal axis (columns): target frames
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ReferenceFrame"] + target_frames)
        for reference_frame, row in zip(reference_frames, matrix):
            writer.writerow([reference_frame] + [f"{value:.6f}" for value in row])


def main():
    """Run the per-frame RMSD workflow from the command line.

    Main workflow:
    1) Read the user's input file.
    2) Split it into frames.
    3) Parse each selected frame with Biopython.
    4) Compare chosen frames against the chosen reference frame.
    5) Write a simple CSV file with one RMSD value per selected frame.
    """
    args = parse_args()

    input_path = args.input.expanduser()
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    output_dir = args.output_dir.expanduser() if args.output_dir else input_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    output_name = args.output_name
    if args.all_vs_all and output_name == DEFAULT_OUTPUT_NAME:
        output_name = "rmsd_all_vs_all.csv"
    output_csv = output_dir / output_name

    # Convert the optional atom list into a set once, then reuse it everywhere.
    selected_atoms = parse_atom_selection(args.atoms)
    residue_numbers = parse_residue_selection(args.residues)

    print(f"Reading PDB file: {input_path}")
    frame_texts = split_frames_from_pdb(input_path)

    if len(frame_texts) < 2:
        raise ValueError(
            "RMSD requires at least two frames in the input PDB. "
            "Provide a multi-frame trajectory with 2 or more frames."
        )

    # Convert the user's frame-selection choices into actual frame numbers.
    # Think of this as the "candidate pool" of frames requested by CLI args.
    selected_frame_numbers = parse_frame_selection(args, len(frame_texts))
    if not selected_frame_numbers:
        raise ValueError("No frames were selected for analysis.")
    print(f"Total frames found in file: {len(frame_texts)}")
    print(f"Atom scope: {args.atom_scope}")
    if selected_atoms is not None:
        print(f"Specific atoms selected: {sorted(selected_atoms)}")
        print("Note: --atoms overrides --atom-scope.")
    if residue_numbers is not None:
        print(f"Residue numbers selected: {sorted(residue_numbers)}")

    parser = PDBParser(QUIET=True)

    if args.all_vs_all:
        if len(selected_frame_numbers) < 2:
            raise ValueError(
                "All-vs-all mode requires at least two selected frames. "
                "Expand --frames/--frame-start/--frame-end or use --compare-all."
            )

        print(f"Mode: all-vs-all matrix")
        print(f"Frames selected for matrix: {selected_frame_numbers}")

        frames_parsed = []
        total = len(selected_frame_numbers)
        for position, frame_number in enumerate(selected_frame_numbers, start=1):
            frame_text = frame_texts[frame_number - 1]
            atom_keys, coords = parse_frame_with_biopython(
                frame_text,
                frame_number,
                parser,
                args.atom_scope,
                selected_atoms=selected_atoms,
                residue_numbers=residue_numbers,
            )
            frames_parsed.append((frame_number, atom_keys, coords))
            if position % 100 == 0 or position == total:
                print(f"Parsed {position}/{total} selected frames")

        frame_numbers, matrix = compute_rmsd_matrix(frames_parsed)
        write_rmsd_matrix_csv(frame_numbers, matrix, output_csv)
        print(f"Saved all-vs-all RMSD matrix CSV to: {output_csv}")
        return

    # NvsM mode: explicit references on rows and explicit targets on columns.
    # This is the easiest way to get a single Excel-friendly table where row and
    # column meaning is fixed and obvious.
    if args.reference_frames:
        if not args.targets:
            raise ValueError(
                "When using --reference-frames, you must also provide --targets. "
                "Example: --reference-frames 2,4,5,9 --targets 45"
            )

        reference_frame_numbers = parse_frame_list_text(args.reference_frames, len(frame_texts))
        target_frame_numbers = parse_frame_list_text(args.targets, len(frame_texts))

        if not reference_frame_numbers:
            raise ValueError("No valid reference frames were provided in --reference-frames.")
        if not target_frame_numbers:
            raise ValueError("No valid target frames were provided in --targets.")

        print("Mode: reference-target matrix (rows=references, columns=targets)")
        print(f"Reference frames (rows): {reference_frame_numbers}")
        print(f"Target frames (columns): {target_frame_numbers}")

        # Parse all frames involved once, then reuse cached coordinates.
        needed_frames = sorted(set(reference_frame_numbers + target_frame_numbers))
        parsed_by_frame = {}
        total = len(needed_frames)

        for position, frame_number in enumerate(needed_frames, start=1):
            frame_text = frame_texts[frame_number - 1]
            atom_keys, coords = parse_frame_with_biopython(
                frame_text,
                frame_number,
                parser,
                args.atom_scope,
                selected_atoms=selected_atoms,
                residue_numbers=residue_numbers,
            )
            parsed_by_frame[frame_number] = (atom_keys, coords)
            if position % 100 == 0 or position == total:
                print(f"Parsed {position}/{total} frames used in matrix")

        matrix = np.zeros((len(reference_frame_numbers), len(target_frame_numbers)), dtype=np.float64)

        for i, ref_frame in enumerate(reference_frame_numbers):
            ref_keys, ref_coords = parsed_by_frame[ref_frame]
            for j, tgt_frame in enumerate(target_frame_numbers):
                tgt_keys, tgt_coords = parsed_by_frame[tgt_frame]
                aligned_target = align_to_reference(tgt_keys, tgt_coords, ref_keys)
                matrix[i, j] = calculate_rmsd(ref_coords, aligned_target)

        write_reference_target_matrix_csv(
            reference_frame_numbers,
            target_frame_numbers,
            matrix,
            output_csv,
        )
        print(f"Saved reference-target RMSD matrix CSV to: {output_csv}")
        return

    if args.reference_frame < 1 or args.reference_frame > len(frame_texts):
        raise ValueError(
            f"Reference frame {args.reference_frame} is out of range for a file with {len(frame_texts)} frames."
        )

    # In Nv1 mode, users can explicitly define targets with --targets.
    # If --targets is not provided, we fall back to the standard selection flags.
    if args.targets:
        explicit_targets = parse_frame_list_text(args.targets, len(frame_texts))
        target_frame_numbers = [n for n in explicit_targets if n != args.reference_frame]
    else:
        # Treat selected frames as comparison targets. Reference can be
        # separate from this set to support 1v1, 2v1, 5v1, etc.
        target_frame_numbers = [n for n in selected_frame_numbers if n != args.reference_frame]

    if not target_frame_numbers:
        raise ValueError(
            "No target frames remain after excluding the reference frame. "
            "Choose at least one frame different from --reference-frame."
        )

    print(f"Mode: target-vs-reference (Nv1)")
    print(f"Reference frame: {args.reference_frame}")
    print(f"Target frames selected for analysis: {target_frame_numbers}")

    total = len(target_frame_numbers)

    # Parse reference frame once.
    reference_text = frame_texts[args.reference_frame - 1]
    reference_keys, reference_coords = parse_frame_with_biopython(
        reference_text,
        args.reference_frame,
        parser,
        args.atom_scope,
        selected_atoms=selected_atoms,
        residue_numbers=residue_numbers,
    )

    rmsd_rows = []

    # Parse only the target frames and compare each one to the reference frame.
    for position, frame_number in enumerate(target_frame_numbers, start=1):
        frame_text = frame_texts[frame_number - 1]
        atom_keys, coords = parse_frame_with_biopython(
            frame_text,
            frame_number,
            parser,
            args.atom_scope,
            selected_atoms=selected_atoms,
            residue_numbers=residue_numbers,
        )
        aligned_coords = align_to_reference(atom_keys, coords, reference_keys)
        rmsd_value = calculate_rmsd(aligned_coords, reference_coords)
        rmsd_rows.append((frame_number, rmsd_value))
        if position % 100 == 0 or position == total:
            print(f"Processed {position}/{total} target frames")

    write_rmsd_series_csv(rmsd_rows, output_csv, args.reference_frame, args.atom_scope)

    print(f"Saved target-vs-reference RMSD CSV to: {output_csv}")


if __name__ == "__main__":
    main()
