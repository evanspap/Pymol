"""Standalone RMSD calculation utilities extracted from the per-frame workflow.
### Version 1.0 (Subject to Updates and Change)
## Author(s): Yagna Devarakonda, Research Volunteer/Assistant, Stony Brook University, NY, USA
## Dr. Evangelos Papadopoulos, Professor, Dana-Farber Cancer Institute
## Date: 7-19-26

**DISCLAIMER**
=============
**Generative AI was used in writing this code, and the authors take
full responsibility for the final content. Please modify code for personal use
if necessary. The authors are not responsible for any errors or issues that
may arise from using this code which originate from team specific
modifications, dependencies, or other factors outside the authors' control.**

What this file contains:
- The direct RMSD formula implementation.
- Alignment helper to keep atom ordering consistent.
- Series helper to compute RMSD for many frames against one reference.
"""

from __future__ import annotations

import argparse
from io import StringIO
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
from Bio.PDB import PDBParser


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}


def align_to_reference(
    frame_keys: Sequence[Tuple],
    frame_coords: np.ndarray,
    reference_keys: Sequence[Tuple],
) -> np.ndarray:
    """Reorder frame coordinates to match the reference-frame atom order."""
    if list(frame_keys) == list(reference_keys):
        return frame_coords

    coord_map = {key: coord for key, coord in zip(frame_keys, frame_coords)}
    try:
        aligned = np.array([coord_map[key] for key in reference_keys], dtype=np.float64)
    except KeyError as exc:
        raise ValueError(
            f"Missing atom name {exc} while aligning a frame to the reference frame."
        ) from exc
    return aligned


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute direct-coordinate RMSD between two same-shaped coordinate arrays.

    RMSD = sqrt(mean(sum((x - y)^2)))
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape.")

    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    return float(np.sqrt(msd))


def compute_rmsd_series(
    frames_parsed: Sequence[Tuple[int, Sequence[Tuple], np.ndarray]],
    reference_frame_index: int,
) -> List[Tuple[int, float]]:
    """Compute one RMSD value per frame against a chosen reference frame.

    Notes:
    - At least two frames are required.
    - Frame IDs are expected to be 1-based IDs carried with each tuple.
    """
    if len(frames_parsed) < 2:
        raise ValueError(
            "RMSD analysis requires at least two parsed frames. "
            "Select two or more frames before running RMSD."
        )

    frame_ids = [frame_number for frame_number, _, _ in frames_parsed]
    if reference_frame_index not in frame_ids:
        raise ValueError(
            f"Reference frame {reference_frame_index} was not included in the selected frame set."
        )

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
    rmsd_rows: List[Tuple[int, float]] = []

    for frame_number, frame_keys, frame_coords in frames_parsed:
        aligned_coords = align_to_reference(frame_keys, frame_coords, reference_keys)
        rmsd_value = calculate_rmsd(aligned_coords, reference_coords)
        rmsd_rows.append((frame_number, rmsd_value))

    return rmsd_rows


def split_frames_from_pdb(pdb_file: Path) -> List[str]:
    """Split a multi-frame PDB file into frame text blocks."""
    frames: List[str] = []
    current: List[str] = []

    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("END"):
                if current:
                    frames.append("".join(current + ["END\n"]))
                    current = []
                continue

            if line.startswith(("ATOM", "HETATM", "TER", "MODEL", "ENDMDL", "CRYST1", "TITLE", "REMARK")):
                current.append(line)

    if current:
        frames.append("".join(current + ["END\n"]))

    return frames


def atom_key(atom) -> Tuple:
    """Build a stable key used to align atoms across frames."""
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
    frame_text: str,
    frame_index: int,
    parser: PDBParser,
    atom_scope: str = "full",
    selected_atoms: set[str] | None = None,
) -> Tuple[List[Tuple], np.ndarray]:
    """Parse one frame and return atom identity keys plus coordinate array."""
    structure = parser.get_structure(f"frame_{frame_index}", StringIO(frame_text))
    model = next(structure.get_models())

    atoms = []
    for atom in model.get_atoms():
        atom_name = atom.get_name().strip().upper()
        if selected_atoms is not None and atom_name not in selected_atoms:
            continue
        if selected_atoms is None and atom_scope == "backbone" and atom_name not in BACKBONE_ATOM_NAMES:
            continue
        atoms.append(atom)

    if not atoms:
        raise ValueError(f"Frame {frame_index} has no atoms after filtering.")

    atoms_sorted = sorted(atoms, key=atom_key)
    atom_keys = [atom_key(atom) for atom in atoms_sorted]
    coords = np.array([atom.get_coord() for atom in atoms_sorted], dtype=np.float64)
    return atom_keys, coords


def calculate_rmsd_from_pdb_frames(
    pdb_path: Path,
    frame_a: int,
    frame_b: int,
    atom_scope: str = "full",
    atoms: Sequence[str] | None = None,
) -> float:
    """Compute one RMSD value between two 1-based frame IDs from a PDB file."""
    frame_texts = split_frames_from_pdb(pdb_path)
    if len(frame_texts) < 2:
        raise ValueError("Input PDB must contain at least two frames.")

    if frame_a < 1 or frame_b < 1 or frame_a > len(frame_texts) or frame_b > len(frame_texts):
        raise ValueError(
            f"Frame IDs must be within 1..{len(frame_texts)}. Received: frame_a={frame_a}, frame_b={frame_b}"
        )

    selected_atoms = {name.strip().upper() for name in atoms} if atoms else None
    parser = PDBParser(QUIET=True)

    keys_a, coords_a = parse_frame_with_biopython(
        frame_texts[frame_a - 1],
        frame_a,
        parser,
        atom_scope=atom_scope,
        selected_atoms=selected_atoms,
    )
    keys_b, coords_b = parse_frame_with_biopython(
        frame_texts[frame_b - 1],
        frame_b,
        parser,
        atom_scope=atom_scope,
        selected_atoms=selected_atoms,
    )

    aligned_b = align_to_reference(keys_b, coords_b, keys_a)
    return calculate_rmsd(coords_a, aligned_b)


def parse_cli_args() -> argparse.Namespace:
    """Parse optional CLI args for direct terminal RMSD output."""
    parser = argparse.ArgumentParser(
        description="Print RMSD values directly from multi-frame PDB data (no CSV output)."
    )
    parser.add_argument("--input", type=Path, required=True, help="Path to multi-frame PDB file.")
    parser.add_argument("--frame-a", type=int, required=True, help="First frame ID (1-based).")
    parser.add_argument("--frame-b", type=int, required=True, help="Second frame ID (1-based).")
    parser.add_argument(
        "--atom-scope",
        choices=["full", "backbone"],
        default="full",
        help="Use all atoms or backbone atoms only.",
    )
    parser.add_argument(
        "--atoms",
        type=str,
        default=None,
        help="Optional comma-separated atom names (example: CA,CB,N). Overrides --atom-scope.",
    )
    return parser.parse_args()


def main() -> None:
    """CLI entrypoint: print one RMSD value for two selected frames."""
    args = parse_cli_args()
    pdb_path = args.input.expanduser()
    if not pdb_path.exists():
        raise FileNotFoundError(f"Input PDB not found: {pdb_path}")

    atoms = [name.strip() for name in args.atoms.split(",") if name.strip()] if args.atoms else None
    rmsd_value = calculate_rmsd_from_pdb_frames(
        pdb_path=pdb_path,
        frame_a=args.frame_a,
        frame_b=args.frame_b,
        atom_scope=args.atom_scope,
        atoms=atoms,
    )
    print(f"RMSD_Angstrom(frame {args.frame_a} vs frame {args.frame_b}): {rmsd_value:.6f}")


if __name__ == "__main__":
    main()
