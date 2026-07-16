"""
True Biopython RMSD matrix calculator.

This version uses Bio.PDB.PDBParser to parse each frame and builds an all-vs-all
RMSD matrix for frames separated by END records in the input PDB trajectory.
"""

import argparse
import csv
import math
import sys
from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute all-vs-all RMSD matrix using Biopython parsing."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=Path(
            r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\apoa__run1_frames_1k.pdb"
        ),
        help="Input PDB file with multiple frames separated by END records.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path(
            r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k"
        ),
        help="Directory for output CSV.",
    )
    parser.add_argument(
        "--output-name",
        default="rmsd_matrix_v2.csv",
        help="Output CSV file name.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional: only process the first N frames.",
    )
    parser.add_argument(
        "--atom-scope",
        choices=["full", "backbone"],
        default="full",
        help="Atom selection scope for RMSD: full or backbone.",
    )
    return parser.parse_args()


def split_frames_from_pdb(pdb_file):
    frames = []
    current = []

    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("END"):
                if current:
                    frames.append("".join(current + ["END\n"]))
                    current = []
                continue

            if line.startswith(("ATOM", "HETATM", "TER", "CRYST1", "HEADER", "REMARK", "TITLE", "MODEL", "ENDMDL", "SEQRES", "COMPND")):
                current.append(line)

    if current:
        frames.append("".join(current + ["END\n"]))

    if len(frames) < 2:
        raise ValueError("Need at least 2 frames in the input PDB.")

    return frames


def atom_key(atom):
    residue = atom.get_parent()
    chain = residue.get_parent() if residue is not None else None
    residue_id = residue.get_id() if residue is not None else ("", 0, "")
    return (
        chain.get_id() if chain is not None else "",
        residue.get_resname() if residue is not None else "",
        residue_id[0],
        residue_id[1],
        residue_id[2],
        atom.get_name(),
        atom.get_altloc() or "",
    )


def parse_frame_with_biopython(frame_text, frame_index, parser, atom_scope):
    structure = parser.get_structure(f"frame_{frame_index}", StringIO(frame_text))
    model = next(structure.get_models())

    atoms = []
    for atom in model.get_atoms():
        if atom_scope == "backbone" and atom.get_name().strip() not in BACKBONE_ATOM_NAMES:
            continue
        atoms.append(atom)
    if not atoms:
        raise ValueError(f"Frame {frame_index} has no atoms.")

    atoms_sorted = sorted(atoms, key=atom_key)
    keys = [atom_key(a) for a in atoms_sorted]
    coords = np.array([a.get_coord() for a in atoms_sorted], dtype=np.float64)
    return keys, coords


def build_coordinate_tensor(frames_parsed):
    ref_keys, _ = frames_parsed[0]
    atom_count = len(ref_keys)
    frame_count = len(frames_parsed)

    vectors = np.empty((frame_count, atom_count * 3), dtype=np.float64)

    for i, (keys, coords) in enumerate(frames_parsed):
        if len(keys) != atom_count:
            raise ValueError(
                f"Frame {i + 1} atom count mismatch: {len(keys)} vs {atom_count}"
            )

        if keys == ref_keys:
            vectors[i, :] = coords.reshape(-1)
            continue

        key_to_coord = {k: c for k, c in zip(keys, coords)}
        try:
            aligned = np.array([key_to_coord[k] for k in ref_keys], dtype=np.float64)
        except KeyError as exc:
            raise ValueError(
                f"Frame {i + 1} missing atom key {exc}."
            ) from exc
        vectors[i, :] = aligned.reshape(-1)

    return vectors, atom_count


def compute_rmsd_matrix(frames_parsed):
    vectors, atom_count = build_coordinate_tensor(frames_parsed)
    sq_norms = np.einsum("ij,ij->i", vectors, vectors)
    gram = vectors @ vectors.T
    msd = (sq_norms[:, None] + sq_norms[None, :] - 2.0 * gram) / float(atom_count)
    msd = np.maximum(msd, 0.0)
    return np.sqrt(msd)


def write_rmsd_csv(rmsd_matrix, output_csv):
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    n = rmsd_matrix.shape[0]
    headers = ["Frame"] + [f"Frame{j}" for j in range(1, n + 1)]

    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for i in range(n):
            row = [str(i + 1)] + [f"{x:.6f}" for x in rmsd_matrix[i]]
            writer.writerow(row)


def main():
    args = parse_args()

    atom_scope = args.atom_scope

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}")
        sys.exit(1)

    print(f"Reading PDB file: {args.input}")
    frames = split_frames_from_pdb(args.input)

    if args.max_frames is not None:
        frames = frames[: args.max_frames]

    print(f"Frames to process: {len(frames)}")
    print(f"Atom scope: {atom_scope}")

    parser = PDBParser(QUIET=True)
    parsed = []
    total = len(frames)
    for idx, frame_text in enumerate(frames, start=1):
        parsed.append(parse_frame_with_biopython(frame_text, idx, parser, atom_scope))
        if idx % 100 == 0 or idx == total:
            print(f"Parsed {idx}/{total} frames")

    print("Computing RMSD matrix...")
    rmsd_matrix = compute_rmsd_matrix(parsed)

    output_csv = args.output_dir / args.output_name
    print(f"Writing RMSD matrix to: {output_csv}")
    write_rmsd_csv(rmsd_matrix, output_csv)
    print("Done.")


if __name__ == "__main__":
    main()
