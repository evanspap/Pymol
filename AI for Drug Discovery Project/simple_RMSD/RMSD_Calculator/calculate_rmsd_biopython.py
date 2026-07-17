"""
Faster RMSD matrix calculator for PDB snapshots separated by END markers.
This script parses the raw ATOM/HETATM coordinate lines directly from each frame,
uses atom identity keys to align the atom order, and writes an all-vs-all RMSD matrix.

Output format:
- CSV matrix with rows/columns labeled by frame number.
- Version 2 label is written in the output filename.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute all-vs-all RMSD matrix from a PDB file using direct coordinate parsing."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        default=Path(r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\apoa__run1_frames_1k.pdb"),
        help="Input PDB file with repeated coordinate snapshots.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path(r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k"),
        help="Directory for the output CSV file.",
    )
    parser.add_argument(
        "--output-name",
        default="rmsd_matrix_v2.csv",
        help="Output CSV filename.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional: only process the first N frames (useful for quick test runs).",
    )
    return parser.parse_args()


def split_frames_from_pdb(pdb_file):
    with open(pdb_file, "r", encoding="utf-8", errors="ignore") as handle:
        lines = handle.readlines()

    frames = []
    current = []
    for line in lines:
        if line.startswith("END"):
            if current:
                frames.append(current)
                current = []
            continue
        if line.startswith(("ATOM", "HETATM")):
            current.append(line)

    if current:
        frames.append(current)

    if len(frames) < 2:
        raise ValueError("The input PDB must contain at least 2 frames/separators.")
    return frames


def atom_key(line):
    chain = line[21].strip()
    resseq = int(line[22:26].strip() or 0)
    icode = line[26:27]
    atom_name = line[12:16].strip()
    resname = line[17:20].strip()
    altloc = line[16:17]
    return (chain, resname, resseq, icode, atom_name, altloc)


def parse_frame_lines(lines):
    atoms = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except ValueError:
            continue
        atoms.append((atom_key(line), np.array([x, y, z], dtype=float)))

    if not atoms:
        raise ValueError("No coordinate atoms found in frame")

    atoms.sort(key=lambda item: item[0])
    keys = [item[0] for item in atoms]
    coords = np.array([item[1] for item in atoms], dtype=float)
    return keys, coords


def build_coordinate_tensor(frames):
    ref_keys, _ = frames[0]
    atom_count = len(ref_keys)
    frame_count = len(frames)
    data = np.empty((frame_count, atom_count * 3), dtype=np.float64)

    for i, (keys_i, coords_i) in enumerate(frames):
        if len(keys_i) != atom_count:
            raise ValueError(
                f"Frame {i + 1} atom count mismatch: {len(keys_i)} vs {atom_count}"
            )

        if keys_i != ref_keys:
            key_to_coord = {k: c for k, c in zip(keys_i, coords_i)}
            try:
                ordered = np.array([key_to_coord[k] for k in ref_keys], dtype=np.float64)
            except KeyError as exc:
                raise ValueError(
                    f"Frame {i + 1} is missing atom key: {exc}"
                ) from exc
            data[i, :] = ordered.reshape(-1)
        else:
            data[i, :] = coords_i.reshape(-1)

    return data, atom_count


def compute_rmsd_matrix(frames):
    vectors, atom_count = build_coordinate_tensor(frames)
    sq_norms = np.einsum("ij,ij->i", vectors, vectors)
    gram = vectors @ vectors.T
    msd = (sq_norms[:, None] + sq_norms[None, :] - 2.0 * gram) / float(atom_count)
    msd = np.maximum(msd, 0.0)
    rmsd = np.sqrt(msd)
    return rmsd


def write_rmsd_matrix(rmsd_matrix, output_csv):
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame_count = rmsd_matrix.shape[0]
    headers = ["Frame"] + [f"Frame{j}" for j in range(1, frame_count + 1)]

    with open(output_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for i in range(frame_count):
            row = [str(i + 1)] + [f"{value:.6f}" for value in rmsd_matrix[i]]
            writer.writerow(row)


def main():
    args = parse_args()
    input_pdb = args.input
    output_dir = args.output_dir
    output_name = args.output_name
    max_frames = args.max_frames

    if not input_pdb.exists():
        print(f"Error: input file not found: {input_pdb}")
        sys.exit(1)

    print(f"Reading PDB file: {input_pdb}")
    frame_lines = split_frames_from_pdb(input_pdb)
    if max_frames is not None:
        frame_lines = frame_lines[:max_frames]
    print(f"Parsed {len(frame_lines)} frames")

    frames = []
    for idx, lines in enumerate(frame_lines, start=1):
        frames.append(parse_frame_lines(lines))
        if idx % 100 == 0:
            print(f"Prepared {idx}/{len(frame_lines)} frames")

    print("Computing RMSD matrix...")
    rmsd_matrix = compute_rmsd_matrix(frames)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / output_name
    print(f"Writing RMSD matrix to: {output_csv}")
    write_rmsd_matrix(rmsd_matrix, output_csv)
    print("Done.")


if __name__ == "__main__":
    main()
