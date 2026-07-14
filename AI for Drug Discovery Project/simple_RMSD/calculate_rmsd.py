"""
RMSD Calculator for Multi-Frame PDB Files
Creates a single CSV where each column contains RMSD values of all frames
against one target frame.
"""

import argparse
import numpy as np
from pathlib import Path
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute RMSD of each frame against every target frame and export to CSV."
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        default=Path(r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\apoa__run1_frames_1k.pdb"),
        help="Input multi-frame PDB file."
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output CSV file path.",
    )
    return parser.parse_args()


def parse_pdb_frames(pdb_file):
    """
    Parse a multi-frame PDB file and extract individual frames.
    Each frame is separated by an 'END' keyword.
    """
    frames = []
    current_frame = []

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        current_frame.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
                elif line.startswith('END'):
                    if current_frame:
                        frames.append(np.array(current_frame, dtype=float))
                        current_frame = []

        if current_frame:
            frames.append(np.array(current_frame, dtype=float))

    except FileNotFoundError:
        print(f"Error: File not found: {pdb_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    return frames


def calculate_rmsd(coords1, coords2):
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    return np.sqrt(msd)


def write_rmsd_matrix(frames, csv_file):
    csv_file.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_file, 'w', newline='') as f:
        header = ["Frame"] + [f"Frame{j}" for j in range(1, len(frames) + 1)]
        f.write(",".join(header) + "\n")

        for i, frame in enumerate(frames, start=1):
            values = [str(i)]
            for target_frame in frames:
                rmsd = calculate_rmsd(frame, target_frame)
                values.append(f"{rmsd:.6f}")
            f.write(",".join(values) + "\n")


def main():
    args = parse_args()
    pdb_file = args.input
    output_csv = args.output

    if not pdb_file.exists():
        print(f"Error: PDB file not found: {pdb_file}")
        sys.exit(1)

    print(f"Reading PDB file: {pdb_file}")
    print("Parsing frames...")
    frames = parse_pdb_frames(pdb_file)
    print(f"Successfully parsed {len(frames)} frames")

    if len(frames) < 2:
        print("Error: Need at least 2 frames to calculate RMSD")
        sys.exit(1)

    print(f"Writing RMSD matrix to: {output_csv}")
    write_rmsd_matrix(frames, output_csv)
    print("Done.")


if __name__ == "__main__":
    main()
