"""Compute and plot comparative RMSD category histograms for backbone and full-atom RMSD.

PowerShell usage (copy one scenario at a time):

Scenario 1: Basic run (output files go to script folder)
    python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb"

Scenario 2: Custom output folder
    python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb" --output-dir "C:\Path\To\results"

Scenario 3: Custom output filenames
    python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb" --output-csv "team_counts.csv" --output-image "team_histogram.png"

Important path note for teams cloning from GitHub:
    Local paths can differ across institutions and machines.
    Replace all example paths with your own system-specific paths.
"""

import argparse
import numpy as np
from collections import Counter
from io import StringIO
from pathlib import Path
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
CATEGORY_ORDER = ["<=2 Å", "2-3 Å", "3-4 Å", "4-5 Å", ">5 Å"]


def parse_args():
    """Parse command-line arguments for portable cross-team execution.

    PowerShell usage examples for any institution (replace placeholders):

    1) Basic run (outputs saved in the same folder as this script):
       python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb"

    2) Custom output directory:
       python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb" --output-dir "C:\Path\To\results"

    3) Custom output file names:
       python .\Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py --input-pdb "C:\Path\To\frames_trajectory.pdb" --output-csv "team_counts.csv" --output-image "team_histogram.png"

    Important:
    - Teams cloning from GitHub will have different local folder paths.
    - Always replace paths in examples with your own machine-specific paths.
    """
    parser = argparse.ArgumentParser(
        description="Compute and plot comparative RMSD category histograms for full-atom vs backbone coordinates from a multi-frame PDB file."
    )
    parser.add_argument(
        "--input-pdb",
        type=Path,
        required=True,
        help="Path to input multi-frame PDB file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for CSV and histogram image (default: script directory).",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default="rmsd_full_vs_backbone_counts.csv",
        help="Output CSV filename for category counts.",
    )
    parser.add_argument(
        "--output-image",
        type=str,
        default="rmsd_comparative_histogram.png",
        help="Output image filename for comparative histogram.",
    )
    return parser.parse_args()


def split_frames_from_pdb(pdb_file):
    """Split multi-frame PDB into per-frame PDB text chunks.

    This mirrors the RMSD calculator approach so both scripts interpret frame
    boundaries consistently.
    """
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

    return frames


def atom_key(atom):
    """Create a stable sortable key for an atom across frames."""
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


def parse_frame_with_biopython(frame_text, frame_index, parser):
    """Parse a single frame with Biopython and return sorted atom arrays."""
    structure = parser.get_structure(f"frame_{frame_index}", StringIO(frame_text))
    model = next(structure.get_models())

    atoms = list(model.get_atoms())
    if not atoms:
        raise ValueError(f"Frame {frame_index} has no atoms.")

    atoms_sorted = sorted(atoms, key=atom_key)
    coords = np.array([a.get_coord() for a in atoms_sorted], dtype=np.float64)
    atom_names = [a.get_name().strip() for a in atoms_sorted]
    return atom_names, coords


def parse_pdb_frames(pdb_file):
    """Read a multi-frame PDB file and return per-frame coordinate bundles.

    Each frame dictionary contains:
    - coords: all atom coordinates for that frame
    - backbone_coords: backbone-only coordinates (N, CA, C, O, OXT)
    - atom_names: atom names in parse order
    - backbone_mask: boolean mask selecting backbone atoms
    """
    frames = []
    frame_texts = split_frames_from_pdb(pdb_file)
    parser = PDBParser(QUIET=True)

    for idx, frame_text in enumerate(frame_texts, start=1):
        atom_names, coords = parse_frame_with_biopython(frame_text, idx, parser)
        backbone_mask = np.array([name in BACKBONE_ATOM_NAMES for name in atom_names])
        backbone_coords = coords[backbone_mask]
        frames.append({
            'coords': coords,
            'backbone_coords': backbone_coords,
            'atom_names': atom_names,
            'backbone_mask': backbone_mask,
        })

    return frames


def calculate_rmsd(coords1, coords2):
    """Compute direct-coordinate RMSD between two same-shaped coordinate arrays."""
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    return np.sqrt(msd)


def get_rmsd_category(rmsd):
    """Map a numeric RMSD value into a categorical bucket."""
    if rmsd <= 2.0:
        return "<=2 Å"
    if rmsd <= 3.0:
        return "2-3 Å"
    if rmsd <= 4.0:
        return "3-4 Å"
    if rmsd <= 5.0:
        return "4-5 Å"
    return ">5 Å"


def generate_histogram_counts(frames):
    """Compute all-vs-all category counts for full-atom and backbone RMSD.

    For each target frame and each compared frame, calculate both RMSDs and
    increment category counters.
    """
    full_counts = Counter()
    backbone_counts = Counter()

    for target in frames:
        for frame in frames:
            # Full-atom comparison.
            full_rmsd = calculate_rmsd(frame['coords'], target['coords'])
            # Backbone-only comparison.
            backbone_rmsd = calculate_rmsd(frame['backbone_coords'], target['backbone_coords'])
            full_counts[get_rmsd_category(full_rmsd)] += 1
            backbone_counts[get_rmsd_category(backbone_rmsd)] += 1

    return full_counts, backbone_counts


def save_counts_csv(full_counts, backbone_counts, output_file):
    """Write category counts to CSV for downstream reporting and reproducibility."""
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        f.write('Category,FullAtomCount,BackboneCount\n')
        for category in CATEGORY_ORDER:
            f.write(f"{category},{full_counts.get(category,0)},{backbone_counts.get(category,0)}\n")


def plot_comparative_histogram(full_counts, backbone_counts, output_image):
    """Create and save side-by-side bar chart for category counts."""
    categories = CATEGORY_ORDER
    full_values = [full_counts.get(cat, 0) for cat in categories]
    backbone_values = [backbone_counts.get(cat, 0) for cat in categories]

    x = np.arange(len(categories))
    width = 0.35

    plt.figure(figsize=(11, 6))
    plt.bar(x - width/2, full_values, width, label='Full-Atom RMSD', color='tab:blue', edgecolor='black')
    plt.bar(x + width/2, backbone_values, width, label='Backbone RMSD', color='tab:orange', edgecolor='black')

    for idx, height in enumerate(full_values):
        plt.annotate(f'{height:,}', xy=(x[idx] - width/2, height), xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=9)
    for idx, height in enumerate(backbone_values):
        plt.annotate(f'{height:,}', xy=(x[idx] + width/2, height), xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=9)

    plt.xticks(x, categories)
    plt.ylabel('Count')
    plt.xlabel('RMSD Category')
    plt.title('Comparative RMSD Category Distribution: Full-Atom vs Backbone')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_image, dpi=200)
    plt.close()


def main():
    """Main workflow: parse args, validate paths, compute counts, export outputs."""
    args = parse_args()

    # Resolve input path and validate it exists.
    pdb_file = args.input_pdb.expanduser()
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    # Default output directory is this script's folder unless user overrides it.
    output_dir = args.output_dir.expanduser() if args.output_dir else Path(__file__).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    frames = parse_pdb_frames(pdb_file)
    if not frames:
        raise ValueError("No frames were parsed from the PDB file.")

    print(f"Parsed {len(frames)} frames from PDB.")
    print(f"Full-atom frame size: {len(frames[0]['coords'])}")
    print(f"Backbone frame size: {len(frames[0]['backbone_coords'])}")

    full_counts, backbone_counts = generate_histogram_counts(frames)

    output_csv = output_dir / args.output_csv
    save_counts_csv(full_counts, backbone_counts, output_csv)
    print(f"Saved comparative counts to {output_csv}")

    output_image = output_dir / args.output_image
    plot_comparative_histogram(full_counts, backbone_counts, output_image)
    print(f"Saved comparative histogram to {output_image}")


if __name__ == '__main__':
    main()
