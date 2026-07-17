"""Compute and plot comparative RMSD category histograms for backbone and full-atom RMSD."""

import numpy as np
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt

BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}
CATEGORY_ORDER = ["<=2 Å", "2-3 Å", "3-4 Å", "4-5 Å", ">5 Å"]


def parse_pdb_frames(pdb_file):
    frames = []
    current_coords = []
    current_names = []

    with open(pdb_file, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except ValueError:
                    continue
                current_coords.append([x, y, z])
                current_names.append(atom_name)
            elif line.startswith('END'):
                if current_coords:
                    coords = np.array(current_coords, dtype=float)
                    atom_names = list(current_names)
                    backbone_mask = np.array([name in BACKBONE_ATOM_NAMES for name in atom_names])
                    backbone_coords = coords[backbone_mask]
                    frames.append({
                        'coords': coords,
                        'backbone_coords': backbone_coords,
                        'atom_names': atom_names,
                        'backbone_mask': backbone_mask,
                    })
                    current_coords = []
                    current_names = []

        if current_coords:
            coords = np.array(current_coords, dtype=float)
            atom_names = list(current_names)
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
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    return np.sqrt(msd)


def get_rmsd_category(rmsd):
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
    full_counts = Counter()
    backbone_counts = Counter()

    for target in frames:
        for frame in frames:
            full_rmsd = calculate_rmsd(frame['coords'], target['coords'])
            backbone_rmsd = calculate_rmsd(frame['backbone_coords'], target['backbone_coords'])
            full_counts[get_rmsd_category(full_rmsd)] += 1
            backbone_counts[get_rmsd_category(backbone_rmsd)] += 1

    return full_counts, backbone_counts


def save_counts_csv(full_counts, backbone_counts, output_file):
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        f.write('Category,FullAtomCount,BackboneCount\n')
        for category in CATEGORY_ORDER:
            f.write(f"{category},{full_counts.get(category,0)},{backbone_counts.get(category,0)}\n")


def plot_comparative_histogram(full_counts, backbone_counts, output_image):
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
    pdb_file = Path(r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\apoa__run1_frames_1k.pdb")
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    frames = parse_pdb_frames(pdb_file)
    if not frames:
        raise ValueError("No frames were parsed from the PDB file.")

    print(f"Parsed {len(frames)} frames from PDB.")
    print(f"Full-atom frame size: {len(frames[0]['coords'])}")
    print(f"Backbone frame size: {len(frames[0]['backbone_coords'])}")

    full_counts, backbone_counts = generate_histogram_counts(frames)

    output_csv = Path(__file__).parent / 'rmsd_full_vs_backbone_counts.csv'
    save_counts_csv(full_counts, backbone_counts, output_csv)
    print(f"Saved comparative counts to {output_csv}")

    output_image = Path(__file__).parent / 'rmsd_comparative_histogram.png'
    plot_comparative_histogram(full_counts, backbone_counts, output_image)
    print(f"Saved comparative histogram to {output_image}")


if __name__ == '__main__':
    main()
