"""Plot comparative RMSD histograms for full-atom vs backbone matrices."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot comparative count and density histograms for full-atom vs backbone RMSD."
    )
    parser.add_argument("--full", type=Path, required=True, help="Full-atom RMSD matrix CSV")
    parser.add_argument("--backbone", type=Path, required=True, help="Backbone RMSD matrix CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for PNGs")
    parser.add_argument("--bin-width", type=float, default=1.0, help="Histogram bin width in Angstrom")
    return parser.parse_args()


def load_matrix_values(csv_path):
    values = []
    with open(csv_path, newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError(f"Invalid RMSD matrix CSV: {csv_path}")
        for row in reader:
            if len(row) < 2:
                continue
            try:
                vals = [float(x) for x in row[1:]]
            except ValueError:
                continue
            values.extend(vals)
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        raise ValueError(f"No RMSD values found in {csv_path}")
    return arr


def build_bins(a, b, bin_width):
    max_val = max(float(np.max(a)), float(np.max(b)), 0.0)
    stop = np.ceil(max_val / bin_width) * bin_width + bin_width
    return np.arange(0.0, stop, bin_width)


def plot_comparison_counts(full_vals, bb_vals, output_path, bin_width):
    bins = build_bins(full_vals, bb_vals, bin_width)
    full_counts, edges = np.histogram(full_vals, bins=bins)
    bb_counts, _ = np.histogram(bb_vals, bins=bins)

    width = bin_width * 0.42
    x = edges[:-1]

    plt.figure(figsize=(12, 6))
    plt.bar(x - width / 2, full_counts, width=width, align="edge", label="Full-atom", alpha=0.85, edgecolor="black")
    plt.bar(x + width / 2, bb_counts, width=width, align="edge", label="Backbone", alpha=0.85, edgecolor="black")
    plt.xlabel("RMSD (A)")
    plt.ylabel("Count")
    plt.title("RMSD Comparison (All Frame Pairs): Full-Atom vs Backbone")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_comparison_density(full_vals, bb_vals, output_path, bin_width):
    bins = build_bins(full_vals, bb_vals, bin_width)
    full_counts, edges = np.histogram(full_vals, bins=bins)
    bb_counts, _ = np.histogram(bb_vals, bins=bins)

    full_density = full_counts.astype(float) / (full_counts.sum() * bin_width)
    bb_density = bb_counts.astype(float) / (bb_counts.sum() * bin_width)
    centers = edges[:-1] + bin_width / 2.0

    plt.figure(figsize=(12, 6))
    plt.plot(centers, full_density, marker="o", linewidth=2, label="Full-atom density")
    plt.plot(centers, bb_density, marker="o", linewidth=2, label="Backbone density")
    plt.xlabel("RMSD (A)")
    plt.ylabel("Probability density")
    plt.title("RMSD Density Comparison (All Frame Pairs): Full-Atom vs Backbone")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    full_vals = load_matrix_values(args.full)
    bb_vals = load_matrix_values(args.backbone)

    count_png = args.output_dir / "rmsd_full_vs_backbone_histogram.png"
    density_png = args.output_dir / "rmsd_full_vs_backbone_density_histogram.png"

    plot_comparison_counts(full_vals, bb_vals, count_png, args.bin_width)
    plot_comparison_density(full_vals, bb_vals, density_png, args.bin_width)

    print(f"Saved comparative count histogram to: {count_png}")
    print(f"Saved comparative density histogram to: {density_png}")


if __name__ == "__main__":
    main()
