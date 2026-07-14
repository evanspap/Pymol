"""Generate RMSD histograms from an RMSD matrix CSV file."""

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot RMSD histograms from a matrix CSV file."
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input RMSD matrix CSV file.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        required=True,
        help="Output directory for histogram images.",
    )
    parser.add_argument(
        "--first-output",
        type=str,
        default="rmsd_frame1_histogram.png",
        help="Filename for the frame-1 histogram image.",
    )
    parser.add_argument(
        "--second-output",
        type=str,
        default="rmsd_frame2_histogram.png",
        help="Filename for the frame-2 histogram image.",
    )
    parser.add_argument(
        "--all-output",
        type=str,
        default="rmsd_all_frames_histogram.png",
        help="Filename for the all-frames histogram image.",
    )
    return parser.parse_args()


def load_rmsd_matrix(csv_path):
    with open(csv_path, newline='', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if not header or len(header) < 2:
            raise ValueError('Input CSV must have at least a header plus one RMSD column.')

        frame1_values = []
        frame2_values = []
        all_values = []
        for row in reader:
            if not row or len(row) < 2:
                continue
            try:
                values = [float(x) for x in row[1:]]
            except ValueError:
                continue
            frame1_values.append(values[0])
            if len(values) > 1:
                frame2_values.append(values[1])
            all_values.extend(values)

    return (
        np.array(frame1_values, dtype=float),
        np.array(frame2_values, dtype=float),
        np.array(all_values, dtype=float),
    )


def resolve_shared_path(path: Path) -> Path:
    path = Path(os.path.normpath(str(path)))
    if path.exists():
        return path

    alternate = Path(str(path).replace('/', os.sep))
    if alternate.exists():
        return alternate

    if path.drive:
        drive_root = Path(f"{path.drive}{os.sep}")
        shortcut_root = drive_root / '.shortcut-targets-by-id'
        if shortcut_root.exists() and shortcut_root.is_dir():
            suffix = Path(*path.parts[-min(6, len(path.parts)):])
            target_name = path.name
            is_file_target = path.suffix != ''

            def matches(candidate: Path) -> bool:
                try:
                    candidate_suffix = Path(*candidate.parts[-len(suffix.parts):])
                except ValueError:
                    return False
                return candidate_suffix == suffix

            for dirpath, dirnames, filenames in os.walk(shortcut_root):
                if is_file_target and target_name in filenames:
                    candidate = Path(dirpath) / target_name
                    if candidate.exists() and matches(candidate):
                        return candidate
                if not is_file_target and target_name in dirnames:
                    candidate = Path(dirpath) / target_name
                    if candidate.exists() and matches(candidate):
                        return candidate

            for dirpath, dirnames, filenames in os.walk(shortcut_root):
                if is_file_target and target_name in filenames:
                    candidate = Path(dirpath) / target_name
                    if candidate.exists():
                        return candidate
                if not is_file_target and target_name in dirnames:
                    candidate = Path(dirpath) / target_name
                    if candidate.exists():
                        return candidate

    return path


def plot_histogram(values, output_path, title, bin_width=1.0, use_counts=False, overlay_normal=False, use_density=False):
    max_val = max(values.max(), 0)
    bins = np.arange(0, np.ceil(max_val / bin_width) * bin_width + bin_width, bin_width)
    counts, edges = np.histogram(values, bins=bins)

    if use_counts:
        plot_values = counts
        ylabel = 'Count'
    elif use_density:
        plot_values = counts.astype(float) / (counts.sum() * bin_width)
        ylabel = 'Probability density'
    else:
        plot_values = counts.astype(float) / counts.sum()
        ylabel = 'Relative frequency'

    mean_val = float(values.mean())
    std_val = float(values.std(ddof=0))

    plt.figure(figsize=(10, 6))
    plt.bar(edges[:-1], plot_values, width=bin_width, align='edge', color='tab:blue', edgecolor='black')
    for idx, height in enumerate(plot_values):
        if height > 0:
            label = f"{int(height)}" if use_counts else f"{height:.3f}"
            plt.text(edges[idx] + bin_width / 2, height, label, ha='center', va='bottom', fontsize=9)

    if overlay_normal and std_val > 0:
        x = np.linspace(values.min(), values.max(), 200)
        normal_pdf = np.exp(-0.5 * ((x - mean_val) / std_val) ** 2) / (std_val * np.sqrt(2 * np.pi))
        if use_counts:
            normal_curve = normal_pdf * values.size * bin_width
        elif use_density:
            normal_curve = normal_pdf
        else:
            normal_curve = normal_pdf * bin_width
        plt.plot(x, normal_curve, color='red', linewidth=2, label='Normal curve')
        plt.legend()

    plt.xlabel('RMSD (Å)')
    plt.ylabel(ylabel)
    plt.title(f"{title}\nmean={mean_val:.3f}, std={std_val:.3f}")
    plt.xticks(edges)
    plt.xlim(edges[0], edges[-1])
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()


def write_histogram_counts_csv(values, output_path, bin_width=1.0):
    max_val = max(values.max(), 0)
    bins = np.arange(0, np.ceil(max_val / bin_width) * bin_width + bin_width, bin_width)
    counts, edges = np.histogram(values, bins=bins)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['bin_start', 'bin_end', 'count'])
        for start, end, count in zip(edges[:-1], edges[1:], counts):
            writer.writerow([float(start), float(end), int(count)])


def main():
    args = parse_args()
    input_csv = resolve_shared_path(args.input)
    output_dir = resolve_shared_path(args.output_dir)

    if not input_csv.exists():
        raise FileNotFoundError(f'Input file not found: {args.input} (resolved to {input_csv})')
    if output_dir.exists() and not output_dir.is_dir():
        raise NotADirectoryError(f'Output path exists and is not a directory: {output_dir}')

    if not output_dir.exists():
        # Try to resolve the output directory from the input path if the user supplied a sibling/shared path.
        if input_csv.parent.exists():
            output_dir = input_csv.parent
        else:
            output_dir.mkdir(parents=True, exist_ok=True)

    print(f'Loading RMSD matrix from: {input_csv}')
    frame1_values, frame2_values, all_values = load_rmsd_matrix(input_csv)
    print(
        f'Got {len(frame1_values)} frame-1 comparisons, '
        f'{len(frame2_values)} frame-2 comparisons, and {len(all_values)} total RMSD values.'
    )

    first_path = output_dir / args.first_output
    second_path = output_dir / args.second_output
    all_path = output_dir / args.all_output
    all_csv = output_dir / 'rmsd_matrix_allframes.csv'

    plot_histogram(frame1_values, first_path, 'Raw Counts: RMSD to Frame 1', bin_width=1.0, use_counts=True)
    print(f'Saved frame-1 histogram to: {first_path}')

    if frame2_values.size > 0:
        plot_histogram(frame2_values, second_path, 'Raw Counts: RMSD to Frame 2', bin_width=1.0, use_counts=True)
        print(f'Saved frame-2 histogram to: {second_path}')
    else:
        print('No frame-2 values found in the input CSV; skipped frame-2 histogram.')

    plot_histogram(all_values, all_path, 'Raw Counts: RMSD across All Frame Comparisons', bin_width=1.0, use_counts=True)
    print(f'Saved all-frames histogram to: {all_path}')

    frame1_normal_path = output_dir / 'rmsd_frame1_normal_overlay.png'
    frame2_normal_path = output_dir / 'rmsd_frame2_normal_overlay.png'
    all_normal_path = output_dir / 'rmsd_all_frames_normal_overlay.png'

    plot_histogram(frame1_values, frame1_normal_path, 'Probability Density with Normal Overlay: RMSD to Frame 1', bin_width=1.0, use_density=True, overlay_normal=True)
    print(f'Saved frame-1 normal-overlay histogram to: {frame1_normal_path}')

    if frame2_values.size > 0:
        plot_histogram(frame2_values, frame2_normal_path, 'Probability Density with Normal Overlay: RMSD to Frame 2', bin_width=1.0, use_density=True, overlay_normal=True)
        print(f'Saved frame-2 normal-overlay histogram to: {frame2_normal_path}')
    else:
        print('No frame-2 values found in the input CSV; skipped frame-2 normal-overlay histogram.')

    plot_histogram(all_values, all_normal_path, 'Probability Density with Normal Overlay: RMSD across All Frame Comparisons', bin_width=1.0, use_density=True, overlay_normal=True)
    print(f'Saved all-frames normal-overlay histogram to: {all_normal_path}')

    write_histogram_counts_csv(all_values, all_csv, bin_width=1.0)
    print(f'Saved all-frames counts CSV to: {all_csv}')


if __name__ == '__main__':
    main()
