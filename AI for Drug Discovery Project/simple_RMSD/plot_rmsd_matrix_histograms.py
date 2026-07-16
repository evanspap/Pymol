"""Generate RMSD histograms from an RMSD matrix CSV file."""

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import table


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
    output_path = resolve_output_path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        plt.savefig(output_path, dpi=200)
    except OSError:
        fallback_path = resolve_output_path(output_path.parent / output_path.name)
        plt.savefig(fallback_path, dpi=200)
        output_path = fallback_path
    plt.close()
    return output_path


def write_histogram_counts_csv(values, output_path, bin_width=1.0):
    max_val = max(values.max(), 0)
    bins = np.arange(0, np.ceil(max_val / bin_width) * bin_width + bin_width, bin_width)
    counts, edges = np.histogram(values, bins=bins)
    output_path = resolve_output_path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['bin_start', 'bin_end', 'count'])
        for start, end, count in zip(edges[:-1], edges[1:], counts):
            writer.writerow([float(start), float(end), int(count)])
    return output_path


def choose_writable_output_dir(requested_dir, input_csv):
    candidates = []
    requested = resolve_shared_path(requested_dir)
    candidates.append(requested)
    if requested != requested_dir:
        candidates.append(requested_dir)
    if input_csv is not None and input_csv.parent.exists():
        candidates.append(input_csv.parent)
    candidates.append(Path.cwd())
    candidates.append(Path(__file__).resolve().parent)

    seen = set()
    for candidate in candidates:
        key = str(candidate.resolve()) if candidate.exists() else str(candidate)
        if key in seen:
            continue
        seen.add(key)
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            probe_path = candidate / '.write_test.tmp'
            with open(probe_path, 'w', encoding='utf-8') as fh:
                fh.write('ok')
            probe_path.unlink(missing_ok=True)
            return candidate
        except OSError:
            continue

    return Path.cwd()


def resolve_output_path(output_path):
    candidate = Path(output_path)
    candidates = [candidate]
    if candidate.parent != candidate:
        candidates.append(candidate.parent / candidate.name)
    candidates.append(Path.cwd() / candidate.name)
    candidates.append(Path(__file__).resolve().parent / candidate.name)
    candidates.append(Path.cwd() / 'rmsd_output' / candidate.name)
    candidates.append(Path(__file__).resolve().parent / 'rmsd_output' / candidate.name)

    seen = set()
    for path in candidates:
        key = str(path.resolve()) if path.exists() else str(path)
        if key in seen:
            continue
        seen.add(key)
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            probe_path = path.parent / '.write_test.tmp'
            with open(probe_path, 'w', encoding='utf-8') as fh:
                fh.write('ok')
            probe_path.unlink(missing_ok=True)
            return path
        except OSError:
            continue

    fallback = Path.cwd() / 'rmsd_output' / candidate.name
    fallback.parent.mkdir(parents=True, exist_ok=True)
    return fallback


def compute_summary_statistics(values):
    if values.size == 0:
        return {
            'n': 0,
            'mean': float('nan'),
            'median': float('nan'),
            'std': float('nan'),
            'min': float('nan'),
            'max': float('nan'),
            'q1': float('nan'),
            'q3': float('nan'),
            'skewness': float('nan'),
            'kurtosis_excess': float('nan'),
        }

    mean_val = float(np.mean(values))
    median_val = float(np.median(values))
    std_val = float(np.std(values, ddof=0))
    q1_val = float(np.quantile(values, 0.25))
    q3_val = float(np.quantile(values, 0.75))

    if std_val > 0:
        skewness = float(np.mean(((values - mean_val) / std_val) ** 3))
        kurtosis_excess = float(np.mean(((values - mean_val) / std_val) ** 4) - 3.0)
    else:
        skewness = 0.0
        kurtosis_excess = 0.0

    return {
        'n': int(values.size),
        'mean': mean_val,
        'median': median_val,
        'std': std_val,
        'min': float(np.min(values)),
        'max': float(np.max(values)),
        'q1': q1_val,
        'q3': q3_val,
        'skewness': skewness,
        'kurtosis_excess': kurtosis_excess,
    }


def write_statistics_csv(summary_rows, output_path):
    output_path = resolve_output_path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['scenario', 'n', 'mean', 'median', 'std', 'min', 'max', 'q1', 'q3', 'skewness', 'kurtosis_excess'])
        for scenario, stats in summary_rows:
            writer.writerow([
                scenario,
                stats['n'],
                stats['mean'],
                stats['median'],
                stats['std'],
                stats['min'],
                stats['max'],
                stats['q1'],
                stats['q3'],
                stats['skewness'],
                stats['kurtosis_excess'],
            ])
    return output_path


def write_statistics_table_png(summary_rows, output_path):
    output_path = resolve_output_path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    headers = ['scenario', 'n', 'mean', 'median', 'std', 'min', 'max', 'q1', 'q3', 'skewness', 'kurtosis_excess']
    rows = []
    for scenario, stats in summary_rows:
        rows.append([
            scenario,
            stats['n'],
            f"{stats['mean']:.3f}",
            f"{stats['median']:.3f}",
            f"{stats['std']:.3f}",
            f"{stats['min']:.3f}",
            f"{stats['max']:.3f}",
            f"{stats['q1']:.3f}",
            f"{stats['q3']:.3f}",
            f"{stats['skewness']:.3f}",
            f"{stats['kurtosis_excess']:.3f}",
        ])

    fig, ax = plt.subplots(figsize=(14, 2.2 + 0.6 * len(rows)))
    ax.axis('off')

    tbl = ax.table(
        cellText=rows,
        colLabels=headers,
        cellLoc='center',
        loc='center',
        colColours=['#d9edf7'] * len(headers),
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.2)

    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold', color='black')
            cell.set_facecolor('#4f81bd')
            cell.set_text_props(color='white')
        else:
            cell.set_facecolor('#f7f7f7' if row % 2 == 0 else '#ffffff')

    plt.tight_layout()
    try:
        plt.savefig(output_path, dpi=200, bbox_inches='tight')
    except OSError:
        fallback_path = resolve_output_path(output_path.parent / output_path.name)
        plt.savefig(fallback_path, dpi=200, bbox_inches='tight')
        output_path = fallback_path
    plt.close(fig)
    return output_path


def main():
    args = parse_args()
    input_csv = resolve_shared_path(args.input)
    output_dir = resolve_shared_path(args.output_dir)

    if not input_csv.exists():
        raise FileNotFoundError(f'Input file not found: {args.input} (resolved to {input_csv})')
    if output_dir.exists() and not output_dir.is_dir():
        output_dir = output_dir.parent

    output_dir = choose_writable_output_dir(output_dir, input_csv)
    print(f'Using output directory: {output_dir}')

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
    stats_csv = output_dir / 'rmsd_statistics.csv'
    stats_png = output_dir / 'rmsd_statistics_table.png'

    first_path = plot_histogram(frame1_values, first_path, 'Raw Counts: RMSD to Frame 1', bin_width=1.0, use_counts=True)
    print(f'Saved frame-1 histogram to: {first_path}')

    if frame2_values.size > 0:
        second_path = plot_histogram(frame2_values, second_path, 'Raw Counts: RMSD to Frame 2', bin_width=1.0, use_counts=True)
        print(f'Saved frame-2 histogram to: {second_path}')
    else:
        print('No frame-2 values found in the input CSV; skipped frame-2 histogram.')

    all_path = plot_histogram(all_values, all_path, 'Raw Counts: RMSD across All Frame Comparisons', bin_width=1.0, use_counts=True)
    print(f'Saved all-frames histogram to: {all_path}')

    frame1_normal_path = output_dir / 'rmsd_frame1_normal_overlay.png'
    frame2_normal_path = output_dir / 'rmsd_frame2_normal_overlay.png'
    all_normal_path = output_dir / 'rmsd_all_frames_normal_overlay.png'

    frame1_normal_path = plot_histogram(frame1_values, frame1_normal_path, 'Probability Density with Normal Overlay: RMSD to Frame 1', bin_width=1.0, use_density=True, overlay_normal=True)
    print(f'Saved frame-1 normal-overlay histogram to: {frame1_normal_path}')

    if frame2_values.size > 0:
        frame2_normal_path = plot_histogram(frame2_values, frame2_normal_path, 'Probability Density with Normal Overlay: RMSD to Frame 2', bin_width=1.0, use_density=True, overlay_normal=True)
        print(f'Saved frame-2 normal-overlay histogram to: {frame2_normal_path}')
    else:
        print('No frame-2 values found in the input CSV; skipped frame-2 normal-overlay histogram.')

    all_normal_path = plot_histogram(all_values, all_normal_path, 'Probability Density with Normal Overlay: RMSD across All Frame Comparisons', bin_width=1.0, use_density=True, overlay_normal=True)
    print(f'Saved all-frames normal-overlay histogram to: {all_normal_path}')

    all_csv = write_histogram_counts_csv(all_values, all_csv, bin_width=1.0)
    print(f'Saved all-frames counts CSV to: {all_csv}')

    summary_rows = [
        ('frame1_vs_all', compute_summary_statistics(frame1_values)),
        ('frame2_vs_all', compute_summary_statistics(frame2_values)),
        ('all_frame_comparisons', compute_summary_statistics(all_values)),
    ]
    stats_csv = write_statistics_csv(summary_rows, stats_csv)
    stats_png = write_statistics_table_png(summary_rows, stats_png)
    print(f'Saved RMSD summary statistics CSV to: {stats_csv}')
    print(f'Saved screenshotable RMSD statistics table image to: {stats_png}')

    for scenario, stats in summary_rows:
        print(
            f'{scenario}: n={stats["n"]}, mean={stats["mean"]:.3f}, median={stats["median"]:.3f}, '
            f'std={stats["std"]:.3f}, q1={stats["q1"]:.3f}, q3={stats["q3"]:.3f}, '
            f'skew={stats["skewness"]:.3f}, kurtosis_excess={stats["kurtosis_excess"]:.3f}'
        )


if __name__ == '__main__':
    main()
