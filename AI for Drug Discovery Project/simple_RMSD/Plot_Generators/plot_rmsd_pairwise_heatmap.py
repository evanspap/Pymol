r"""Generate a 2D heatmap from an all-vs-all RMSD matrix CSV.

Version: 1.0.0
Author: Yagnavalkya Devarakonda, Stony Brook University
Dr. Evangelos Papadopoulos, Dana-Farber Cancer Institute
with help from GitHub Copilot
Date: 2026-07-25

Disclaimer:
- Generative AI was used to assist in drafting this script.
- Users should still review plots and validate that the input matrix layout is
  correct for their workflow before interpreting the output scientifically.

Expected input:
- Row 1 contains frame labels for the matrix columns.
- Column 1 contains frame labels for the matrix rows.
- Remaining cells contain pairwise RMSD values.

Typical PowerShell usage:
py -3 .\plot_rmsd_pairwise_heatmap.py \
  --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" \
  --output-dir "C:\path\to\your\plot_output"

Beginner note:
- The horizontal axis and vertical axis both represent frame position in the
  pairwise matrix.
- The color of each square shows the RMSD value for the row-frame versus the
  column-frame comparison.
- Tick marks can be shown at a readable interval so the viewer understands what
  region of the matrix they are looking at, even for very large matrices.

High-level workflow:
1. Read the RMSD matrix from the CSV file.
2. Convert the numeric part of the matrix into a NumPy array.
3. Create a heatmap using matplotlib.
4. Label the axes so the viewer understands what row and column positions mean.
5. Add tick marks and grid lines so large matrices remain readable.
6. Save the result as a PNG image.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    """Define the command-line interface for the heatmap generator.

    The options here are intentionally explicit so users can control:
    - which matrix is plotted,
    - where the image is written,
    - how the plot is colored,
    - whether frame labels are shown,
    - and how densely the axis tick marks are drawn.
    """
    parser = argparse.ArgumentParser(
        description="Generate a 2D heatmap image from an all-vs-all RMSD matrix CSV."
    )
    # The input matrix is required because this plotting script should never be
    # tied to one lab's data path. The user must explicitly provide the matrix.
    parser.add_argument(
        "--input",
        "-i",
        type=Path,
        required=True,
        help="Path to the all-vs-all RMSD matrix CSV.",
    )
    # The output directory is also required so the user decides where plots are
    # written instead of the script assuming a machine-specific location.
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        required=True,
        help="Directory where the heatmap image will be written.",
    )
    # Output filename can be customized in case a team wants naming conventions
    # that match a paper, presentation, or shared project folder.
    parser.add_argument(
        "--output-name",
        default="pairwise_rmsd_heatmap.png",
        help="Filename for the heatmap image. Default: pairwise_rmsd_heatmap.png",
    )
    # Users can override the title when they want a more descriptive label than
    # the default title derived from the input filename.
    parser.add_argument(
        "--title",
        default=None,
        help="Optional custom title for the heatmap. Default: derived from input filename.",
    )
    # The colormap controls how low and high RMSD values are colored.
    parser.add_argument(
        "--colormap",
        default="viridis",
        help="Matplotlib colormap to use. Default: viridis",
    )
    # DPI controls output sharpness. Higher DPI usually means a clearer image.
    parser.add_argument(
        "--dpi",
        type=int,
        default=250,
        help="Output image DPI. Default: 250",
    )
    # Some users want every frame label printed; others want a cleaner plot.
    parser.add_argument(
        "--show-frame-labels",
        choices=("Y", "N"),
        default="N",
        help="Show every frame label on both axes? Use Y or N. Default: N",
    )
    # Tick spacing helps prevent the axis from becoming unreadable on large
    # matrices with hundreds or thousands of frames.
    parser.add_argument(
        "--tick-step",
        type=int,
        default=0,
        help=(
            "Axis tick spacing for large matrices when labels are shown. "
            "Use 0 for automatic spacing. Default: 0"
        ),
    )
    # The user can decide whether the axis text should show simple matrix index
    # positions or the actual frame labels from the CSV file.
    parser.add_argument(
        "--axis-label-mode",
        choices=("index", "frame"),
        default="index",
        help=(
            "Label the axes by numeric position in the matrix ('index') or by "
            "actual frame labels from the CSV ('frame'). Default: index"
        ),
    )
    return parser.parse_args()


def load_pairwise_matrix(csv_path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read the matrix CSV and return row labels, column labels, and numeric values.

    This loader is intentionally strict about the matrix body being numeric so a
    malformed RMSD CSV fails early with a clear error.
    """
    # `utf-8-sig` is used so files with a UTF-8 BOM still read cleanly.
    # `errors="replace"` makes the read more tolerant if the CSV was produced
    # by software that inserted odd characters.
    with open(csv_path, "r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.reader(handle)
        rows = [row for row in reader]

    if not rows:
        raise ValueError("Input CSV is empty.")
    if len(rows[0]) < 2:
        raise ValueError("Input CSV must contain a header row and at least one RMSD column.")

    # The first row defines which frame each matrix column corresponds to.
    column_labels = [cell.strip() for cell in rows[0][1:] if cell.strip()]
    # `row_labels` stores the frame label for each matrix row.
    row_labels: list[str] = []
    # `numeric_rows` stores the numeric RMSD values that will become the heatmap.
    numeric_rows: list[list[float]] = []

    for row_index, raw_row in enumerate(rows[1:], start=2):
        # Skip empty lines completely.
        if not raw_row:
            continue
        frame_label = raw_row[0].strip()
        # Skip rows that do not have a frame label in column 1.
        if not frame_label:
            continue

        # Read exactly as many RMSD cells as the header says should exist.
        values_text = raw_row[1 : 1 + len(column_labels)]
        if len(values_text) != len(column_labels):
            raise ValueError(
                f"Row {row_index} has {len(values_text)} RMSD values but the header expects {len(column_labels)}."
            )

        # Every matrix value must be numeric. If not, the plot should fail with
        # a clear message rather than silently building a wrong heatmap.
        try:
            values = [float(item.strip()) for item in values_text]
        except ValueError as exc:
            raise ValueError(f"Row {row_index} contains a non-numeric RMSD value.") from exc

        row_labels.append(frame_label)
        numeric_rows.append(values)

    # Convert the nested Python list into a NumPy array so matplotlib can plot it.
    matrix = np.array(numeric_rows, dtype=float)
    if matrix.size == 0:
        raise ValueError("No numeric RMSD rows were found in the input CSV.")

    return row_labels, column_labels, matrix


def choose_tick_positions(size: int, requested_step: int) -> np.ndarray:
    """Choose readable tick positions for large matrices.

    Small matrices can show every label. Large matrices need spacing so the plot
    remains legible.
    """
    # For small matrices, every tick can be shown without making the figure too busy.
    if size <= 25:
        return np.arange(size)
    # If the user explicitly requests a step, obey that request.
    if requested_step > 0:
        return np.arange(0, size, requested_step)
    # Default to 10-frame spacing for large matrices so the user can quickly
    # judge where they are in the plot without being overwhelmed by labels.
    auto_step = 10
    return np.arange(0, size, auto_step)


def build_axis_tick_labels(
    labels: list[str],
    tick_positions: np.ndarray,
    axis_label_mode: str,
) -> list[str]:
    """Create the text that will appear at each selected tick position.

    Two display modes are supported:
    - 'index': show 1-based matrix positions so viewers can quickly orient themselves
    - 'frame': show the actual frame labels taken from the CSV
    """
    # When `frame` mode is requested, show the literal frame labels from the CSV.
    if axis_label_mode == "frame":
        return [labels[index] for index in tick_positions]
    # Otherwise, show simple 1-based matrix positions: 1, 2, 3, ...
    return [str(index + 1) for index in tick_positions]


def plot_heatmap(
    row_labels: list[str],
    column_labels: list[str],
    matrix: np.ndarray,
    output_path: Path,
    title: str,
    colormap: str,
    dpi: int,
    show_frame_labels: bool,
    tick_step: int,
    axis_label_mode: str,
) -> Path:
    """Render the RMSD matrix as a 2D heatmap image.

    This function controls the actual appearance of the figure:
    - figure size
    - heatmap colors
    - colorbar label
    - x and y axis titles
    - tick density and tick labels
    """
    # Ensure the destination folder exists before we try to save the PNG.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n_rows, n_cols = matrix.shape
    # Scale the figure roughly with matrix size, while keeping it within a
    # practical range so it does not become absurdly small or large.
    figure_width = max(8.0, min(20.0, n_cols / 45.0))
    figure_height = max(6.0, min(18.0, n_rows / 45.0))

    # Create the matplotlib figure using the size chosen above.
    plt.figure(figsize=(figure_width, figure_height))
    # `imshow` converts the numeric matrix into a color-coded image.
    image = plt.imshow(matrix, cmap=colormap, aspect="auto", interpolation="nearest")
    colorbar = plt.colorbar(image)
    colorbar.set_label("RMSD (Å)")

    # The title explains what file or dataset the viewer is looking at.
    plt.title(title)
    # The horizontal axis tells the viewer which column-frame is being compared.
    plt.xlabel("Column frame index in the pairwise RMSD matrix (the frame each row-frame is compared against)")
    # The vertical axis tells the viewer which row-frame is currently being evaluated.
    plt.ylabel("Row frame index in the pairwise RMSD matrix")

    if show_frame_labels:
        # Select readable tick positions first, then create labels from those positions.
        x_ticks = choose_tick_positions(len(column_labels), tick_step)
        y_ticks = choose_tick_positions(len(row_labels), tick_step)
        x_tick_labels = build_axis_tick_labels(column_labels, x_ticks, axis_label_mode)
        y_tick_labels = build_axis_tick_labels(row_labels, y_ticks, axis_label_mode)
        plt.xticks(x_ticks, x_tick_labels, rotation=90, fontsize=7)
        plt.yticks(y_ticks, y_tick_labels, fontsize=7)
    else:
        # Even when frame labels are hidden, we still show numeric index ticks so
        # the viewer can tell which region of the matrix is being inspected.
        x_ticks = choose_tick_positions(len(column_labels), tick_step)
        y_ticks = choose_tick_positions(len(row_labels), tick_step)
        plt.xticks(x_ticks, [str(index + 1) for index in x_ticks], fontsize=8)
        plt.yticks(y_ticks, [str(index + 1) for index in y_ticks], fontsize=8)

    # Add visible demarcation lines so viewers can track frame position across
    # the matrix more easily. We place them every major tick interval.
    x_tick_step = tick_step if tick_step > 0 else 10
    y_tick_step = tick_step if tick_step > 0 else 10
    x_minor_ticks = np.arange(-0.5, n_cols, x_tick_step)
    y_minor_ticks = np.arange(-0.5, n_rows, y_tick_step)
    plt.gca().set_xticks(x_minor_ticks, minor=True)
    plt.gca().set_yticks(y_minor_ticks, minor=True)
    # White minor-grid lines overlay the heatmap without visually overpowering it.
    plt.grid(which="minor", color="white", linestyle="-", linewidth=0.4, alpha=0.8)
    plt.tick_params(which="minor", bottom=False, left=False)

    # Tight layout reduces the chance that labels or titles are clipped.
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return output_path


def main() -> int:
    """Read the input matrix, create the heatmap, and report the output path.

    This is the end-to-end driver that beginners can mentally follow as:
    - parse inputs
    - load matrix
    - decide plot title/output path
    - draw heatmap
    - print saved file information
    """
    args = parse_args()
    input_csv = Path(args.input)
    output_dir = Path(args.output_dir)

    # Fail early if the user gave a non-existent file path.
    if not input_csv.is_file():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    # Read and validate the matrix before any plotting is attempted.
    row_labels, column_labels, matrix = load_pairwise_matrix(input_csv)

    # If the user did not provide a title, derive a readable one from the file name.
    title = args.title or f"Pairwise RMSD Heatmap\n{input_csv.stem}"
    output_path = output_dir / args.output_name
    # Hand everything to the plotting function that renders and saves the figure.
    saved_path = plot_heatmap(
        row_labels=row_labels,
        column_labels=column_labels,
        matrix=matrix,
        output_path=output_path,
        title=title,
        colormap=args.colormap,
        dpi=args.dpi,
        show_frame_labels=args.show_frame_labels == "Y",
        tick_step=args.tick_step,
        axis_label_mode=args.axis_label_mode,
    )

    # Print a short success summary so a terminal user knows where the output went.
    print(f"Saved heatmap to: {saved_path}")
    print(f"Matrix size: {matrix.shape[0]} rows x {matrix.shape[1]} columns")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())