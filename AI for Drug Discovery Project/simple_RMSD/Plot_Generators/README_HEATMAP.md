# RMSD Pairwise Heatmap Generator

## Overview

This guide explains how to use:

- `plot_rmsd_pairwise_heatmap.py`

This script reads an all-vs-all RMSD matrix CSV and creates a 2D heatmap image.
It is intended to be reusable for any team working with RMSD matrix data.

## Important Disclaimer

Generative AI was used to assist in drafting this plotting script and its documentation.
Users should still review the input data, inspect the output image carefully, and validate that the plotted matrix is the correct one before drawing scientific conclusions.

## What This Script Does

The script takes a pairwise RMSD matrix and converts it into a heatmap image.

In the heatmap:

- the horizontal axis represents the matrix columns
- the vertical axis represents the matrix rows
- each colored square represents the RMSD value for one row-frame versus one column-frame comparison
- the color bar shows the RMSD scale in angstroms

## What Kind Of Input File Is Required

The script expects a matrix-style CSV file where:

- row 1 contains the frame labels for the matrix columns
- column 1 contains the frame labels for the matrix rows
- the rest of the cells contain numeric RMSD values

This is the standard output format for an all-vs-all RMSD matrix.

## What The Axes Mean

### Horizontal axis

The horizontal axis represents the column frame position in the pairwise matrix.
This is the frame that the current row frame is being compared against.

### Vertical axis

The vertical axis represents the row frame position in the pairwise matrix.
This is the frame whose RMSD is being measured against the frame on the horizontal axis.

### Color

The color of each square tells you the magnitude of the RMSD value.
Smaller RMSD values and larger RMSD values are shown as different colors according to the selected colormap.

## Why Tick Marks And Grid Lines Matter

For large matrices, it becomes hard to see where you are in the plot.
This script includes:

- visible axis tick marks
- readable axis labels
- minor grid lines
- regular demarcation spacing

By default, large matrices use 10-frame spacing so a user can visually track position across the heatmap.

## Script Location

The script lives here:

- `Plot_Generators\plot_rmsd_pairwise_heatmap.py`

## Quick Start

Open PowerShell in the `Plot_Generators` folder and run:

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output"
```

## Most Important Inputs

### `--input`

Path to the all-vs-all RMSD matrix CSV.

### `--output-dir`

Folder where the heatmap PNG file will be written.

## Optional Inputs

### `--output-name`

Lets you choose a custom output filename.

Example:

```powershell
--output-name "apoa_pairwise_heatmap.png"
```

### `--title`

Lets you set a custom plot title.

Example:

```powershell
--title "ApoA Run1 Pairwise RMSD Heatmap"
```

### `--colormap`

Lets you choose a Matplotlib colormap.

Examples:

- `viridis`
- `plasma`
- `inferno`
- `magma`
- `coolwarm`

### `--dpi`

Controls the output image quality.
Higher values make a sharper image but increase file size.

### `--show-frame-labels Y|N`

Controls whether the axis labels show the actual frame labels from the CSV.

- `Y`: show labels
- `N`: do not show labels

For very large matrices, `N` is usually easier to read.

### `--tick-step`

Controls how often tick labels appear.

Examples:

- `10` means every 10 frames
- `50` means every 50 frames
- `0` means let the script decide automatically

### `--axis-label-mode index|frame`

Controls whether shown labels represent:

- `index`: 1-based numeric position in the matrix
- `frame`: actual frame labels from the CSV

For large matrices, `index` is usually clearer.

## Cookbook

### Scenario 1: Simplest possible run

Use this when you just want a heatmap from your matrix.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output"
```

### Scenario 2: Custom output filename

Use this when you want the PNG to have a meaningful name.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --output-name "apoa_run1_heatmap.png"
```

### Scenario 3: Custom plot title

Use this when you want the plot title to describe the dataset directly.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --title "ApoA Run1 Pairwise RMSD Heatmap"
```

### Scenario 4: Show frame labels on the axes

Use this when the matrix is small enough that full labels will still be readable.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --show-frame-labels Y --axis-label-mode frame
```

### Scenario 5: Keep labels simple by using matrix index positions

Use this when you want labels but the actual frame IDs are too crowded.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --show-frame-labels Y --axis-label-mode index --tick-step 25
```

### Scenario 6: Large matrix with readable demarcations

Use this for large pairwise matrices where you want the plot to stay readable.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --show-frame-labels N --tick-step 10
```

### Scenario 7: Use a different colormap

Use this when you want a different visual style for presentations or comparison plots.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --colormap plasma
```

### Scenario 8: Create a higher-resolution image

Use this when you want sharper output for slides or figures.

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --dpi 400
```

## Example Using The ApoA Matrix

If you want to test it with the ApoA matrix you used earlier, this is the style of command:

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "G:\path\to\your\all 1000 frames all vs all RMSD 7-22-26.csv" --output-dir "C:\path\to\your\plot_output"
```

This example is still intentionally generic in spirit. Users should replace the paths with whatever locations exist on their own system.

## Output Produced

The script writes a PNG file.

Default filename:

- `pairwise_rmsd_heatmap.png`

The image is saved in the folder given by `--output-dir`.

## Typical Interpretation

When looking at the heatmap:

- diagonal values should generally be near zero
- blocks of similar color may indicate clusters of similar frames
- sharp color changes may indicate transitions between structural groups
- the axis tick spacing helps the viewer identify where in the matrix these patterns occur

## Troubleshooting

### Problem: "Input CSV not found"

Cause:
The path provided to `--input` is incorrect.

Fix:
Double-check the file path and ensure the file exists.

### Problem: Heatmap looks too crowded

Fix:
- use `--show-frame-labels N`
- increase `--tick-step`
- use index labels instead of frame labels

Example:

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output" --show-frame-labels N --tick-step 25
```

### Problem: Labels are too small

Fix:
Hide full frame labels and rely on numeric tick marks.

### Problem: The matrix does not load

Cause:
The CSV may not be formatted like a pairwise RMSD matrix.

Fix:
Check that:

- row 1 is the header row
- column 1 is the row-label column
- matrix cells are numeric

## Help Command

Use this to see the script options in PowerShell:

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py -h
```
