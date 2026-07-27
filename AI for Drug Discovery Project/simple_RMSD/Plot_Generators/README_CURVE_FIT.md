# Curve_Fit Cookbook (Beginner Friendly)

## Overview

This guide explains how to use:

- Curve_Fit .py

This script reads a multi-frame PDB file, computes RMSD distributions, builds probability histograms, and overlays Gaussian curve fits.

The script now uses adaptive model selection so the fit can follow:

- unimodal symmetric shapes
- unimodal skewed shapes
- bimodal shapes

It is designed to be reusable for any team and is not hardcoded to lab-specific file paths.

## Important Disclaimer

Generative AI may have been used to assist in writing this script and guide.
Always validate results before scientific interpretation.

## What This Script Does

From a PDB trajectory file, this script can generate:

- full-atom all-frames-vs-all-frames RMSD histogram
- backbone all-frames-vs-all-frames RMSD histogram
- custom atom-selection all-frames-vs-all-frames RMSD histogram
- adaptive fit overlay for each histogram
- combined figure when multiple comparisons are selected
- optional pairwise RMSD CSV tables

## Key Requirement Covered

- Histogram x-axis bin width defaults to 0.2 Angstrom.
- Y-axis is probability per bin.
- X-axis range automatically extends to include the full observed RMSD range across selected comparisons.

## Beginner Concept: How The Curve Is Chosen

The script does not force one single curve type.
Instead, it tries multiple fit families and chooses the best one.

Models considered:

- gaussian
- skew_gaussian
- gaussian_mixture_2

You can also force behavior with --fit-mode:

- adaptive (default)
- gaussian (force Gaussian-only)
- none (disable fit overlay)

Selection method:

- the script compares model quality using BIC (Bayesian Information Criterion)
- lower BIC is better
- extra guardrails reduce false bimodal overfitting

In plain language:

- if your histogram is roughly bell-shaped and symmetric, you will usually get gaussian
- if your histogram is one peak with a long tail, you will usually get skew_gaussian
- if your histogram has two peaks, you can get gaussian_mixture_2

## Script Location

Plot_Generators/Curve_Fit .py

## Quick Start

Open PowerShell in the Plot_Generators folder and run:

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full backbone --bin-width 0.2
```

## Input And Output

### Required input

- --pdb : path to a multi-frame PDB file

### Typical output files

- rmsd_hist_gaussian_full.png
- rmsd_hist_gaussian_backbone.png
- rmsd_hist_gaussian_<custom_label>.png
- rmsd_hist_gaussian_combined.png
- rmsd_gaussian_fit_summary.csv
- rmsd_pairs_<comparison>.csv (only if requested)

## Main CLI Options

### --pdb

Path to multi-frame PDB input.

### --output-dir

Where all output files are saved.

### --comparisons

Pick one or more:

- full
- backbone
- custom

### --custom-selection

Defines custom atom groups by atom names.
Format:

- Label=ATOM1,ATOM2,ATOM3

Example:

- SideChain=CB,CG,CD

You can repeat --custom-selection multiple times.

### --bin-width

Histogram bin width in Angstrom.
Default is 0.2.

### --x-max-angstrom

Optional fixed x-axis maximum (Angstrom).

- If omitted, the script auto-uses the maximum RMSD observed in the run.
- If provided, the script uses at least that value (useful for standardizing plots across runs).

### Fit summary CSV columns (important)

The summary CSV includes these fit-specific columns:

- FitModel
- FitBIC
- FitAmplitude
- FitMu
- FitSigma

Interpretation tips:

- FitModel tells you which distribution family was selected
- FitBIC is useful when comparing model quality (lower is better)
- FitAmplitude/FitMu/FitSigma are standardized descriptors of the selected curve

### --include-diagonal Y|N

- Y includes frame-vs-same-frame zeros in the distribution.
- N excludes same-frame pairs.

### --write-pairwise-csv Y|N

- Y writes per-comparison pair tables.
- N skips pair-table CSV files.

### --fit-mode adaptive|gaussian|none

Controls whether and how a fit curve is drawn.

- adaptive: model is selected automatically
- gaussian: force Gaussian-only fit
- none: do not draw fit curve

## Cookbook Scenarios

### Scenario 1: Full only

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full --bin-width 0.2
```

### Scenario 2: Backbone only

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons backbone --bin-width 0.2
```

### Scenario 3: Full + backbone together

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full backbone --bin-width 0.2
```

### Scenario 4: Custom atom group only

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons custom --custom-selection SideChain=CB,CG,CD --bin-width 0.2
```

### Scenario 5: Multiple custom groups

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons custom --custom-selection SideChain=CB,CG,CD --custom-selection Polar=N,O,S --bin-width 0.2
```

### Scenario 6: Full + backbone + custom in one run

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full backbone custom --custom-selection SideChain=CB,CG,CD --bin-width 0.2
```

### Scenario 7: Exclude diagonal and save pairwise tables

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full backbone --bin-width 0.2 --include-diagonal N --write-pairwise-csv Y
```

### Scenario 8: Force a standard x-axis maximum for cross-run comparison

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full --bin-width 0.2 --x-max-angstrom 12.0
```

### Scenario 9: Force Gaussian-only fit

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full --bin-width 0.2 --fit-mode gaussian
```

### Scenario 10: Disable fit line (histogram only)

```powershell
py -3 ".\Curve_Fit .py" --pdb "C:\path\to\your\frames.pdb" --output-dir "C:\path\to\your\curve_fit_output" --comparisons full --bin-width 0.2 --fit-mode none
```

## How To Read The Output

- Bars: probability mass in each RMSD bin.
- Red fit line: adaptive model approximation of the distribution.
- Fit summary CSV: numerical fit parameters for each comparison.

### Reading fit model names in the legend

- gaussian fit: symmetric one-peak fit
- skew_gaussian fit: one-peak but skewed fit
- gaussian_mixture_2 fit: two-component mixture fit

If --fit-mode none is used, no fit line or fit legend entries are drawn.

## Troubleshooting

### Problem: PDB file not found

Cause:
Path in --pdb is wrong.

Fix:
Double-check path spelling and quotes.

### Problem: Not enough frames

Cause:
PDB has fewer than 2 parsed frames.

Fix:
Use a multi-frame PDB with frame separators (for example MODEL/ENDMDL or END).

### Problem: Custom selection fails

Cause:
Format is incorrect or atom names are missing in frames.

Fix:
Use Label=ATOM1,ATOM2 and verify those atom names exist in your PDB.

### Problem: Script runs but curve fit looks rough

Cause:
Very small datasets can produce unstable fit estimation.

Fix:
Use more frames or compare a broader structural region.

### Problem: Model selected is not what you expected

Cause:
Your histogram shape may not strongly support the model you expected, or differences are subtle.

Fix:

- inspect the histogram bars and fitted line visually
- compare FitBIC values across repeated runs or parameter choices
- use more frames when possible for more stable distribution estimation

## Help Command

```powershell
py -3 ".\Curve_Fit .py" --help
```
