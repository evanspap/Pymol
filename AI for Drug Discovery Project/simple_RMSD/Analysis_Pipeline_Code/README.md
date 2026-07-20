# Analysis Pipeline (Very Beginner Guide)

This folder helps you do RMSD analysis in 3 simple steps.

## Before you run anything
1. Open PowerShell.
2. Move into this folder.

```powershell
Set-Location "C:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol\AI for Drug Discovery Project\simple_RMSD\Analysis_Pipeline_Code"
```

3. Replace these placeholders in commands below:
- <path-to-your-pdb-file>
- <path-to-output-folder>

---

## Step 1: Calculate RMSD
Main script:
- RMSD_per_frame_biopython_7-19-26.py

This script supports two modes:
1. Target-vs-reference (Nv1): compare one reference frame against one or many target frames.
2. All-vs-all matrix: compare all selected frames against each other.

### Easiest copy-paste commands

```powershell
# 1v1 example: compare frame 1 vs frame 25
python .\RMSD_per_frame_biopython_7-19-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 25
```

```powershell
# 1vN example: compare frame 1 vs frames 2,4,5,9
python .\RMSD_per_frame_biopython_7-19-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 2,4,5,9
```

```powershell
# Nv1 with backbone atoms only
python .\RMSD_per_frame_biopython_7-19-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 2,4,5,9 --atom-scope backbone
```

```powershell
# All-vs-all matrix for selected frames
python .\RMSD_per_frame_biopython_7-19-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --frames 1-100 --all-vs-all
```

```powershell
# All-vs-all matrix for every frame in the file
python .\RMSD_per_frame_biopython_7-19-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all --all-vs-all
```

### Output files from Step 1
1. Nv1 mode output:
- CSV with columns: Frame, RMSD_Angstrom

2. All-vs-all mode output:
- Square RMSD matrix CSV (frame IDs as row/column labels)

---

## Step 2: Classify frames using RMSD threshold
Script:
- Classification_of_Frames_7-19-26.py

Goal:
1. Read RMSD matrix CSV from Step 1.
2. Group frames into classes by threshold.
3. Save class CSV outputs.

---

## Step 3 (Optional): Plot backbone vs full-atom behavior
Script:
- Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py

Goal:
1. Create comparison plots.
2. Save histogram outputs.

---

## If you only want one RMSD value in terminal (no CSV)
Script:
- RMSD_calculation_function_7-19-26.py

```powershell
python .\RMSD_calculation_function_7-19-26.py --input "<path-to-your-pdb-file>" --frame-a 5 --frame-b 12
```

---

## Recommended run order
1. Run Step 1 (RMSD generation)
2. Run Step 2 (classification)
3. Run Step 3 only if you want plots

---

## Troubleshooting in plain English
1. Error says file not found:
- Check your input path text first.

2. Error says no frames or bad frame number:
- Make sure the PDB has multiple frames.
- Make sure frame numbers exist in your file.

3. No output created:
- Confirm output folder path exists (or create it first).
