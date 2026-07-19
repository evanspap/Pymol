# Analysis Pipeline (Beginner Guide)

This folder contains a simple 3-step workflow for RMSD analysis.

## What this pipeline does
1. Calculates RMSD values from a PDB trajectory.
2. Groups frames into classes using RMSD thresholds.
3. Optionally creates comparison plots for backbone vs full-atom RMSD.

---

## Step 1: Calculate RMSD (Biopython)
Use this file:
- RMSD_calculation_biopython_7-19-26.py

Goal:
- Start with your multi-frame PDB file.
- Compute frame-vs-frame RMSD values.
- Save the result as an all-vs-all RMSD matrix CSV.

Output from this step:
- A CSV matrix of RMSD values (this is the input for Step 2).

---

## Step 2: Classify frames by RMSD threshold
Use this file:
- Classification_of_Frames_7-19-26.py

Goal:
- Read the RMSD matrix CSV from Step 1.
- Group frames into classes based on RMSD cutoff rules.
- Export useful CSV outputs (for example, Class 1 members and related submatrices).

Output from this step:
- Class membership files and class-specific RMSD CSV outputs.

---

## Step 3 (Optional): Visualize backbone vs full-atom differences
Use this file:
- Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py

Goal:
- Compare how RMSD distributions look for:
  - Full-atom coordinates
  - Backbone-only coordinates
- Save plot images and/or count summaries.

Output from this step:
- Histogram image(s) and comparison summary files.

---

## Recommended order to run
1. Run RMSD_calculation_biopython_7-19-26.py
2. Run Classification_of_Frames_7-19-26.py
3. (Optional) Run Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py

---

## Beginner tips
- Keep all input and output paths clear and organized.
- Use descriptive output filenames (for example: full vs backbone).
- If a script says a file is missing, check the path first.
- Save the command you used in your notes for reproducibility.

---

## In one sentence
This pipeline takes a PDB trajectory, turns it into an RMSD matrix, classifies frames by similarity, and optionally visualizes backbone-vs-full behavior.
