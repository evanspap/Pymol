# RMSD Calculator for Multi-Frame PDB Files

## Overview
This script calculates the RMSD (Root Mean Square Deviation) between Frame 1 and every other frame in a multi-frame PDB file. Each frame should be separated by an `END` keyword in the PDB file.

## Files
- `calculate_rmsd.py` - Main RMSD calculation script
- `rmsd_results.csv` - Output file with RMSD values for each frame

## Requirements
- Python 3.7+
- NumPy

## Installation

### Install NumPy
```bash
pip install numpy
```

## Usage

### Run the script
```bash
python calculate_rmsd.py
```

### Output
The script will:
1. Parse the PDB file and extract all frames
2. Print RMSD values to the console
3. Save results to `rmsd_results.csv` with columns:
   - `Frame`: Frame number (1-indexed)
   - `RMSD_vs_Frame1`: RMSD relative to Frame 1 (in Ångströms)
   - `RMSD_vs_Frame2`: RMSD relative to Frame 2 (in Ångströms)

## Example Output
```
Reading PDB file: G:\...apoa__run1_frames_1k.pdb
Parsing frames...
Successfully parsed 1000 frames

Frame 1 has 2500 atoms
Frame 2 has 2500 atoms

============================================================
RMSD Calculation: Frame 1 vs Frame 2
============================================================
Frame    RMSD (Å)       Description
------------------------------------------------------------
1        0.0000         vs Frame 1
2        1.2345         vs Frame 1
3        1.5432         vs Frame 1
...
```

## How it Works
1. **Parse PDB**: Reads ATOM and HETATM records from the PDB file
2. **Extract Coordinates**: Extracts X, Y, Z coordinates from each atom
3. **Calculate RMSD**: For each frame, computes: RMSD = √(mean((Δx² + Δy² + Δz²)))
4. **Save Results**: Exports results to CSV for further analysis

## Notes
- The script uses the standard RMSD formula without alignment (translation/rotation)
- For aligned RMSD, consider using PyMOL or MDAnalysis libraries
- The input file path is hardcoded in the script; modify it to analyze different PDB files
