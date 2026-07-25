# Plot Generators

This folder contains plotting scripts for visualizing RMSD outputs.

## Files
- plot_rmsd_backbone_vs_full_histograms.py
- plot_rmsd_category_histogram.py
- plot_rmsd_full_vs_backbone_histograms.py
- plot_rmsd_matrix_histograms.py
- plot_rmsd_pairwise_heatmap.py

## PowerShell Quick Start
1. Move to this folder:

```powershell
Set-Location "C:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol\AI for Drug Discovery Project\simple_RMSD\Plot_Generators"
```

2. Example run (matrix histogram):

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\plot_rmsd_matrix_histograms.py
```

3. Example run (category histogram):

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\plot_rmsd_category_histogram.py
```

4. Example run (full vs backbone histogram):

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\plot_rmsd_full_vs_backbone_histograms.py
```

5. Example run (pairwise RMSD heatmap):

```powershell
py -3 .\plot_rmsd_pairwise_heatmap.py --input "C:\path\to\your\all-vs-all_rmsd_matrix.csv" --output-dir "C:\path\to\your\plot_output"
```

## Notes
- Scripts may require CSV/RMSD outputs generated beforehand; check each script's expected input path.
- If this repository is cloned from GitHub on another machine, local paths will differ.
- Replace the paths in these commands with paths on your own system.
- If Python is already on PATH, you can replace the long interpreter path with python.
