# RMSD Calculator

This folder contains scripts to compute all-vs-all RMSD matrices from multi-frame PDB trajectory files.

## Files
- calculate_rmsd.py
- calculate_rmsd_biopython.py
- calculate_rmsd_biopython_true.py

## Recommended Script
Use calculate_rmsd_biopython_true.py for the most robust parser-based workflow.

## PowerShell Quick Start
1. Move to this folder:

```powershell
Set-Location "C:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol\AI for Drug Discovery Project\simple_RMSD\RMSD_Calculator"
```

2. Run with your PDB file:

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\calculate_rmsd_biopython_true.py -i "C:\Path\To\your_frames.pdb"
```

3. Optional custom output location and name:

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\calculate_rmsd_biopython_true.py -i "C:\Path\To\your_frames.pdb" -o "C:\Path\To\results" --output-name "team_rmsd.csv"
```

4. Optional backbone-only mode and frame limit:

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\calculate_rmsd_biopython_true.py -i "C:\Path\To\your_frames.pdb" --atom-scope backbone --max-frames 500
```

5. See all options:

```powershell
& "c:/Users/yagd9/Documents/Stony Brook AI Drug Discovery/Github/Pymol/AI for Drug Discovery Project/.venv/Scripts/python.exe" .\calculate_rmsd_biopython_true.py --help
```

## Notes
- If you cloned this repository from GitHub, your local path will likely be different.
- Replace script and data paths in the commands above with your local paths.
- If your environment already has Python on PATH, you can replace the long interpreter path with python.
