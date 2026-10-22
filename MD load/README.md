# MD Analysis

A PyMOL plugin that loads a protein and its trajectory as frames (PyMOL states)
in one object. Version 1.0.0 focuses only on loading and viewing.

## Open in PyMOL

Paste these commands into the PyMOL command line:

```pymol
run H:/My Drive/VSCode_Github/MD_classification/Pymol_plugin/__init__.py
md_analysis_gui
```

The window is titled **MD Analysis**. For permanent installation, zip the
`Pymol_plugin` directory (with `__init__.py` directly inside it), install that
ZIP through **Plugin > Plugin Manager > Install New Plugin**, then restart.
Open **Plugin > MD Analysis**.

## Load frames

1. Click **Browse...** beside **Topology (.gro or .pdb)** and select a GRO or PDB file from your filesystem. The selected full path appears in the field.
2. For a multi-model PDB, leave **Trajectory** empty: its models become frames.
3. Otherwise browse to a matching trajectory, such as XTC, TRR or DCD.
4. Optionally choose the first/last trajectory frame and interval. Frame numbers
   start at 1; last frame -1 means all.
5. Click **Load**, then use PyMOL's frame slider or Play.

The protein structure and separate trajectory must have identical atom counts
and atom order. A protein-only structure requires a protein-only trajectory.
Supported topology and trajectory formats depend on the installed PyMOL build.
For a separate trajectory, use a single-state topology. Trajectory frame 1
replaces the topology coordinates, so the topology is not an extra frame.

Repeated loads choose unique object names and preserve existing objects.
A failed load removes the newly created object. No files are written.
The entire selected trajectory is loaded into memory. Large files take time
and the loader window can pause during loading; use an interval or frame range
to load fewer states. This version does not perform alignment or analysis.

**Dry run** validates paths and options without loading coordinates.
It cannot verify atom compatibility or whether PyMOL supports a file format.

## PyMOL commands

```pymol
md_load protein.pdb, trajectory.xtc
md_load ensemble.pdb
md_load protein.pdb, trajectory.dcd, md_protein, 1, 100, 5
```

Syntax:

```text
md_load structure, trajectory="", name=md_protein, start=1, stop=-1, interval=1, dry_run=0
```

For paths containing commas, use Python syntax in PyMOL:

```python
python
md_load("D:/data, project/protein.pdb", "D:/data, project/trajectory.xtc")
python end
```

## Standalone invocation

Use a Python interpreter that provides PyMOL. A real load launches its GUI;
a dry run needs only standard Python.

```powershell
python .\__init__.py --structure "D:\data\protein.pdb" --trajectory "D:\data\trajectory.xtc" --dry-run
```

Run integration checks with PyMOL's Python:

```powershell
& 'C:\Users\geras\AppData\Local\Schrodinger\PyMOL2\python.exe' -B .\test_plugin.py
```

Plugin author: OpenAI Codex, for the MD_classification project.

