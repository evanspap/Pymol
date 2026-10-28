# README: `Pymol_plugin/__init__.py` - MD Analysis Plugin

This guide describes the plugin implementation in [`__init__.py`](__init__.py).
The plugin loads a protein structure and, optionally, its molecular-dynamics
trajectory into PyMOL. Trajectory frames are represented as states of one PyMOL
object so they can be viewed with the frame controls or animation playback.

**Version:** 1.0.2

**Author:** OpenAI Codex, for the MD_classification project.

## Scope and Limitations

The plugin loads and views structures and trajectories, with optional
per-frame principal-component alignment. This is not reference-structure global
alignment, local alignment, or RMSD analysis.

The plugin does not write or modify input files. PyMOL reads the selected data
into memory. Large trajectories may take time to load and may use substantial
memory; use the frame range and interval controls to reduce the number of
loaded frames when appropriate.

## Requirements

- PyMOL with the `pymol` and `pymol.Qt` modules available.
- The file formats selected must be supported by the installed PyMOL build.
- When a separate trajectory is provided, the structure/topology and trajectory
   must have the same atom count and atom order. For example, a protein-only
   topology must be paired with a trajectory containing the same protein atoms.
- A separate trajectory requires a single-state structure/topology. A
   multi-model PDB is loaded directly and should be supplied without a separate
   trajectory.

## Open the Plugin in PyMOL

For a one-time session, run the plugin file and open its loader window from the
PyMOL command line:

```pymol
run /path/to/MD_classification/Pymol_plugin/__init__.py
md_analysis_gui
```

Replace the example path with the location of this repository on your computer.
The window is titled **MD Analysis**.

For permanent installation, create a ZIP archive of the `Pymol_plugin` folder
with `__init__.py` directly inside the archive's plugin folder. In PyMOL, use
**Plugin > Plugin Manager > Install New Plugin**, select the ZIP, and restart
PyMOL if requested. Open the installed plugin from **Plugin > MD Analysis**.

## Load Data with the GUI

1. In **Topology (.gro or .pdb)**, browse to a structure/topology file.
2. For a multi-model PDB, leave **Trajectory** blank. Each `MODEL` becomes a
    state in the PyMOL object.
3. For a separate trajectory, browse to the matching trajectory file. Examples
    include XTC, TRR, and DCD; the exact supported formats depend on PyMOL.
4. Optionally enable **Principal-axis alignment** to center each frame's protein
   at `(0, 0, 0)` and align its major, middle, and minor axes to X, Y, and Z.
5. Set the object name if desired. If that name is already in use, the plugin
    chooses a suffixed unique name rather than replacing the existing object.
6. Optionally set the first frame, last frame, and frame interval. These controls
    are enabled only when a separate trajectory is selected.
7. Use **Dry run** to validate paths and options without loading coordinates, or
    use **Load** to load the selected data.
8. Use PyMOL's frame slider or **Play** to view the states.

Frame numbers are one-based. A last-frame value of `-1` means all remaining
frames. An interval of `5`, for example, loads every fifth frame in the selected
range.

## PyMOL Commands

The registered command signature is:

```text
md_load structure, trajectory="", name=md_protein, start=1, stop=-1, interval=1, dry_run=0, align_principal_axes=0
```

Examples:

```pymol
md_load protein.pdb, trajectory.xtc
md_load protein.pdb, trajectory.xtc, align_principal_axes=1
md_load ensemble.pdb
md_load protein.pdb, trajectory.dcd, md_protein, 1, 100, 5
```

The first example loads a structure and a matching trajectory. The second loads
all models from a multi-model PDB. The third loads trajectory frames 1 through
100, taking every fifth frame.

Principal-axis alignment is off by default. When enabled, each state's protein
coordinates are centered at the origin and its principal axes are mapped from
largest to smallest variance onto X, Y, and Z. Axis signs are kept consistent
between sequential states where possible. Non-protein atoms in the same object
receive the same translation and rotation as the protein.

The PyMOL command parser treats commas as argument separators. For paths that
contain commas, call the function using Python syntax instead:

```python
python
md_load("D:/data, project/protein.pdb", "D:/data, project/trajectory.xtc")
python end
```

### Load and Cleanup Behavior

- The structure is loaded as one object. PDB model records are retained as
   states rather than split into separate objects.
- A separate trajectory is attached to the single-state structure object. The
   structure coordinates are replaced by the trajectory's first selected frame;
   the topology snapshot is not an extra trajectory frame.
- Existing objects are not overwritten. The requested object name is given a
   numeric suffix when that name is already present.
- If loading raises an exception after creating a new object, that new object is
   deleted. Objects that existed before the load are retained.
- After loading, protein atoms are shown as a cartoon and non-protein atoms as
   lines. The first state is selected and the camera is zoomed to the object.

## Dry Run

Dry run validates that the required input files exist and that the object name
and frame options are valid. It prints the command and reports successful
validation without importing PyMOL in the standalone CLI path or changing the
PyMOL scene.

Dry run cannot verify that the files can be decoded by the installed PyMOL, or
that a topology and trajectory actually have compatible atoms. Those checks
require a real load.

## Standalone Command-Line Use

Run `__init__.py` with a Python interpreter that can import PyMOL. A real load
initializes PyMOL; a dry run only checks paths and options.

```powershell
python .\__init__.py --structure "D:\data\protein.pdb" --trajectory "D:\data\trajectory.xtc" --dry-run
python .\__init__.py --structure "D:\data\protein.pdb" --trajectory "D:\data\trajectory.xtc" --align-principal-axes
python .\__init__.py --structure "D:\data\ensemble.pdb" --name ensemble
python .\__init__.py --help
python .\__init__.py --version
```

Available options:

| Option | Meaning | Default |
| --- | --- | --- |
| `--structure` | Required structure/topology path | None |
| `--trajectory` | Optional matching trajectory path | Empty |
| `--name` | Requested PyMOL object name | `md_protein` |
| `--start` | First trajectory frame, one-based | `1` |
| `--stop` | Last frame; `-1` means all remaining frames | `-1` |
| `--interval` | Load every Nth frame | `1` |
| `--align-principal-axes` | Center and align each state's protein using PCA | Off |
| `--dry-run` | Validate without loading coordinates | Off |
| `--version` | Print plugin version and exit | - |

## Implementation Map

- `_command` formats a command line for status output using the current
   operating system's quoting convention.
- `_options` resolves paths and validates the name and frame options. The GUI,
   CLI, and loader share this function.
- `md_load` is the central loader used by both the GUI and PyMOL command. It
   performs validation, dry-run handling, unique-name selection, PyMOL loading,
   display setup, and cleanup after a failed load.
- `md_analysis_gui` builds the file-picker window and delegates actual work to
   `md_load`.
- `_register` adds `md_load` and `md_analysis_gui` to PyMOL's command namespace
   when PyMOL is available.
- `__init_plugin__` is the PyMOL Plugin Manager hook and adds the menu item.
- `main` implements standalone command-line use.

The integration checks for this plugin are documented separately in
[`README_test_plugin.md`](README_test_plugin.md) and implemented in
[`test_plugin.py`](test_plugin.py).

