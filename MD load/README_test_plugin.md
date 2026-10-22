# README: `Pymol_plugin/test_plugin.py` - Integration Tests

This guide describes the integration test script in [`test_plugin.py`](test_plugin.py).
The tests exercise the real MD Analysis plugin in [`__init__.py`](__init__.py)
through PyMOL's command API. They use small generated input files rather than
the research project's full structures or trajectories.

## Purpose

The test suite checks that the current loader:

- Loads a multi-model PDB as multiple states of one PyMOL object.
- Loads a separate DCD trajectory onto a matching single-state structure.
- Selects the requested range of trajectory frames.
- Centers and aligns every loaded state to its protein's principal axes when requested.
- Chooses a unique PyMOL object name when the requested name is already in use.
- Removes an object created during a failed load without deleting older objects.
- Validates dry-run behavior and invalid options.
- Prints the expected command-line help/header and exits with the expected status.

These are integration tests, not mocked unit tests. They start an embedded
PyMOL instance and load actual structures and trajectory data. They test
principal-axis alignment, not alignment to a separate reference structure.

## Requirements

- A Python interpreter that can import both `pymol` and `pymol2`. Use the Python
  environment supplied with your PyMOL installation, or an environment where
  those modules are installed and configured.
- PyMOL must be able to create a `SingletonPyMOL` instance and read the test's
  PDB and DCD files.
- The Python standard library modules used by the script (`unittest`,
  `tempfile`, `struct`, and others) are included with Python.

A regular Python installation without PyMOL will stop at the top-level PyMOL
imports before any test runs. The CLI dry-run case does not need PyMOL on its
own, but the complete test module does because the other tests use the embedded
PyMOL API.

## Run the Tests

Run the command from the `MD_classification` repository root. Replace the
example interpreter path with the Python executable used by your PyMOL
installation:

```powershell
& "C:\path\to\pymol\python.exe" -B -m unittest Pymol_plugin.test_plugin
```

Alternatively, change into the plugin directory and run the script directly:

```powershell
Set-Location "C:\path\to\MD_classification\Pymol_plugin"
& "C:\path\to\pymol\python.exe" -B .\test_plugin.py
```

The `-B` option tells Python not to write `.pyc` bytecode cache files. The
script uses `unittest`, so a successful run prints a summary of tests and ends
with `OK`. A failing test prints a traceback and a failure/error summary.

If PowerShell reports that `python` cannot be found, use the full path to the
PyMOL Python executable as shown above, or activate the environment that
contains PyMOL before running the command.

## Generated Test Data

Each test gets a new temporary directory created by `tempfile.TemporaryDirectory`.
The script generates three small inputs there:

- `protein with spaces.pdb`: a single-state structure with three alanine CA
  atoms. The path deliberately contains spaces to exercise path handling.
- `ensemble.pdb`: a three-model PDB. The first CA atom has x coordinates 0, 5,
  and 10 in models 1, 2, and 3.
- `trajectory.dcd`: a minimal little-endian CHARMM DCD with three atoms and
  three frames. Its corresponding first-atom x coordinates are 0, 5, and 10.

`pdb_frame` formats the synthetic PDB records. `write_dcd` writes the binary DCD
records and their byte-count markers. No research data is read or changed.

## Test Cases

| Test method | Behavior checked |
| --- | --- |
| `test_multimodel_and_unique_names` | Multi-model PDB states, unique repeated names, and final-state coordinates. |
| `test_separate_trajectory` | Three DCD frames load as states and retain the expected final coordinate. |
| `test_principal_axis_alignment_centers_and_orders_each_state` | Each trajectory state is centered and has major-to-minor variance along X, Y, and Z. |
| `test_trajectory_frame_selection` | Selecting frames 2-3 yields two states, with original frame 2 first. |
| `test_failed_load_removes_only_new_object` | A failed load preserves an existing object and cleans up the failed object's addition. |
| `test_dry_run` | Dry run prints status and does not change the PyMOL object set. |
| `test_bad_options` | Invalid start, stop, interval, and object name raise `ValueError`. |
| `test_cli_header_and_dry_run` | Missing CLI arguments show documentation and exit 2; a valid CLI dry run exits 0. |

## Isolation and Cleanup

`LoaderTests.setUp` creates the temporary fixtures and records the PyMOL object
names that existed before a test. `tearDown` deletes only objects added by that
test and then removes the temporary directory. This allows the tests to share
the embedded PyMOL instance without intentionally deleting objects that were
already present before a test began.

The suite uses `unittest` assertions and launches the plugin CLI checks with
`subprocess.run`. CLI subprocesses use the same Python executable as the test
process, which is why the test process must use a PyMOL-capable interpreter.

## Interpreting Failures

- Import errors for `pymol` or `pymol2` usually mean the wrong Python
  interpreter is being used.
- DCD loading failures may indicate that the active PyMOL build does not support
  this DCD variant or cannot initialize in the current environment.
- Coordinate/state assertion failures indicate that the loaded object differs
  from the known three-frame fixtures or expected frame mapping.
- CLI return-code failures include captured standard error in the assertion
  message when available.

The test suite verifies the current structure and trajectory loading contract.
When a new behavior such as global alignment is implemented, add focused tests
for its reference state, atom selection, coordinate changes, and unchanged
loading behavior here.