# Save Selections PyMOL Plugin

This tool is part of the `Protein_Tools` package. Install the parent folder
once through PyMOL's Plugin Manager; do not install this subfolder separately.

Run `save_selections` to open a checklist of all named selections. The list is
independent of enabled state and of row highlights in PyMOL's object panel.
Check one or more names, optionally use **Select All** or **Select None**, choose
an output folder, and click **Save CSVs**. From the command line:

```pymol
save_selections /path/to/output
```

Omit the directory to use PyMOL's current working directory. Each selected name
is written to its own `<selection>_residues.csv` file. The CSV columns are
`Chain,Res_Num,Res_Name`, with one row per CA atom. Selections without CA atoms
produce a header-only CSV. Unsafe filename characters are replaced with
underscores; collisions receive numeric suffixes.

The dialog is also available under **Plugin > Protein Tools > Choose Selections
to Save...**. Enabled molecular objects are exported separately by `save_objects`.

Run `python3 -m unittest discover -s Protein_Tools/save_selection -p
'test_plugin.py'` from the repository root to execute its mocked tests.