# Save Enabled Objects PyMOL Plugin

Part of the `Protein_Tools` package. Install the parent folder once through
PyMOL's Plugin Manager; do not install this subfolder separately.

Run `save_objects` to export each enabled molecular object to its own CSV in
the current working directory. Disabled objects and non-molecular groups are
skipped. To choose another folder, run:

```pymol
save_objects /path/to/output
```

Omit the directory to use PyMOL's current working directory. Each
`<object>_residues.csv` contains `Chain,Res_Num,Res_Name` columns and one row per
CA atom. Objects without CA atoms produce a header-only CSV. The action is also
available under **Plugin > Protein Tools > Save Enabled Objects (CSV)**.

Run `python3 -m unittest discover -s Protein_Tools/save_objects -p
'test_plugin.py'` from the repository root to execute the mocked tests.