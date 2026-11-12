# Protein Tools PyMOL Plugin

## Install

Install the `Protein_Tools` folder once through PyMOL's Plugin Manager, then
restart PyMOL. The package loader initializes all child plugins and places their
menu actions under **Plugin > Protein Tools**. Do not install the child folders
individually.

## Commands

The package provides four PyMOL commands:

```pymol
label_around [selection]
show_around [selection [, distance [, output_name]]]
save_selections [directory]
save_objects [directory]
```

`label_around [selection]` labels CA atoms with one-letter amino-acid codes and
residue identifiers. The default source is `sele`.

`show_around [selection [, distance [, output_name]]]` selects complete protein
residues within the distance cutoff. Defaults are `sele`, 6 angstroms, and
`nearby_residues`.

`save_selections [directory]` opens a checklist of named selections. Select the
rows in this dialog, then choose **Save CSVs**. All named selections appear in
the checklist, whether enabled or disabled; this is separate from row
highlighting in PyMOL's object panel. The optional directory sets the initial
output folder; otherwise PyMOL's current working directory is used.

`save_objects [directory]` exports enabled molecular objects without opening a
checklist. Groups and disabled objects are skipped. Its optional directory also
defaults to PyMOL's current working directory.

Each exporter writes one `<name>_residues.csv` file per chosen selection or
object. Files contain one row per CA atom with the header `Chain,Res_Num,Res_Name`.
Entries with no CA atoms produce a header-only CSV. Unsafe filename characters
are replaced with underscores; duplicate output names get a numeric suffix.

The checklist is also available from **Plugin > Protein Tools > Choose
Selections to Save...**. The object exporter is at **Plugin > Protein Tools >
Save Enabled Objects (CSV)**.

Run the test suite from the repository root with:

```shell
python3 -m unittest discover -s Protein_Tools -t . -p 'test_plugin.py'
```