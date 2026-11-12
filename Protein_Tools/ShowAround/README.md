# ShowAround PyMOL Plugin

Create a selection containing complete protein residues with at least one atom
within a chosen distance of a source selection. The defaults use PyMOL's main
selection, `sele`, a 6 angstrom radius, and the output selection name
`nearby_residues`.

## Use

Install the parent `Protein_Tools` package through PyMOL's Plugin Manager. The
parent package loads this plugin and adds it to the shared menu. Do not install
this subfolder separately.

For a one-session manual load, run the source file:

```pymol
run /path/to/ShowAround/__init__.py
```

Then run the command with up to three comma-separated values:

```pymol
show_around
show_around ligand, 8, ligand_neighborhood
```

The syntax is `show_around selection, distance, output_name`:

- `selection` is the existing selection to measure from. It defaults to `sele`,
	PyMOL's current selection. For example, `ligand` uses the atoms in the selection
	named `ligand` as the reference.
- `distance` is the maximum atom-to-atom distance in angstroms. It defaults to `6`.
- `output_name` names the new selection that PyMOL creates. It defaults to
	`nearby_residues`. In `show_around ligand, 8, nearby`, `nearby` is simply the name you
	chose for the result; it has no special meaning and could instead be
	`ligand_neighborhood` or another valid selection name.

The menu item is **Plugin > Protein Tools > Show Around (6 A)**. The command's
detailed help is available with `help show_around`.

The generated selection filters to `polymer.protein` before expanding matching
atoms with `byres`, so nearby residues are selected in full and non-protein atoms
are excluded.

## Tests

Run `python -m unittest` from this folder. The tests mock PyMOL's command API and
cover default/custom selections, radius validation, and output-name validation.