# Label Around PyMOL Plugin

Label the alpha-carbon (`CA`) atoms in a PyMOL selection with a one-letter amino
acid code followed by the residue identifier. Labels stay at their CA atom
positions; the command clears any per-atom label offsets left by earlier runs.

## Use

Install the parent `Protein_Tools` package through PyMOL's Plugin Manager. The
parent package loads this plugin and adds it to the shared menu. Do not install
this subfolder separately.

For a one-session manual load, run the source file:

```pymol
run /path/to/label_around/__init__.py
```

Then label the current selection (`sele`):

```pymol
label_around
```

Or provide a selection:

```pymol
label_around chain A
```

For example, alanine residue 42 is labeled `A42`. The residue identifier uses
PyMOL's `resi` value, preserving insertion codes. Standard amino acids use their
one-letter codes; selenocysteine (`SEC`) maps to `U`, pyrrolysine (`PYL`) to `O`,
and unrecognized residue names to `X`.

The command is also available as **Plugin > Protein Tools > Label CA Residues**.
Command help is available with `help label_around`.

## Tests

Run `python -m unittest` from this folder. Tests mock PyMOL's command API and
check CA-only labeling and amino-acid labels.