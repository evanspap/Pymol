# UniProt Structures PyMOL Plugin

Load UniProt-linked experimental structures and the matching AlphaFold model into PyMOL.
The plugin accepts either a UniProt accession, for example `P06730`, or a UniProt entry
name, for example `IF4E_HUMAN`.

For each UniProt entry it will:

- resolve the UniProt accession and entry name;
- read PDB cross-references from UniProt;
- download the AlphaFold model;
- download experimental coordinate files from RCSB, trying `.pdb` first and `.cif` if needed;
- keep only the chains mapped to that UniProt entry;
- organize objects by experimental method;
- align each PDB object to the AlphaFold object using chain A when present.

## Quick Use

```pymol
run H:/My Drive/VSCode_Github/Pymol/UNIPROT_plugin/__init__.py
uniprot_structures P06730
```

Equivalent entry-name example:

```pymol
uniprot_structures IF4E_HUMAN
```

For a safer memory-friendly test:

```pymol
delete all
uniprot_structures P06730, 1, 1, 1, align, 10
```

## PyMOL Layout

For `P06730` / `IF4E_HUMAN`, the root group is:

```text
IF4E_HUMAN
```

Inside it:

- `AF_IF4E_HUMAN` stays directly in the root group.
- `NMR` contains NMR structures.
- `X-ray` contains X-ray crystallography structures.
- `EM` contains electron microscopy structures.
- `Other_methods` is created only if a PDB entry uses another method.

Each experimental structure is kept as one PyMOL object named with only its PDB code,
for example `7XTP` or `3U7X`. If a PDB contains multiple UniProt-mapped chains, they
stay together inside that single PDB object.

PyMOL object and group names are global. If you load multiple UniProt entries in the
same PyMOL session, PyMOL may add suffixes such as `NMR_2` to avoid name collisions,
but the groups are still placed under the relevant UniProt root group.

## Alignment Behavior

The AlphaFold object is the reference. For each PDB object, the plugin aligns using
chain A if chain A exists:

```pymol
align 7XTP & c. A, AF_IF4E_HUMAN
```

In code this is equivalent to using the mobile selection:

```text
(7XTP) and chain A
```

The transform is applied to the PDB object, so the object remains a single object. If
chain A is not present, the plugin falls back to aligning the whole PDB object.

Available alignment methods:

- `align`, default;
- `super`;
- `cealign`.

## Install

In PyMOL:

1. Open `Plugin > Plugin Manager`.
2. Choose `Install New Plugin`.
3. Select this file:

```text
H:\My Drive\VSCode_Github\Pymol\UNIPROT_plugin\__init__.py
```

You can also load it for the current PyMOL session:

```pymol
run H:/My Drive/VSCode_Github/Pymol/UNIPROT_plugin/__init__.py
```

## GUI

After installation, open:

```text
Plugin > UniProt Structures
```

If the menu does not open the window, reload the source file and run:

```pymol
uniprot_structures_gui
```

## Command

```pymol
uniprot_structures query [, include_pdb [, include_alphafold [, align [, alignment_method [, limit]]]]]
```

Examples:

```pymol
uniprot_structures P06730
uniprot_structures IF4E_HUMAN, 1, 1, 1, super
uniprot_structures P06730, 1, 1, 1, align, 10
```

Arguments:

- `query`: UniProt accession or entry name.
- `include_pdb`: `1` or `0`, default `1`.
- `include_alphafold`: `1` or `0`, default `1`.
- `align`: `1` or `0`, default `1`.
- `alignment_method`: `align`, `super`, or `cealign`, default `align`.
- `limit`: maximum number of PDB cross-references to load; `0` means all.

## Downloads And Cache

The plugin downloads real coordinate files to your filesystem and then loads them into
PyMOL. By default, files are cached in your system temp directory under:

```text
pymol_uniprot_structures\<UniProt entry name>
```

The exact folder is printed in the PyMOL console after loading, for example:

```text
[UniProt plugin] downloaded files are in: C:\Users\...\AppData\Local\Temp\pymol_uniprot_structures\IF4E_HUMAN
```

Failed coordinate downloads are skipped, so one missing PDB/mmCIF file does not stop
loading or alignment of the remaining structures.

## Memory Notes

Some UniProt entries have many PDB cross-references, and EM/mmCIF structures can be
large. Loading all structures for several UniProt entries in one PyMOL session can hit
memory limits. Use the `limit` argument for testing or large proteins:

```pymol
uniprot_structures P06730, 1, 1, 1, align, 5
```

The command currently expects one UniProt query at a time. A comma-separated or
space-separated list is not treated as a batch list.

## Data Sources

- UniProt REST API resolves the accession / entry name and reads PDB cross-references.
- AlphaFold DB API provides the AlphaFold model download URL.
- RCSB downloads provide the experimental coordinate files.