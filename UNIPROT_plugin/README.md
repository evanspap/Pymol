# UniProt Structures PyMOL Plugin

Download UniProt-linked PDB/mmCIF structures and AlphaFold models, organize them by
experimental method, and align one mapped protein chain per structure.

## Activate the updated source

In the PyMOL command line:

```pymol
run H:/My Drive/VSCode_Github/Pymol/UNIPROT_plugin/__init__.py
uniprot_structures_gui
```

For permanent installation, reinstall this source through Plugin > Plugin Manager >
Install New Plugin, then restart PyMOL. An already installed copy will not update
just because the source on H: changed.

## GUI

Open Plugin > UniProt Structures, or use `uniprot_structures_gui`.

- Enter one or more UniProt accessions (`P06730`) or entry names (`IF4E_HUMAN`).
  Pasted names separated by spaces, newlines, tabs, commas or semicolons are
  processed as separate proteins. Exact duplicates are processed once. For example:
  `BCL2_HUMAN CAB39_HUMAN PYR1_HUMAN CARM1_HUMAN PEBB_HUMAN`.
- Alternatively, **Browse...** beside **List file** selects a TXT, CSV or TSV file.
  A selected file takes priority over the single-entry field; clear it to use that field.
- Choose a **Download folder** for persistent files, or leave blank for the system
  temporary cache. Each protein has its own subfolder.
- Choose PDB categories independently: **NMR**, **X-ray**, **EM**, and **Other**.
  Uncheck **Download AlphaFold model** for experimental-only downloads. For example,
  NMR-only means keeping only NMR checked among the PDB categories and unchecking AlphaFold.
  For **AlphaFold only**, uncheck **Download PDB structures** and keep AlphaFold checked.
  Category filtering happens before the PDB limit and before any coordinate downloads.
- Choose an alignment method (`align`, `super`, `cealign`). Optionally check
  **Fit only to AlphaFold regions with pLDDT (B-factor) >**; its default threshold is **50**.
- Set **PDB limit per protein**; zero means all available PDB references.
- For many proteins/structures, use **Download only** to save files without adding
  objects to PyMOL or consuming memory for the full collection.
- Click **Start**. The window stays responsive, displays elapsed time, a progress
  bar, and an activity log. The PyMOL console receives the same activity messages.
- **Cancel** stops cooperatively. Finished files and objects are retained. A running
  network read may take up to its timeout (60 seconds); a PyMOL load or alignment
  must finish before cancellation takes effect. Closing a busy dialog requests cancellation.

Progress is reported separately for each protein: structures processed / total,
percentage, remaining, freshly downloaded, reused from cache, download failures,
and loaded objects. Transfer messages also show MB received and, when the server
supplies a size, that file's transfer percentage. **Processed** includes failed jobs;
100% processed does not mean every structure downloaded or aligned successfully.
The final batch message includes the error count. Detailed failures remain in the log.

## Show one loaded protein

The Qt dialog has a **View loaded proteins** selector. Choose a protein and click
**Show only this protein**. This enables its molecular objects and their parent
folders, even if subfolders were individually disabled, and disables other scene
objects. **Show** adds the selected protein to the current view; **Hide** switches
off its parent group. Click **Refresh** after loading a saved session externally.
The selector also recognizes repeated groups such as `BCL2_HUMAN_2`.

These commands work with already loaded sessions, including through the Tk fallback:

```pymol
uniprot_show BCL2_HUMAN
uniprot_show CAB39_HUMAN
uniprot_show BCL2_HUMAN, 0
uniprot_hide BCL2_HUMAN
```

The first two commands isolate the chosen protein; `, 0` shows it alongside the
current view. They change object switches and the camera, preserving atom
coordinates and existing representations (including hidden waters). They do not
restore representations previously removed with H/hide. Empty method folders may
remain disabled because they contain no structures. No files are downloaded.

`disable ! BCL2_HUMAN` excludes only the parent's name from disabling; it also
switches off differently named children such as `X-ray` and PDB-code objects.
Use the controls above to restore and isolate the complete protein hierarchy.

## Input files

A plain TXT list can contain one entry per line or whitespace-separated entries:

```text
# UniProt accessions or entry names
P06730
IF4E_HUMAN
P04637
```

CSV/TSV can contain a recognized UniProt column, for example:

```csv
accession,note
P06730,eIF4E
P04637,p53
```

Accepted column headings include `uniprot`, `uniprot_id`, `uniprot_accession`,
`accession`, `entry`, and `entry_name` (spaces or hyphens also accepted).
Without a recognized heading, the first CSV/TSV column is used. UTF-8 files,
including UTF-8 BOM, are supported. Blank lines, lines beginning with `#`, and
exact duplicate entries are ignored. Protein names as free text and gene symbols
are not resolved: use UniProt accessions or UniProt entry names.

## Visualization and alignment

Each PDB object keeps **one UniProt-mapped protein chain**: chain A when A is mapped
and present, otherwise the first available mapped chain. If UniProt provides no
chain mapping, the first protein chain is used, preferring A. If a supplied mapping
matches no protein chain, the structure is skipped and the reason is logged.
Other protein chains are removed from the PyMOL object. **All non-protein atoms
are retained**, including ligands, cofactors, ions, waters, and nucleic acids,
regardless of their chain identifier or distance from the selected protein chain.
Cached coordinate files remain complete. All model states, including NMR ensembles, are
retained. The selected chain is reported in the activity log.

Structures are aligned as they load, using the protein chain against a single
protein chain in the AlphaFold reference. If AlphaFold is absent, disabled, or fails
to load, the first successfully loaded PDB becomes the reference (unless the
AlphaFold confidence filter is enabled; see below). Alignment failures
are logged and the affected objects are hidden and disabled, so they do not appear
as if they were successfully fitted. Successful protein chains are shown as cartoons, non-water ligands and other
non-protein components as sticks, and isolated atoms (including waters and ions)
as nonbonded markers. Fitting uses only protein atoms and transforms the entire
object together, preserving ligand positions relative to the protein.
With alignment switched off, selected chains are displayed without fitting.

The root group is the UniProt entry name, with `AF_<entry_name>` directly inside it.
Experimental PDB objects use their PDB codes. Method groups are created only when
at least one structure loads successfully, using the entry-name prefix, method,
and number of loaded PDB objects. For example:

```text
TFE2_HUMAN
  AF_TFE2_HUMAN
  TFE2_NMR_1
  TFE2_X-ray_4
```

No EM folder is created in this example. Counts reflect loaded PDB objects after
category filtering, the PDB limit, and download/load failures; NMR model states are
not counted as separate structures. Objects whose alignment failed are still loaded
and included in the count, although hidden. Counts update during loading, including
when a job is cancelled. Repeated runs add a collision suffix when needed to preserve
existing groups. These are PyMOL groups; coordinate files stay in the existing
per-protein cache folder, avoiding duplicate downloads.
Only the objects created by the current job have their representations changed.

## AlphaFold confidence-filtered fitting

The optional filter restricts only the **reference atoms used for fitting** to
protein atoms in the AlphaFold chain whose B-factor field is **strictly greater
than 50**, or your chosen threshold. The field stores pLDDT in AlphaFold files:
[AlphaFold FAQ](https://www.alphafold.ebi.ac.uk/faq). Experimental PDB B-factors
are not filtered or interpreted as confidence. No atoms are deleted by this option;
the full protein and retained ligands are transformed together and remain available.

The option is off by default and requires an AlphaFold reference. If AlphaFold is
missing, restricted fitting is not silently replaced with an experimental B-factor
filter: the affected PDB objects are retained hidden and an error is logged. The
same happens when too few reference atoms pass the cutoff or fitting fails. Without
this option, the existing first-PDB reference fallback still applies.

## Commands

Existing positional syntax is preserved, with new trailing options:

```pymol
uniprot_structures query, include_pdb, include_alphafold, align, alignment_method, limit, download_only, cache_dir, pdb_methods, af_confidence_only, af_plddt_min
```

Defaults: both sources enabled, alignment enabled, method `align`, limit `0`,
`download_only=0`, and the system temporary cache.

```pymol
uniprot_structures P06730
uniprot_structures P06730, 1, 1, 1, align, 10
uniprot_structures P06730, 1, 1, 1, align, 0, 1, D:/UniProt_cache/P06730
```

The last example downloads all selected files without loading them into PyMOL.

New trailing arguments: `pdb_methods=all`, `af_confidence_only=0`, and
`af_plddt_min=50`. Use `NMR`, `X-ray`, `EM`, `Other`, or combinations such as
`NMR+EM`; `none` excludes all experimental categories. Examples:

```pymol
uniprot_structures TFE2_HUMAN, pdb_methods=NMR, include_alphafold=0
uniprot_structures TFE2_HUMAN, pdb_methods=X-ray, include_alphafold=0
uniprot_structures TFE2_HUMAN, pdb_methods=EM, include_alphafold=0
uniprot_structures TFE2_HUMAN, include_pdb=0
uniprot_structures TFE2_HUMAN, pdb_methods=NMR+EM, af_confidence_only=1, af_plddt_min=50
```

Batch file command (limit applies to each protein):

```pymol
uniprot_structures_file H:/My Drive/uniprot_list.txt, 10
uniprot_structures_file H:/My Drive/uniprot_list.txt, 0, 1, D:/UniProt_cache
```

`uniprot_structures_file` also accepts trailing keyword options `pdb_methods`,
`include_alphafold`, `af_confidence_only`, and `af_plddt_min`, with the same defaults.

GUI jobs run in a background worker. Command-line calls run synchronously and print
progress to the console. Only one job can run at a time. Batch entries run sequentially;
an entry failure is logged and does not prevent later entries from being attempted.

## Downloads and memory

The default cache is `<system temp>/pymol_uniprot_structures/<entry_name>`.
The exact path is reported. RCSB downloads try PDB first, then mmCIF on a 404.
Downloads use fixed-size buffers and atomic file replacement, so interrupted or
failed transfers do not leave a partial file that is reused as a completed download.
Existing nonempty files are reused; older cache files are not structurally validated.

Pruning chains in place avoids a second full PyMOL object, but PyMOL must initially
load each complete coordinate file. Large assemblies and NMR ensembles can still
require substantial memory. Loading many structures retains them in the scene;
the worker and progress bar do not impose a memory ceiling. Use a modest PDB limit
for visualization and **Download only** for bulk acquisition.

## Verification

`test_plugin.py` covers list parsing, progress and cache accounting, atomic-download
failure/cancellation handling, batch failure recovery, and actual headless PyMOL
chain selection, fitting, and failed-alignment visibility. It also checks method
filters, nonempty counted groups, and confidence-filtered reference selections. Run with a Python
interpreter that provides PyMOL (for example PyMOL's bundled Python):

```powershell
& 'C:\Users\geras\AppData\Local\Schrodinger\PyMOL2\python.exe' -B 'H:\My Drive\VSCode_Github\Pymol\UNIPROT_plugin\test_plugin.py'
```

## Data and API sources

- [UniProt REST API](https://www.uniprot.org/help/api_queries)
- [AlphaFold DB](https://alphafold.ebi.ac.uk/)
- [RCSB coordinate downloads](https://files.rcsb.org/download/)
- [PyMOL command implementation](https://github.com/schrodinger/pymol-open-source/blob/master/modules/pymol/importing.py)
