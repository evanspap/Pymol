# PyMOL Contact Map Plugin

**Version v5**
**Purpose:** inspect residue contacts and residue-pair distances for structures and trajectories loaded in PyMOL.

The plugin opens a separate window inside your PyMOL session. It reads the selected molecular object, updates as you change frames, and can save the displayed plot as an image. It does not alter the molecule's coordinates or write changes to your structure or trajectory files.

![Illustration of the current-frame contact map and split mean-distance/standard-deviation view](contact_map_modes.svg)

*Illustration only, not a screenshot. In the split view, mean distances occupy the upper triangle and distance standard deviations occupy the lower triangle.*

## What You Can Plot

Choose a view from the **Plot** menu in the plugin window:

| Plot option | What is displayed | Atoms and frames |
| --- | --- | --- |
| **Current frame cutoff** | Binary contact map: a pair is colored as a contact when its CA-to-CA distance is at or below the cutoff. | CA atoms in the current PyMOL frame. |
| **Current frame residue-center distances** | Pairwise Euclidean distances between residue geometric centers. Each center averages all loaded protein atoms belonging to that residue, including any hydrogens present. | Current PyMOL frame. |
| **Mean + standard deviation** | One statistics view with **Show mean** and **Show standard deviation** checkboxes. Display either statistic alone or both side by side. | All object states, using the selected distance definition. |
| **Split mean / SD across diagonal** | Mean distance in the upper triangle; standard deviation in the lower triangle. The diagonal shows the zero mean self-distance. Each triangle has its own colorbar when **Auto max** is selected. | All states/frames in the selected object. |

The separate Mean and SD dropdown entries were removed in version 1.7.0. Select **Mean + standard deviation**, then choose the statistic(s) to display within that view. In split mode both triangles remain present; their profiles and color ranges can still be controlled independently.

For both statistics layouts, the **Distance definition** dropdown selects what is measured in each frame:

| Definition | Per-frame residue-pair measurement |
| --- | --- |
| **C-alpha to C-alpha** | Distance between the two residues' CA positions; the default. |
| **Residue center to center (all atoms)** | Distance between geometric, equally weighted centers of all loaded protein atoms in each residue, including loaded hydrogens. Not a mass-weighted center. |

Both mean and population SD use this same selected measurement across all states. Changing the definition recalculates both statistics. Current-frame cutoff contacts continue to use CA distances, independently of this selector.

The closest-residue-distance option is temporarily unavailable in the dropdown. Residue centers now use all loaded atoms rather than a heavy-atom-only selection. The plugin does not add missing atoms or hydrogens; when no hydrogens are present in the input, the center uses the atoms that are available. Current-frame center plots, temporal mean/SD, and frame exports use the same all-atom definition.

Every heatmap shows its column residue labels along the **top** edge and row residue labels along the **left** edge in the same order. Matching residues therefore start together near the upper-left corner. This applies to cutoff, center-distance, combined, split, and exported frame heatmaps. Line plots keep their existing residue-axis orientations.

The dedicated **Line plots** section contains **Show mean line plot** and **Show SD line plot**. Both are optional and initially unchecked. In the combined view each line plot appears below its own heatmap, with residues on the x-axis. The mean line's y-axis is **Mean of means (A)**: the average of eligible pairwise mean distances for that residue. The SD line's y-axis is **Mean of SDs (A)**: the average of eligible pairwise distance SDs, not an SD calculated from the mean line.

In split mode the SD line is below the matrix, with residues on the x-axis; the mean line is to its right, with residue positions aligned vertically to the heatmap and mean distance on its horizontal axis. Profiles omit self-pairs and neighbor-filtered pairs. Both profiles honor the Mean pair-selection range. **SD's own cutoffs change its heatmap, not its line profile or profile CSV.** A residue with no eligible partners has a gap rather than an artificial zero. Lines have no circle markers and use larger plotting areas; split profiles stay clear of the colorbar column. Color scaling and gray exclusions never modify the raw matrices.

Every plot mode now uses the single **Export…** button in the upper-right corner. Its popup offers the current view with any selected CSV datasets, or a range of loaded frames. There are no separate visible Save plot, Export CSV, or Export frames buttons. Matrix CSVs retain raw values; profile CSVs contain the range-filtered column averages used by the lines. See [Export Any Plot as CSV](#export-any-plot-as-csv) for the data formats.

## What's New in v5

- Temporarily removed the closest-residue-distance dropdown option.
- Residue centers average all loaded protein atoms, including any hydrogens present, in current-frame calculations, temporal statistics, and frame exports.
- Column residue labels appear at the top of every heatmap, matching the left-axis order. Line-plot axes are unchanged.
- Mean and SD colorbars each have draggable upper/lower cutoff handles, precise numeric fields, and a third value at the midpoint of the full scale. Moving the upper handle turns off Auto max; lower handles cannot cross upper handles. The full scale stays anchored to the data maximum, or a larger manually entered upper value.
- SD line profiles respond only to Mean cutoffs, not SD cutoffs. Lines are marker-free, wider, and taller, with scrolling when needed.
- Image exports use an off-screen snapshot instead of rendering live controls. Typed destination folders are created as needed; incomplete current-view bundles are removed on failure. PNG/JPEG and frame-series exports remain available.

## What's New in v4

- Mean and SD lower/upper numeric limits sit directly beside the corresponding viridis bar: upper bound at the top, lower bound at the bottom. Log/Linear and Auto max are located there too.
- Mean limits act as a shared residue-pair filter. The corresponding cells in SD turn dark gray below Mean's lower bound, or light gray above Mean's upper bound. These pairs no longer contribute to the SD line profile.
- For pairs inside the Mean range, SD's own limits and color scale still apply independently. Mean exclusions take precedence when both filters exclude a pair in different directions.
- All changes update the heatmaps and line profiles immediately using cached raw statistics; no manual Refresh click is needed.
- Cutoff contacts have a genuinely binary, two-color legend. Adjust the cutoff with the slider or numeric field at that legend; it no longer occupies a top-row control.
- Export actions share one popup, top rows use compact spacing, and smaller plot areas can scroll instead of cropping line plots. Exported images include the complete canvas, not just its visible scroll viewport, and omit interactive control widgets.

## Install and Open

### Requirements

- PyMOL with its bundled Qt interface and NumPy. Matplotlib is not needed.
- A protein object loaded in the PyMOL session. The plugin does not load files itself.
- For trajectory analysis, one PyMOL object with multiple states. A multi-model PDB naturally loads as multiple states; a separate trajectory must be attached to a compatible topology.

### One-time use in a PyMOL session

1. Start PyMOL and open the structure or trajectory you want to analyze.
2. Locate PyMOL's command input area. These are **PyMOL commands**, not PowerShell commands.
3. Load the plugin by executing the following block. `runpy.run_path` is used because the Windows folder path contains spaces:

   ```pymol
   python
   import runpy
   runpy.run_path(r"C:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol-Public\contact map\__init__.py")
   python end
   ```

4. Open the window:

   ```pymol
   contact_map_gui
   ```

The window is called **Protein Contact Map**. Running `contact_map_gui` again raises the existing window if it is already open. The plugin also defines a **Contact Map** menu item when installed through a PyMOL plugin manager that supports this plugin format.

You need to load the plugin once per PyMOL session unless you install it permanently through your PyMOL build's Plugin Manager. Reopening PyMOL clears commands loaded with the one-time method.

### Load a protein and trajectory

For a separate topology and trajectory, first load the topology and then attach its matching trajectory to the same object:

```pymol
load "D:/data/protein.gro", protein
load_traj "D:/data/trajectory.xtc", protein
```

Use the actual paths to your files. Forward slashes are convenient in PyMOL commands on Windows; quote paths that contain spaces. Your PyMOL build must support the formats, and the topology and trajectory must have matching atom counts and ordering.

For an ensemble PDB that contains `MODEL`/`ENDMDL` records, load it as one object; each model becomes a state:

```pymol
load "D:/data/ensemble.pdb", ensemble
```

You can open the contact-map window before or after loading data. If you load an object after opening it, click **Refresh objects** to update the object list.

## First-use Walkthrough

Follow these steps once the data is loaded in PyMOL and the **Protein Contact Map** window is open.

1. **Select the object.** Open the **Object** dropdown and choose the protein or trajectory object. CA plots require CA atoms; center statistics use all loaded protein atoms. If the object is missing from the list, click **Refresh objects**.
2. **Check the frame count.** The label below the controls shows the current frame and total states for a current-frame plot. A single-state structure can show one frame; aggregate mean and SD then have only that one observation.
3. **Start with a current-frame view.** Choose **Current frame cutoff** to display a binary contact map for the frame PyMOL is currently showing. The status text reports the number of CA residues and contacts.
4. **Adjust the cutoff at its legend.** The cutoff starts at `8.0 A`. Use the vertical slider or numeric field beside/below the binary colorbar. Lowering it requires residues to be closer to count as a contact; raising it includes more distant pairs. The legend has only **1: contact** and **0: no contact**, not a continuous distance gradient.
5. **Choose the neighbor exclusion.** **Exclude same-chain neighbors** now starts at `0`, so no off-diagonal sequence neighbors are excluded by default. Increase it only if you want to omit nearby same-chain entries in the selected residue order. Self-distances stay zero in distance matrices and are excluded from profile averages.
6. **Inspect residue-center distances if needed.** Choose **Current frame residue-center distances**. For each residue, the plugin averages the coordinates of all loaded protein atoms, including any hydrogens, then measures center-to-center distances in the current frame. The cutoff is disabled because this mode does not threshold distances.
7. **Move through the trajectory.** Drag PyMOL's frame slider or press Play. Current-frame plots refresh from the visible PyMOL frame about every 350 ms. Pause playback before saving if you want a specific frame image.
8. **Open statistics.** Choose **Mean + standard deviation**. Initially both **Show mean** and **Show standard deviation** are checked. Uncheck either to enlarge the other plot; if both are unchecked, the window asks you to select a statistic.
9. **Choose the distance definition.** Select CA-to-CA or **Residue center to center (all atoms)**. Mean and SD both use that definition. The status line reports frame progress during the first calculation.
10. **Choose ranges on the colorbars.** Drag either horizontal cutoff handle directly along the Mean or SD colorbar. The Upper and Lower fields show precise cutoffs; a third value marks the full scale's midpoint (logarithmic midpoint in Log mode). Dragging the upper handle switches off **Auto max**. You can also uncheck Auto max and type the upper bound; the lower field stays editable. For example, set Mean to `5-20 A` and SD to `0.5-3 A`. Boundaries are inclusive. Mean limits filter pairs in both maps; SD limits additionally filter only SD heatmap colors. Handles cannot cross. Typed limits must satisfy lower <= upper; invalid input retains the last valid view and reports a status error.
11. **Show line plots.** Check **Show mean line plot** and/or **Show SD line plot**. Both start off and can be toggled independently. In the combined view, keep the corresponding heatmap visible to see its line beneath it. The Mean line averages in-range pairwise means; the SD line averages pairwise SDs for pairs passing **Mean's range only**. Changing SD cutoffs does not change the SD line. Exclusions do not contribute zeros. Mean cutoff changes update both lines immediately; Log/Linear changes only heatmap colors. Lines have no point markers and use expanded plotting areas. If a small window cannot show the full canvas, use the plot-area scrollbars.
12. **Compare both statistics in one matrix.** Choose **Split mean / SD across diagonal**. Read the upper triangle as mean distance and the lower triangle as distance SD.
13. **Choose scales independently.** Choose Log or Linear beside each viridis bar. Mean may use Log while SD uses Linear. The same colorbar controls work in the combined and split layouts. In split mode enable the Mean line for the right-side plot and the SD line for the bottom plot.
14. **Read the gray exclusions.** Below the Mean lower bound a pair is dark gray in both maps; above Mean's upper bound it is light gray in both. For other pairs, SD's own range determines its viridis/gray color. Structural exclusions remain white. Mean-excluded pairs do not contribute to either line; pairs gray only because of SD cutoffs still contribute to the SD line. Raw matrix CSV values are unchanged.
15. **Export images and data.** Click the single upper-right **Export…** button. Keep **Current view and selected CSV datasets** for the displayed aggregate/frame plot. Choose a destination folder, optionally save a PNG/JPEG image, and check any combination of matrix/profile CSVs. Uncheck all CSVs for image-only output, or disable the image for CSV-only output. Click the dialog's **Export** button; it creates a uniquely named subfolder and reports its path.
16. **Export frame numbers.** In that same popup, select **Selected frames**, enter First, Last, and Interval, and choose PNGs and/or matching frame CSVs. A statistics view exports per-frame distances under its selected distance definition, not repeated temporal Mean/SD images.

The plugin does not change coordinates or input files. Statistics use all states equally and are cached by object, state count, distance definition, and neighbor exclusion. Display ranges and visibility use these cached results. Click **Refresh objects** after modifying coordinates in an existing object to clear cached statistics; coordinate edits do not automatically invalidate the cache.

## Save Plots Frame by Frame

Open the upper-right **Export…** button and choose **Selected frames**. This saves PNG plots, their matching raw matrix CSVs, or both for a selected range of loaded states. Choose **Current view and selected CSV datasets** in the same popup when you need the displayed plot or aggregate statistics instead.

### Choose what each frame represents

| Starting plot mode | What Export frames saves for each state |
| --- | --- |
| Current frame cutoff | Binary CA-contact map using the current cutoff and neighbor filter. |
| Current frame residue-center distances | Geometric residue-center distance map using all loaded protein atoms. |
| Mean + standard deviation or Split mean / SD | Individual-frame distance map using the selected CA-to-CA or all-atom residue-center **Distance definition**. |

Frame exports from statistics views are **not repeated copies of the temporal mean or SD plots**. Each map contains distances measured in that one frame. Use the popup's **Current view** scope for temporal Mean/SD and their line profiles. No aggregate line plots are placed on frame-export images.

### Export steps

1. Select the molecular **Object**. Only frames/states already loaded in PyMOL are available.
2. Select the starting plot mode from the table above. For statistics views, choose the desired distance definition first. For example, select **Mean + standard deviation**, then **Residue center to center (all atoms)**, to export residue-center distance maps for each frame.
3. Configure the neighbor exclusion and color scale. Frame export from a current-frame mode uses that mode's global scale settings. Export from a statistics view uses the **Mean** scale, Min, and Auto/fixed Max settings; SD settings do not describe an individual-frame distance map.
4. Click **Export…**, then choose **Selected frames** in the popup's Export dropdown.
5. Click **Browse…** beside Destination and choose a writable folder, or enter a path. Missing destination folders are created automatically when permissions allow.
6. Enter **First**, **Last (inclusive)**, and **Interval**. Numbering is one-based. For example, `1`, `10`, `2` exports `1, 3, 5, 7, 9`. These numbers refer to loaded PyMOL states, not original simulation timestamps.
7. Leave both **Save PNG plot per frame** and **Save matching frame matrix CSVs** checked for paired output, or uncheck either for images-only/CSV-only output. At least one type is required. Frame-series images are PNG; the JPEG choice applies to Current view only.
8. Click **Export**. A progress dialog shows the current frame and completed work; the PyMOL console also prints per-frame progress. Settings are temporarily disabled and playback is paused so each output corresponds to a stable calculation.
9. Wait for completion or press **Cancel**. Cancellation is checked between frames; the current frame's calculation may finish before cancellation takes effect.
10. Read the status line for the output directory. Each invocation creates a new, uniquely named subfolder inside your chosen destination, so an earlier export is not overwritten. Completed frames remain saved if you cancel. The original PyMOL frame, view mode, timer, and playback state are restored afterward.

### Output files and interpretation

Each completed frame has the same base filename for its PNG and CSV, for example:

```text
protein_centers_frame_000001.png
protein_centers_frame_000001.csv
protein_centers_frame_000002.png
protein_centers_frame_000002.csv
```

The CSV contains the raw residue-by-residue matrix for that exact frame, with labels on both axes. Contact CSVs contain `0`/`1`; distance CSVs contain Angstroms with `NA` for structurally excluded neighbors. Color transforms and display limits affect PNGs but do not clip the CSV numbers. Files are saved after each frame, not held until the whole batch finishes. If writing the current frame fails, its incomplete PNG/CSV pair is removed; earlier completed frames remain.

With **Auto max**, each image uses its own frame's maximum. Use a fixed Max if colors must be comparable between images. A positive Min cannot exceed a frame's effective maximum; if that happens, the export reports a range error.

## Export Any Plot as CSV

### Choose the data and destination

1. Select the **Object** and the **Plot** mode. Wait until the calculation completes without an error.
2. Set the current frame for a frame-based plot, or the distance definition and neighbor filter for a statistics plot.
3. If exporting a line profile, first set the **Mean** cutoff range. That range determines which partners contribute to both exported profiles. SD cutoffs affect heatmap colors only.
4. Click upper-right **Export…** and keep **Current view and selected CSV datasets**. The popup shows a checkbox list of datasets available for the current plot mode:

| Plot mode | CSV choices |
| --- | --- |
| Current frame cutoff | **Cutoff contacts matrix (0/1)** |
| Current frame residue-center distances | **Residue-center distance matrix (A)** |
| Mean + standard deviation | **Mean CSV**, **SD CSV**, **Mean residue profile**, **SD residue profile** |
| Split mean / SD across diagonal | All four statistics choices plus **Split matrix (upper mean / lower SD)** |

5. Choose the destination folder with Browse. Check precisely the CSV datasets you want; multiple selections are allowed. Clear **Save plot image** for CSV-only output.
6. Click the dialog's **Export** button. The plugin creates a unique subfolder in your chosen destination and writes all selected files there. The status line reports the output path. Canceling the popup writes nothing.
7. You can save the current PNG/JPEG in the same operation. The old separate on-window CSV selector and export buttons are no longer needed.

Profiles are available for export even when their display checkboxes are off, and Mean/SD CSVs remain available when either heatmap is hidden. In split mode the right-hand Mean profile exports the same residue/value data as the Mean profile in the combined layout; its visual orientation does not change the CSV format. No data are available to export when there is no valid object or the plot has a calculation/range error.

### Heatmap CSV format

Heatmaps export square matrices. The first row contains column residue labels; the first column contains row residue labels. Both axes follow the same residue order as the plot. Numeric values are saved with up to 12 significant digits.

```csv
residue,A:ALA1,A:ALA2,A:ALA3
A:ALA1,0,12,20
A:ALA2,12,0,8
A:ALA3,20,8,0
```

This example is a mean-distance matrix in Angstroms. A cutoff-contact matrix instead contains binary `0`/`1` values. In cutoff mode, excluded neighbors and self-pairs are represented as `0`, matching the binary map. For distance matrices, unavailable neighbor-filtered pairs are `NA`, while self-distances are `0`.

The split matrix exports the raw mean values above the diagonal, raw SD values below it, and zero mean self-distances on the diagonal. Both statistics have units of Angstroms, but the split matrix is not symmetric and does not represent a single symmetric distance measure. Export Mean and SD separately if you need two complete symmetric matrices.

Heatmap exports are **not** transformed by Log/Linear scaling or clipped to display limits. Dark/light gray cells still retain their raw numeric values. Matrix exports do not contain image colors, per-frame coordinates, or an extra line-profile column.

### Line-plot CSV format

Mean and SD profiles export one row per residue and two columns: residue label and the plotted column-average value.

```csv
residue,column_mean_distance_angstrom
A:ALA1,12
A:ALA2,12
A:ALA3,NA
```

The Mean profile's value column is `column_mean_distance_angstrom`; the SD profile uses `column_mean_distance_sd_angstrom`. The SD profile is the average of eligible pairwise distance SDs, not the SD of the mean profile.

Self-pairs and neighbor exclusions do not contribute. Both profiles exclude pairs outside Mean's range. SD's own range does not filter its profile. `NA` means no eligible partners remain, not zero. Mean range changes affect exported profiles; SD range, Log/Linear coloring, and profile visibility do not.

### Filenames and reproducibility

Output filenames include the object and dataset. Current-frame filenames include the loaded frame number; statistics filenames include `ca` or `centers`. Each export uses a new subfolder so previous files are not overwritten. Record the distance definition, loaded states, neighbor exclusion, and ranges alongside CSVs; these are not embedded as metadata rows. Playback and automatic frame refresh are paused while the export popup is open.

## How the Measurements Work

### Current-frame cutoff contacts

For CA coordinates $\mathbf{x}_i$ and $\mathbf{x}_j$, the contact is present when

$$
\lVert \mathbf{x}_i - \mathbf{x}_j \rVert_2 \leq c,
$$

where $c$ is the cutoff in Angstroms. The map is symmetric. Same-chain neighbors within the selected sequence-order exclusion and the diagonal are omitted.

### Current-frame residue-center distances

For each residue, the plugin averages the coordinates of all its loaded protein atoms, including hydrogens when present, to get a geometric center:

$$
\mathbf{r}_i = \frac{1}{N_i}\sum_{a=1}^{N_i}\mathbf{x}_{ia}.
$$

The plotted value for a residue pair is $\lVert\mathbf{r}_i-\mathbf{r}_j\rVert_2$ in the current frame. This is a geometric center, not a mass-weighted center and not the closest distance between atoms.

### Mean and standard-deviation distances

For residue pair $i,j$ over $F$ object states, let $d_{ij}^{(f)}$ be the selected CA or all-atom center-to-center distance in frame $f$. The mean distance is

$$
\bar d_{ij} = \frac{1}{F}\sum_{f=1}^{F} d_{ij}^{(f)}.
$$

The plugin displays the **population** standard deviation:

$$
\sigma_{ij} = \sqrt{\frac{1}{F}\sum_{f=1}^{F}\left(d_{ij}^{(f)}-\bar d_{ij}\right)^2}.
$$

The side summaries report the overall mean of the eligible pair distances and overall population SD, or summaries of per-pair SDs. Diagonal self-pairs and excluded same-chain neighbors are omitted from these summaries.

Side summaries describe raw eligible distances, not only colored pairs. Both line profiles honor the shared Mean selection; SD's own range changes its colors only. Filtering occurs after per-pair temporal statistics; it does not discard individual frames before mean/SD computation.

## Color Scales and Range

In statistics views drag the two cutoff handles on each colorbar or edit its Lower/Upper fields. Each bar displays the upper cutoff, lower cutoff, and midpoint of its full scale. The full scale extends from zero to the raw maximum (or a larger manual upper bound), so handle movement does not continually rescale its track. The two statistics retain independent viridis scales, while Mean's range filters pair eligibility in both SD colors and the SD profile. With Auto max checked, the upper field is read-only and shows the raw maximum; dragging the upper handle disables Auto max. The lower field remains editable.

- **Log** maps in-range values using $[\log(1+x)-\log(1+L)]/[\log(1+U)-\log(1+L)]$, retaining valid zero values when $L=0$.
- **Linear** maps in-range values using $(x-L)/(U-L)$.
- Below Mean's lower bound: dark gray `(80, 80, 80)` in both maps; above Mean's upper bound: light gray `(210, 210, 210)` in both maps. Mean exclusions take precedence over SD's own range colors. For pairs passing Mean's filter, SD applies its own bounds normally. Boundaries are inclusive.
- Mean and SD each retain their own settings when switching between the combined and split layouts. A degenerate all-zero SD range renders valid zeros with the first viridis color.
- The cutoff-contact legend is binary: purple means `0: no contact`, yellow means `1: contact`. The adjacent slider and Cutoff field adjust the distance threshold in Angstroms. There is no continuous viridis scale or Log/Linear control in this binary mode.

White denotes structurally unavailable pairs in distance maps. True self-distance zeros use viridis if zero is in-range, or dark gray if Min is positive. Cutoff contact maps remain binary, with zero meaning no included contact.

## Controls at a Glance

| Control | Default | Applies to |
| --- | --- | --- |
| Object | First available molecular object | All views |
| Cutoff slider / field | `8.0 A`, on binary legend | Current frame cutoff only |
| Exclude same-chain neighbors | `0` sequence-order entries | All views |
| Plot | Current frame cutoff | All views |
| Color scale | Log, at each continuous bar | Distance/statistics views; cutoff mode is binary |
| Auto max | On | Distance plots; cutoff contacts always use 0-1 |
| Max (A) | `40.00` when enabled | Distance plots with Auto max off |
| Show mean / Show standard deviation | Both on | Combined statistics view |
| Distance definition | C-alpha to C-alpha | Both statistics layouts |
| Show mean line plot / Show SD line plot | Both off, under Line plots | Combined: bottom of own map; split: mean right, SD bottom |
| Mean colorbar limits | Log, lower `0`, Auto max | Mean map / upper triangle; shared selection in SD |
| SD colorbar limits | Log, lower `0`, Auto max | SD map / lower triangle, in addition to Mean selection |
| Export… | Upper-right, enabled after valid calculation | One dialog for current PNG/JPEG, selected CSVs, or frame-series PNG/CSV |

## Requirements and Limitations

- PyMOL with its bundled Qt interface and NumPy. Matplotlib is not required.
- Contact and aggregate CA-distance views require protein CA atoms in the selected object.
- Residue-center distance view requires at least one eligible residue pair; centers average the protein atoms actually loaded for each residue. Missing hydrogens are not generated automatically.
- Aggregate views require matching selected residue identities across states. Segment, chain, residue number/insertion code, and residue name distinguish residues.
- The plugin analyzes coordinates already loaded into PyMOL. It does not align, transform, or modify those coordinates.
- Large trajectories can take time to calculate the first aggregate map. Progress is shown in the status line; subsequent requests with the same object and neighbor setting use cached statistics.
- An empty map can result when the object has no protein residues, or when the same-chain neighbor exclusion removes every eligible pair.
- Closest-pair distance is not currently offered in the dropdown. No external analysis package is required for the available metrics.
- Distances are ordinary Cartesian Euclidean distances, without periodic minimum-image correction. Unwrap/repair broken molecules before loading if appropriate for the trajectory.
- Profiles use a linear numeric distance axis and are not changed by the heatmap's Log/Linear transform. A missing point means no eligible partners remain inside Mean's range, not a distance of zero.

## Version History

| Version | Incremental changes |
| --- | --- |
| **1.0.0** | Initial live CA contact map with object selection, distance cutoff, same-chain neighbor exclusion, frame refresh, and a separate PyMOL Qt window. |
| **1.1.0** | Added plot saving, viridis coloring, and aggregate plot controls. The first aggregate implementation used contact occupancy; this was replaced in 1.2.0. |
| **1.2.0** | Corrected aggregate plots to mean CA-pair distance and distance SD; added the optional split-diagonal view and zero-based logarithmic color mapping. |
| **1.3.0** | Added square Log/Linear controls, Auto max, and a user-set maximum range with clipping. |
| **1.4.0** | Separated the current-frame cutoff contact map from a new current-frame residue-center distance map based on geometric centers of non-hydrogen protein atoms. |
| **1.5.0** | Added optional per-residue column-average profiles below mean and SD heatmaps, plus independent scale mode and maximum controls for each split-diagonal half. |
| **1.6.0** | Added optional labeled CSV export for Mean distance and Distance standard deviation matrices. Excluded/filtered matrix entries are represented as `NA`. |
| **1.7.0** | Consolidated mean/SD into one selectable statistics view; added CA, heavy-atom-center and closest-heavy-atom metrics; independent Min/Max and Log/Linear controls in both statistics layouts; gray below/above-range cells; range-filtered profiles; split SD profile below and mean profile right; zero default neighbor exclusion; tightened Object spacing and retained raw CSV exports. |
| **1.8.0** | Added destination-selectable CSV exports in every plot mode: current-frame cutoff and residue-center matrices, raw mean/SD matrices, the combined split matrix, and range-filtered Mean/SD residue profiles. |
| **v2** | Moved the optional line-plot checkboxes to a dedicated Line plots section, labeled the mean-of-means and mean-of-SDs axes explicitly, added a Residue label beneath each line plot, and tightened the Cutoff label/input spacing. Release renamed to v2; functionality unchanged. |
| **v3** | Renamed the third distance definition to closest inter residue atomic distance. Added optional frame-series PNG/CSV export with a destination chooser, frame range/interval, progress and cancellation, unique output folders, and restoration of the original session state. |
| **v4** | Consolidated export into one upper-right dialog with multiple CSV selection and frame ranges; moved live range/scale controls onto the colorbars; linked Mean pair-selection bounds to SD heatmap/profile exclusions; added a binary cutoff legend with an on-bar slider/field; compacted the header and made the full plot area scrollable. |
| **v5** | Temporarily removed the closest-distance option; changed residue centers to include all loaded atoms, including hydrogens; moved column residue labels to the top of every heatmap. Added draggable Mean/SD cutoff handles, Mean-only SD profile filtering, expanded marker-free line plots, and isolated off-screen export rendering with failed-bundle cleanup. |

## Troubleshooting

**I cannot see the line plots.** Select **Mean + standard deviation**, keep the desired heatmap(s) checked under Show mean / Show standard deviation, and enable the corresponding **Show mean line plot** / **Show SD line plot** checkboxes under **Line plots**. These options start unchecked and are separate from the CSV selector. If the new section is absent, restart PyMOL and load the updated plugin file; an already-open window from a previously loaded version does not acquire new controls automatically.

**The object menu is empty.** Load a molecular object in PyMOL and press **Refresh objects**.

**The plot says no protein CA atoms were found.** This view requires protein residues with CA atoms. Check the selected object and loaded state.

**The residue-center view has too few pairs.** Reduce **Exclude same-chain neighbors** or select an object with more protein residues. The view needs at least one eligible residue pair.

**A map has gray cells.** Dark gray is below the lower cutoff; light gray is above the upper cutoff. Values remain in raw matrix CSVs. Mean exclusions remove pairs from both profiles; SD-only exclusions do not remove pairs from the SD profile.

**A profile has gaps or is empty.** No off-diagonal partners survive the Mean range and structural/neighbor exclusions. Check Mean cutoffs for either line. Changing SD cutoffs cannot restore excluded pairs or alter its line. Gaps are not zero.

**The range is invalid.** Ensure Min does not exceed Max, including the effective data maximum when Auto max is checked. A static trajectory's SD may be zero everywhere; Min must be zero to use its zero-valued auto range.

**Both statistics disappeared.** Check at least one of Show mean / Show standard deviation in the combined view.

**CSV export is disabled.** Load/select a valid molecular object and fix any calculation or range errors. The selector lists only datasets for the current plot mode. Profile exports do not require enabling their display checkboxes.

**The CSV contains values above the displayed maximum.** This is intentional for heatmap matrices: exports preserve raw numbers. Select a residue-profile export to obtain the range-filtered averages used in a line plot.

**Frame-export images do not contain mean/SD statistics.** This is intentional: they contain individual-state measurements. In the same Export dialog, switch to Current view to save aggregate statistics.

**Line plots do not all fit on the screen.** The smaller header leaves more plot space. If necessary use the plot area's scrollbars; exports still save the entire canvas. No line-plot data are dropped because a window is small.

**I cannot find the frame exports.** Look at the final status message for the new subfolder created inside your chosen destination. Completed frames are written one by one. An empty subfolder can remain if cancellation or an error occurs before the first frame is saved.

## Run the Tests

From the repository root in PowerShell, use a Python interpreter bundled with PyMOL:

```powershell
& "C:\Users\yagd9\AppData\Local\Schrodinger\PyMOL2_new\python.exe" -B -m unittest discover -s "contact map" -p "test_contact_map.py" -v
```

The tests use temporary synthetic structures rather than research trajectories. They check all distance definitions, temporal means/SDs, gray range colors, filtered profiles, independent scales, visibility controls, image export, and every CSV export option, including destination selection, cancellation, frame updates, and raw-versus-profile filtering.

**The plot does not change while playing.** Confirm the object has multiple PyMOL states and that the global frame slider is moving. Current-frame views refresh from PyMOL's current frame.

**The plugin command is not found.** Load the plugin file first in the PyMOL session, then run `contact_map_gui`.