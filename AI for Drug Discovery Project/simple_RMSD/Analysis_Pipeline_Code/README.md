# Analysis Pipeline (Very Beginner Guide)

This folder is for a simple RMSD class project.

## What RMSD means here
RMSD means how different two frames are from each other.

Mathematical formula used in these scripts:
- Let N be the number of aligned atoms.
- Let r_i and s_i be the 3D coordinates of atom i in the two frames being compared.
- RMSD is:
	RMSD = sqrt((1/N) * sum_{i=1..N} ||r_i - s_i||^2)
- Practical meaning:
	First compute coordinate differences atom-by-atom, square them, average over all atoms, then take the square root.

In this project:
- A reference frame is the frame you compare against.
- A target frame is the frame you compare to the reference.
- If you compare one reference frame to one target frame, that is 1v1.
- If you compare one reference frame to many target frames, that is 1vN.
- If you compare many reference frames to one target frame, that is Nvs1.
- If you compare every selected frame to every other selected frame, that is all-vs-all.

If you only want a single RMSD number, use the barebones helper file.
If you want a table, use the main RMSD script.

## What to do first
1. Open PowerShell.
2. Move into this folder.

```powershell
Set-Location "C:\Users\yagd9\Documents\Stony Brook AI Drug Discovery\Github\Pymol\AI for Drug Discovery Project\simple_RMSD\Analysis_Pipeline_Code"
```

3. Replace these two placeholders in the commands below:
- `<path-to-your-pdb-file>`
- `<path-to-output-folder>`

That is all you need to change in most cases.

---

## 1. RMSD tools

There are two different RMSD files in this folder.

### A. Barebones helper file
Use this file if you only want one RMSD number printed in the terminal:
- RMSD_calculation_function_7-21-26.py

This file is the simplest one.
It does not make tables.
It does not save a CSV.
It only prints one RMSD value for one pair of frames.

What this file is best for:
- A quick spot-check of two specific frames.
- A sanity check before running the larger CSV pipeline.
- Teaching/demo use when you want one number and no extra files.

Important behavior to know first:
- Frame numbering is 1-based.
- You must provide exactly two frame numbers: `--frame-a` and `--frame-b`.
- `--atoms` overrides `--atom-scope` when both are present.
- Output is terminal-only by default (no file is created).

#### Barebones command pattern
```powershell
python .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a <first-frame-id> --frame-b <second-frame-id>
```

How to read this pattern:
- `--input` is the path to your multi-frame PDB file.
- `--frame-a` is the first frame (reference side for this one comparison).
- `--frame-b` is the second frame (target side for this one comparison).
- The script aligns matching atoms, then applies:
	RMSD = sqrt((1/N) * sum_{i=1..N} ||r_i - s_i||^2)

Expected terminal output format:
- `RMSD_Angstrom(frame X vs frame Y): <value>`

#### Windows-friendly invocation note
If `python` is not recognized on your machine, use the launcher form:

```powershell
py -3.13 .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a <first-frame-id> --frame-b <second-frame-id>
```

#### Mini cookbook for barebones scenarios (zero doubt)

Scenario A1: basic full-atom 1v1 comparison

```powershell
python .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a 5 --frame-b 10
```

What this does:
- Uses all available atoms shared by both frames after alignment.
- Prints one RMSD value for frame 5 vs frame 10.

Use it when:
- You want the fastest full-structure difference check.

Scenario A2: backbone-only 1v1 comparison

```powershell
python .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a 5 --frame-b 10 --atom-scope backbone
```

What this does:
- Uses only `N`, `CA`, `C`, `O`, `OXT` (when present).
- Ignores side-chain atoms.
- Prints one backbone RMSD value.

Use it when:
- You want global fold/motion signal with less side-chain noise.

Scenario A3: custom atom subset comparison

```powershell
python .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a 5 --frame-b 10 --atoms CA,CB,N
```

What this does:
- Uses only `CA`, `CB`, and `N` atoms.
- Overrides `--atom-scope` if both are provided.

Use it when:
- You need a focused comparison on a specific atom set.

Scenario A4: far-apart frame comparison

```powershell
python .\RMSD_calculation_function_7-21-26.py --input "<path-to-your-pdb-file>" --frame-a 50 --frame-b 500
```

What this does:
- Compares frames farther apart in trajectory time.
- Often useful for seeing larger conformational drift.

Use it when:
- You want early-vs-late trajectory difference quickly.

Scenario A5: run from anywhere with full script path

```powershell
python "<path-to-this-script>\RMSD_calculation_function_7-21-26.py" --input "<path-to-your-pdb-file>" --frame-a 5 --frame-b 10
```

What this does:
- Runs correctly even if your current folder is not the script folder.

Use it when:
- You prefer absolute paths or are working from another directory.

#### Exactly what does NOT happen in barebones mode
- No CSV is produced.
- No matrix is produced.
- No multi-scenario batch output is produced.
- No plotting is produced.

#### Common barebones mistakes and fixes
1. `python was not found`
- Fix: use `py -3.13` instead of `python`.

2. `ModuleNotFoundError: No module named 'Bio'`
- Fix: install dependencies in that interpreter:

```powershell
py -3.13 -m pip install biopython numpy
```

3. `Frame IDs must be within 1..N`
- Fix: choose frame numbers inside the actual frame count.

4. `Frame X has no atoms after filtering`
- Fix: relax your atom filter or use `--atom-scope full`.

5. Path errors
- Fix: use quotes around paths and verify the file exists.

When to use it:
- Use it when you want one quick RMSD number.
- Use it when you do not need a table.

### B. Full-featured RMSD script
Use this file if you want a table or a matrix:
- RMSD_per_frame_biopython_7-21-26.py

This file is the main RMSD script.
It can make CSV tables.
It can make a reference-vs-target table.
It can make an all-vs-all matrix.
It can also make a reference-frames-vs-targets matrix.

What "backbone-only comparison" means:
- The script uses only backbone atom names: `N`, `CA`, `C`, `O`, and `OXT`.
- It ignores side-chain atoms.
- This usually gives a cleaner view of overall protein fold/motion.
- Use it when side-chain flexibility is distracting your trend.
- Use full-atom mode when side-chain behavior is important to your question.

Command-order convention for this file:
1. Always write reference first.
2. Always write targets second.

Use this pattern to avoid confusion:
- `--reference-frame ... --targets ...`
- `--reference-frames ... --targets ...`

Very important output rule:
- `--reference-frame ... --targets ...` gives a simple list CSV.
- `--reference-frames ... --targets ...` gives a matrix CSV.
- If you want reference frames on the vertical axis and target frames on the horizontal axis, use `--reference-frames ... --targets ...`.
- This is true even if you only have one reference frame.

Zero-doubt shortcut:
- Want a simple two-column list of RMSD values? Use `--reference-frame`.
- Want an Excel-style matrix with row labels and column labels? Use `--reference-frames`.
- One reference frame can still be written as a matrix row. In that case, still use `--reference-frames 3`, not `--reference-frame 3`.

#### Scenario B1: one reference frame and many target frames
Use this when you want one reference frame on the left and several target frames across the top.

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 1 --targets 25,45,90
```

What this command means:
- `--reference-frames 1` means frame 1 is the reference row.
- `--targets 25,45,90` means frames 25, 45, and 90 are the target columns.
- The output is one CSV table.
- The row labels are the reference frames.
- The column labels are the target frames.
- This is matrix mode, not simple list mode.

When to use it:
- Use it when you want a table where the rows and columns are not the same.
- This is the easiest way to read a reference-vs-target comparison.
- Use this when you want one reference on the vertical axis and many targets on the horizontal axis.

#### Scenario B2: all-vs-all matrix
Use this when you want every chosen frame compared to every other chosen frame.

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all --all-vs-all
```

What this command means:
- `--compare-all` tells the script to use every frame in the PDB file.
- `--all-vs-all` tells the script to compare all selected frames to each other.
- The output is a square CSV matrix.
- The same frame numbers appear on both the rows and the columns.

When to use it:
- Use it when you want the most complete comparison table.
- Use it when you want to compare every frame to every other frame.

#### Mini cookbook for sample commands (step-by-step)

Scenario B3: default run (all frames vs reference frame 1, full atoms)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>"
```

How to read it:
- `--input` points to your multi-frame PDB.
- `--output-dir` is where the CSV is saved.
- No reference is given, so default reference is frame 1.
- No targets are given, so all available frames are used as targets.
- No atom filter is given, so full-atom mode is used.

What you get:
- One CSV with frame-by-frame RMSD against reference frame 1.
- This is usually the fastest first check of trajectory drift.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B4: backbone-only run (all frames vs reference frame 1)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --atom-scope backbone
```

How to read it:
- Same as Scenario B3, but adds `--atom-scope backbone`.
- Only `N`, `CA`, `C`, `O`, `OXT` are used.
- Side-chain atoms are excluded.

What you get:
- One CSV with backbone RMSD values.
- Often smoother than full-atom RMSD.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

When this is useful:
- You want global structure motion, not side-chain wobble.

Scenario B5: one custom reference frame vs specific target frames (1vN)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 5 --targets 12,48
```

How to read it:
- `--reference-frame 5` sets the single reference frame.
- `--targets 12,48` compares only frames 12 and 48 to frame 5.

What you get:
- A small CSV containing RMSD values for those target frames only.
- This is a simple list CSV, not a matrix.

Important note:
- Even though this is conceptually "1 vs N", it is not written as a row-and-column matrix.
- If you want frame 5 on the vertical axis and target frames across the top, use `--reference-frames 5 --targets ...` instead.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

When this is useful:
- You already know the frames you care about.

Scenario B6: one reference vs many targets (explicit 1vN)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 2,4,5,9
```

How to read it:
- `--reference-frame 1` fixes the reference.
- `--targets 2,4,5,9` lists only those targets.

What you get:
- A compact table of RMSD values for a specific subset of frames.
- This is a simple list CSV, not a matrix.

Important note:
- This output will not place target frame numbers as horizontal column headers.
- To force that orientation, switch to matrix mode with `--reference-frames 1 --targets 2,4,5,9`.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B7: one reference vs target range with custom atom list

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 10-25 --atoms CA,CB,N
```

How to read it:
- `--targets 10-25` includes every frame from 10 to 25.
- `--atoms CA,CB,N` uses only those atom names.
- `--atoms` overrides `--atom-scope` if both are provided.

What you get:
- RMSD values focused on a custom set of atoms.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

When this is useful:
- You want a feature-focused RMSD (for example backbone+CB behavior).

Scenario B8: explicit all-frames selection with default comparison style

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all
```

How to read it:
- `--compare-all` means all frames are selected from the file.
- Comparison is still default Nv1 style unless `--all-vs-all` is also used.

What you get:
- Same shape of output as Scenario B3, but explicitly declared.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B9: reference vs target range

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 10-25
```

How to read it:
- Reference is frame 1.
- Targets are all frames in range 10 through 25.

What you get:
- A focused RMSD table for a continuous segment of your trajectory.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B10: reference vs discrete targets

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 5,12,48,77
```

How to read it:
- Compares frame 1 to only the listed targets.

What you get:
- Sparse targeted comparison values.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B11: reference vs mixed target list and ranges

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 1-10,20,25-30
```

How to read it:
- Combines ranges and single frames in one target selector.

What you get:
- Flexible selection without multiple commands.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B12: different reference frame and custom results folder

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --reference-frame 10 --output-dir "<path-to-results-folder>"
```

How to read it:
- Uses frame 10 as reference.
- Saves output to a folder of your choice.

What you get:
- RMSD profile relative to frame 10 instead of frame 1.

Expected output filename (default):
- `<path-to-results-folder>\\rmsd_per_frame.csv`

Scenario B13: custom atoms plus explicit targets

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frame 1 --targets 1-10,20 --atoms CA,CB,N
```

How to read it:
- Reference first, targets second, then custom atoms.

What you get:
- Feature-specific RMSD over your selected target set.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B14: all-vs-all matrix for every frame in file

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --compare-all --all-vs-all
```

How to read it:
- `--compare-all` selects every frame.
- `--all-vs-all` switches to square matrix mode.

What you get:
- A square matrix CSV where row and column frame sets are identical.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_all_vs_all.csv`

Scenario B15: many references vs one target (Nvs1 matrix)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 2,4,5,9,12 --targets 45
```

How to read it:
- `--reference-frames` defines multiple row frames.
- `--targets 45` defines one column frame.

What you get:
- A matrix-style CSV with references on rows and target(s) on columns.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Scenario B16: one reference vs many targets (same structure as Scenario B1)

```powershell
python .\RMSD_per_frame_biopython_7-21-26.py --input "<path-to-your-pdb-file>" --output-dir "<path-to-output-folder>" --reference-frames 1 --targets 25,45,90
```

How to read it:
- Uses matrix mode spelling (`--reference-frames`) for a single reference.
- Keeps orientation fixed: reference rows, target columns.

What you get:
- A clean table ready for Excel with reference vertical and targets horizontal.

Why this scenario matters:
- Some users think `--reference-frame 3 --targets ...` will make a horizontal matrix.
- It will not.
- To force matrix orientation, use `--reference-frames 3 --targets ...` instead.

Simple memory rule:
- Singular `--reference-frame` = list output.
- Plural `--reference-frames` = matrix output.

Expected output filename (default):
- `<path-to-output-folder>\\rmsd_per_frame.csv`

Beginner checklist before running any command:
- Confirm the PDB file path exists.
- Confirm frame IDs you typed are within the file range.
- Confirm output folder exists and is writable.
- If using Windows and `python` is not found, use your explicit interpreter path.

Filename override note:
- You can change the saved filename in any command by adding `--output-name your_file_name.csv`.
- If you do that, the expected default names above no longer apply.

#### What the output looks like
1. Reference-vs-target table:
- Reference frames are on the left side as rows.
- Target frames are across the top as columns.
- Each cell in the table is one RMSD value.

Example:
- Row 1 = reference frame 1
- Columns = target frames 25, 45, 90

This means the table shows how frame 1 compares to each of those target frames.

2. All-vs-all table:
- Every frame is compared to every other frame.
- The table is square because the row frames and column frames are the same set.
- This is the most complete comparison option.

3. Barebones terminal output:
- One RMSD value is printed directly in the terminal.
- No table is made.
- No CSV is created.

---

## 2. Classify frames
Use this file:
- Classification_of_Frames_7-19-26.py

It groups frames into classes after you make the RMSD table.

What this step does:
- It reads the RMSD matrix or table output from Step 1.
- It groups similar frames together based on a threshold.
- It writes new CSV outputs that show which frames belong together.

---

## 3. Make plots (optional)
Use this file:
- Backbone_vs_Full_Atom_RMSD_Histogram_7-19-26.py

It makes simple comparison plots.

What this step does:
- It compares backbone-only behavior to full-atom behavior.
- It makes a plot or histogram so you can see the difference visually.
- This step is optional.

---

## Recommended order
1. Make the RMSD table
2. Classify the frames
3. Make plots if you need them

If you are brand new, just start with Step 1 and only do Step 2 or Step 3 after Step 1 works.

---

## Common problems
1. File not found
- Check the input path.
- Make sure the PDB file path is typed exactly right.

2. Wrong frame number
- Make sure the number exists in the PDB file.
- Make sure the frame number is not bigger than the number of frames in the file.

3. No output file
- Make sure the output folder exists.
- Make sure you typed the output folder path correctly.
