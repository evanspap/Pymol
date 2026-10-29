"""MD Analysis — protein and trajectory loader for PyMOL.
Version: 1.0.2
Author: OpenAI Codex, for the MD_classification project.

Summary:
Load a protein structure and optional matching trajectory as states (frames)
of one PyMOL object. Multi-model PDB files work without a separate trajectory.
Optional principal-axis alignment centers each state's protein at the origin
and maps its major, middle and minor axes to X, Y and Z.

Usage in PyMOL:
    run /path/to/Pymol_plugin/__init__.py
    md_analysis_gui
    md_load protein.pdb, trajectory.xtc
    md_load protein.pdb, trajectory.xtc, align_principal_axes=1
    md_load ensemble.pdb

Standalone usage:
    python __init__.py --structure protein.pdb --trajectory trajectory.xtc
    python __init__.py --structure ensemble.pdb --dry-run

Files and frame options come from arguments. Dry run checks paths and options
without importing PyMOL or loading coordinates. Missing arguments show this
header, including usage, examples, author and summary. No files are written.
Topology and trajectory must have the same atoms in the same order.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

__version__ = "1.0.2"
_DIALOGS = []


def _command(parts):
    return subprocess.list2cmdline(parts) if os.name == "nt" else shlex.join(parts)


def _options(structure, trajectory="", name="md_protein", start=1, stop=-1, interval=1):
    structure = Path(structure).expanduser().resolve()
    trajectory = Path(trajectory).expanduser().resolve() if trajectory else None
    for path in (structure, trajectory):
        if path is not None and not path.is_file():
            raise ValueError(f"File not found: {path}")
    name = name.strip() or "md_protein"
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name):
        raise ValueError("Object name must start with a letter and contain letters, numbers or underscores.")
    start, stop, interval = int(start), int(stop), int(interval)
    if start < 1 or interval < 1 or (stop != -1 and stop < start):
        raise ValueError("Start and interval must be positive; stop must be -1 (all) or at least start.")
    if trajectory is None and (start, stop, interval) != (1, -1, 1):
        raise ValueError("Frame range options apply to separate trajectory files. Use defaults for a multi-model PDB.")
    return structure, trajectory, name, start, stop, interval


def _align_principal_axes(name):
    """Center each state on its protein and rotate its PCA axes to XYZ."""
    import numpy as np
    from pymol import cmd

    previous_axes = None
    state_count = cmd.count_states(name)
    protein_selection = f"({name}) and polymer.protein"
    for state in range(1, state_count + 1):
        coordinates = np.asarray(cmd.get_coords(name, state=state), dtype=float)
        protein_coordinates = cmd.get_coords(protein_selection, state=state)
        if protein_coordinates is None or len(protein_coordinates) == 0:
            protein_coordinates = coordinates
        protein_coordinates = np.asarray(protein_coordinates, dtype=float)
        center = protein_coordinates.mean(axis=0)
        centered_protein = protein_coordinates - center
        covariance = centered_protein.T @ centered_protein
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        axes = eigenvectors[:, np.argsort(eigenvalues)[::-1]]

        if np.linalg.det(axes) < 0:
            axes[:, 2] *= -1
        if previous_axes is None:
            for axis in (0, 1):
                dominant_component = np.argmax(np.abs(axes[:, axis]))
                if axes[dominant_component, axis] < 0:
                    axes[:, axis] *= -1
            if np.linalg.det(axes) < 0:
                axes[:, 2] *= -1
        else:
            sign_options = np.asarray(((1, 1, 1), (1, -1, -1),
                                       (-1, 1, -1), (-1, -1, 1)))
            candidates = axes[:, :, np.newaxis] * sign_options.T[np.newaxis, :, :]
            scores = np.einsum("ij,ijk->k", previous_axes, candidates)
            axes = candidates[:, :, np.argmax(scores)]

        aligned_coordinates = (coordinates - center) @ axes
        cmd.load_coords(aligned_coordinates, name, state=state)
        previous_axes = axes


def md_load(structure, trajectory="", name="md_protein", start=1, stop=-1,
            interval=1, dry_run=0, align_principal_axes=0):
    """Load structure plus optional trajectory; return the unique object name."""
    structure, trajectory, name, start, stop, interval = _options(
        structure, trajectory, name, start, stop, interval)
    dry_run = bool(int(dry_run))
    parts = ["md_load", str(structure)]
    if trajectory:
        parts += ["trajectory=" + str(trajectory)]
    align_principal_axes = bool(int(align_principal_axes))
    parts += [f"name={name}", f"start={start}", f"stop={stop}",
              f"interval={interval}", f"dry_run={int(dry_run)}",
              f"align_principal_axes={int(align_principal_axes)}"]
    print("Running: " + _command(parts), flush=True)
    if dry_run:
        print("Dry run: paths and options valid. No coordinates loaded or files written.")
        print("Atom order, atom count and trajectory format are checked during loading.")
        return None
    from pymol import cmd
    existing = set(cmd.get_names("all"))
    base = name
    number = 2
    while name in existing:
        name = f"{base}_{number}"
        number += 1
    try:
        # multiplex=0 keeps MODEL records together as states of a single object.
        cmd.load(str(structure), name, multiplex=0, discrete=0, zoom=0)
        if not cmd.count_atoms(name):
            raise ValueError("The structure contains no atoms.")
        if trajectory:
            if cmd.count_states(name) != 1:
                raise ValueError("For a separate trajectory, choose a single-state topology/structure.")
            # Populate from state 1: the topology snapshot is not an extra frame.
            cmd.load_traj(str(trajectory), name, state=1, start=start,
                          stop=stop, interval=interval)
        if align_principal_axes:
            _align_principal_axes(name)
        cmd.hide("everything", name)
        cmd.show("cartoon", f"({name}) and polymer.protein")
        cmd.show("lines", f"({name}) and not polymer.protein")
        cmd.set("all_states", 0, name)
        cmd.set("static_singletons", 0, name)
        cmd.frame(1)
        cmd.zoom(name)
        count = cmd.count_states(name)
        print(f"Loaded {name}: {count} frame(s). Use the PyMOL frame slider or Play to view the trajectory.", flush=True)
        return name
    except Exception:
        # Only the new object is removed; objects from earlier loads are retained.
        cmd.delete(name)
        raise


def md_analysis_gui():
    """Open the MD Analysis loader window inside PyMOL."""
    from pymol.Qt import QtWidgets
    from pymol import plugins
    try:
        parent = plugins.get_qtwindow()
    except (AttributeError, RuntimeError):
        parent = None
    dialog = QtWidgets.QDialog(parent)
    dialog.setWindowTitle("MD Analysis")
    dialog.setMinimumWidth(650)
    layout = QtWidgets.QVBoxLayout(dialog)
    form = QtWidgets.QFormLayout()
    structure = QtWidgets.QLineEdit(dialog)
    structure.setPlaceholderText("Select a topology file (.gro or .pdb)")
    trajectory = QtWidgets.QLineEdit(dialog)
    trajectory.setPlaceholderText("Optional: XTC, TRR, DCD or another PyMOL-supported trajectory")
    for label, edit, filters in (
        ("Topology (.gro or .pdb)", structure, "Topology files (*.gro *.pdb);;GROMACS topology (*.gro);;PDB structure (*.pdb);;Other structures (*.cif *.mmcif *.psf *.mol2);;All files (*)"),
        ("Trajectory (optional)", trajectory, "Trajectories (*.xtc *.trr *.dcd *.trj *.nc *.netcdf);;All files (*)"),
    ):
        row = QtWidgets.QHBoxLayout()
        row.addWidget(edit)
        button = QtWidgets.QPushButton("Browse...", dialog)
        def browse(_checked=False, field=edit, title=label, file_filter=filters):
            path, _ = QtWidgets.QFileDialog.getOpenFileName(dialog, title, "", file_filter)
            if path:
                field.setText(path)
        button.clicked.connect(browse)
        row.addWidget(button)
        form.addRow(label, row)
    name = QtWidgets.QLineEdit("md_protein", dialog)
    form.addRow("Object name", name)
    align_axes = QtWidgets.QCheckBox(
        "Align each frame to protein principal axes and center at (0, 0, 0)", dialog)
    form.addRow("Principal-axis alignment", align_axes)
    spins = []
    for label, minimum, maximum, default in (
        ("First trajectory frame (1-based)", 1, 2147483647, 1),
        ("Last trajectory frame (-1 = all)", -1, 2147483647, -1),
        ("Load every Nth frame", 1, 2147483647, 1),
    ):
        spin = QtWidgets.QSpinBox(dialog)
        spin.setRange(minimum, maximum)
        spin.setValue(default)
        form.addRow(label, spin)
        spins.append(spin)
    def update_frames():
        for spin in spins:
            spin.setEnabled(bool(trajectory.text().strip()))
    trajectory.textChanged.connect(update_frames)
    update_frames()
    layout.addLayout(form)
    note = QtWidgets.QLabel(
        "For a multi-model PDB, leave Trajectory empty. For a separate trajectory, "
        "the topology must match its atom count and order. Each frame becomes a PyMOL state.", dialog)
    note.setWordWrap(True)
    layout.addWidget(note)
    log = QtWidgets.QPlainTextEdit(dialog)
    log.setReadOnly(True)
    log.setMaximumHeight(130)
    layout.addWidget(log)
    buttons = QtWidgets.QHBoxLayout()
    load_button = QtWidgets.QPushButton("Load", dialog)
    dry_button = QtWidgets.QPushButton("Dry run", dialog)
    close_button = QtWidgets.QPushButton("Close", dialog)
    def load(dry):
        try:
            values = [spin.value() for spin in spins] if trajectory.text().strip() else [1, -1, 1]
            _options(structure.text().strip(), trajectory.text().strip(), name.text(), *values)
            log.appendPlainText("Checking paths..." if dry else "Loading... Large trajectories may take time; see the PyMOL console.")
            for button in (load_button, dry_button, close_button):
                button.setEnabled(False)
            dialog.repaint()
            result = md_load(structure.text().strip(), trajectory.text().strip(),
                             name.text(), *values, dry_run=int(dry),
                             align_principal_axes=int(align_axes.isChecked()))
            if dry:
                log.appendPlainText("Dry run passed. No coordinates loaded. Atom compatibility is checked on Load.")
            else:
                from pymol import cmd
                log.appendPlainText(f"Loaded {result}: {cmd.count_states(result)} frame(s). Use PyMOL's frame slider or Play.")
        except Exception as exc:
            log.appendPlainText(f"Error: {exc}")
            QtWidgets.QMessageBox.warning(dialog, "MD Analysis", str(exc))
        finally:
            for button in (load_button, dry_button, close_button):
                button.setEnabled(True)
    load_button.clicked.connect(lambda: load(False))
    dry_button.clicked.connect(lambda: load(True))
    close_button.clicked.connect(dialog.close)
    for button in (load_button, dry_button, close_button):
        buttons.addWidget(button)
    layout.addLayout(buttons)
    dialog.finished.connect(lambda _code: _DIALOGS.remove(dialog) if dialog in _DIALOGS else None)
    _DIALOGS.append(dialog)
    dialog.show()
    return dialog


def _register():
    try:
        from pymol import cmd
    except ImportError:
        return
    cmd.extend("md_load", md_load)
    cmd.extend("md_analysis_gui", md_analysis_gui)


def __init_plugin__(app=None):
    _register()
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("MD Analysis", md_analysis_gui)


class _Parser(argparse.ArgumentParser):
    def error(self, message):
        print(__doc__, flush=True)
        super().error(message)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    print("Running: " + _command([sys.executable, __file__, *argv]), flush=True)
    parser = _Parser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--trajectory", default="")
    parser.add_argument("--name", default="md_protein")
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--stop", type=int, default=-1)
    parser.add_argument("--interval", type=int, default=1)
    parser.add_argument("--align-principal-axes", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--version", action="version", version="MD Analysis " + __version__)
    args = parser.parse_args(argv)
    try:
        _options(args.structure, args.trajectory, args.name, args.start, args.stop, args.interval)
        if not args.dry_run:
            import pymol
            pymol.finish_launching()
        md_load(args.structure, args.trajectory, args.name, args.start,
                args.stop, args.interval, int(args.dry_run),
                int(args.align_principal_axes))
    except (ValueError, OSError, ImportError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
_register()

