# Name: Save Selections
# Version: 1.2
# Date: 2026-11-12
# Author: William Sandor, Evangelos Papadopoulos, Gerasimos Papadopoulos
# Summary: Choose named selections from a checklist and export each to a CA residue CSV file.
# Usage: save_selections [directory] (opens a checklist; defaults to PyMOL's current working directory)

"""PyMOL plugin for exporting named selections as residue CSV files."""

from __future__ import annotations

import csv
import re
from pathlib import Path

try:
    from pymol import cmd
except Exception:  # Allow importing and testing outside PyMOL.
    cmd = None


CSV_HEADER = ("Chain", "Res_Num", "Res_Name")


def _safe_stem(name):
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return stem or "selection"


def _unique_filename(name, used):
    stem = _safe_stem(name)
    candidate = f"{stem}_residues.csv"
    suffix = 2
    while candidate in used:
        candidate = f"{stem}_{suffix}_residues.csv"
        suffix += 1
    used.add(candidate)
    return candidate


def export_selections(selection_names, directory=""):
    """Write selected named selections to separate CSV files."""
    if cmd is None:
        raise RuntimeError("save_selections must be run inside PyMOL.")

    output_dir = Path(directory).expanduser() if str(directory).strip() else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)
    names = list(dict.fromkeys(selection_names))
    used_filenames = set()
    written_files = []
    for selection_name in names:
        rows = []
        ca_selection = f"({selection_name}) and name CA"
        cmd.iterate(
            ca_selection,
            "rows.append((chain, resi, resn))",
            space={"rows": rows},
            quiet=1,
        )

        output_path = output_dir / _unique_filename(selection_name, used_filenames)
        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(CSV_HEADER)
            writer.writerows(rows)
        written_files.append(str(output_path))
        print(f"[save_selections] Wrote {len(rows)} CA residues from '{selection_name}' to {output_path}.")

    return written_files


def _open_selection_dialog(selection_names, directory=""):
    try:
        from pymol import plugins
        from pymol.Qt import QtCore, QtWidgets
    except Exception as exc:
        raise RuntimeError(f"The selection checklist requires PyMOL's Qt interface: {exc}")

    parent = None
    try:
        parent = plugins.get_qtwindow()
    except Exception:
        pass

    dialog = QtWidgets.QDialog(parent)
    dialog.setWindowTitle("Save Selected Selections")
    dialog.setMinimumWidth(420)
    layout = QtWidgets.QVBoxLayout(dialog)
    layout.addWidget(QtWidgets.QLabel("Choose named selections to export:", dialog))

    list_widget = QtWidgets.QListWidget(dialog)
    for name in selection_names:
        item = QtWidgets.QListWidgetItem(name)
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(QtCore.Qt.Unchecked)
        list_widget.addItem(item)
    layout.addWidget(list_widget)

    list_buttons = QtWidgets.QHBoxLayout()
    select_all_button = QtWidgets.QPushButton("Select All", dialog)
    select_none_button = QtWidgets.QPushButton("Select None", dialog)
    list_buttons.addWidget(select_all_button)
    list_buttons.addWidget(select_none_button)
    layout.addLayout(list_buttons)
    select_all_button.clicked.connect(
        lambda: [list_widget.item(i).setCheckState(QtCore.Qt.Checked) for i in range(list_widget.count())]
    )
    select_none_button.clicked.connect(
        lambda: [list_widget.item(i).setCheckState(QtCore.Qt.Unchecked) for i in range(list_widget.count())]
    )

    folder_row = QtWidgets.QHBoxLayout()
    folder_edit = QtWidgets.QLineEdit(str(Path(directory).expanduser() if directory else Path.cwd()), dialog)
    browse_button = QtWidgets.QPushButton("Browse...", dialog)
    folder_row.addWidget(folder_edit)
    folder_row.addWidget(browse_button)
    layout.addLayout(folder_row)

    def browse_folder():
        path = QtWidgets.QFileDialog.getExistingDirectory(dialog, "CSV output folder", folder_edit.text())
        if path:
            folder_edit.setText(path)

    browse_button.clicked.connect(browse_folder)

    action_row = QtWidgets.QHBoxLayout()
    save_button = QtWidgets.QPushButton("Save CSVs", dialog)
    cancel_button = QtWidgets.QPushButton("Cancel", dialog)
    action_row.addStretch(1)
    action_row.addWidget(save_button)
    action_row.addWidget(cancel_button)
    layout.addLayout(action_row)
    cancel_button.clicked.connect(dialog.reject)

    def save_checked():
        selected_names = [
            list_widget.item(i).text()
            for i in range(list_widget.count())
            if list_widget.item(i).checkState() == QtCore.Qt.Checked
        ]
        if not selected_names:
            QtWidgets.QMessageBox.warning(dialog, "Save Selections", "Select at least one selection.")
            return
        try:
            files = export_selections(selected_names, folder_edit.text())
        except Exception as exc:
            QtWidgets.QMessageBox.critical(dialog, "Save Selections", str(exc))
            return
        QtWidgets.QMessageBox.information(dialog, "Save Selections", f"Saved {len(files)} CSV file(s).")
        dialog.accept()

    save_button.clicked.connect(save_checked)
    execute = getattr(dialog, "exec", None) or dialog.exec_
    execute()


def save_selections(directory=""):
    """Open a checklist to choose named selections for CA-residue CSV export.

    Usage: save_selections [directory]

    The checklist includes enabled and disabled named selections. Select the
    entries to export; the CSV columns are Chain, Res_Num, and Res_Name. The
    default output directory is PyMOL's current working directory.
    """
    if cmd is None:
        raise RuntimeError("save_selections must be run inside PyMOL.")

    selection_names = cmd.get_names("selections", enabled_only=0) or []
    if not selection_names:
        print("[save_selections] No named selections found.")
        return []
    _open_selection_dialog(selection_names, directory)
    return []


def _run_default():
    try:
        save_selections()
    except Exception as exc:
        print(f"[save_selections] {exc}")


def _register_commands():
    if cmd is not None:
        cmd.extend("save_selections", save_selections)


def __init_plugin__(app=None):
    _register_commands()
    try:
        from pymol.plugins import addmenuitem

        addmenuitem("Protein Tools|Choose Selections to Save...", _run_default)
    except Exception:
        try:
            app.rootmenu.addmenuitem(
                "Plugin",
                "command",
                "Choose Selections to Save...",
                label="Choose Selections to Save...",
                command=_run_default,
            )
        except Exception:
            pass


_register_commands()