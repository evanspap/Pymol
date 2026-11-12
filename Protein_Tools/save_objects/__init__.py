# Name: Save Enabled Objects
# Version: 1.0
# Date: 2026-11-12
# Author: William Sandor, Evangelos Papadopoulos, Gerasimos Papadopoulos
# Summary: Export enabled molecular objects to individual CA residue CSV files.
# Usage: save_objects [directory] (defaults to PyMOL's current working directory)

"""PyMOL plugin for exporting enabled molecular objects as residue CSV files."""

from __future__ import annotations

try:
    from pymol import cmd
except Exception:  # Allow importing and testing outside PyMOL.
    cmd = None

from ..save_selection import export_selections


def save_objects(directory=""):
    """Export each enabled molecular object as a separate CA residue CSV."""
    if cmd is None:
        raise RuntimeError("save_objects must be run inside PyMOL.")

    enabled_names = cmd.get_names("objects", enabled_only=1) or []
    object_names = [name for name in enabled_names if cmd.get_type(name) == "object:molecule"]
    if not object_names:
        print("[save_objects] No enabled molecular objects found.")
        return []
    return export_selections(object_names, directory)


def _run_default():
    try:
        save_objects()
    except Exception as exc:
        print(f"[save_objects] {exc}")


def _register_commands():
    if cmd is not None:
        cmd.extend("save_objects", save_objects)


def __init_plugin__(app=None):
    _register_commands()
    try:
        from pymol.plugins import addmenuitem

        addmenuitem("Protein Tools|Save Enabled Objects (CSV)", _run_default)
    except Exception:
        try:
            app.rootmenu.addmenuitem(
                "Plugin",
                "command",
                "Save Enabled Objects (CSV)",
                label="Save Enabled Objects (CSV)",
                command=_run_default,
            )
        except Exception:
            pass


_register_commands()