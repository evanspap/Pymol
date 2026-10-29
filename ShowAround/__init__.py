# Name: Show Around
# Version: 1.0
# Date: 2026-10-29
# Author: William Sandor, Evangelos Papadopoulos, Gerasimos Papadopoulos
# Summary: Select complete protein residues within a distance of a source selection.
# Usage: show_around [selection [, distance [, output_name]]] (defaults to sele, 6 A, nearby_residues)

"""PyMOL plugin for selecting protein residues near a source selection."""

from __future__ import annotations

import math
import re

try:
    from pymol import cmd
except Exception:  # Allow importing and testing outside PyMOL.
    cmd = None


def show_around(selection="sele", distance=6.0, name="nearby_residues"):
    """Create a selection of whole protein residues near another selection.

    Usage: show_around [selection [, distance [, name]]]

    Arguments:
      selection: Existing PyMOL selection used as the distance reference.
                 Defaults to 'sele', PyMOL's current selection.
      distance:  Maximum atom-to-atom distance in Angstroms. Defaults to 6.0.
    name:      Name for the new output selection. Defaults to
             'nearby_residues'. This
                 is an ordinary selection name, not a special PyMOL keyword.

    Any protein residue with at least one atom within the distance is included
    in full. For example, 'show_around ligand, 8, nearby' measures from ligand,
    uses an 8 Angstrom cutoff, and stores the result as the selection 'nearby'.
    The result is returned as the number of selected atoms.
    """
    if cmd is None:
        raise RuntimeError("ShowAround must be run inside PyMOL.")

    selection = str(selection).strip()
    if not selection:
        raise ValueError("Source selection cannot be empty.")

    try:
        distance = float(distance)
    except (TypeError, ValueError):
        raise ValueError("Distance must be a positive number in Angstroms.")
    if not math.isfinite(distance) or distance <= 0:
        raise ValueError("Distance must be a positive number in Angstroms.")

    name = str(name).strip()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name):
        raise ValueError("Selection name must start with a letter or underscore and contain only letters, numbers, or underscores.")

    if cmd.count_atoms(selection) == 0:
        raise ValueError(f"Source selection '{selection}' is empty or does not exist.")

    expression = f"byres ((polymer.protein) within {distance:g} of ({selection}))"
    count = cmd.select(name, expression, quiet=1)
    print(f"[ShowAround] Created '{name}' with {count} protein atoms within {distance:g} A of '{selection}'.")
    return count


def _run_default():
    try:
        show_around()
    except Exception as exc:
        print(f"[ShowAround] {exc}")


def _register_commands():
    if cmd is not None:
        cmd.extend("show_around", show_around)


def __init_plugin__(app=None):
    _register_commands()
    try:
        from pymol.plugins import addmenuitem

        addmenuitem("Protein Tools|Show Around (6 A)", _run_default)
    except Exception:
        try:
            app.rootmenu.addmenuitem(
                "Plugin",
                "command",
                "Show Around (6 A)",
                label="Show Around (6 A)",
                command=_run_default,
            )
        except Exception:
            pass


_register_commands()