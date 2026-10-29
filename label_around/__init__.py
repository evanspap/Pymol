"""PyMOL plugin for labeling alpha-carbon atoms with residue identifiers."""

from __future__ import annotations

try:
    from pymol import cmd
except Exception:  # Allow importing and testing outside PyMOL.
    cmd = None


_AMINO_ACID_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "ASX": "B", "GLX": "Z", "SEC": "U", "PYL": "O", "MSE": "M",
}


def label_around(selection="sele"):
    """Label selected CA atoms as one-letter residue code plus residue number.

    Usage: label_around [selection]

    Args:
      selection: PyMOL selection to label; defaults to 'sele'. Only atoms named
                 CA are labeled.
    Labels use the residue identifier (resi), including insertion codes. For
    example, an alanine with resi '42' is labeled 'A42'. Unknown residue names
    are labeled with 'X'. Labels remain at their CA atom positions.
    """
    if cmd is None:
        raise RuntimeError("label_around must be run inside PyMOL.")

    selection = str(selection).strip()
    if not selection:
        raise ValueError("Selection cannot be empty.")

    ca_selection = f"({selection}) and name CA"
    count = cmd.count_atoms(ca_selection)
    if count == 0:
        raise ValueError(f"Selection '{selection}' contains no CA atoms.")

    expression = f"{_AMINO_ACID_CODES!r}.get(resn, 'X') + str(resi)"
    cmd.label(ca_selection, expression, quiet=1)
    cmd.unset("label_position", ca_selection, quiet=1)
    cmd.unset("label_placement_offset", ca_selection, quiet=1)
    print(f"[label_around] Labeled {count} CA atoms in '{selection}' at their atom positions.")
    return count


def _run_default():
    try:
        label_around()
    except Exception as exc:
        print(f"[label_around] {exc}")


def _register_commands():
    if cmd is not None:
        cmd.extend("label_around", label_around)


def __init_plugin__(app=None):
    _register_commands()
    try:
        from pymol.plugins import addmenuitem

        addmenuitem("Protein Tools|Label CA Residues", _run_default)
    except Exception:
        try:
            app.rootmenu.addmenuitem(
                "Plugin",
                "command",
                "Label CA Residues",
                label="Label CA Residues",
                command=_run_default,
            )
        except Exception:
            pass


_register_commands()