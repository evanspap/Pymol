# Name: Protein Tools
# Version: 1.1
# Date: 2026-11-12
# Author: William Sandor, Evangelos Papadopoulos, Gerasimos Papadopoulos
# Summary: Load selection, neighborhood, and CSV export tools as one PyMOL plugin.
# Usage: Install Protein_Tools as a plugin; commands include label_around, show_around, save_selections, and save_objects.

"""PyMOL plugin package that loads all Protein Tools subplugins."""

from __future__ import annotations

from . import ShowAround, label_around, save_objects, save_selection

_PLUGINS = (label_around, ShowAround, save_selection, save_objects)


def __init_plugin__(app=None):
    """Initialize each child plugin and register its menu item."""
    for plugin in _PLUGINS:
        plugin.__init_plugin__(app)