"""Draw a fixed, labeled XYZ coordinate triad in the PyMOL scene.

Version: 1.0.0

PyMOL usage:
    run /path/to/Pymol_plugins/axes/__init__.py
    draw_xyz_axes

The Plugin menu also provides ``Draw XYZ Axes`` after installation.
"""


def draw_xyz_axes():
    """Draw red X, green Y, and blue Z cylinders from the origin."""
    from pymol import cmd
    from pymol.cgo import CYLINDER, cyl_text
    from pymol.vfont import plain

    axes_object = [
        CYLINDER, 0., 0., 0., 10., 0., 0., 0.2, 1.0, 0.0, 0.0, 1.0, 0.0, 0.,
        CYLINDER, 0., 0., 0., 0., 10., 0., 0.2, 0.0, 1.0, 0.0, 0., 1.0, 0.,
        CYLINDER, 0., 0., 0., 0., 0., 10., 0.2, 0.0, 0.0, 1.0, 0., 0.0, 1.0,
    ]
    label_axes = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
    cyl_text(axes_object, plain, [11., 0., 0.], "X", 0.25, axes=label_axes)
    cyl_text(axes_object, plain, [0., 11., 0.], "Y", 0.25, axes=label_axes)
    cyl_text(axes_object, plain, [0., 0., 11.], "Z", 0.25, axes=label_axes)

    cmd.delete("axes")
    cmd.load_cgo(axes_object, "axes")
    cmd.orient()
    return "axes"


def _register():
    try:
        from pymol import cmd
    except ImportError:
        return
    cmd.extend("draw_xyz_axes", draw_xyz_axes)


def __init_plugin__(app=None):
    _register()
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Draw XYZ Axes", draw_xyz_axes)


_register()