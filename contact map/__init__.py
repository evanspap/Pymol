"""Contact Map: display a live CA-residue contact map for a PyMOL object.

Version: 1.6.0
Purpose:
    Open a separate Qt window showing current-frame cutoff contacts or
    residue-center distances, plus mean pairwise distances and standard
    deviations across frames. Mean and SD plots can include per-residue column
    profiles. Split-diagonal plots have independent scales for both halves.

Usage in PyMOL:
    run /path/to/Pymol_plugins/contact_map/__init__.py
    contact_map_gui

Install the folder as a PyMOL plugin to add Contact Map to the Plugin menu.
Requires PyMOL's bundled Qt and NumPy; it does not require Matplotlib.
"""
from __future__ import annotations

from pathlib import Path
import re
import csv

import numpy as np

__version__ = "1.6.0"
_WINDOWS = []
_VIRIDIS_STOPS = np.asarray((
    (0.267, 0.005, 0.329),
    (0.283, 0.141, 0.458),
    (0.254, 0.265, 0.530),
    (0.207, 0.372, 0.553),
    (0.164, 0.471, 0.558),
    (0.128, 0.567, 0.551),
    (0.135, 0.659, 0.518),
    (0.267, 0.749, 0.441),
    (0.478, 0.821, 0.318),
    (0.741, 0.873, 0.150),
    (0.993, 0.906, 0.144),
))


def _contact_matrix(coordinates, chains, cutoff, minimum_sequence_separation):
    """Return the binary contact map for one frame."""
    distances = _distance_matrix(coordinates)
    contacts = distances <= float(cutoff)
    contacts[_excluded_pairs(chains, minimum_sequence_separation)] = False
    return contacts


def _distance_matrix(coordinates):
    """Return all pairwise Euclidean distances between CA coordinates."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("Coordinates must have shape (residues, 3).")
    differences = coordinates[:, np.newaxis, :] - coordinates[np.newaxis, :, :]
    return np.linalg.norm(differences, axis=2)


def _column_profile(matrix):
    """Average eligible off-diagonal values for each residue column."""
    values = np.asarray(matrix, dtype=float).copy()
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("A per-residue profile requires a square matrix.")
    np.fill_diagonal(values, np.nan)
    valid_counts = np.sum(np.isfinite(values), axis=0)
    totals = np.nansum(values, axis=0)
    means = np.full(values.shape[1], np.nan, dtype=float)
    np.divide(totals, valid_counts, out=means, where=valid_counts > 0)
    return means


def _write_matrix_csv(path, matrix, labels):
    """Write a labeled raw matrix, representing filtered entries as NA."""
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (len(labels), len(labels)):
        raise ValueError("The distance matrix and residue labels do not match.")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["residue", *labels])
        for label, row in zip(labels, matrix):
            writer.writerow([
                label,
                *[f"{value:.12g}" if np.isfinite(value) else "NA" for value in row],
            ])


def _excluded_pairs(chains, minimum_sequence_separation):
    """Mark diagonal and nearby same-chain residue pairs as unavailable."""
    chains = np.asarray(chains, dtype=object)
    positions = np.arange(len(chains))
    nearby_in_same_chain = (
        (np.abs(positions[:, np.newaxis] - positions[np.newaxis, :])
         <= int(minimum_sequence_separation))
        & (chains[:, np.newaxis] == chains[np.newaxis, :])
    )
    return nearby_in_same_chain | np.eye(len(chains), dtype=bool)


def _viridis_pixels(matrix, scale_mode="log", max_value=None):
    """Map values from zero to a selected maximum using log1p or linear scale."""
    values = np.asarray(matrix, dtype=float)
    pixels = np.full((*values.shape, 3), 255, dtype=np.uint8)
    finite = np.isfinite(values)
    if not np.any(finite):
        return pixels, None, None

    maximum = max(0.0, float(values[finite].max())) if max_value is None else float(max_value)
    if maximum < 0 or (maximum == 0 and max_value is not None):
        raise ValueError("The color-scale maximum must be greater than zero.")
    clipped = np.clip(values[finite], 0.0, maximum)
    normalized = np.zeros(values.shape, dtype=float)
    if maximum == 0:
        normalized[finite] = 0.0
    elif scale_mode == "log":
        normalized[finite] = np.log1p(clipped) / np.log1p(maximum)
    elif scale_mode == "linear":
        normalized[finite] = clipped / maximum
    else:
        raise ValueError("Scale mode must be 'log' or 'linear'.")

    palette_position = np.clip(normalized, 0.0, 1.0) * (len(_VIRIDIS_STOPS) - 1)
    lower = np.floor(palette_position).astype(int)
    upper = np.minimum(lower + 1, len(_VIRIDIS_STOPS) - 1)
    fraction = (palette_position - lower)[..., np.newaxis]
    colors = _VIRIDIS_STOPS[lower] * (1.0 - fraction) + _VIRIDIS_STOPS[upper] * fraction
    pixels[finite] = np.rint(colors[finite] * 255).astype(np.uint8)
    return pixels, 0.0, maximum


def _scale_ticks(maximum, scale_mode):
    if maximum is None:
        return []
    if maximum <= 0:
        return [0.0]
    if scale_mode == "linear":
        return [maximum * fraction for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)]
    ticks = [0.0, maximum]
    first_power = int(np.floor(np.log10(maximum)))
    for power in range(first_power, -1, -1):
        tick = 10.0 ** power
        if 0 < tick < maximum:
            ticks.append(tick)
    return sorted(set(ticks))


def _format_scale_value(value):
    return f"{value:.0e}" if 0 < value < 0.01 else f"{value:.3g}"


def _load_qt():
    from pymol.Qt import QtCore, QtGui, QtWidgets
    return QtCore, QtGui, QtWidgets


class ContactMapCanvas:
    """Create a Qt widget that paints a residue contact matrix."""

    def __new__(cls, *args, **kwargs):
        QtCore, QtGui, QtWidgets = _load_qt()

        class CanvasWidget(QtWidgets.QWidget):
            def __init__(self, parent=None):
                super().__init__(parent)
                self._panels = []
                self._message = "No CA residues to plot."
                self.setMinimumSize(760, 540)

            def set_panels(self, panels, scale_mode, max_value,
                           split_scale_settings=None, show_profile=False):
                self._panels = []
                for panel in panels:
                    if panel.get("split"):
                        mean_matrix = np.asarray(panel["mean_matrix"], dtype=float)
                        standard_deviation = np.asarray(panel["std_matrix"], dtype=float)
                        if split_scale_settings is None:
                            mean_scale_mode, mean_max_value = scale_mode, max_value
                            std_scale_mode, std_max_value = scale_mode, max_value
                        else:
                            mean_scale_mode, mean_max_value = split_scale_settings["mean"]
                            std_scale_mode, std_max_value = split_scale_settings["std"]
                        mean_pixels, mean_minimum, mean_maximum = _viridis_pixels(
                            mean_matrix, mean_scale_mode, mean_max_value)
                        std_pixels, std_minimum, std_maximum = _viridis_pixels(
                            standard_deviation, std_scale_mode, std_max_value)
                        upper = np.triu_indices(mean_matrix.shape[0], k=0)
                        lower = np.tril_indices(mean_matrix.shape[0], k=-1)
                        pixels = np.full((*mean_matrix.shape, 3), 255, dtype=np.uint8)
                        pixels[upper] = mean_pixels[upper]
                        pixels[lower] = std_pixels[lower]
                        matrix = mean_matrix
                        minimum, maximum = mean_minimum, mean_maximum
                        panel_mean_scale_mode = mean_scale_mode
                        panel_std_scale_mode = std_scale_mode
                    else:
                        matrix = np.asarray(panel["matrix"], dtype=float)
                        pixels, minimum, maximum = _viridis_pixels(
                            matrix, scale_mode, max_value)
                        std_minimum = std_maximum = None
                    image_data = pixels.tobytes()
                    height, width = matrix.shape
                    image = QtGui.QImage(
                        image_data, width, height, width * 3,
                        QtGui.QImage.Format_RGB888,
                    )
                    self._panels.append({
                        **panel,
                        "matrix": matrix,
                        "labels": list(panel["labels"]),
                        "image_data": image_data,
                        "image": image,
                        "minimum": minimum,
                        "maximum": maximum,
                        "std_minimum": std_minimum,
                        "std_maximum": std_maximum,
                        "scale_mode": scale_mode,
                        "mean_scale_mode": panel_mean_scale_mode if panel.get("split") else None,
                        "std_scale_mode": panel_std_scale_mode if panel.get("split") else None,
                        "show_profile": bool(show_profile and panel.get("profile") is not None),
                    })
                self.update()

            def set_message(self, message):
                self._panels = []
                self._message = message
                self.update()

            def paintEvent(self, _event):
                painter = QtGui.QPainter(self)
                painter.setRenderHint(QtGui.QPainter.Antialiasing, False)
                painter.fillRect(self.rect(), QtGui.QColor("white"))
                if not self._panels:
                    painter.setPen(QtGui.QColor(90, 98, 100))
                    painter.drawText(self.rect(), QtCore.Qt.AlignCenter,
                                     self._message)
                    painter.end()
                    return

                panel_width = self.width() / len(self._panels)
                painter.setFont(QtGui.QFont("Segoe UI", 8))
                for panel_index, panel in enumerate(self._panels):
                    panel_left = panel_index * panel_width
                    profile_visible = panel["show_profile"]
                    profile_reserve = 235 if profile_visible else 0
                    left, top, right, bottom = 66, 62, 154, 72 + profile_reserve
                    map_size = min(panel_width - left - right,
                                   self.height() - top - bottom)
                    if map_size <= 0:
                        continue
                    plot = QtCore.QRectF(panel_left + left, top, map_size, map_size)
                    painter.drawImage(plot, panel["image"])
                    painter.setPen(QtGui.QPen(QtGui.QColor(65, 74, 76), 1))
                    painter.drawRect(plot)
                    painter.setPen(QtGui.QColor(35, 44, 46))
                    painter.setFont(QtGui.QFont("Segoe UI", 10, QtGui.QFont.DemiBold))
                    painter.drawText(QtCore.QRectF(panel_left + left, 12,
                                                    panel_width - left - right, 38),
                                     QtCore.Qt.AlignCenter | QtCore.Qt.TextWordWrap,
                                     panel["title"])
                    painter.setFont(QtGui.QFont("Segoe UI", 8))

                    matrix = panel["matrix"]
                    residue_count = matrix.shape[0]
                    tick_count = min(7, residue_count)
                    ticks = np.unique(np.linspace(0, residue_count - 1,
                                                  tick_count, dtype=int))
                    painter.setPen(QtGui.QColor(55, 65, 67))
                    for index in ticks:
                        x = plot.left() + (index + 0.5) * map_size / residue_count
                        y = plot.top() + (index + 0.5) * map_size / residue_count
                        label = panel["labels"][index]
                        painter.drawLine(QtCore.QPointF(x, plot.bottom()),
                                         QtCore.QPointF(x, plot.bottom() + 4))
                        painter.drawLine(QtCore.QPointF(plot.left() - 4, y),
                                         QtCore.QPointF(plot.left(), y))
                        painter.drawText(QtCore.QRectF(x - 30, plot.bottom() + 6, 60, 28),
                                         QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop,
                                         label)
                        painter.drawText(QtCore.QRectF(panel_left + 1, y - 10,
                                                        left - 10, 20),
                                         QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter,
                                         label)

                    painter.save()
                    painter.translate(panel_left + 14, plot.center().y())
                    painter.rotate(-90)
                    painter.drawText(QtCore.QRectF(-map_size / 2, -10, map_size, 20),
                                     QtCore.Qt.AlignCenter, "Residue")
                    painter.restore()
                    x_label_y = plot.bottom() + 42 if profile_visible else self.height() - 30
                    painter.drawText(QtCore.QRectF(plot.left(), x_label_y, map_size, 20),
                                     QtCore.Qt.AlignCenter, "Residue")

                    self._paint_colorbar(painter, panel, plot, top)
                    if profile_visible:
                        self._paint_profile(painter, panel, plot)
                painter.end()

            def _paint_profile(self, painter, panel, plot):
                values = np.asarray(panel["profile"], dtype=float)
                finite = np.isfinite(values)
                chart_top = plot.bottom() + 70
                chart_bottom = self.height() - 42
                chart = QtCore.QRectF(plot.left(), chart_top, plot.width(),
                                      chart_bottom - chart_top)
                if chart.height() < 45:
                    return
                painter.setPen(QtGui.QColor(35, 44, 46))
                painter.setFont(QtGui.QFont("Segoe UI", 9, QtGui.QFont.DemiBold))
                painter.drawText(QtCore.QRectF(chart.left(), chart_top - 23,
                                                chart.width(), 20),
                                 QtCore.Qt.AlignCenter, panel["profile_title"])
                painter.setFont(QtGui.QFont("Segoe UI", 8))
                if not np.any(finite):
                    painter.drawText(chart, QtCore.Qt.AlignCenter,
                                     "No eligible residue pairs")
                    return

                maximum = float(np.nanmax(values))
                y_maximum = maximum if maximum > 0 else 1.0
                painter.setPen(QtGui.QPen(QtGui.QColor(215, 222, 221), 1))
                for fraction in (0.0, 0.5, 1.0):
                    y = chart.bottom() - fraction * chart.height()
                    painter.drawLine(QtCore.QPointF(chart.left(), y),
                                     QtCore.QPointF(chart.right(), y))
                    painter.setPen(QtGui.QColor(75, 84, 86))
                    painter.drawText(QtCore.QRectF(chart.left() - 54, y - 9, 48, 18),
                                     QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter,
                                     f"{y_maximum * fraction:.3g}")
                    painter.setPen(QtGui.QPen(QtGui.QColor(215, 222, 221), 1))

                painter.setPen(QtGui.QPen(QtGui.QColor(30, 112, 124), 2))
                path = QtGui.QPainterPath()
                drawing = False
                count = len(values)
                for index, value in enumerate(values):
                    if not np.isfinite(value):
                        drawing = False
                        continue
                    x = (chart.left() + index * chart.width() / max(1, count - 1))
                    y = chart.bottom() - (value / y_maximum) * chart.height()
                    point = QtCore.QPointF(x, y)
                    if drawing:
                        path.lineTo(point)
                    else:
                        path.moveTo(point)
                        drawing = True
                painter.drawPath(path)

                tick_count = min(6, count)
                ticks = np.unique(np.linspace(0, count - 1, tick_count, dtype=int))
                painter.setPen(QtGui.QColor(70, 79, 81))
                for index in ticks:
                    x = (chart.left() + index * chart.width() / max(1, count - 1))
                    painter.drawLine(QtCore.QPointF(x, chart.bottom()),
                                     QtCore.QPointF(x, chart.bottom() + 4))
                    painter.drawText(QtCore.QRectF(x - 34, chart.bottom() + 5, 68, 18),
                                     QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop,
                                     panel["labels"][index])

                painter.save()
                painter.translate(chart.left() - 48, chart.center().y())
                painter.rotate(-90)
                painter.drawText(QtCore.QRectF(-chart.height() / 2, -10,
                                                chart.height(), 18),
                                 QtCore.Qt.AlignCenter, panel["profile_ylabel"])
                painter.restore()

            def _paint_colorbar(self, painter, panel, plot, top):
                x = plot.right() + 14
                if panel.get("split"):
                    half_height = plot.height() / 2
                    scale_height = max(70, half_height - 86)
                    self._paint_scale(
                        painter, x, top + 8, scale_height,
                        panel["maximum"], panel["mean_scale_mode"], "Mean distance (A)",
                        panel["mean_summary"],
                    )
                    self._paint_scale(
                        painter, x, top + half_height + 12, scale_height,
                        panel["std_maximum"], panel["std_scale_mode"], "Distance SD (A)",
                        panel["std_summary"],
                    )
                else:
                    self._paint_scale(
                        painter, x, top + 8, max(70, plot.height() - 156),
                        panel["maximum"], panel["scale_mode"],
                        panel["legend"], panel["summary"],
                    )

            def _paint_scale(self, painter, x, y, bar_height, maximum,
                             scale_mode, title, summary):
                bar_width = 16
                painter.setPen(QtGui.QColor(55, 65, 67))
                scale_label = f"{title}\n{scale_mode.capitalize()} scale"
                painter.drawText(QtCore.QRectF(x - 5, y - 4, 136, 42),
                                 QtCore.Qt.AlignLeft | QtCore.Qt.TextWordWrap, scale_label)
                bar_top = y + 42
                if maximum is None:
                    painter.drawText(QtCore.QRectF(x, bar_top, 124, 24),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter,
                                     "No valid distances")
                    return

                if maximum == 0:
                    color = QtGui.QColor(*np.rint(_VIRIDIS_STOPS[0] * 255).astype(int))
                    painter.fillRect(QtCore.QRectF(x, bar_top, bar_width, bar_height), color)
                else:
                    gradient_steps = min(160, max(40, int(bar_height)))
                    for step in range(gradient_steps):
                        fraction = 1.0 - step / max(1, gradient_steps - 1)
                        palette_position = fraction * (len(_VIRIDIS_STOPS) - 1)
                        lower = int(np.floor(palette_position))
                        upper = min(lower + 1, len(_VIRIDIS_STOPS) - 1)
                        blend = palette_position - lower
                        color_values = (_VIRIDIS_STOPS[lower] * (1 - blend)
                                        + _VIRIDIS_STOPS[upper] * blend)
                        color = QtGui.QColor(*(np.rint(color_values * 255).astype(int)))
                        y_position = bar_top + step * bar_height / gradient_steps
                        painter.fillRect(QtCore.QRectF(x, y_position, bar_width,
                                                       bar_height / gradient_steps + 1), color)
                painter.setPen(QtGui.QPen(QtGui.QColor(65, 74, 76), 1))
                painter.drawRect(QtCore.QRectF(x, bar_top, bar_width, bar_height))
                painter.setPen(QtGui.QColor(55, 65, 67))
                ticks = _scale_ticks(maximum, scale_mode)
                for value in ticks:
                    if maximum > 0 and scale_mode == "log":
                        position = np.log1p(value) / np.log1p(maximum)
                    else:
                        position = value / maximum if maximum > 0 else 0.0
                    y_position = bar_top + (1.0 - position) * bar_height
                    painter.drawLine(QtCore.QPointF(x + bar_width, y_position),
                                     QtCore.QPointF(x + bar_width + 4, y_position))
                    painter.drawText(QtCore.QRectF(x + bar_width + 7, y_position - 9,
                                                    67, 18),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter,
                                     _format_scale_value(value))
                painter.setPen(QtGui.QColor(45, 54, 56))
                painter.drawText(QtCore.QRectF(x - 3, bar_top + bar_height + 8,
                                                138, 62),
                                 QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop
                                 | QtCore.Qt.TextWordWrap,
                                 "\n".join(summary))

        return CanvasWidget(*args, **kwargs)


class ContactMapWindow:
    """Build a separate PyMOL Qt window that tracks the current frame."""

    def __new__(cls):
        QtCore, QtGui, QtWidgets = _load_qt()
        from pymol import cmd

        class Window(QtWidgets.QDialog):
            def __init__(self):
                super().__init__(None, QtCore.Qt.Window)
                self.setWindowTitle("Protein Contact Map")
                self.resize(1180, 860)
                self._last_key = None
                self._aggregate_cache = {}

                layout = QtWidgets.QVBoxLayout(self)
                controls = QtWidgets.QHBoxLayout()
                controls.addWidget(QtWidgets.QLabel("Object", self))
                self.object_combo = QtWidgets.QComboBox(self)
                self.object_combo.setMinimumWidth(190)
                controls.addWidget(self.object_combo)

                controls.addWidget(QtWidgets.QLabel("Cutoff", self))
                self.cutoff_spin = QtWidgets.QDoubleSpinBox(self)
                self.cutoff_spin.setRange(1.0, 30.0)
                self.cutoff_spin.setSingleStep(0.5)
                self.cutoff_spin.setValue(8.0)
                self.cutoff_spin.setSuffix(" A")
                controls.addWidget(self.cutoff_spin)

                controls.addWidget(QtWidgets.QLabel("Exclude same-chain neighbors <=", self))
                self.separation_spin = QtWidgets.QSpinBox(self)
                self.separation_spin.setRange(0, 20)
                self.separation_spin.setValue(3)
                self.separation_spin.setSuffix(" residues")
                controls.addWidget(self.separation_spin)

                self.refresh_button = QtWidgets.QPushButton("Refresh objects", self)
                controls.addWidget(self.refresh_button)
                self.save_button = QtWidgets.QPushButton("Save plot…", self)
                controls.addWidget(self.save_button)
                layout.addLayout(controls)

                view_controls = QtWidgets.QHBoxLayout()
                view_controls.addWidget(QtWidgets.QLabel("Plot", self))
                self.plot_mode = QtWidgets.QComboBox(self)
                self.plot_mode.addItem("Current frame cutoff", "contacts")
                self.plot_mode.addItem("Current frame residue-center distances", "frame_centers")
                self.plot_mode.addItem("Mean distance", "mean")
                self.plot_mode.addItem("Distance standard deviation", "std")
                self.plot_mode.addItem("Mean + standard deviation", "both")
                self.plot_mode.addItem("Split mean / SD across diagonal", "split")
                view_controls.addWidget(self.plot_mode)
                view_controls.addSpacing(18)
                self.global_scale_label = QtWidgets.QLabel("Color scale", self)
                view_controls.addWidget(self.global_scale_label)
                self.log_scale_button = QtWidgets.QToolButton(self)
                self.log_scale_button.setText("Log")
                self.log_scale_button.setCheckable(True)
                self.log_scale_button.setChecked(True)
                self.log_scale_button.setFixedSize(48, 28)
                self.linear_scale_button = QtWidgets.QToolButton(self)
                self.linear_scale_button.setText("Linear")
                self.linear_scale_button.setCheckable(True)
                self.linear_scale_button.setFixedSize(58, 28)
                self.scale_button_group = QtWidgets.QButtonGroup(self)
                self.scale_button_group.setExclusive(True)
                self.scale_button_group.addButton(self.log_scale_button)
                self.scale_button_group.addButton(self.linear_scale_button)
                view_controls.addWidget(self.log_scale_button)
                view_controls.addWidget(self.linear_scale_button)
                self.auto_max_check = QtWidgets.QCheckBox("Auto max", self)
                self.auto_max_check.setChecked(True)
                view_controls.addWidget(self.auto_max_check)
                self.global_max_label = QtWidgets.QLabel("Max (A)", self)
                view_controls.addWidget(self.global_max_label)
                self.max_range_spin = QtWidgets.QDoubleSpinBox(self)
                self.max_range_spin.setRange(0.1, 10000.0)
                self.max_range_spin.setSingleStep(1.0)
                self.max_range_spin.setDecimals(2)
                self.max_range_spin.setValue(40.0)
                self.max_range_spin.setEnabled(False)
                self.max_range_spin.setFixedWidth(86)
                view_controls.addWidget(self.max_range_spin)
                view_controls.addStretch(1)
                layout.addLayout(view_controls)

                profile_controls = QtWidgets.QHBoxLayout()
                self.profile_check = QtWidgets.QCheckBox(
                    "Show per-residue column average below heatmap", self)
                self.profile_check.setEnabled(False)
                profile_controls.addWidget(self.profile_check)
                self.export_csv_button = QtWidgets.QPushButton("Export raw CSV…", self)
                self.export_csv_button.setEnabled(False)
                profile_controls.addWidget(self.export_csv_button)
                profile_controls.addStretch(1)
                layout.addLayout(profile_controls)

                self.split_scale_row = QtWidgets.QWidget(self)
                split_scale_layout = QtWidgets.QHBoxLayout(self.split_scale_row)
                split_scale_layout.setContentsMargins(0, 0, 0, 0)
                self.mean_split_controls = self._create_split_scale_controls(
                    "Mean distance (upper triangle)")
                self.std_split_controls = self._create_split_scale_controls(
                    "Standard deviation (lower triangle)")
                split_scale_layout.addWidget(self.mean_split_controls["group"])
                split_scale_layout.addWidget(self.std_split_controls["group"])
                self.split_scale_row.hide()
                layout.addWidget(self.split_scale_row)

                self.frame_label = QtWidgets.QLabel("Frame: -", self)
                layout.addWidget(self.frame_label)
                self.canvas = ContactMapCanvas(self)
                layout.addWidget(self.canvas, 1)
                self.status_label = QtWidgets.QLabel(
                    "Cutoff contacts use CA atoms; frame distances use residue heavy-atom centers.", self)
                layout.addWidget(self.status_label)

                self.object_combo.currentTextChanged.connect(self._selection_changed)
                self.cutoff_spin.valueChanged.connect(self._settings_changed)
                self.separation_spin.valueChanged.connect(self._settings_changed)
                self.refresh_button.clicked.connect(self.refresh_objects)
                self.plot_mode.currentIndexChanged.connect(self._settings_changed)
                self.log_scale_button.toggled.connect(self._scale_button_changed)
                self.linear_scale_button.toggled.connect(self._scale_button_changed)
                self.auto_max_check.toggled.connect(self._toggle_auto_max)
                self.max_range_spin.valueChanged.connect(self._settings_changed)
                self.profile_check.toggled.connect(self._settings_changed)
                self.export_csv_button.clicked.connect(self.export_raw_csv)
                for scale_controls in (self.mean_split_controls, self.std_split_controls):
                    scale_controls["log"].toggled.connect(self._split_scale_changed)
                    scale_controls["linear"].toggled.connect(self._split_scale_changed)
                    scale_controls["auto"].toggled.connect(
                        lambda checked, controls=scale_controls: self._toggle_split_auto(
                            controls, checked))
                    scale_controls["maximum"].valueChanged.connect(self._settings_changed)
                self.save_button.clicked.connect(self.save_plot)

                self.timer = QtCore.QTimer(self)
                self.timer.setInterval(350)
                self.timer.timeout.connect(self.refresh_if_needed)
                self.refresh_objects()
                self.timer.start()

            def _create_split_scale_controls(self, title):
                group = QtWidgets.QGroupBox(title, self)
                controls = QtWidgets.QHBoxLayout(group)
                log_button = QtWidgets.QRadioButton("Log", group)
                log_button.setChecked(True)
                linear_button = QtWidgets.QRadioButton("Linear", group)
                scale_group = QtWidgets.QButtonGroup(group)
                scale_group.setExclusive(True)
                scale_group.addButton(log_button)
                scale_group.addButton(linear_button)
                auto_check = QtWidgets.QCheckBox("Auto max", group)
                auto_check.setChecked(True)
                maximum_spin = QtWidgets.QDoubleSpinBox(group)
                maximum_spin.setRange(0.1, 10000.0)
                maximum_spin.setSingleStep(1.0)
                maximum_spin.setDecimals(2)
                maximum_spin.setValue(40.0)
                maximum_spin.setSuffix(" A")
                maximum_spin.setEnabled(False)
                maximum_spin.setFixedWidth(90)
                controls.addWidget(log_button)
                controls.addWidget(linear_button)
                controls.addWidget(auto_check)
                controls.addWidget(maximum_spin)
                return {
                    "group": group,
                    "log": log_button,
                    "linear": linear_button,
                    "auto": auto_check,
                    "maximum": maximum_spin,
                }

            def refresh_objects(self):
                previous = self.object_combo.currentText()
                names = []
                for object_name in cmd.get_names("objects"):
                    try:
                        if cmd.get_type(object_name) == "object:molecule":
                            names.append(object_name)
                    except Exception:
                        continue
                self.object_combo.blockSignals(True)
                self.object_combo.clear()
                self.object_combo.addItems(names)
                if previous in names:
                    self.object_combo.setCurrentText(previous)
                self.object_combo.blockSignals(False)
                self._aggregate_cache.clear()
                self._last_key = None
                self.update_contact_map()

            def _selection_changed(self, _text):
                self._last_key = None
                self.update_contact_map()

            def _settings_changed(self, _value):
                self._last_key = None
                self.update_contact_map()

            def _scale_button_changed(self, checked):
                if checked:
                    self._last_key = None
                    self.update_contact_map()

            def _split_scale_changed(self, checked):
                if checked:
                    self._last_key = None
                    self.update_contact_map()

            def _toggle_auto_max(self, auto_max):
                self.max_range_spin.setEnabled(not auto_max)
                self._last_key = None
                self.update_contact_map()

            def _toggle_split_auto(self, controls, auto_max):
                controls["maximum"].setEnabled(not auto_max)
                self._last_key = None
                self.update_contact_map()

            def _scale_settings(self):
                scale_mode = "log" if self.log_scale_button.isChecked() else "linear"
                if self.plot_mode.currentData() == "contacts":
                    return scale_mode, 1.0
                max_value = None if self.auto_max_check.isChecked() else self.max_range_spin.value()
                return scale_mode, max_value

            def _split_scale_settings(self):
                settings = {}
                for key, controls in (("mean", self.mean_split_controls),
                                      ("std", self.std_split_controls)):
                    mode = "log" if controls["log"].isChecked() else "linear"
                    maximum = (None if controls["auto"].isChecked()
                               else controls["maximum"].value())
                    settings[key] = (mode, maximum)
                return settings

            def _update_control_states(self, mode):
                contact_mode = mode == "contacts"
                split_mode = mode == "split"
                for widget in (self.global_scale_label, self.log_scale_button,
                               self.linear_scale_button, self.auto_max_check,
                               self.global_max_label, self.max_range_spin):
                    widget.setVisible(not split_mode)
                self.split_scale_row.setVisible(split_mode)
                self.cutoff_spin.setEnabled(contact_mode)
                distance_mode = not contact_mode
                self.global_scale_label.setEnabled(not split_mode)
                self.auto_max_check.setEnabled(distance_mode and not split_mode)
                self.max_range_spin.setEnabled(
                    distance_mode and not split_mode and not self.auto_max_check.isChecked())
                self.profile_check.setEnabled(mode in {"mean", "std"})
                self.export_csv_button.setEnabled(mode in {"mean", "std"})

            def _frame_data(self, object_name, state):
                model = cmd.get_model(
                    f"({object_name}) and polymer.protein and name CA", state=state)
                atoms = list(model.atom)
                if not atoms:
                    raise ValueError(f"No protein CA atoms found in frame {state}.")
                coordinates = np.asarray([atom.coord for atom in atoms], dtype=float)
                chains = [atom.chain or "" for atom in atoms]
                labels = [f"{atom.chain or '-'}:{atom.resn}{atom.resi}" for atom in atoms]
                return coordinates, chains, labels

            def _frame_residue_center_data(self, object_name, state):
                model = cmd.get_model(
                    f"({object_name}) and polymer.protein and not hydro", state=state)
                residue_atoms = {}
                for atom in model.atom:
                    key = (atom.chain or "", atom.resi, atom.resn)
                    residue_atoms.setdefault(key, []).append(atom.coord)
                if not residue_atoms:
                    raise ValueError(f"No protein heavy atoms found in frame {state}.")
                keys = list(residue_atoms)
                centers = np.asarray([
                    np.mean(residue_atoms[key], axis=0) for key in keys
                ], dtype=float)
                chains = [key[0] for key in keys]
                labels = [f"{chain or '-'}:{resn}{resi}"
                          for chain, resi, resn in keys]
                return centers, chains, labels

            def _aggregate_maps(self, object_name, state_count, separation):
                cache_key = (object_name, state_count, separation)
                if cache_key in self._aggregate_cache:
                    return self._aggregate_cache[cache_key]

                mean = None
                sum_squared_deviations = None
                expected_labels = None
                valid_pairs = None
                total_pair_distance = 0.0
                total_pair_distance_squared = 0.0
                total_pair_observations = 0
                progress_step = max(1, state_count // 20)
                for state in range(1, state_count + 1):
                    coordinates, chains, labels = self._frame_data(object_name, state)
                    if expected_labels is None:
                        expected_labels = labels
                        mean = np.zeros((len(labels), len(labels)), dtype=np.float64)
                        sum_squared_deviations = np.zeros_like(mean)
                        excluded = _excluded_pairs(chains, separation)
                        np.fill_diagonal(excluded, False)
                        valid_pairs = ~excluded
                        if not np.any(np.triu(valid_pairs, k=1)):
                            raise ValueError(
                                "No residue pairs remain after the same-chain neighbor filter."
                            )
                    elif labels != expected_labels:
                        raise ValueError(
                            f"CA residue selection changes in frame {state}; "
                            "distance statistics require matching residues in every frame."
                        )

                    distances = _distance_matrix(coordinates)
                    delta = np.zeros_like(mean)
                    delta[valid_pairs] = distances[valid_pairs] - mean[valid_pairs]
                    mean[valid_pairs] += delta[valid_pairs] / state
                    sum_squared_deviations[valid_pairs] += (
                        delta[valid_pairs] * (distances[valid_pairs] - mean[valid_pairs])
                    )

                    upper_valid = np.triu(valid_pairs, k=1)
                    pair_distances = distances[upper_valid]
                    total_pair_distance += float(pair_distances.sum())
                    total_pair_distance_squared += float(np.square(pair_distances).sum())
                    total_pair_observations += int(pair_distances.size)
                    if state % progress_step == 0 or state == state_count:
                        self.status_label.setText(
                            f"Calculating distance statistics: frame {state}/{state_count}…")
                        QtWidgets.QApplication.processEvents()

                standard_deviation = np.zeros_like(sum_squared_deviations)
                standard_deviation[valid_pairs] = np.sqrt(
                    sum_squared_deviations[valid_pairs] / state_count
                )
                mean[~valid_pairs] = np.nan
                standard_deviation[~valid_pairs] = np.nan
                overall_mean = (total_pair_distance / total_pair_observations
                                if total_pair_observations else 0.0)
                overall_variance = (
                    total_pair_distance_squared / total_pair_observations - overall_mean ** 2
                    if total_pair_observations else 0.0
                )
                overall_standard_deviation = float(np.sqrt(max(0.0, overall_variance)))
                result = (mean, standard_deviation, expected_labels,
                          overall_mean, overall_standard_deviation)
                self._aggregate_cache[cache_key] = result
                return result

            def refresh_if_needed(self):
                object_name = self.object_combo.currentText()
                if not object_name:
                    return
                frame = max(1, int(cmd.get_frame()))
                mode = self.plot_mode.currentData()
                state_count = max(1, int(cmd.count_states(object_name)))
                frame_key = frame if mode in {"contacts", "frame_centers"} else None
                scale_mode, max_value = self._scale_settings()
                split_scale_settings = self._split_scale_settings() if mode == "split" else None
                profile_key = self.profile_check.isChecked() if mode in {"mean", "std"} else False
                key = (object_name, frame_key, state_count, self.cutoff_spin.value(),
                      self.separation_spin.value(), mode, scale_mode, max_value,
                      split_scale_settings, profile_key)
                if key != self._last_key:
                    self.update_contact_map()

            def closeEvent(self, event):
                self.timer.stop()
                super().closeEvent(event)

            def update_contact_map(self):
                object_name = self.object_combo.currentText()
                if not object_name:
                    self.frame_label.setText("Frame: -")
                    self.status_label.setText("Load a protein object in PyMOL to begin.")
                    self.canvas.set_message("No molecular objects available.")
                    return

                try:
                    frame = max(1, int(cmd.get_frame()))
                    state_count = max(1, int(cmd.count_states(object_name)))
                    state = min(frame, state_count)
                    cutoff = self.cutoff_spin.value()
                    separation = self.separation_spin.value()
                    mode = self.plot_mode.currentData()
                    self._update_control_states(mode)
                    scale_mode, max_value = self._scale_settings()
                    if mode == "contacts":
                        coordinates, chains, labels = self._frame_data(object_name, state)
                        matrix = _contact_matrix(coordinates, chains, cutoff, separation).astype(float)
                        pair_values = matrix[np.triu_indices(len(labels), k=1)]
                        contact_count = int(pair_values.sum())
                        panels = [{
                            "matrix": matrix,
                            "labels": labels,
                            "title": f"Contacts in frame {state}",
                            "legend": "Contact (log scale)",
                            "summary": [f"Frame: {state}/{state_count}",
                                        f"Contacts: {contact_count}"],
                            "has_zero": bool(np.any(matrix == 0)),
                        }]
                        self.frame_label.setText(f"Frame: {state}/{state_count}")
                        status = (f"{len(labels)} CA residues | {contact_count} contacts | "
                                  f"cutoff {cutoff:.1f} A")
                    elif mode == "frame_centers":
                        centers, chains, labels = self._frame_residue_center_data(
                            object_name, state)
                        matrix = _distance_matrix(centers)
                        excluded = _excluded_pairs(chains, separation)
                        np.fill_diagonal(excluded, False)
                        matrix[excluded] = np.nan
                        pair_distances = matrix[np.triu_indices(len(labels), k=1)]
                        pair_distances = pair_distances[np.isfinite(pair_distances)]
                        if not pair_distances.size:
                            raise ValueError(
                                "No residue-center pairs remain after the neighbor filter."
                            )
                        panels = [{
                            "matrix": matrix,
                            "labels": labels,
                            "title": f"Residue-center distances in frame {state}",
                            "legend": "Center distance (A)",
                            "summary": [
                                f"Frame: {state}/{state_count}",
                                f"Mean pair distance: {np.mean(pair_distances):.3f} A",
                                f"Maximum distance: {np.max(pair_distances):.3f} A",
                            ],
                        }]
                        self.frame_label.setText(f"Frame: {state}/{state_count}")
                        status = (f"{len(labels)} residue centers | "
                                  f"mean pair distance {np.mean(pair_distances):.3f} A")
                    else:
                        (mean, standard_deviation, labels, overall_mean,
                         overall_standard_deviation) = self._aggregate_maps(
                            object_name, state_count, separation)
                        upper_indices = np.triu_indices(len(labels), k=1)
                        mean_values = mean[upper_indices]
                        sd_values = standard_deviation[upper_indices]
                        mean_summary = [
                            f"Frames: {state_count}",
                            f"Mean distance: {overall_mean:.3f} A",
                            f"Distance SD: {overall_standard_deviation:.3f} A",
                        ]
                        sd_summary = [
                            f"Frames: {state_count}",
                            f"Mean pair SD: {np.nanmean(sd_values):.3f} A",
                            f"Max pair SD: {np.nanmax(sd_values):.3f} A",
                        ]
                        mean_panel = {
                            "matrix": mean,
                            "labels": labels,
                            "title": "Mean CA-pair distance",
                            "legend": "Mean distance (A)",
                            "summary": mean_summary,
                            "profile": _column_profile(mean),
                            "profile_title": "Mean distance by residue",
                            "profile_ylabel": "Mean CA distance (A)",
                        }
                        sd_panel = {
                            "matrix": standard_deviation,
                            "labels": labels,
                            "title": "Distance standard deviation",
                            "legend": "Distance SD (A)",
                            "summary": sd_summary,
                            "profile": _column_profile(standard_deviation),
                            "profile_title": "Mean pair SD by residue",
                            "profile_ylabel": "Mean pair SD (A)",
                        }
                        if mode == "mean":
                            panels = [mean_panel]
                            status = (f"Mean CA-pair distance across {state_count} frames | "
                                      f"overall mean {overall_mean:.3f} A; "
                                      f"overall SD {overall_standard_deviation:.3f} A")
                        elif mode == "std":
                            panels = [sd_panel]
                            status = (f"Population SD of CA-pair distances across {state_count} frames | "
                                      f"mean pair SD {np.nanmean(sd_values):.3f} A")
                        elif mode == "split":
                            panels = [{
                                "split": True,
                                "matrix": mean,
                                "mean_matrix": mean,
                                "std_matrix": standard_deviation,
                                "labels": labels,
                                "title": "Mean distance (upper) / SD (lower)",
                                "mean_summary": mean_summary,
                                "std_summary": sd_summary,
                            }]
                            status = ("Upper triangle: mean CA-pair distance; "
                                      "lower triangle: distance SD")
                        else:
                            panels = [mean_panel, sd_panel]
                            status = f"Mean distance and population SD across {state_count} frames"
                        self.frame_label.setText(f"Statistics: {state_count} frames")

                    split_scale_settings = (
                        self._split_scale_settings() if mode == "split" else None
                    )
                    show_profile = (
                        self.profile_check.isChecked() and mode in {"mean", "std"}
                    )
                    self.canvas.set_panels(
                        panels, scale_mode, max_value,
                        split_scale_settings=split_scale_settings,
                        show_profile=show_profile,
                    )
                    self.status_label.setText(status)
                    frame_key = frame if mode in {"contacts", "frame_centers"} else None
                    self._last_key = (object_name, frame_key, state_count, cutoff,
                                      separation, mode, scale_mode, max_value,
                                      split_scale_settings,
                                      self.profile_check.isChecked()
                                      if mode in {"mean", "std"} else False)
                except Exception as exc:
                    self.status_label.setText(f"Could not update contact map: {exc}")
                    self.canvas.set_message("Contact map update failed; see status below.")
                    self._last_key = None

            def save_plot(self):
                if not self.canvas._panels:
                    self.status_label.setText("There is no contact map to save yet.")
                    return
                object_name = self.object_combo.currentText() or "contact_map"
                safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", object_name)
                mode = self.plot_mode.currentData()
                if mode in {"contacts", "frame_centers"}:
                    suffix = f"frame_{max(1, int(cmd.get_frame())):04d}"
                else:
                    suffix = mode
                default_path = f"{safe_name}_{suffix}_contact_map.png"
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self, "Save contact-map plot", default_path,
                    "PNG image (*.png);;JPEG image (*.jpg *.jpeg)",
                )
                if not path:
                    return
                if not Path(path).suffix:
                    path += ".png"
                scale = 2
                image = QtGui.QImage(
                    self.canvas.width() * scale,
                    self.canvas.height() * scale,
                    QtGui.QImage.Format_ARGB32,
                )
                image.setDevicePixelRatio(scale)
                image.fill(QtGui.QColor("white"))
                painter = QtGui.QPainter(image)
                self.canvas.render(painter)
                painter.end()
                if not image.save(path):
                    self.status_label.setText(f"Could not save plot: {path}")
                    return
                self.status_label.setText(f"Saved plot: {path}")

            def export_raw_csv(self):
                mode = self.plot_mode.currentData()
                if mode not in {"mean", "std"} or not self.canvas._panels:
                    self.status_label.setText(
                        "Raw matrix CSV export is available for Mean distance and "
                        "Distance standard deviation."
                    )
                    return
                object_name = self.object_combo.currentText() or "protein"
                suffix = "mean_distance" if mode == "mean" else "distance_standard_deviation"
                default_path = f"{re.sub(r'[^A-Za-z0-9_-]', '_', object_name)}_{suffix}.csv"
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self, "Export raw distance matrix", default_path,
                    "CSV file (*.csv)",
                )
                if not path:
                    return
                if not Path(path).suffix:
                    path += ".csv"
                try:
                    panel = self.canvas._panels[0]
                    _write_matrix_csv(path, panel["matrix"], panel["labels"])
                except Exception as exc:
                    self.status_label.setText(f"Could not export CSV: {exc}")
                    return
                self.status_label.setText(f"Saved raw matrix CSV: {path}")

        return Window()


def contact_map_gui():
    """Open or raise the live contact-map window."""
    global _WINDOWS
    _WINDOWS = [window for window in _WINDOWS if window.isVisible()]
    if _WINDOWS:
        _WINDOWS[0].raise_()
        _WINDOWS[0].activateWindow()
        return _WINDOWS[0]
    window = ContactMapWindow()
    _WINDOWS.append(window)
    window.show()
    return window


def _register():
    try:
        from pymol import cmd
    except ImportError:
        return
    cmd.extend("contact_map_gui", contact_map_gui)


def __init_plugin__(app=None):
    _register()
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("Contact Map", contact_map_gui)


_register()