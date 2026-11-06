"""Contact Map: display a live CA-residue contact map for a PyMOL object.

Version: v6
Purpose:
    Open a separate Qt window showing current-frame cutoff contacts or
    residue-center distances, plus mean pairwise distances and standard
    deviations across frames with CA or all-atom residue-center distances.
    Statistics share one configurable view with on-colorbar range
    controls, linked Mean-to-SD pair filtering, and range-filtered profiles.
    One Export dialog handles current plots, selected CSVs, and frame series.
    Split profiles appear below (SD) and
    to the right (mean). Every plot mode offers destination-selectable CSV
    export of its matrices and, for statistics, its residue profiles. Optional
    frame-series export saves one PNG and/or matching CSV per selected state.

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
import tempfile
import shutil

import numpy as np

__version__ = "v5"
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


def _column_profile(matrix, minimum=0.0, maximum=None):
    """Average eligible off-diagonal values for each residue column."""
    values = np.asarray(matrix, dtype=float).copy()
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("A per-residue profile requires a square matrix.")
    np.fill_diagonal(values, np.nan)
    values[values < minimum] = np.nan
    if maximum is not None:
        values[values > maximum] = np.nan
    valid_counts = np.sum(np.isfinite(values), axis=0)
    totals = np.nansum(values, axis=0)
    means = np.full(values.shape[1], np.nan, dtype=float)
    np.divide(totals, valid_counts, out=means, where=valid_counts > 0)
    return means


def _residue_distances(atoms, metric):
    """Return distances, labels, and stable identities for one protein state."""
    residues = {}
    for atom in atoms:
        if metric == "ca" and atom.name != "CA":
            continue
        key = (atom.segi, atom.chain, atom.resi, atom.resn)
        residues.setdefault(key, []).append(atom.coord)
    keys = list(residues)
    if not keys:
        raise ValueError("No protein atoms are available for this distance definition.")
    coordinates = [np.asarray(residues[key], dtype=float) for key in keys]
    if metric in {"ca", "centers"}:
        matrix = _distance_matrix([points.mean(axis=0) for points in coordinates])
    elif metric == "closest":
        matrix = np.zeros((len(keys), len(keys)), dtype=float)
        for first in range(len(keys)):
            for second in range(first + 1, len(keys)):
                squared_minimum = np.inf
                for block_start in range(0, len(coordinates[first]), 64):
                    block = coordinates[first][block_start:block_start + 64]
                    delta = block[:, None, :] - coordinates[second][None, :, :]
                    squared_minimum = min(squared_minimum, float(np.square(delta).sum(axis=2).min()))
                matrix[first, second] = matrix[second, first] = np.sqrt(squared_minimum)
    else:
        raise ValueError("Unknown residue distance definition.")
    chains = [(segment, chain) for segment, chain, _, _ in keys]
    labels = [f"{segment + '/' if segment else ''}{chain or '-'}:{name}{residue}"
              for segment, chain, residue, name in keys]
    return matrix, chains, labels, keys


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


def _write_profile_csv(path, values, labels, value_label):
    """Export each residue's plotted column average, retaining gaps as NA."""
    values = np.asarray(values, dtype=float)
    if values.shape != (len(labels),):
        raise ValueError("The profile and residue labels do not match.")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["residue", value_label])
        for label, value in zip(labels, values):
            writer.writerow([label, f"{value:.12g}" if np.isfinite(value) else "NA"])


def _excluded_pairs(chains, minimum_sequence_separation):
    """Mark diagonal and nearby same-chain residue pairs as unavailable."""
    chain_values = list(chains)
    chains = np.empty(len(chain_values), dtype=object)
    chains[:] = chain_values
    positions = np.arange(len(chains))
    nearby_in_same_chain = (
        (np.abs(positions[:, np.newaxis] - positions[np.newaxis, :])
         <= int(minimum_sequence_separation))
        & (chains[:, np.newaxis] == chains[np.newaxis, :])
    )
    return nearby_in_same_chain | np.eye(len(chains), dtype=bool)


def _viridis_pixels(matrix, scale_mode="log", max_value=None, min_value=0.0):
    """Map values from zero to a selected maximum using log1p or linear scale."""
    values = np.asarray(matrix, dtype=float)
    pixels = np.full((*values.shape, 3), 255, dtype=np.uint8)
    finite = np.isfinite(values)
    if not np.any(finite):
        return pixels, None, None

    maximum = max(0.0, float(values[finite].max())) if max_value is None else float(max_value)
    minimum = float(min_value)
    if minimum < 0 or maximum < minimum:
        raise ValueError("The range must satisfy 0 <= minimum <= maximum.")
    clipped = np.clip(values[finite], minimum, maximum)
    normalized = np.zeros(values.shape, dtype=float)
    if maximum == minimum:
        normalized[finite] = 0.0
    elif scale_mode == "log":
        normalized[finite] = ((np.log1p(clipped) - np.log1p(minimum))
                              / (np.log1p(maximum) - np.log1p(minimum)))
    elif scale_mode == "linear":
        normalized[finite] = (clipped - minimum) / (maximum - minimum)
    else:
        raise ValueError("Scale mode must be 'log' or 'linear'.")

    palette_position = np.clip(normalized, 0.0, 1.0) * (len(_VIRIDIS_STOPS) - 1)
    lower = np.floor(palette_position).astype(int)
    upper = np.minimum(lower + 1, len(_VIRIDIS_STOPS) - 1)
    fraction = (palette_position - lower)[..., np.newaxis]
    colors = _VIRIDIS_STOPS[lower] * (1.0 - fraction) + _VIRIDIS_STOPS[upper] * fraction
    pixels[finite] = np.rint(colors[finite] * 255).astype(np.uint8)
    pixels[finite & (values < minimum)] = (80, 80, 80)
    pixels[finite & (values > maximum)] = (210, 210, 210)
    return pixels, minimum, maximum


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


def _statistics_display(mean, standard_deviation, settings):
    """Share Mean pair-selection limits with SD colors and profile averages."""
    mean_mode, mean_minimum, mean_maximum = settings["mean"]
    sd_mode, sd_minimum, sd_maximum = settings["std"]
    mean_pixels, mean_minimum, mean_maximum = _viridis_pixels(
        mean, mean_mode, mean_maximum, mean_minimum)
    sd_pixels, sd_minimum, sd_maximum = _viridis_pixels(
        standard_deviation, sd_mode, sd_maximum, sd_minimum)
    finite_mean = np.isfinite(mean)
    below = finite_mean & (mean < mean_minimum)
    above = finite_mean & (mean > mean_maximum)
    finite_sd = np.isfinite(standard_deviation)
    sd_pixels[finite_sd & below] = (80, 80, 80)
    sd_pixels[finite_sd & above] = (210, 210, 210)
    filtered_sd = np.asarray(standard_deviation, dtype=float).copy()
    filtered_sd[below | above | ~finite_mean] = np.nan
    return {
        "mean": {"pixels": mean_pixels, "minimum": mean_minimum,
                 "maximum": mean_maximum, "scale_mode": mean_mode,
                 "profile": _column_profile(mean, mean_minimum, mean_maximum)},
        "std": {"pixels": sd_pixels, "minimum": sd_minimum,
                "maximum": sd_maximum, "scale_mode": sd_mode,
                "profile": _column_profile(filtered_sd)},
    }


def _format_scale_value(value):
    return f"{value:.0e}" if 0 < value < 0.01 else f"{value:.3g}"


def _load_qt():
    from pymol.Qt import QtCore, QtGui, QtWidgets
    return QtCore, QtGui, QtWidgets


class ContactMapRangeSlider:
    def __new__(cls, *args, **kwargs):
        QtCore, QtGui, QtWidgets = _load_qt()

        class RangeSlider(QtWidgets.QWidget):
            def __init__(self, parent=None):
                super().__init__(parent)
                self.minimum = 0.0
                self.maximum = 1.0
                self.domain = 1.0
                self.scale_mode = "linear"
                self.active_handle = None
                self.on_change = None
                self.setCursor(QtCore.Qt.SizeVerCursor)
                self.setToolTip("Drag the upper or lower cutoff handle along the colorbar")

            def set_range(self, minimum, maximum, domain, scale_mode):
                self.minimum, self.maximum = minimum, maximum
                self.domain = max(domain, maximum, minimum, 0.01)
                self.scale_mode = scale_mode
                self.update()

            def position(self, value):
                fraction = (np.log1p(value) / np.log1p(self.domain)
                            if self.scale_mode == "log" else value / self.domain)
                return (1.0 - fraction) * max(1, self.height() - 1)

            def paintEvent(self, event):
                painter = QtGui.QPainter(self)
                for value in (self.minimum, self.maximum):
                    position = self.position(value)
                    painter.setPen(QtGui.QPen(QtGui.QColor("white"), 5))
                    painter.drawLine(QtCore.QPointF(0, position), QtCore.QPointF(self.width(), position))
                    painter.setPen(QtGui.QPen(QtGui.QColor(30, 40, 42), 2))
                    painter.drawLine(QtCore.QPointF(0, position), QtCore.QPointF(self.width(), position))
                painter.end()

            def mousePressEvent(self, event):
                if event.button() == QtCore.Qt.LeftButton:
                    self.active_handle = min(("maximum", "minimum"),
                        key=lambda name: abs(event.pos().y() - self.position(getattr(self, name))))
                    self._drag(event.pos().y())

            def mouseMoveEvent(self, event):
                if self.active_handle:
                    self._drag(event.pos().y())

            def mouseReleaseEvent(self, event):
                if self.active_handle:
                    self._drag(event.pos().y())
                    self.active_handle = None

            def _drag(self, position):
                fraction = float(np.clip(1.0 - position / max(1, self.height() - 1), 0, 1))
                value = (float(np.expm1(fraction * np.log1p(self.domain)))
                         if self.scale_mode == "log" else fraction * self.domain)
                value = min(value, self.maximum) if self.active_handle == "minimum" else max(value, self.minimum)
                value = round(value, 2)
                if self.active_handle == "minimum":
                    value = min(value, float(np.floor(self.maximum * 100) / 100))
                if self.on_change:
                    self.on_change(self.active_handle, value)

        return RangeSlider(*args, **kwargs)


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
                self.scale_controls = {}
                self.export_render = False

            def set_panels(self, panels, scale_mode, max_value,
                           split_scale_settings=None, show_profile=False):
                rendered_panels = []
                for panel in panels:
                    if panel.get("split"):
                        mean_matrix = np.asarray(panel["mean_matrix"], dtype=float)
                        standard_deviation = np.asarray(panel["std_matrix"], dtype=float)
                        if split_scale_settings is None:
                            mean_scale_mode, mean_min_value, mean_max_value = scale_mode, 0.0, max_value
                            std_scale_mode, std_min_value, std_max_value = scale_mode, 0.0, max_value
                        else:
                            mean_scale_mode, mean_min_value, mean_max_value = split_scale_settings["mean"]
                            std_scale_mode, std_min_value, std_max_value = split_scale_settings["std"]
                        display = _statistics_display(mean_matrix, standard_deviation, {
                            "mean": (mean_scale_mode, mean_min_value, mean_max_value),
                            "std": (std_scale_mode, std_min_value, std_max_value),
                        })
                        mean_pixels = display["mean"]["pixels"]
                        mean_minimum, mean_maximum = display["mean"]["minimum"], display["mean"]["maximum"]
                        std_pixels = display["std"]["pixels"]
                        std_minimum, std_maximum = display["std"]["minimum"], display["std"]["maximum"]
                        upper = np.triu_indices(mean_matrix.shape[0], k=0)
                        lower = np.tril_indices(mean_matrix.shape[0], k=-1)
                        pixels = np.full((*mean_matrix.shape, 3), 255, dtype=np.uint8)
                        pixels[upper] = mean_pixels[upper]
                        pixels[lower] = std_pixels[lower]
                        matrix = mean_matrix
                        minimum, maximum = mean_minimum, mean_maximum
                        panel_mean_scale_mode = mean_scale_mode
                        panel_std_scale_mode = std_scale_mode
                        profile = display["std"]["profile"]
                        right_profile = display["mean"]["profile"]
                    else:
                        matrix = np.asarray(panel["matrix"], dtype=float)
                        panel_scale_mode, panel_minimum, panel_maximum = scale_mode, 0.0, max_value
                        if panel.get("statistic") and split_scale_settings:
                            panel_scale_mode, panel_minimum, panel_maximum = split_scale_settings[panel["statistic"]]
                        pixels, minimum, maximum = _viridis_pixels(
                            matrix, panel_scale_mode, panel_maximum, panel_minimum)
                        std_minimum = std_maximum = None
                        profile = _column_profile(matrix, minimum, maximum) if panel.get("statistic") else None
                        if panel.get("statistic") and panel.get("mean_matrix") is not None:
                            display = _statistics_display(
                                panel["mean_matrix"], panel["std_matrix"], split_scale_settings)
                            result = display[panel["statistic"]]
                            pixels, minimum, maximum = result["pixels"], result["minimum"], result["maximum"]
                            profile = result["profile"]
                        right_profile = None
                    image_data = pixels.tobytes()
                    height, width = matrix.shape
                    image = QtGui.QImage(
                        image_data, width, height, width * 3,
                        QtGui.QImage.Format_RGB888,
                    )
                    rendered_panels.append({
                        **panel,
                        "matrix": matrix,
                        "labels": list(panel["labels"]),
                        "image_data": image_data,
                        "image": image,
                        "minimum": minimum,
                        "maximum": maximum,
                        "std_minimum": std_minimum,
                        "std_maximum": std_maximum,
                        "scale_mode": scale_mode if panel.get("split") else panel_scale_mode,
                        "mean_scale_mode": panel_mean_scale_mode if panel.get("split") else None,
                        "std_scale_mode": panel_std_scale_mode if panel.get("split") else None,
                        "profile": profile,
                        "right_profile": right_profile,
                        "show_profile": bool(panel.get("show_profile", show_profile) and profile is not None),
                        "show_right_profile": bool(panel.get("show_right_profile", False)),
                    })
                self._panels = rendered_panels
                profiles_visible = any(panel["show_profile"] or panel["show_right_profile"]
                                       for panel in rendered_panels)
                self.setMinimumHeight(800 if profiles_visible else 480)
                self.setMinimumWidth(1200 if len(rendered_panels) > 1 else
                                     1000 if any(panel["show_right_profile"] for panel in rendered_panels) else 760)
                self.update()

            def set_message(self, message):
                self._panels = []
                self._message = message
                self.update()

            def paintEvent(self, _event):
                self._used_scale_controls = set()
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
                    profile_reserve = 330 if profile_visible else 0
                    right_reserve = 320 if panel["show_right_profile"] else 0
                    left, top, right, bottom = 76, 88, 170 + right_reserve, 72 + profile_reserve
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
                    painter.drawText(QtCore.QRectF(panel_left + left, 6,
                                                    map_size, 28),
                                     QtCore.Qt.AlignCenter | QtCore.Qt.TextWordWrap,
                                     panel["title"])
                    painter.setFont(QtGui.QFont("Segoe UI", 8))

                    matrix = panel["matrix"]
                    residue_count = matrix.shape[0]
                    tick_count = min(7, residue_count, max(2, int(map_size / 85) + 1))
                    ticks = np.unique(np.linspace(0, residue_count - 1,
                                                  tick_count, dtype=int))
                    painter.setPen(QtGui.QColor(55, 65, 67))
                    for index in ticks:
                        x = plot.left() + (index + 0.5) * map_size / residue_count
                        y = plot.top() + (index + 0.5) * map_size / residue_count
                        label = panel["labels"][index]
                        painter.drawLine(QtCore.QPointF(x, plot.top()),
                                         QtCore.QPointF(x, plot.top() - 4))
                        painter.drawLine(QtCore.QPointF(plot.left() - 4, y),
                                         QtCore.QPointF(plot.left(), y))
                        painter.drawText(QtCore.QRectF(x - 30, plot.top() - 26, 60, 20),
                                         QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter,
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
                    x_label_y = plot.top() - 50
                    painter.drawText(QtCore.QRectF(plot.left(), x_label_y, map_size, 20),
                                     QtCore.Qt.AlignCenter, "Residue")

                    self._paint_colorbar(painter, panel, plot, top)
                    if profile_visible:
                        self._paint_profile(painter, panel, plot)
                    if panel["show_right_profile"]:
                        self._paint_right_profile(painter, panel, plot)
                painter.end()
                for key, controls in self.scale_controls.items():
                    controls["group"].setVisible(key in self._used_scale_controls and not self.export_render)

            def _profile_rect(self, panel, plot):
                chart_top = plot.bottom() + 70
                chart_bottom = self.height() - 60
                chart_width = (plot.width() + (300 if panel["show_right_profile"] else 0)
                               if panel.get("split") else self.width() / len(self._panels) - 116)
                return QtCore.QRectF(plot.left(), chart_top, chart_width, chart_bottom - chart_top)

            def _paint_profile(self, painter, panel, plot):
                values = np.asarray(panel["profile"], dtype=float)
                finite = np.isfinite(values)
                chart = self._profile_rect(panel, plot)
                chart_top = chart.top()
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

                tick_count = min(6, count, max(2, int(chart.width() / 85) + 1))
                ticks = np.unique(np.linspace(0, count - 1, tick_count, dtype=int))
                painter.setPen(QtGui.QColor(70, 79, 81))
                for index in ticks:
                    x = (chart.left() + index * chart.width() / max(1, count - 1))
                    painter.drawLine(QtCore.QPointF(x, chart.bottom()),
                                     QtCore.QPointF(x, chart.bottom() + 4))
                    painter.drawText(QtCore.QRectF(x - 34, chart.bottom() + 5, 68, 18),
                                     QtCore.Qt.AlignHCenter | QtCore.Qt.AlignTop,
                                     panel["labels"][index])

                painter.drawText(QtCore.QRectF(chart.left(), self.height() - 26,
                                                chart.width(), 20),
                                 QtCore.Qt.AlignCenter, "Residue")

                painter.save()
                painter.translate(chart.left() - 48, chart.center().y())
                painter.rotate(-90)
                painter.drawText(QtCore.QRectF(-chart.height() / 2, -10,
                                                chart.height(), 18),
                                 QtCore.Qt.AlignCenter, panel["profile_ylabel"])
                painter.restore()

            def _paint_right_profile(self, painter, panel, plot):
                values = np.asarray(panel["right_profile"], dtype=float)
                chart = QtCore.QRectF(plot.right() + 24, plot.top(), 230, plot.height())
                finite = np.isfinite(values)
                painter.setPen(QtGui.QColor(45, 54, 56))
                painter.drawText(QtCore.QRectF(chart.left() - 3, chart.top() - 40,
                                                chart.width() + 6, 36),
                                 QtCore.Qt.AlignCenter | QtCore.Qt.TextWordWrap,
                                 "Mean distance by residue")
                if not finite.any():
                    painter.drawText(chart, QtCore.Qt.AlignCenter | QtCore.Qt.TextWordWrap,
                                     "No pairs in range")
                    return
                maximum = max(float(values[finite].max()), 1e-12)
                for fraction in (0.0, 0.5, 1.0):
                    x = chart.left() + fraction * chart.width()
                    painter.setPen(QtGui.QColor(215, 222, 221))
                    painter.drawLine(QtCore.QPointF(x, chart.top()), QtCore.QPointF(x, chart.bottom()))
                    painter.setPen(QtGui.QColor(55, 65, 67))
                    painter.drawText(QtCore.QRectF(x - 24, chart.bottom() + 5, 48, 20),
                                     QtCore.Qt.AlignCenter, f"{fraction * maximum:.3g}")
                path = QtGui.QPainterPath()
                drawing = False
                count = len(values)
                painter.setPen(QtGui.QPen(QtGui.QColor(30, 112, 124), 2))
                for index, value in enumerate(values):
                    if not np.isfinite(value):
                        drawing = False
                        continue
                    point = QtCore.QPointF(chart.left() + value / maximum * chart.width(),
                                          chart.top() + (index + 0.5) * chart.height() / count)
                    if drawing:
                        path.lineTo(point)
                    else:
                        path.moveTo(point)
                    drawing = True
                painter.drawPath(path)
                painter.setPen(QtGui.QColor(55, 65, 67))
                painter.drawText(QtCore.QRectF(chart.left(), chart.bottom() + 28, chart.width(), 20),
                                 QtCore.Qt.AlignCenter, "Mean distance (A)")
                for index in np.unique(np.linspace(0, count - 1, min(6, count), dtype=int)):
                    y = chart.top() + (index + 0.5) * chart.height() / count
                    painter.drawText(QtCore.QRectF(chart.right() + 3, y - 10, 64, 20),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter, panel["labels"][index])

            def _paint_colorbar(self, painter, panel, plot, top):
                x = plot.right() + (335 if panel["show_right_profile"] else 14)
                def domain_for(key):
                    matrix = np.asarray(panel.get("std_matrix" if key == "std" else "mean_matrix",
                                                  panel["matrix"]), dtype=float)
                    finite = matrix[np.isfinite(matrix)]
                    return float(finite.max()) if finite.size else 0.0
                if panel.get("binary"):
                    self._paint_binary_scale(painter, x, top + 8, max(100, plot.height() - 150), panel)
                elif panel.get("split"):
                    half_height = (self.height() - top - 20) / 2
                    scale_height = max(35, half_height - 205)
                    self._paint_scale(
                        painter, x, top + 8, scale_height,
                        panel["minimum"], panel["maximum"], panel["mean_scale_mode"], "Mean distance (A)",
                        panel["mean_summary"], "mean", domain_for("mean"),
                    )
                    self._paint_scale(
                        painter, x, top + half_height + 12, scale_height,
                        panel["std_minimum"], panel["std_maximum"], panel["std_scale_mode"], "Distance SD (A)",
                        panel["std_summary"], "std", domain_for("std"),
                    )
                else:
                    self._paint_scale(
                        painter, x, top + 8, max(75, plot.height() - 200),
                        panel["minimum"], panel["maximum"], panel["scale_mode"],
                            panel["legend"], panel["summary"], panel.get("statistic", "frame"),
                            domain_for(panel.get("statistic")) if panel.get("statistic") else None,
                    )

            def _paint_scale(self, painter, x, y, bar_height, minimum, maximum,
                                scale_mode, title, summary, control_key, domain=None):
                bar_width = 16
                painter.setPen(QtGui.QColor(55, 65, 67))
                scale_label = f"{title}\n{scale_mode.capitalize()} scale"
                painter.drawText(QtCore.QRectF(x - 5, y - 4, 150, 30),
                                 QtCore.Qt.AlignLeft | QtCore.Qt.TextWordWrap, scale_label)
                bar_top = y + 85
                if maximum is None:
                    painter.drawText(QtCore.QRectF(x, bar_top, 124, 24),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter,
                                     "No valid distances")
                    return

                full_domain = max(domain, maximum, minimum) if domain is not None else None

                if maximum == minimum and (full_domain is None or full_domain == 0):
                    color = QtGui.QColor(*np.rint(_VIRIDIS_STOPS[0] * 255).astype(int))
                    painter.fillRect(QtCore.QRectF(x, bar_top, bar_width, bar_height), color)
                else:
                    gradient_steps = min(160, max(40, int(bar_height)))
                    for step in range(gradient_steps):
                        fraction = 1.0 - step / max(1, gradient_steps - 1)
                        if full_domain is not None:
                            value = (np.expm1(fraction * np.log1p(full_domain))
                                     if scale_mode == "log" else fraction * full_domain)
                            fraction = ((np.log1p(value) - np.log1p(minimum)) /
                                        max(np.log1p(maximum) - np.log1p(minimum), 1e-12)
                                        if scale_mode == "log" else
                                        (value - minimum) / max(maximum - minimum, 1e-12))
                        else:
                            value = minimum + fraction * (maximum - minimum)
                        fraction = float(np.clip(fraction, 0, 1))
                        palette_position = fraction * (len(_VIRIDIS_STOPS) - 1)
                        lower = int(np.floor(palette_position))
                        upper = min(lower + 1, len(_VIRIDIS_STOPS) - 1)
                        blend = palette_position - lower
                        color_values = (_VIRIDIS_STOPS[lower] * (1 - blend)
                                        + _VIRIDIS_STOPS[upper] * blend)
                        color = QtGui.QColor(*(np.rint(color_values * 255).astype(int)))
                        if value < minimum:
                            color = QtGui.QColor(80, 80, 80)
                        elif value > maximum:
                            color = QtGui.QColor(210, 210, 210)
                        y_position = bar_top + step * bar_height / gradient_steps
                        painter.fillRect(QtCore.QRectF(x, y_position, bar_width,
                                                       bar_height / gradient_steps + 1), color)
                painter.setPen(QtGui.QPen(QtGui.QColor(65, 74, 76), 1))
                painter.drawRect(QtCore.QRectF(x, bar_top, bar_width, bar_height))
                painter.setPen(QtGui.QColor(55, 65, 67))
                ticks = ([] if minimum is None else list(np.expm1(np.linspace(
                    np.log1p(minimum), np.log1p(maximum), 3))) if scale_mode == "log"
                    else list(np.linspace(minimum, maximum, 3)))
                if minimum == maximum:
                    ticks = [minimum]
                if full_domain is not None:
                    ticks = [float(np.expm1(0.5 * np.log1p(full_domain)))
                             if scale_mode == "log" else full_domain / 2]
                    if self.export_render:
                        for value, label, offset in ((maximum, "Upper", -12),
                                                     (minimum, "Lower", bar_height - 12)):
                            painter.drawText(QtCore.QRectF(x + bar_width + 7, bar_top + offset, 120, 24),
                                             QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter,
                                             f"{label}: {_format_scale_value(value)} A")
                for value in ticks:
                    if (full_domain is None and not self.export_render and control_key in self.scale_controls
                            and (np.isclose(value, minimum) or np.isclose(value, maximum))):
                        continue
                    if full_domain is not None:
                        position = ((np.log1p(value) / np.log1p(full_domain)
                                     if scale_mode == "log" else value / full_domain)
                                    if full_domain > 0 else 0.5)
                    elif maximum > minimum and scale_mode == "log":
                        position = ((np.log1p(value) - np.log1p(minimum))
                                    / (np.log1p(maximum) - np.log1p(minimum)))
                    else:
                        position = (value - minimum) / (maximum - minimum) if maximum > minimum else 0.0
                    y_position = bar_top + (1.0 - position) * bar_height
                    painter.drawLine(QtCore.QPointF(x + bar_width, y_position),
                                     QtCore.QPointF(x + bar_width + 4, y_position))
                    painter.drawText(QtCore.QRectF(x + bar_width + 7, y_position - 9,
                                                    67, 18),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter,
                                     _format_scale_value(value))
                painter.setPen(QtGui.QColor(45, 54, 56))
                for offset, color, label in ((0, (80, 80, 80), "Below range"),
                                             (18, (210, 210, 210), "Above range")):
                    painter.fillRect(QtCore.QRectF(x, bar_top + bar_height + 12 + offset, 12, 12),
                                     QtGui.QColor(*color))
                    painter.drawText(QtCore.QRectF(x + 17, bar_top + bar_height + 9 + offset, 110, 18),
                                     QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter, label)
                painter.drawText(QtCore.QRectF(x - 3, bar_top + bar_height + 50,
                                                138, 62),
                                 QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop
                                 | QtCore.Qt.TextWordWrap,
                                 "\n".join(summary))

                self._place_scale_controls(control_key, x, y, bar_height, maximum,
                                           minimum, full_domain, scale_mode)

            def _place_scale_controls(self, key, x, y, bar_height, maximum,
                                      minimum, domain, scale_mode):
                controls = self.scale_controls.get(key)
                if controls is None or self.export_render:
                    return
                self._used_scale_controls.add(key)
                controls["group"].setGeometry(int(x), int(y), 155, int(bar_height + 210))
                controls["log"].setGeometry(0, 29, 52, 23)
                controls["linear"].setGeometry(55, 29, 65, 23)
                controls["auto"].setGeometry(0, 55, 115, 22)
                for name in ("log", "linear", "auto", "maximum"):
                    controls[name].show()
                controls["maximum"].setGeometry(24, 73, controls["maximum"].width(), 25)
                if controls["auto"].isChecked() and maximum is not None:
                    controls["maximum"].blockSignals(True)
                    controls["maximum"].setValue(maximum)
                    controls["maximum"].blockSignals(False)
                if "minimum" in controls:
                    controls["minimum"].setGeometry(24, int(73 + bar_height), controls["minimum"].width(), 25)
                    controls["minimum"].show()
                    controls["range_slider"].setGeometry(-3, 85, 25, max(30, int(bar_height)))
                    controls["range_slider"].set_range(minimum, maximum, domain, scale_mode)
                    controls["range_slider"].show()
                if key == "frame":
                    for name in ("cutoff", "cutoff_label", "slider"):
                        controls[name].hide()

            def _paint_binary_scale(self, painter, x, y, height, panel):
                painter.setPen(QtGui.QColor(45, 54, 56))
                painter.drawText(QtCore.QRectF(x, y, 150, 30), QtCore.Qt.AlignLeft, "Cutoff contacts")
                bar_top = y + 38
                painter.fillRect(QtCore.QRectF(x, bar_top, 18, height / 2), QtGui.QColor(253, 231, 37))
                painter.fillRect(QtCore.QRectF(x, bar_top + height / 2, 18, height / 2), QtGui.QColor(68, 1, 84))
                painter.drawText(QtCore.QRectF(x + 25, bar_top, 120, 25), QtCore.Qt.AlignLeft, "1: contact")
                painter.drawText(QtCore.QRectF(x + 25, bar_top + height - 23, 120, 25),
                                 QtCore.Qt.AlignLeft, "0: no contact")
                painter.drawText(QtCore.QRectF(x, bar_top + height + 45, 145, 55),
                                 QtCore.Qt.TextWordWrap, "\n".join(panel["summary"]))
                if self.export_render:
                    painter.drawText(QtCore.QRectF(x, bar_top + height + 8, 145, 25),
                                     QtCore.Qt.AlignLeft, f"Cutoff: {panel['cutoff']:.1f} A")
                    return
                controls = self.scale_controls.get("frame")
                if controls is None:
                    return
                self._used_scale_controls.add("frame")
                controls["group"].setGeometry(int(x), int(y), 155, int(height + 145))
                for name in ("log", "linear", "auto", "maximum"):
                    controls[name].hide()
                controls["slider"].setGeometry(23, 65, 20, max(30, int(height - 54)))
                controls["cutoff_label"].setFixedWidth(43)
                controls["cutoff_label"].setGeometry(0, int(height + 48), 43, 25)
                controls["cutoff"].setGeometry(49, int(height + 48), 100, 25)
                for name in ("slider", "cutoff", "cutoff_label"):
                    controls[name].show()

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
                self._raw_matrices = None
                self._csv_datasets = {}
                self._updating = False
                self._exporting_frames = False

                layout = QtWidgets.QVBoxLayout(self)
                layout.setContentsMargins(6, 6, 6, 6)
                layout.setSpacing(3)
                controls = QtWidgets.QHBoxLayout()
                controls.setSpacing(6)
                object_label = QtWidgets.QLabel("Object", self)
                object_label.setFixedWidth(object_label.fontMetrics().horizontalAdvance("Object") + 4)
                controls.addWidget(object_label)
                self.object_combo = QtWidgets.QComboBox(self)
                self.object_combo.setFixedWidth(190)
                controls.addWidget(self.object_combo)

                cutoff_label = QtWidgets.QLabel("Cutoff", self)
                cutoff_label.setFixedWidth(cutoff_label.fontMetrics().horizontalAdvance("Cutoff") + 4)
                self.cutoff_spin = QtWidgets.QDoubleSpinBox(self)
                self.cutoff_spin.setRange(1.0, 30.0)
                self.cutoff_spin.setSingleStep(0.5)
                self.cutoff_spin.setValue(8.0)
                self.cutoff_spin.setSuffix(" A")
                self.cutoff_spin.setFixedWidth(90)

                controls.addWidget(QtWidgets.QLabel("Exclude same-chain neighbors <=", self))
                self.separation_spin = QtWidgets.QSpinBox(self)
                self.separation_spin.setRange(0, 20)
                self.separation_spin.setValue(0)
                self.separation_spin.setSuffix(" residues")
                controls.addWidget(self.separation_spin)

                self.refresh_button = QtWidgets.QPushButton("Refresh objects", self)
                controls.addWidget(self.refresh_button)
                self.save_button = QtWidgets.QPushButton("Save plot…", self)
                self.save_button.hide()
                self.export_frames_button = QtWidgets.QPushButton("Export frames…", self)
                self.export_frames_button.setEnabled(False)
                self.export_frames_button.hide()
                controls.addStretch(1)
                self.export_button = QtWidgets.QPushButton("Export…", self)
                self.export_button.setEnabled(False)
                controls.addWidget(self.export_button)
                layout.addLayout(controls)

                view_controls = QtWidgets.QHBoxLayout()
                view_controls.addWidget(QtWidgets.QLabel("Plot", self))
                self.plot_mode = QtWidgets.QComboBox(self)
                self.plot_mode.addItem("Current frame cutoff", "contacts")
                self.plot_mode.addItem("Current frame residue-center distances", "frame_centers")
                self.plot_mode.addItem("Mean + standard deviation", "both")
                self.plot_mode.addItem("Split mean / SD across diagonal", "split")
                view_controls.addWidget(self.plot_mode)
                view_controls.addSpacing(18)
                self.global_scale_label = QtWidgets.QLabel("Color scale", self)
                self.global_scale_label.hide()
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
                self.auto_max_check = QtWidgets.QCheckBox("Auto max", self)
                self.auto_max_check.setChecked(True)
                self.auto_max_check.hide()
                self.global_max_label = QtWidgets.QLabel("Max (A)", self)
                self.global_max_label.hide()
                self.max_range_spin = QtWidgets.QDoubleSpinBox(self)
                self.max_range_spin.setRange(0.0, 10000.0)
                self.max_range_spin.setSingleStep(1.0)
                self.max_range_spin.setDecimals(2)
                self.max_range_spin.setValue(40.0)
                self.max_range_spin.setEnabled(False)
                self.max_range_spin.setFixedWidth(86)
                view_controls.addStretch(1)
                layout.addLayout(view_controls)

                self.statistics_row = QtWidgets.QWidget(self)
                statistics_controls = QtWidgets.QHBoxLayout(self.statistics_row)
                statistics_controls.setContentsMargins(0, 0, 0, 0)
                statistics_controls.addWidget(QtWidgets.QLabel("Distance definition", self))
                self.distance_metric = QtWidgets.QComboBox(self)
                self.distance_metric.addItem("C-alpha to C-alpha", "ca")
                self.distance_metric.addItem("Residue center to center (all atoms)", "centers")
                self.distance_metric.setToolTip(
                    "Residue centers are geometric averages of all loaded protein atoms, including any hydrogens present.")
                statistics_controls.addWidget(self.distance_metric)
                self.show_mean_check = QtWidgets.QCheckBox("Show mean", self)
                self.show_std_check = QtWidgets.QCheckBox("Show standard deviation", self)
                self.show_mean_check.setChecked(True)
                self.show_std_check.setChecked(True)
                statistics_controls.addWidget(self.show_mean_check)
                statistics_controls.addWidget(self.show_std_check)
                statistics_controls.addStretch(1)
                layout.addWidget(self.statistics_row)

                self.line_plot_controls = QtWidgets.QWidget(self)
                profile_controls = QtWidgets.QHBoxLayout(self.line_plot_controls)
                profile_controls.setContentsMargins(0, 0, 0, 0)
                profile_controls.addWidget(QtWidgets.QLabel("Line plots", self.line_plot_controls))
                self.profile_check = QtWidgets.QCheckBox(
                    "Show mean line plot", self.line_plot_controls)
                self.profile_check.setToolTip(
                    "Residues on the x-axis; mean of eligible pairwise mean distances on the y-axis.")
                self.profile_check.setEnabled(False)
                profile_controls.addWidget(self.profile_check)
                self.std_profile_check = QtWidgets.QCheckBox(
                    "Show SD line plot", self.line_plot_controls)
                self.std_profile_check.setToolTip(
                    "Residues on the x-axis; mean of eligible pairwise distance SDs on the y-axis.")
                self.std_profile_check.setEnabled(False)
                profile_controls.addWidget(self.std_profile_check)
                profile_controls.addStretch(1)
                self.line_plot_controls.hide()
                layout.addWidget(self.line_plot_controls)

                self.csv_statistic = QtWidgets.QComboBox(self)
                self.csv_statistic.setMinimumWidth(210)
                self.csv_statistic.setToolTip("Select the heatmap matrix or residue profile to export")
                self.csv_statistic.hide()
                self.export_csv_button = QtWidgets.QPushButton("Export CSV…", self)
                self.export_csv_button.setEnabled(False)
                self.export_csv_button.hide()

                self.split_scale_row = QtWidgets.QWidget(self)
                split_scale_layout = QtWidgets.QHBoxLayout(self.split_scale_row)
                split_scale_layout.setContentsMargins(0, 0, 0, 0)
                self.mean_split_controls = self._create_split_scale_controls(
                    "Mean distance / upper triangle")
                self.std_split_controls = self._create_split_scale_controls(
                    "Standard deviation / lower triangle")
                self.split_scale_row.hide()

                self.frame_label = QtWidgets.QLabel("Frame: -", self)
                layout.addWidget(self.frame_label)
                self.canvas = ContactMapCanvas(self)
                self.plot_scroll = QtWidgets.QScrollArea(self)
                self.plot_scroll.setWidgetResizable(True)
                self.plot_scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
                self.plot_scroll.setWidget(self.canvas)
                layout.addWidget(self.plot_scroll, 1)
                frame_group = QtWidgets.QWidget(self.canvas)
                for widget in (self.log_scale_button, self.linear_scale_button,
                               self.auto_max_check, self.max_range_spin,
                               self.cutoff_spin, cutoff_label):
                    widget.setParent(frame_group)
                cutoff_slider = QtWidgets.QSlider(QtCore.Qt.Vertical, frame_group)
                cutoff_slider.setRange(10, 300)
                cutoff_slider.setValue(int(self.cutoff_spin.value() * 10))
                cutoff_slider.setToolTip("Adjust the contact-distance cutoff in Angstroms")
                cutoff_slider.valueChanged.connect(lambda value: self.cutoff_spin.setValue(value / 10.0))
                self.cutoff_spin.valueChanged.connect(lambda value: cutoff_slider.setValue(round(value * 10)))
                self.canvas.scale_controls = {
                    "mean": self.mean_split_controls,
                    "std": self.std_split_controls,
                    "frame": {"group": frame_group, "log": self.log_scale_button,
                              "linear": self.linear_scale_button, "auto": self.auto_max_check,
                              "maximum": self.max_range_spin, "cutoff": self.cutoff_spin,
                              "cutoff_label": cutoff_label, "slider": cutoff_slider},
                }
                for control in self.canvas.scale_controls.values():
                    control["group"].setParent(self.canvas)
                    control["group"].hide()
                self.status_label = QtWidgets.QLabel(
                    "Cutoff contacts use CA atoms; residue centers use all loaded protein atoms.", self)
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
                self.std_profile_check.toggled.connect(self._settings_changed)
                self.distance_metric.currentIndexChanged.connect(self._settings_changed)
                self.show_mean_check.toggled.connect(self._settings_changed)
                self.show_std_check.toggled.connect(self._settings_changed)
                self.export_csv_button.clicked.connect(self.export_raw_csv)
                for scale_controls in (self.mean_split_controls, self.std_split_controls):
                    scale_controls["log"].toggled.connect(self._split_scale_changed)
                    scale_controls["linear"].toggled.connect(self._split_scale_changed)
                    scale_controls["auto"].toggled.connect(
                        lambda checked, controls=scale_controls: self._toggle_split_auto(
                            controls, checked))
                    scale_controls["maximum"].valueChanged.connect(self._settings_changed)
                    scale_controls["minimum"].valueChanged.connect(self._settings_changed)
                self.save_button.clicked.connect(self.save_plot)
                self.export_frames_button.clicked.connect(self.export_frames_dialog)
                self.export_button.clicked.connect(self.export_dialog)

                self.timer = QtCore.QTimer(self)
                self.timer.setInterval(350)
                self.timer.timeout.connect(self.refresh_if_needed)
                self.refresh_objects()
                self.timer.start()

            def _create_split_scale_controls(self, title):
                group = QtWidgets.QWidget(self)
                group.setToolTip(title)
                log_button = QtWidgets.QRadioButton("Log", group)
                log_button.setChecked(True)
                linear_button = QtWidgets.QRadioButton("Linear", group)
                scale_group = QtWidgets.QButtonGroup(group)
                scale_group.setExclusive(True)
                scale_group.addButton(log_button)
                scale_group.addButton(linear_button)
                auto_check = QtWidgets.QCheckBox("Auto max", group)
                auto_check.setChecked(True)
                minimum_spin = QtWidgets.QDoubleSpinBox(group)
                minimum_spin.setRange(0.0, 9999.99)
                minimum_spin.setDecimals(2)
                minimum_spin.setSuffix(" A")
                minimum_spin.setFixedWidth(90)
                maximum_spin = QtWidgets.QDoubleSpinBox(group)
                maximum_spin.setRange(0.0, 10000.0)
                maximum_spin.setSingleStep(1.0)
                maximum_spin.setDecimals(2)
                maximum_spin.setValue(40.0)
                maximum_spin.setSuffix(" A")
                maximum_spin.setEnabled(False)
                maximum_spin.setFixedWidth(90)
                minimum_spin.setToolTip("Lower colorbar bound in Angstroms")
                maximum_spin.setToolTip("Upper colorbar bound in Angstroms")
                minimum_spin.setFixedWidth(126)
                maximum_spin.setFixedWidth(126)
                minimum_spin.setPrefix("Lower: ")
                maximum_spin.setPrefix("Upper: ")
                range_slider = ContactMapRangeSlider(group)
                controls = {"minimum": minimum_spin, "maximum": maximum_spin, "auto": auto_check}
                range_slider.on_change = lambda name, value: self._range_handle_changed(controls, name, value)
                return {
                    "group": group,
                    "log": log_button,
                    "linear": linear_button,
                    "auto": auto_check,
                    "minimum": minimum_spin,
                    "maximum": maximum_spin,
                    "range_slider": range_slider,
                }

            def _range_handle_changed(self, controls, name, value):
                if name == "maximum":
                    controls["auto"].blockSignals(True)
                    controls["auto"].setChecked(False)
                    controls["auto"].blockSignals(False)
                    controls["maximum"].setEnabled(True)
                controls[name].blockSignals(True)
                controls[name].setValue(value)
                controls[name].blockSignals(False)
                self._settings_changed(value)

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
                    settings[key] = (mode, controls["minimum"].value(), maximum)
                return settings

            def _update_control_states(self, mode):
                contact_mode = mode == "contacts"
                split_mode = mode in {"both", "split"}
                self.statistics_row.setVisible(split_mode)
                self.line_plot_controls.setVisible(split_mode)
                self.show_mean_check.setEnabled(mode == "both")
                self.show_std_check.setEnabled(mode == "both")
                self.cutoff_spin.setEnabled(contact_mode)
                distance_mode = not contact_mode
                self.global_scale_label.setEnabled(not split_mode)
                self.auto_max_check.setEnabled(distance_mode and not split_mode)
                self.max_range_spin.setEnabled(
                    distance_mode and not split_mode and not self.auto_max_check.isChecked())
                self.profile_check.setEnabled(split_mode)
                self.std_profile_check.setEnabled(split_mode)
                self.csv_statistic.setEnabled(bool(self._csv_datasets))
                self.export_csv_button.setEnabled(bool(self._csv_datasets))

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
                    f"({object_name}) and polymer.protein", state=state)
                residue_atoms = {}
                for atom in model.atom:
                    key = (atom.chain or "", atom.resi, atom.resn)
                    residue_atoms.setdefault(key, []).append(atom.coord)
                if not residue_atoms:
                    raise ValueError(f"No protein atoms found in frame {state}.")
                keys = list(residue_atoms)
                centers = np.asarray([
                    np.mean(residue_atoms[key], axis=0) for key in keys
                ], dtype=float)
                chains = [key[0] for key in keys]
                labels = [f"{chain or '-'}:{resn}{resi}"
                          for chain, resi, resn in keys]
                return centers, chains, labels

            def _aggregate_maps(self, object_name, state_count, separation):
                metric = self.distance_metric.currentData()
                cache_key = (object_name, state_count, separation, metric)
                if cache_key in self._aggregate_cache:
                    return self._aggregate_cache[cache_key]

                mean = None
                sum_squared_deviations = None
                expected_labels = None
                expected_identities = None
                valid_pairs = None
                total_pair_distance = 0.0
                total_pair_distance_squared = 0.0
                total_pair_observations = 0
                progress_step = max(1, state_count // 20)
                for state in range(1, state_count + 1):
                    selection = "polymer.protein" + (" and name CA" if metric == "ca" else "")
                    model = cmd.get_model(f"({object_name}) and ({selection})", state=state)
                    distances, chains, labels, identities = _residue_distances(model.atom, metric)
                    if expected_labels is None:
                        expected_labels = labels
                        expected_identities = identities
                        mean = np.zeros((len(labels), len(labels)), dtype=np.float64)
                        sum_squared_deviations = np.zeros_like(mean)
                        excluded = _excluded_pairs(chains, separation)
                        np.fill_diagonal(excluded, False)
                        valid_pairs = ~excluded
                        if not np.any(np.triu(valid_pairs, k=1)):
                            raise ValueError(
                                "No residue pairs remain after the same-chain neighbor filter."
                            )
                    elif identities != expected_identities:
                        raise ValueError(
                            f"Residue selection changes in frame {state}; "
                            "distance statistics require matching residues in every frame."
                        )

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

            def _plot_key(self):
                object_name = self.object_combo.currentText()
                mode = self.plot_mode.currentData()
                return (
                    object_name, cmd.get_frame() if mode in {"contacts", "frame_centers"} else None,
                    cmd.count_states(object_name) if object_name else 0,
                    self.cutoff_spin.value(), self.separation_spin.value(), mode,
                    self._scale_settings(), self._split_scale_settings(),
                    self.distance_metric.currentData(), self.show_mean_check.isChecked(),
                    self.show_std_check.isChecked(), self.profile_check.isChecked(),
                    self.std_profile_check.isChecked(),
                )

            def refresh_if_needed(self):
                if self._updating or self._exporting_frames or not self.object_combo.currentText():
                    return
                if self._plot_key() != self._last_key:
                    self.update_contact_map()

            def closeEvent(self, event):
                if self._exporting_frames:
                    event.ignore()
                    return
                self.timer.stop()
                super().closeEvent(event)

            def update_contact_map(self):
                if self._updating or self._exporting_frames:
                    return
                self._updating = True
                self.export_button.setEnabled(False)
                self.export_frames_button.setEnabled(False)
                self.export_csv_button.setEnabled(False)
                self._set_csv_datasets({})
                self._raw_matrices = None
                object_name = self.object_combo.currentText()
                if not object_name:
                    self.frame_label.setText("Frame: -")
                    self.status_label.setText("Load a protein object in PyMOL to begin.")
                    self.canvas.set_message("No molecular objects available.")
                    self._updating = False
                    return

                try:
                    plot_key = self._plot_key()
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
                            "binary": True, "cutoff": cutoff,
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
                        self._raw_matrices = {
                            "mean": mean, "std": standard_deviation,
                            "labels": labels, "object": object_name,
                            "metric": self.distance_metric.currentData(),
                        }
                        metric_name = self.distance_metric.currentText()
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
                            "statistic": "mean",
                            "mean_matrix": mean, "std_matrix": standard_deviation,
                            "matrix": mean,
                            "labels": labels,
                            "title": f"Mean distance: {metric_name}",
                            "legend": "Mean distance (A)",
                            "summary": mean_summary,
                            "profile": _column_profile(mean),
                            "profile_title": "Mean of means by residue",
                            "profile_ylabel": "Mean of means (A)",
                            "show_profile": self.profile_check.isChecked(),
                        }
                        sd_panel = {
                            "statistic": "std",
                            "mean_matrix": mean, "std_matrix": standard_deviation,
                            "matrix": standard_deviation,
                            "labels": labels,
                            "title": f"Distance SD: {metric_name}",
                            "legend": "Distance SD (A)",
                            "summary": sd_summary,
                            "profile": _column_profile(standard_deviation),
                            "profile_title": "Mean of SDs by residue",
                            "profile_ylabel": "Mean of SDs (A)",
                            "show_profile": self.std_profile_check.isChecked(),
                        }
                        if mode == "split":
                            panels = [{
                                "split": True,
                                "matrix": mean,
                                "mean_matrix": mean,
                                "std_matrix": standard_deviation,
                                "labels": labels,
                                "title": f"Mean (upper) / SD (lower): {metric_name}",
                                "mean_summary": mean_summary,
                                "std_summary": sd_summary,
                                "show_profile": self.std_profile_check.isChecked(),
                                "show_right_profile": self.profile_check.isChecked(),
                                "profile_title": "Mean of SDs by residue",
                                "profile_ylabel": "Mean of SDs (A)",
                            }]
                            status = f"{metric_name} | Upper: mean; lower: population SD"
                        else:
                            panels = []
                            if self.show_mean_check.isChecked():
                                panels.append(mean_panel)
                            if self.show_std_check.isChecked():
                                panels.append(sd_panel)
                            status = f"Mean distance and population SD across {state_count} frames"
                        self.frame_label.setText(f"Statistics: {state_count} frames")

                    split_scale_settings = (
                        self._split_scale_settings() if mode in {"both", "split"} else None
                    )
                    self.canvas.set_panels(
                        panels, scale_mode, max_value,
                        split_scale_settings=split_scale_settings,
                    )
                    if not panels:
                        self.canvas.set_message("Select Show mean and/or Show standard deviation.")
                    self.status_label.setText(status)
                    self._prepare_csv_datasets(mode, object_name, state)
                    self.export_frames_button.setEnabled(True)
                    self.export_button.setEnabled(True)
                    self._last_key = plot_key
                except Exception as exc:
                    self.status_label.setText(f"Could not update contact map: {exc}")
                    if not self.canvas._panels:
                        self.canvas.set_message("Contact map update failed; see status below.")
                    self._raw_matrices = None
                    self._set_csv_datasets({})
                    self._last_key = self._plot_key()
                finally:
                    self._updating = False

            def _frame_export_panel(self, object_name, state, mode, metric, cutoff, separation):
                if mode == "contacts":
                    coordinates, chains, labels = self._frame_data(object_name, state)
                    matrix = _contact_matrix(coordinates, chains, cutoff, separation).astype(float)
                    title = "Cutoff contacts"
                    legend = "Contact (0/1)"
                else:
                    frame_metric = "centers" if mode == "frame_centers" else metric
                    selection = "polymer.protein" + (" and name CA" if frame_metric == "ca" else "")
                    model = cmd.get_model(f"({object_name}) and ({selection})", state=state)
                    matrix, chains, labels, _ = _residue_distances(model.atom, frame_metric)
                    excluded = _excluded_pairs(chains, separation)
                    np.fill_diagonal(excluded, False)
                    matrix[excluded] = np.nan
                    title = {
                        "ca": "C-alpha to C-alpha distance",
                        "centers": "Residue-center distance",
                        "closest": "closest inter residue atomic distance",
                    }[frame_metric]
                    legend = "Distance (A)"
                return {
                    "matrix": matrix, "labels": labels,
                    "binary": mode == "contacts", "cutoff": cutoff,
                    "title": f"{object_name}: {title}, frame {state}",
                    "legend": legend,
                    "summary": [f"Frame: {state}", f"Residues: {len(labels)}"],
                }

            def export_dialog(self):
                if self._updating or self._exporting_frames or not self._csv_datasets:
                    return
                dialog = QtWidgets.QDialog(self)
                dialog.setWindowTitle("Export plots and data")
                form = QtWidgets.QFormLayout(dialog)
                scope = QtWidgets.QComboBox(dialog)
                scope.setObjectName("export_scope")
                scope.addItem("Current view and selected CSV datasets", "current")
                scope.addItem("Selected frames", "frames")
                form.addRow("Export", scope)
                destination = QtWidgets.QLineEdit(str(Path.home()), dialog)
                destination.setObjectName("export_destination")
                destination_row = QtWidgets.QHBoxLayout()
                destination_row.addWidget(destination)
                browse = QtWidgets.QPushButton("Browse…", dialog)
                destination_row.addWidget(browse)

                def choose_folder():
                    folder = QtWidgets.QFileDialog.getExistingDirectory(dialog, "Export destination", destination.text())
                    if folder:
                        destination.setText(folder)

                browse.clicked.connect(choose_folder)
                form.addRow("Destination", destination_row)
                save_image = QtWidgets.QCheckBox("Save PNG plot", dialog)
                save_image.setObjectName("export_png")
                save_image.setChecked(True)
                form.addRow(save_image)
                image_format = QtWidgets.QComboBox(dialog)
                image_format.addItem("PNG", "png")
                image_format.addItem("JPEG", "jpg")
                form.addRow("Current-view image format", image_format)
                dataset_group = QtWidgets.QGroupBox("CSV datasets (select any)", dialog)
                dataset_layout = QtWidgets.QVBoxLayout(dataset_group)
                dataset_checks = {}
                snapshot = {key: {**dataset, "values": dataset["values"].copy(),
                                   "labels": list(dataset["labels"])}
                            for key, dataset in self._csv_datasets.items()}
                for key, dataset in snapshot.items():
                    checkbox = QtWidgets.QCheckBox(dataset["title"], dataset_group)
                    checkbox.setObjectName("export_csv_" + key)
                    checkbox.setChecked(True)
                    dataset_layout.addWidget(checkbox)
                    dataset_checks[key] = checkbox
                form.addRow(dataset_group)
                frame_group = QtWidgets.QGroupBox("Frame range (loaded PyMOL states)", dialog)
                frame_form = QtWidgets.QFormLayout(frame_group)
                first, last, interval = (QtWidgets.QSpinBox(frame_group) for _ in range(3))
                for name, control in (("first", first), ("last", last), ("interval", interval)):
                    control.setObjectName("export_frame_" + name)
                count = int(cmd.count_states(self.object_combo.currentText()))
                for control in (first, last, interval):
                    control.setRange(1, count)
                first.setValue(1)
                last.setValue(count)
                interval.setValue(1)
                frame_form.addRow("First", first)
                frame_form.addRow("Last (inclusive)", last)
                frame_form.addRow("Interval", interval)
                frame_csv = QtWidgets.QCheckBox("Save matching frame matrix CSVs", frame_group)
                frame_csv.setObjectName("export_frame_csv")
                frame_csv.setChecked(True)
                frame_form.addRow(frame_csv)
                form.addRow(frame_group)
                frame_group.hide()
                scope.currentIndexChanged.connect(
                    lambda: (dataset_group.setVisible(scope.currentData() == "current"),
                             frame_group.setVisible(scope.currentData() == "frames"),
                             image_format.setEnabled(scope.currentData() == "current"),
                             save_image.setText("Save PNG plot per frame" if scope.currentData() == "frames"
                                                else "Save plot image")))
                save_image.setText("Save plot image")
                buttons = QtWidgets.QDialogButtonBox(
                    QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, dialog)
                buttons.button(QtWidgets.QDialogButtonBox.Ok).setText("Export")
                buttons.accepted.connect(dialog.accept)
                buttons.rejected.connect(dialog.reject)
                form.addRow(buttons)
                timer_active = self.timer.isActive()
                playing = bool(cmd.get_movie_playing())
                self.timer.stop()
                cmd.mstop()
                try:
                    if dialog.exec_() != QtWidgets.QDialog.Accepted:
                        return
                    if not destination.text().strip():
                        raise ValueError("Choose an export destination folder.")
                    folder = Path(destination.text().strip()).expanduser()
                    if scope.currentData() == "frames":
                        result = self._export_frame_series(
                            folder, first.value(), last.value(), interval.value(),
                            save_image.isChecked(), frame_csv.isChecked())
                        status = "Canceled" if result["canceled"] else "Completed"
                        self.status_label.setText(
                            f"{status}: {result['completed']} frames in {result['directory']}")
                    else:
                        keys = [key for key, checkbox in dataset_checks.items() if checkbox.isChecked()]
                        output = self._export_current_bundle(
                            folder, save_image.isChecked(), keys, snapshot, image_format.currentData())
                        self.status_label.setText(f"Saved export in {output}")
                except Exception as exc:
                    self.status_label.setText(f"Export failed: {exc}")
                    QtWidgets.QMessageBox.warning(self, "Export failed", str(exc))
                finally:
                    if timer_active:
                        self.timer.start()
                    if playing:
                        cmd.mplay()

            def _export_current_bundle(self, destination, save_image, keys, datasets=None,
                                       image_format="png"):
                datasets = self._csv_datasets if datasets is None else datasets
                if not save_image and not keys:
                    raise ValueError("Select a PNG plot and/or one or more CSV datasets.")
                destination = Path(destination).expanduser()
                destination.mkdir(parents=True, exist_ok=True)
                if any(key not in datasets for key in keys):
                    raise ValueError("A selected CSV dataset is unavailable.")
                if image_format not in {"png", "jpg"}:
                    raise ValueError("Choose PNG or JPEG output.")
                if save_image and not self.canvas._panels:
                    raise ValueError("Select a visible plot before saving an image.")
                safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", self.object_combo.currentText())
                output = Path(tempfile.mkdtemp(prefix=safe_name + "_export_", dir=str(destination)))
                try:
                    if save_image:
                        self._save_canvas_image(output / (safe_name + "_plot." + image_format))
                    for key in keys:
                        dataset = datasets[key]
                        path = output / f"{safe_name}_{dataset['suffix']}.csv"
                        if dataset.get("value_label"):
                            _write_profile_csv(path, dataset["values"], dataset["labels"], dataset["value_label"])
                        else:
                            _write_matrix_csv(path, dataset["values"], dataset["labels"])
                except Exception:
                    shutil.rmtree(output, ignore_errors=True)
                    raise
                return output

            def export_frames_dialog(self):
                object_name = self.object_combo.currentText()
                if not object_name or self._updating or self._exporting_frames:
                    return
                count = int(cmd.count_states(object_name))
                if count < 1:
                    self.status_label.setText("No object states are available to export.")
                    return
                dialog = QtWidgets.QDialog(self)
                dialog.setWindowTitle("Export frame plots and CSVs")
                form = QtWidgets.QFormLayout(dialog)
                mode = self.plot_mode.currentData()
                source = ("Current-frame cutoff contacts" if mode == "contacts" else
                          "Current-frame residue-center distances" if mode == "frame_centers" else
                          "Per-frame distance: " + self.distance_metric.currentText())
                source_label = QtWidgets.QLabel(source, dialog)
                source_label.setWordWrap(True)
                form.addRow("Data", source_label)
                destination = QtWidgets.QLineEdit(str(Path.home()), dialog)
                destination_row = QtWidgets.QHBoxLayout()
                destination_row.addWidget(destination)
                browse = QtWidgets.QPushButton("Browse…", dialog)
                destination_row.addWidget(browse)

                def choose_destination():
                    folder = QtWidgets.QFileDialog.getExistingDirectory(
                        dialog, "Choose frame-export destination", destination.text())
                    if folder:
                        destination.setText(folder)

                browse.clicked.connect(choose_destination)
                form.addRow("Destination folder", destination_row)
                first = QtWidgets.QSpinBox(dialog)
                last = QtWidgets.QSpinBox(dialog)
                interval = QtWidgets.QSpinBox(dialog)
                for control in (first, last, interval):
                    control.setRange(1, count)
                first.setValue(1)
                last.setValue(count)
                interval.setValue(1)
                form.addRow("First frame (1-based)", first)
                form.addRow("Last frame (inclusive)", last)
                form.addRow("Frame interval", interval)
                save_images = QtWidgets.QCheckBox("Save a PNG plot for each frame", dialog)
                save_csv = QtWidgets.QCheckBox("Save the matching raw matrix CSV", dialog)
                save_images.setChecked(True)
                save_csv.setChecked(True)
                form.addRow(save_images)
                form.addRow(save_csv)
                buttons = QtWidgets.QDialogButtonBox(
                    QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, dialog)
                buttons.button(QtWidgets.QDialogButtonBox.Ok).setText("Export")
                buttons.rejected.connect(dialog.reject)
                buttons.accepted.connect(dialog.accept)
                form.addRow(buttons)
                if dialog.exec_() != QtWidgets.QDialog.Accepted:
                    return
                try:
                    result = self._export_frame_series(
                        Path(destination.text()).expanduser(), first.value(), last.value(),
                        interval.value(), save_images.isChecked(), save_csv.isChecked())
                except Exception as exc:
                    self.status_label.setText(f"Frame export failed: {exc}")
                    QtWidgets.QMessageBox.warning(self, "Frame export failed", str(exc))
                    return
                outcome = "Canceled" if result["canceled"] else "Completed"
                self.status_label.setText(
                    f"{outcome}: {result['completed']} frame(s) saved in {result['directory']}")

            def _export_frame_series(self, destination, first, last, interval,
                                     save_images=True, save_csv=True):
                if self._updating or self._exporting_frames:
                    raise ValueError("An update or export is already in progress.")
                object_name = self.object_combo.currentText()
                count = int(cmd.count_states(object_name)) if object_name else 0
                if not 1 <= first <= last <= count or interval < 1:
                    raise ValueError(f"Choose a range between 1 and {count} and a positive interval.")
                if not (save_images or save_csv):
                    raise ValueError("Select PNG plots, CSV matrices, or both.")
                destination = Path(destination).expanduser()
                destination.mkdir(parents=True, exist_ok=True)
                mode = self.plot_mode.currentData()
                metric = self.distance_metric.currentData()
                cutoff, separation = self.cutoff_spin.value(), self.separation_spin.value()
                if mode in {"both", "split"}:
                    scale_mode, minimum, maximum = self._split_scale_settings()["mean"]
                else:
                    scale_mode, maximum = self._scale_settings()
                    minimum = 0.0
                source = "contacts" if mode == "contacts" else "centers" if mode == "frame_centers" else metric
                safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", object_name)
                folder = Path(tempfile.mkdtemp(
                    prefix=f"{safe_name}_{source}_frames_", dir=str(destination)))
                states = range(first, last + 1, interval)
                original_frame = int(cmd.get_frame())
                movie_playing = bool(cmd.get_movie_playing())
                timer_active = self.timer.isActive()
                self._exporting_frames = True
                self.timer.stop()
                locked = [widget for widget in self.findChildren(QtWidgets.QWidget)
                          if isinstance(widget, (QtWidgets.QAbstractButton, QtWidgets.QComboBox,
                                                 QtWidgets.QAbstractSpinBox, QtWidgets.QLineEdit,
                                                 QtWidgets.QSlider)) or hasattr(widget, "on_change")]
                enabled_states = [(widget, widget.isEnabled()) for widget in locked]
                progress = QtWidgets.QProgressDialog("Preparing export", "Cancel", 0, len(states), self)
                progress.setWindowTitle("Exporting frame plots and CSVs")
                progress.setWindowModality(QtCore.Qt.WindowModal)
                progress.setMinimumDuration(0)
                progress.setAutoClose(False)
                completed = 0
                canceled = False
                try:
                    cmd.mstop()
                    for widget, _ in enabled_states:
                        widget.setEnabled(False)
                    progress.show()
                    for index, state in enumerate(states):
                        progress.setLabelText(f"{object_name}: frame {state} ({index + 1}/{len(states)})")
                        progress.setValue(index)
                        QtWidgets.QApplication.processEvents()
                        if progress.wasCanceled():
                            canceled = True
                            break
                        panel = self._frame_export_panel(
                            object_name, state, mode, metric, cutoff, separation)
                        state_stem = f"{safe_name}_{source}_frame_{state:06d}"
                        png_path = folder / (state_stem + ".png")
                        csv_path = folder / (state_stem + ".csv")
                        try:
                            if save_images:
                                self.canvas.set_panels([panel], scale_mode, maximum)
                                if minimum:
                                    self.canvas.set_panels(
                                        [{**panel, "statistic": "mean"}], scale_mode, maximum,
                                        {"mean": (scale_mode, minimum, maximum)})
                                self._save_canvas_image(png_path)
                            if save_csv:
                                _write_matrix_csv(csv_path, panel["matrix"], panel["labels"])
                        except Exception:
                            for output in (png_path, csv_path):
                                if output.exists():
                                    output.unlink()
                            raise
                        completed += 1
                        print(f"Frame export: {completed}/{len(states)} - frame {state}", flush=True)
                    progress.setValue(completed)
                finally:
                    progress.close()
                    self._exporting_frames = False
                    for widget, enabled in enabled_states:
                        widget.setEnabled(enabled)
                    cmd.frame(original_frame)
                    self._last_key = None
                    self.update_contact_map()
                    if timer_active:
                        self.timer.start()
                    if movie_playing:
                        cmd.mplay()
                return {"directory": folder, "completed": completed, "canceled": canceled}

            def _save_canvas_image(self, path):
                scale = 2
                snapshot = type(self.canvas)()
                snapshot.setMinimumSize(0, 0)
                snapshot.resize(self.canvas.size())
                snapshot._panels = list(self.canvas._panels)
                snapshot.export_render = True
                image = QtGui.QImage(
                    snapshot.width() * scale, snapshot.height() * scale,
                    QtGui.QImage.Format_ARGB32)
                if image.isNull():
                    snapshot.deleteLater()
                    raise MemoryError("Not enough memory to render this plot. Reduce the plot window size and retry.")
                image.setDevicePixelRatio(scale)
                image.fill(QtGui.QColor("white"))
                painter = QtGui.QPainter(image)
                try:
                    if not painter.isActive():
                        raise RuntimeError("Could not initialize the plot image renderer.")
                    snapshot.render(painter, QtCore.QPoint(), QtGui.QRegion(),
                                    QtWidgets.QWidget.DrawWindowBackground)
                finally:
                    painter.end()
                    snapshot.deleteLater()
                if not image.save(str(path)):
                    raise OSError(f"Could not save plot: {path}")

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
                try:
                    self._save_canvas_image(path)
                except OSError as exc:
                    self.status_label.setText(str(exc))
                    return
                self.status_label.setText(f"Saved plot: {path}")

            def _set_csv_datasets(self, datasets):
                previous = self.csv_statistic.currentData()
                if previous is not None:
                    self._csv_selection = previous
                self._csv_datasets = datasets
                self.csv_statistic.blockSignals(True)
                self.csv_statistic.clear()
                for key, dataset in datasets.items():
                    self.csv_statistic.addItem(dataset["title"], key)
                previous_index = self.csv_statistic.findData(
                    getattr(self, "_csv_selection", previous))
                if previous_index >= 0:
                    self.csv_statistic.setCurrentIndex(previous_index)
                self.csv_statistic.blockSignals(False)
                self.csv_statistic.setEnabled(bool(datasets))
                self.export_csv_button.setEnabled(bool(datasets))
                self.export_button.setEnabled(bool(datasets))

            def _prepare_csv_datasets(self, mode, object_name, state):
                datasets = {}
                if mode in {"contacts", "frame_centers"}:
                    panel = self.canvas._panels[0]
                    suffix = f"{mode}_frame_{state:04d}"
                    datasets[mode] = {
                        "title": ("Cutoff contacts matrix (0/1)" if mode == "contacts"
                                  else "Residue-center distance matrix (A)"),
                        "values": panel["matrix"], "labels": panel["labels"],
                        "suffix": suffix,
                    }
                elif self._raw_matrices is not None:
                    raw = self._raw_matrices
                    settings = self._split_scale_settings()
                    for statistic, title, suffix in (
                        ("mean", "Mean CSV", "mean_distance"),
                        ("std", "SD CSV", "distance_standard_deviation"),
                    ):
                        datasets[statistic] = {
                            "title": title, "values": raw[statistic], "labels": raw["labels"],
                            "suffix": raw["metric"] + "_" + suffix,
                        }
                    if mode == "split":
                        combined = np.triu(raw["mean"]) + np.tril(raw["std"], k=-1)
                        datasets["split"] = {
                            "title": "Split matrix (upper mean / lower SD)",
                            "values": combined, "labels": raw["labels"],
                            "suffix": raw["metric"] + "_split_mean_sd",
                        }
                    for statistic, title, value_label in (
                        ("mean", "Mean residue profile", "column_mean_distance_angstrom"),
                        ("std", "SD residue profile", "column_mean_distance_sd_angstrom"),
                    ):
                        display = _statistics_display(raw["mean"], raw["std"], settings)
                        datasets[statistic + "_profile"] = {
                            "title": title,
                            "values": display[statistic]["profile"],
                            "labels": raw["labels"], "value_label": value_label,
                            "suffix": raw["metric"] + "_" + statistic + "_residue_profile",
                        }
                for dataset in datasets.values():
                    dataset["object"] = object_name
                self._set_csv_datasets(datasets)

            def export_raw_csv(self):
                dataset = self._csv_datasets.get(self.csv_statistic.currentData())
                if dataset is None:
                    self.status_label.setText("Calculate a valid plot before exporting CSV.")
                    return
                values, labels = dataset["values"].copy(), list(dataset["labels"])
                value_label = dataset.get("value_label")
                safe_name = re.sub(r"[^A-Za-z0-9_-]", "_", dataset["object"])
                default_path = f"{safe_name}_{dataset['suffix']}.csv"
                path, _ = QtWidgets.QFileDialog.getSaveFileName(
                    self, "Export plot data as CSV", default_path,
                    "CSV file (*.csv)",
                )
                if not path:
                    return
                if not Path(path).suffix:
                    path += ".csv"
                try:
                    if value_label is None:
                        _write_matrix_csv(path, values, labels)
                    else:
                        _write_profile_csv(path, values, labels, value_label)
                except Exception as exc:
                    self.status_label.setText(f"Could not export CSV: {exc}")
                    return
                self.status_label.setText(f"Saved CSV: {path}")

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