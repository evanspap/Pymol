"""Run with PyMOL's Python: python -B -m unittest discover -s 'contact map'."""
import csv
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
from pymol2 import SingletonPyMOL
from pymol.Qt import QtCore, QtGui, QtWidgets


def pdb_frame(second_ca, second_cb):
    points = [(1, "CA", 0.), (1, "CB", 2.), (2, "CA", second_ca),
              (2, "CB", second_cb), (3, "CA", 20.), (3, "CB", 22.)]
    return "".join(
        f"ATOM  {serial:5d} {name:^4s} ALA A{residue:4d}    "
        f"{coordinate:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C\n"
        for serial, (residue, name, coordinate) in enumerate(points, 1)
    )


class ContactMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
        cls.instance = SingletonPyMOL()
        cls.instance.start()
        spec = importlib.util.spec_from_file_location(
            "contact_map_under_test", Path(__file__).with_name("__init__.py"))
        cls.plugin = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.plugin)

    @classmethod
    def tearDownClass(cls):
        cls.instance.stop()

    def setUp(self):
        self.cmd = self.instance.cmd
        self.cmd.read_pdbstr(
            "MODEL        1\n" + pdb_frame(10., 6.) + "ENDMDL\n"
            "MODEL        2\n" + pdb_frame(14., 10.) + "ENDMDL\nEND\n",
            "distance_test")
        self.cmd.frame(1)
        self.window = self.plugin.ContactMapWindow()

    def tearDown(self):
        self.window.close()
        self.application.processEvents()
        self.cmd.delete("distance_test")

    def select_mode(self, mode):
        self.window.plot_mode.setCurrentIndex(self.window.plot_mode.findData(mode))

    def test_distance_definitions_and_temporal_statistics(self):
        self.assertEqual(self.plugin.__version__, "v9")
        self.assertEqual([self.window.distance_metric.itemData(index)
                  for index in range(self.window.distance_metric.count())], ["ca", "centers"])
        self.select_mode("both")
        for metric, expected_mean in (("ca", 12.), ("centers", 9.)):
            with self.subTest(metric=metric):
                self.window.distance_metric.setCurrentIndex(
                    self.window.distance_metric.findData(metric))
                self.assertIsNotNone(self.window._raw_matrices, self.window.status_label.text())
                self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], expected_mean)
                self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 2.)
        self.assertEqual(len(self.window._aggregate_cache), 2)

    def test_residue_centers_include_loaded_hydrogens_in_all_paths(self):
        def atom_frame(first_hydrogen):
            atoms = [(1, "CA", "C", 0.), (1, "H", "H", first_hydrogen),
                     (2, "CA", "C", 10.), (2, "H", "H", 12.)]
            return "".join(
                f"ATOM  {serial:5d} {name:^4s} ALA A{residue:4d}    "
                f"{coordinate:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00          {element:>2s}\n"
                for serial, (residue, name, element, coordinate) in enumerate(atoms, 1))
        self.cmd.read_pdbstr(
            "MODEL        1\n" + atom_frame(6.) + "ENDMDL\n"
            "MODEL        2\n" + atom_frame(12.) + "ENDMDL\nEND\n", "hydrogen_test")
        try:
            self.window.refresh_objects()
            self.window.object_combo.setCurrentText("hydrogen_test")
            self.select_mode("frame_centers")
            self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 8.)
            self.select_mode("both")
            self.window.distance_metric.setCurrentIndex(self.window.distance_metric.findData("centers"))
            self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 6.5)
            self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 1.5)
            for state, expected in ((1, 8.), (2, 5.)):
                panel = self.window._frame_export_panel(
                    "hydrogen_test", state, "both", "centers", 8., 0)
                self.assertAlmostEqual(panel["matrix"][0, 1], expected)
        finally:
            self.cmd.delete("hydrogen_test")

    def test_range_colors_and_column_profiles(self):
        matrix = np.array([[0., 2., 8.], [2., 0., 5.], [8., 5., 0.]])
        for mode in ("log", "linear"):
            pixels, minimum, maximum = self.plugin._viridis_pixels(matrix, mode, 6., 3.)
            self.assertEqual(tuple(pixels[0, 1]), (80, 80, 80))
            self.assertEqual(tuple(pixels[0, 2]), (210, 210, 210))
            self.assertEqual((minimum, maximum), (3., 6.))
        np.testing.assert_allclose(
            self.plugin._column_profile(matrix, 3., 6.), [np.nan, 5., 5.], equal_nan=True)
        pixels, minimum, maximum = self.plugin._viridis_pixels(np.zeros((2, 2)))
        self.assertEqual((minimum, maximum), (0., 0.))
        with self.assertRaises(ValueError):
            self.plugin._viridis_pixels(matrix, "linear", 3., 6.)

    def test_unified_visibility_default_and_controls(self):
        self.assertEqual(self.window.separation_spin.value(), 0)
        self.assertEqual([self.window.plot_mode.itemData(index)
                          for index in range(self.window.plot_mode.count())],
                         ["contacts", "frame_centers", "both", "split"])
        self.select_mode("both")
        self.assertEqual(len(self.window.canvas._panels), 2)
        self.window.show_std_check.setChecked(False)
        self.assertEqual([panel["statistic"] for panel in self.window.canvas._panels], ["mean"])
        self.window.show_mean_check.setChecked(False)
        self.assertEqual(self.window.canvas._panels, [])
        self.window.show_std_check.setChecked(True)
        self.assertEqual([panel["statistic"] for panel in self.window.canvas._panels], ["std"])

    def test_independent_ranges_profiles_and_split_images(self):
        self.select_mode("both")
        controls = self.window.mean_split_controls
        controls["auto"].setChecked(False)
        controls["maximum"].setValue(15.)
        controls["minimum"].setValue(10.)
        self.window.std_split_controls["linear"].setChecked(True)
        self.window.std_split_controls["auto"].setChecked(False)
        self.window.std_split_controls["maximum"].setValue(1.)
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        panels = self.window.canvas._panels
        np.testing.assert_allclose(panels[0]["profile"], [12., 12., np.nan], equal_nan=True)
        np.testing.assert_allclose(panels[1]["profile"], [np.nan, np.nan, np.nan], equal_nan=True)
        self.select_mode("split")
        panel = self.window.canvas._panels[0]
        self.assertTrue(panel["show_profile"])
        self.assertTrue(panel["show_right_profile"])
        self.assertEqual(panel["mean_scale_mode"], "log")
        self.assertEqual(panel["std_scale_mode"], "linear")
        self.assertEqual((panel["minimum"], panel["maximum"], panel["std_maximum"]),
                         (10., 15., 1.))
        self.assertEqual(panel["image"].pixelColor(2, 0).getRgb()[:3], (210, 210, 210))
        self.assertEqual(panel["image"].pixelColor(0, 1).getRgb()[:3], (210, 210, 210))
        pixmap = QtGui.QPixmap(self.window.canvas.size())
        self.window.canvas.render(pixmap)
        self.assertFalse(pixmap.isNull())

    def test_csv_stays_raw_and_image_includes_profiles(self):
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.mean_split_controls["auto"].setChecked(False)
        self.window.mean_split_controls["maximum"].setValue(1.)
        with tempfile.TemporaryDirectory() as directory:
            output = str(Path(directory) / "matrix.csv")
            with patch.object(QtWidgets.QFileDialog, "getSaveFileName", return_value=(output, "CSV")):
                self.window.export_raw_csv()
            with open(output, newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(rows[0][0], "residue")
            self.assertEqual(float(rows[1][2]), 12.)
            self.window.csv_statistic.setCurrentIndex(1)
            with patch.object(QtWidgets.QFileDialog, "getSaveFileName", return_value=(output, "CSV")):
                self.window.export_raw_csv()
            with open(output, newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            self.assertEqual(float(rows[1][2]), 2.)
            image_path = str(Path(directory) / "profiles.png")
            with patch.object(QtWidgets.QFileDialog, "getSaveFileName", return_value=(image_path, "PNG")):
                self.window.save_plot()
            self.assertFalse(QtGui.QImage(image_path).isNull())

    def test_line_plot_controls_visible_and_cutoff_spacing(self):
        self.window.show()
        self.select_mode("both")
        self.application.processEvents()
        self.assertTrue(self.window.line_plot_controls.isVisible())
        for checkbox in (self.window.profile_check, self.window.std_profile_check):
            self.assertTrue(checkbox.isVisible())
            self.assertTrue(checkbox.isEnabled())
            self.assertFalse(checkbox.isChecked())
            checkbox.click()
        panels = self.window.canvas._panels
        self.assertTrue(panels[0]["show_profile"])
        self.assertTrue(panels[1]["show_profile"])
        np.testing.assert_allclose(panels[0]["profile"], [16., 10., 14.])
        np.testing.assert_allclose(panels[1]["profile"], [1., 2., 1.])
        self.assertEqual(panels[0]["profile_ylabel"], "Mean of means (A)")
        self.assertEqual(panels[1]["profile_ylabel"], "Mean of SDs (A)")
        self.window.profile_check.click()
        self.assertFalse(self.window.canvas._panels[0]["show_profile"])
        self.assertTrue(self.window.canvas._panels[1]["show_profile"])
        cutoff_label = next(label for label in self.window.findChildren(QtWidgets.QLabel)
                            if label.text() == "Cutoff")
        gap = (self.window.cutoff_spin.geometry().left()
               - cutoff_label.geometry().right() - 1)
        self.assertGreaterEqual(gap, 0)
        self.assertLessEqual(gap, 10)
        self.select_mode("contacts")
        self.assertFalse(self.window.line_plot_controls.isVisible())

    def test_empty_profile_and_segment_identity(self):
        np.testing.assert_allclose(
            self.plugin._column_profile(np.zeros((1, 1))), [np.nan], equal_nan=True)
        atoms = [SimpleNamespace(segi=segment, chain="A", resi="1", resn="ALA",
                                 name="CA", coord=[coordinate, 0., 0.])
                 for segment, coordinate in (("first", 0.), ("second", 10.))]
        matrix, chains, labels, identities = self.plugin._residue_distances(atoms, "centers")
        self.assertEqual(matrix.shape, (2, 2))
        self.assertNotEqual(labels[0], labels[1])
        self.assertFalse(self.plugin._excluded_pairs(chains, 3)[0, 1])

    def test_csv_export_for_every_mode_and_dataset(self):
        expected_choices = {
            "contacts": {"contacts"},
            "frame_centers": {"frame_centers"},
            "both": {"mean", "std", "mean_profile", "std_profile"},
            "split": {"mean", "std", "split", "mean_profile", "std_profile"},
        }
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "chosen destination"
            destination.mkdir()
            for mode, keys in expected_choices.items():
                self.select_mode(mode)
                self.assertTrue(self.window.export_csv_button.isEnabled())
                self.assertTrue(self.window.csv_statistic.isEnabled())
                self.assertEqual(set(self.window._csv_datasets), keys)
                for key in keys:
                    with self.subTest(mode=mode, dataset=key):
                        self.window.csv_statistic.setCurrentIndex(
                            self.window.csv_statistic.findData(key))
                        target = destination / f"{mode}_{key}.csv"
                        with patch.object(QtWidgets.QFileDialog, "getSaveFileName",
                                          return_value=(str(target), "CSV")):
                            self.window.export_raw_csv()
                        with target.open(newline="", encoding="utf-8") as handle:
                            rows = list(csv.reader(handle))
                        expected = self.window._csv_datasets[key]["values"]
                        self.assertEqual(rows[0][0], "residue")
                        self.assertEqual(len(rows), len(expected) + 1)
                        if key.endswith("_profile"):
                            self.assertEqual(len(rows[0]), 2)
                            np.testing.assert_allclose(
                                [float(row[1]) if row[1] != "NA" else np.nan
                                 for row in rows[1:]], expected, equal_nan=True)
                        else:
                            np.testing.assert_allclose(
                                [[float(value) if value != "NA" else np.nan
                                  for value in row[1:]] for row in rows[1:]],
                                expected, equal_nan=True)

    def test_exported_profiles_honor_ranges_and_preserve_selection(self):
        self.select_mode("both")
        controls = self.window.mean_split_controls
        controls["auto"].setChecked(False)
        controls["maximum"].setValue(15.)
        controls["minimum"].setValue(10.)
        self.window.csv_statistic.setCurrentIndex(
            self.window.csv_statistic.findData("mean_profile"))
        controls["linear"].setChecked(True)
        self.assertEqual(self.window.csv_statistic.currentData(), "mean_profile")
        np.testing.assert_allclose(self.window._csv_datasets["mean_profile"]["values"],
                                   [12., 12., np.nan], equal_nan=True)
        self.assertEqual(self.window._csv_datasets["mean"]["values"][0, 2], 20.)
        self.select_mode("split")
        self.assertEqual(self.window.csv_statistic.currentData(), "mean_profile")
        np.testing.assert_allclose(self.window._csv_datasets["split"]["values"],
                                   [[0., 12., 20.], [2., 0., 8.], [0., 2., 0.]])

    def test_csv_cancellation_current_frame_and_invalid_plot(self):
        self.select_mode("contacts")
        self.cmd.frame(2)
        self.window.refresh_if_needed()
        self.assertIn("frame_0002", self.window._csv_datasets["contacts"]["suffix"])
        self.assertEqual(self.window._csv_datasets["contacts"]["values"][1, 2], 1.)
        with patch.object(QtWidgets.QFileDialog, "getSaveFileName", return_value=("", "")), \
                patch.object(self.plugin, "_write_matrix_csv") as writer:
            self.window.export_raw_csv()
            writer.assert_not_called()
        self.select_mode("both")
        self.window.mean_split_controls["auto"].setChecked(False)
        self.window.mean_split_controls["maximum"].setValue(5.)
        self.window.mean_split_controls["minimum"].setValue(10.)
        self.assertEqual(self.window.mean_split_controls["minimum"].value(), 5.)
        self.assertTrue(self.window._csv_datasets)
        self.assertTrue(self.window.export_csv_button.isEnabled())

    def test_frame_series_exports_matching_png_and_csv(self):
        for mode, expected in (("contacts", (0., 0.)),
                               ("frame_centers", (7., 11.)),
                               ("both", (7., 11.))):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as destination:
                self.select_mode(mode)
                self.window.distance_metric.setCurrentIndex(
                    self.window.distance_metric.findData("centers"))
                self.cmd.frame(2)
                self.window.refresh_if_needed()
                original_mode = self.window.plot_mode.currentData()
                original_timer = self.window.timer.isActive()
                result = self.window._export_frame_series(destination, 1, 2, 1)
                self.assertEqual(result["completed"], 2)
                self.assertFalse(result["canceled"])
                self.assertEqual(len(list(result["directory"].glob("*.png"))), 2)
                csv_paths = sorted(result["directory"].glob("*.csv"))
                self.assertEqual(len(csv_paths), 2)
                for index, path in enumerate(csv_paths):
                    with path.open(newline="", encoding="utf-8") as handle:
                        rows = list(csv.reader(handle))
                    self.assertAlmostEqual(float(rows[1][2]), expected[index])
                    self.assertFalse(QtGui.QImage(str(path.with_suffix(".png"))).isNull())
                self.assertEqual(self.cmd.get_frame(), 2)
                self.assertEqual(self.window.plot_mode.currentData(), original_mode)
                self.assertEqual(self.window.timer.isActive(), original_timer)
                self.assertFalse(self.window._exporting_frames)

    def test_frame_series_csv_only_stride_and_unique_folder(self):
        self.select_mode("frame_centers")
        with tempfile.TemporaryDirectory() as destination:
            first = self.window._export_frame_series(destination, 1, 2, 2, False, True)
            second = self.window._export_frame_series(destination, 1, 1, 1, False, True)
            self.assertNotEqual(first["directory"], second["directory"])
            self.assertEqual(first["completed"], 1)
            self.assertEqual(len(list(first["directory"].glob("*.csv"))), 1)
            self.assertEqual(len(list(first["directory"].glob("*.png"))), 0)
            with self.assertRaises(ValueError):
                self.window._export_frame_series(destination, 2, 1, 1)
            with self.assertRaises(ValueError):
                self.window._export_frame_series(destination, 1, 2, 1, False, False)

    def test_frame_series_cancellation_restores_session(self):
        self.select_mode("contacts")
        self.cmd.frame(2)
        with tempfile.TemporaryDirectory() as destination, \
                patch.object(QtWidgets.QProgressDialog, "wasCanceled", return_value=True):
            result = self.window._export_frame_series(destination, 1, 2, 1)
            self.assertTrue(result["canceled"])
            self.assertEqual(result["completed"], 0)
            self.assertEqual(list(result["directory"].iterdir()), [])
        self.assertEqual(self.cmd.get_frame(), 2)
        self.assertTrue(self.window.timer.isActive())
        self.assertFalse(self.window._exporting_frames)
        self.assertTrue(self.window.export_frames_button.isEnabled())

    def test_frame_series_error_removes_partial_pair_and_restores_session(self):
        self.select_mode("frame_centers")
        self.cmd.frame(2)
        with tempfile.TemporaryDirectory() as destination:
            with patch.object(self.plugin, "_write_matrix_csv", side_effect=OSError("test failure")):
                with self.assertRaisesRegex(OSError, "test failure"):
                    self.window._export_frame_series(destination, 1, 2, 1)
            folders = list(Path(destination).iterdir())
            self.assertEqual(len(folders), 1)
            self.assertEqual(list(folders[0].iterdir()), [])
        self.assertEqual(self.cmd.get_frame(), 2)
        self.assertTrue(self.window.timer.isActive())
        self.assertFalse(self.window._exporting_frames)

    def test_v4_shared_mean_filter_updates_sd_colors_and_profiles_live(self):
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        mean = self.window.mean_split_controls
        mean["auto"].setChecked(False)
        mean["maximum"].setValue(15.)
        mean["minimum"].setValue(10.)
        sd_panel = self.window.canvas._panels[1]
        self.assertEqual(sd_panel["image"].pixelColor(2, 0).getRgb()[:3], (210, 210, 210))
        self.assertEqual(sd_panel["image"].pixelColor(2, 1).getRgb()[:3], (80, 80, 80))
        np.testing.assert_allclose(sd_panel["profile"], [2., 2., np.nan], equal_nan=True)
        np.testing.assert_allclose(self.window._csv_datasets["std_profile"]["values"],
                                   sd_panel["profile"], equal_nan=True)
        mean["minimum"].setValue(0.)
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [2., 2., 2.])
        mean["maximum"].setValue(25.)
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [1., 2., 1.])
        sd = self.window.std_split_controls
        sd["auto"].setChecked(False)
        sd["maximum"].setValue(1.)
        sd["minimum"].setValue(0.5)
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [np.nan] * 3, equal_nan=True)
        np.testing.assert_allclose(self.window._csv_datasets["std_profile"]["values"], [np.nan] * 3, equal_nan=True)
        self.assertEqual(self.window._raw_matrices["mean"][0, 2], 20.)

    def test_v4_colorbar_widgets_single_export_and_binary_cutoff(self):
        self.window.show()
        self.application.processEvents()
        visible_exports = [button for button in self.window.findChildren(QtWidgets.QPushButton)
                           if button.isVisible() and "Export" in button.text()]
        self.assertEqual(visible_exports, [self.window.export_button])
        self.assertFalse(self.window.save_button.isVisible())
        self.assertTrue(self.window.cutoff_spin.isVisible())
        self.assertEqual(self.window.cutoff_spin.parent(), self.window.canvas.scale_controls["frame"]["group"])
        slider = self.window.canvas.scale_controls["frame"]["slider"]
        slider.setValue(110)
        self.assertEqual(self.window.cutoff_spin.value(), 11.)
        self.assertEqual(self.window.canvas._panels[0]["matrix"][0, 1], 1.)
        image = self.window.canvas._panels[0]["image"]
        self.assertEqual({image.pixelColor(column, row).getRgb()[:3]
                          for row in range(image.height()) for column in range(image.width())},
                         {(68, 1, 84), (253, 231, 37)})
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        self.application.processEvents()
        for controls in (self.window.mean_split_controls, self.window.std_split_controls):
            self.assertTrue(controls["minimum"].isVisible())
            self.assertTrue(controls["maximum"].isVisible())
            self.assertEqual(controls["group"].parent(), self.window.canvas)
        self.assertGreaterEqual(self.window.canvas.minimumHeight(), 540)
        self.assertIs(self.window.plot_scroll.widget(), self.window.canvas)

    def test_v4_export_dialog_selected_csvs_and_frames(self):
        self.select_mode("both")
        problems = []
        with tempfile.TemporaryDirectory() as destination:
            def choose_csvs():
                dialog = self.application.activeModalWidget()
                try:
                    dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                    dialog.findChild(QtWidgets.QCheckBox, "export_png").setChecked(False)
                    for checkbox in dialog.findChildren(QtWidgets.QCheckBox):
                        if checkbox.objectName().startswith("export_csv_"):
                            checkbox.setChecked(checkbox.objectName() in {"export_csv_mean", "export_csv_std_profile"})
                    dialog.accept()
                except Exception as exc:
                    problems.append(exc)
                    dialog.reject()
            QtCore.QTimer.singleShot(0, choose_csvs)
            self.window.export_dialog()
            self.assertFalse(problems)
            paths = list(Path(destination).glob("*/*.csv"))
            self.assertEqual(len(paths), 2)
            self.assertTrue(any("mean_distance" in path.name for path in paths))
            self.assertTrue(any("std_residue_profile" in path.name for path in paths))
            def choose_frames():
                dialog = self.application.activeModalWidget()
                try:
                    dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                    dialog.findChild(QtWidgets.QComboBox, "export_scope").setCurrentIndex(1)
                    dialog.findChild(QtWidgets.QCheckBox, "export_png").setChecked(False)
                    dialog.findChild(QtWidgets.QSpinBox, "export_frame_first").setValue(2)
                    dialog.findChild(QtWidgets.QSpinBox, "export_frame_last").setValue(2)
                    dialog.accept()
                except Exception as exc:
                    problems.append(exc)
                    dialog.reject()
            QtCore.QTimer.singleShot(0, choose_frames)
            self.window.export_dialog()
            self.assertFalse(problems)
            mean_paths = list(Path(destination).glob("*/*mean_distance.csv"))
            values = []
            for path in mean_paths:
                with path.open(newline="", encoding="utf-8") as handle:
                    values.append(float(list(csv.reader(handle))[1][2]))
            self.assertIn(14., values)

    def test_colorbar_handle_drags_filter_sd_profile(self):
        standalone = self.plugin.ContactMapRangeSlider()
        standalone.resize(25, 100)
        standalone.set_range(0., 1.275, 1.275, "linear")
        changes = []
        standalone.on_change = lambda name, value: changes.append((name, value))
        standalone.active_handle = "minimum"
        standalone._drag(0)
        self.assertEqual(changes, [("minimum", 1.27)])
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        self.window.show()
        self.application.processEvents()
        for controls in (self.window.mean_split_controls, self.window.std_split_controls):
            controls["linear"].setChecked(True)
        self.application.processEvents()

        def drag(controls, name, fraction):
            slider = controls["range_slider"]
            self.assertTrue(slider.isVisible())
            start = QtCore.QPointF(10, slider.position(getattr(slider, name)))
            end = QtCore.QPointF(10, (1 - fraction) * (slider.height() - 1))
            for kind, position, button, buttons in (
                    (QtCore.QEvent.MouseButtonPress, start, QtCore.Qt.LeftButton, QtCore.Qt.LeftButton),
                    (QtCore.QEvent.MouseMove, end, QtCore.Qt.NoButton, QtCore.Qt.LeftButton),
                    (QtCore.QEvent.MouseButtonRelease, end, QtCore.Qt.LeftButton, QtCore.Qt.NoButton)):
                self.application.sendEvent(slider, QtGui.QMouseEvent(
                    kind, position, button, buttons, QtCore.Qt.NoModifier))
                self.application.processEvents()

        mean = self.window.mean_split_controls
        drag(mean, "maximum", 0.75)
        drag(mean, "minimum", 0.5)
        self.assertTrue(mean["auto"].isChecked())
        self.assertAlmostEqual(mean["maximum"].value(), 15., delta=0.2)
        self.assertAlmostEqual(mean["minimum"].value(), 10., delta=0.2)
        original = self.window.canvas._panels[1]["profile"].copy()
        sd = self.window.std_split_controls
        drag(sd, "maximum", 0.5)
        drag(sd, "minimum", 0.25)
        self.assertTrue(np.isfinite(original).any())
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [np.nan] * 3, equal_nan=True)
        np.testing.assert_allclose(self.window._csv_datasets["std_profile"]["values"], [np.nan] * 3, equal_nan=True)
        self.assertEqual(self.window.canvas._panels[1]["image"].pixelColor(1, 0).getRgb()[:3], (210, 210, 210))
        self.assertGreaterEqual(self.window.canvas.width(), 1200)
        self.assertGreaterEqual(self.window.canvas.height(), 800)
        self.select_mode("split")
        self.application.processEvents()
        panel = self.window.canvas._panels[0]
        plot = QtCore.QRectF(76, 88, 310, 310)
        chart = self.window.canvas._profile_rect(panel, plot)
        self.assertGreater(chart.width(), plot.width())
        self.assertGreaterEqual(chart.height(), 200)
        self.assertGreater(chart.top(), plot.bottom())
        self.window.profile_check.setChecked(False)
        chart = self.window.canvas._profile_rect(self.window.canvas._panels[0], plot)
        self.assertGreater(chart.top(), plot.bottom())
        self.assertGreaterEqual(chart.height(), 200)

    def test_bundle_export_nested_destination_and_failed_write_cleanup(self):
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        self.window.show()
        self.application.processEvents()
        controls = self.window.mean_split_controls
        geometry = controls["group"].geometry()
        panels = self.window.canvas._panels
        with tempfile.TemporaryDirectory() as destination:
            nested = Path(destination) / "chosen folder" / "nested output"
            with patch.object(self.window.canvas, "render", side_effect=AssertionError("live canvas rendered")):
                output = self.window._export_current_bundle(nested, True, list(self.window._csv_datasets))
            self.assertEqual(output.parent, nested)
            image = QtGui.QImage(str(next(output.glob("*.png"))))
            self.assertFalse(image.isNull())
            self.assertEqual(image.width(), self.window.canvas.width() * 2)
            self.assertEqual(image.height(), self.window.canvas.height() * 2)
            self.assertEqual(len(list(output.glob("*.csv"))), len(self.window._csv_datasets))
            self.assertIs(self.window.canvas._panels, panels)
            self.assertEqual(controls["group"].geometry(), geometry)
            self.assertTrue(controls["group"].isVisible())
            with patch.object(self.plugin, "_write_matrix_csv", side_effect=OSError("disk write failure")):
                with self.assertRaisesRegex(OSError, "disk write failure"):
                    self.window._export_current_bundle(nested, True, ["mean"])
            self.assertEqual(list(nested.iterdir()), [output])

    def test_v8_split_controls_clear_profiles(self):
        self.select_mode("split")
        self.window.show()
        for mean_visible, sd_visible in ((False, False), (True, False), (False, True), (True, True)):
            self.window.profile_check.setChecked(mean_visible)
            self.window.std_profile_check.setChecked(sd_visible)
            self.application.processEvents()
            canvas = self.window.canvas
            mean_controls = self.window.mean_split_controls["group"].geometry()
            sd_controls = self.window.std_split_controls["group"].geometry()
            self.assertFalse(mean_controls.intersects(sd_controls))
            self.assertEqual(mean_controls.top(), sd_controls.top())
            self.assertTrue(canvas.rect().contains(mean_controls))
            self.assertTrue(canvas.rect().contains(sd_controls))
            self.assertEqual(mean_controls.top(), 24)
            right = 434
            profile_count = int(mean_visible) + int(sd_visible)
            bottom = 72 + (300 if profile_count else 0)
            map_size = min(canvas.width() - 76 - right, canvas.height() - 88 - bottom, 340)
            self.assertLessEqual(mean_controls.height(), 370)
            self.assertLessEqual(sd_controls.height(), 370)
            plot = QtCore.QRectF(76, 88, map_size, map_size)
            charts = []
            for profile_index in range(profile_count):
                chart = canvas._profile_rect(canvas._panels[0], plot, profile_index)
                charts.append(chart)
                self.assertGreaterEqual(chart.height(), 200)
                self.assertGreaterEqual(chart.width(), 400 if profile_count == 2 else 700)
                self.assertGreater(chart.top(), plot.bottom())
                for controls in (mean_controls, sd_controls):
                    self.assertFalse(chart.intersects(QtCore.QRectF(controls)))
            if len(charts) == 2:
                self.assertEqual(charts[0].top(), charts[1].top())
                self.assertEqual(charts[0].height(), charts[1].height())
                self.assertGreaterEqual(charts[1].left() - charts[0].right(), 100)
                self.assertFalse(charts[0].intersects(charts[1]))
            with patch.object(canvas, "_paint_profile", wraps=canvas._paint_profile) as horizontal, \
                    patch.object(canvas, "_paint_right_profile", wraps=canvas._paint_right_profile) as sideways:
                canvas.render(QtGui.QPixmap(canvas.size()))
                self.assertEqual(horizontal.call_count, profile_count)
                self.assertFalse(sideways.called)
                expected_titles = (["Mean of means by residue"] if mean_visible else []) + \
                                  (["Mean of SDs by residue"] if sd_visible else [])
                self.assertEqual([call.args[1]["profile_title"] for call in horizontal.call_args_list], expected_titles)
        with tempfile.TemporaryDirectory() as destination:
            path = Path(destination) / "split_v8.png"
            self.window._save_canvas_image(path)
            image = QtGui.QImage(str(path))
            self.assertFalse(image.isNull())
            self.assertEqual(image.width(), self.window.canvas.width() * 2)
            self.assertEqual(image.height(), self.window.canvas.height() * 2)

    def test_cutoffs_preserve_data_range_and_gray_direction(self):
        matrix = np.array([[0., 5., 10.], [5., 0., 15.], [10., 15., 0.]])
        settings = {"mean": ("linear", 5., 12.), "std": ("linear", 1., 3.)}
        display = self.plugin._statistics_display(matrix, matrix / 5, settings)
        self.assertEqual(display["mean"]["range"], 15.)
        self.assertEqual(display["std"]["range"], 3.)
        self.assertEqual(tuple(display["mean"]["pixels"][0, 0]), (80, 80, 80))
        self.assertEqual(tuple(display["mean"]["pixels"][1, 2]), (210, 210, 210))
        settings["mean"] = ("linear", 4., 14.)
        self.assertEqual(self.plugin._statistics_display(matrix, matrix / 5, settings)["mean"]["range"], 15.)

    def test_v7_current_distance_metric_and_cutoffs(self):
        self.select_mode("frame_centers")
        self.assertIn("residue distances", self.window.plot_mode.currentText())
        self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 7.)
        self.window.frame_distance_metric.setCurrentIndex(self.window.frame_distance_metric.findData("ca"))
        self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 10.)
        controls = self.window.distance_controls
        controls["auto"].setChecked(False)
        controls["maximum"].setValue(15.)
        controls["minimum"].setValue(5.)
        self.assertEqual(self.window.canvas._panels[0]["color_range"], 20.)
        self.assertEqual(self.window.canvas._panels[0]["image"].pixelColor(2, 0).getRgb()[:3], (210, 210, 210))
        self.assertEqual(self.window.canvas._panels[0]["image"].pixelColor(0, 0).getRgb()[:3], (80, 80, 80))
        self.assertFalse(controls["range"].isReadOnly())
        with tempfile.TemporaryDirectory() as destination:
            result = self.window._export_frame_series(destination, 1, 1, 1)
            with next(result["directory"].glob("*.csv")).open(newline="", encoding="utf-8") as handle:
                self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), 10.)

    def test_v7_plugin_frame_scrubbing_and_statistics_range(self):
        self.select_mode("frame_centers")
        self.window.frame_slider.setValue(2)
        self.assertEqual(self.cmd.get_frame(), 2)
        self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 11.)
        self.window.current_frame_spin.setValue(1)
        self.assertEqual(self.window.frame_slider.value(), 1)
        self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 7.)
        self.select_mode("both")
        self.window.first_frame_spin.setValue(2)
        self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 14.)
        self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 0.)
        self.window.first_frame_spin.setValue(1)
        self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 12.)
        self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 2.)

    def test_v7_command_retains_statistics_frame_bounds(self):
        self.select_mode("both")
        self.window.first_frame_spin.setValue(2)
        settings = self.window._export_settings()
        with tempfile.TemporaryDirectory() as destination:
            results = self.plugin.contact_map_export(
                destination, objects=["distance_test"], mode="both", metric="ca", save_image=False,
                datasets=["mean", "std"], settings=settings)
            for suffix, expected in (("mean_distance", 14.), ("distance_standard_deviation", 0.)):
                with next(results[0]["directory"].glob("*" + suffix + ".csv")).open(newline="", encoding="utf-8") as handle:
                    self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), expected)
        self.assertEqual((self.window.first_frame_spin.value(), self.window.last_frame_spin.value()), (2, 2))

    def test_v7_export_range_prefill_and_checkbox_toggle(self):
        self.window.first_frame_spin.setValue(2)
        problems = []
        with tempfile.TemporaryDirectory() as destination:
            def choose_range():
                dialog = self.application.activeModalWidget()
                try:
                    self.assertEqual(dialog.findChild(QtWidgets.QComboBox, "export_scope").currentData(), "frames")
                    first = dialog.findChild(QtWidgets.QSpinBox, "export_frame_first")
                    last = dialog.findChild(QtWidgets.QSpinBox, "export_frame_last")
                    self.assertEqual((first.value(), last.value()), (2, 2))
                    first.setValue(1)
                    last.setValue(1)
                    self.assertEqual((first.value(), last.value()), (1, 1))
                    first.setValue(2)
                    self.assertEqual((first.value(), last.value()), (2, 2))
                    frame_csv = dialog.findChild(QtWidgets.QCheckBox, "export_frame_csv")
                    self.assertTrue(frame_csv.isEnabled())
                    frame_csv.click()
                    self.assertFalse(frame_csv.isChecked())
                    frame_csv.click()
                    self.assertTrue(frame_csv.isChecked())
                    objects = dialog.findChild(QtWidgets.QListWidget, "export_objects")
                    point = objects.visualItemRect(objects.item(0)).center()
                    for expected in (QtCore.Qt.Unchecked, QtCore.Qt.Checked):
                        for kind, buttons in ((QtCore.QEvent.MouseButtonPress, QtCore.Qt.LeftButton),
                                              (QtCore.QEvent.MouseButtonRelease, QtCore.Qt.NoButton)):
                            self.application.sendEvent(objects.viewport(), QtGui.QMouseEvent(
                                kind, QtCore.QPointF(point), QtCore.Qt.LeftButton, buttons, QtCore.Qt.NoModifier))
                        self.assertEqual(objects.item(0).checkState(), expected)
                    dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                    dialog.findChild(QtWidgets.QCheckBox, "export_png").setChecked(False)
                    dialog.accept()
                except Exception as exc:
                    problems.append(exc)
                    dialog.reject()
            QtCore.QTimer.singleShot(0, choose_range)
            self.window.export_range_button.click()
            self.assertFalse(problems)
            self.assertEqual(len(list(Path(destination).glob("*/*000002.csv"))), 1)
            self.assertEqual((self.window.first_frame_spin.value(), self.window.last_frame_spin.value()), (2, 2))

    def test_fixed_zero_to_range_colors_and_sd_profile_filters(self):
        mean = np.array([[0., 5., 10.], [5., 0., 15.], [10., 15., 0.]])
        sd = np.array([[0., 1., 2.], [1., 0., 3.], [2., 3., 0.]])
        settings = {"mean": ("linear", 0., 20.), "std": ("linear", 0., 4.),
                    "mean_range": 20., "std_range": 4.}
        before = self.plugin._statistics_display(mean, sd, settings)
        settings.update(mean=("linear", 3., 12.), std=("linear", 1.2, 2.5))
        after = self.plugin._statistics_display(mean, sd, settings)
        np.testing.assert_array_equal(before["mean"]["pixels"][0, 1], after["mean"]["pixels"][0, 1])
        expected = self.plugin._viridis_pixels(mean, "linear", 12., 3., 20.)[0]
        np.testing.assert_array_equal(expected[0, 1], after["mean"]["pixels"][0, 1])
        np.testing.assert_allclose(after["std"]["profile"], [2., np.nan, 2.], equal_nan=True)
        np.testing.assert_allclose((after["mean"]["range"], after["std"]["range"]), (20., 4.))
        self.select_mode("both")
        controls = self.window.mean_split_controls
        controls["auto"].setChecked(False)
        controls["maximum"].setValue(15.)
        controls["minimum"].setValue(5.)
        self.assertEqual(self.window.canvas._panels[0]["mean_range"], 20.)
        self.assertEqual(self.window.canvas._panels[0]["maximum"], 15.)
        self.window.show()
        self.application.processEvents()
        self.assertFalse(controls["range"].isReadOnly())
        self.assertEqual(controls["range"].value(), 20.)

    def test_frame_selection_ranges_exact_frames_and_validation(self):
        parse = self.plugin._parse_frame_selection
        self.assertEqual(parse("10-100", 250), list(range(10, 101)))
        self.assertEqual(parse("50,200", 250), [50, 200])
        self.assertEqual(parse("50", 250), [50])
        self.assertEqual(parse("10-100:10,50,200", 250), list(range(10, 101, 10)) + [200])
        self.assertEqual(parse([2, 1, 2], 2), [2, 1])
        for invalid in ("", "0", "251", "100-10", "10-100:0", "x", [1.5]):
            with self.subTest(selection=invalid), self.assertRaises(ValueError):
                parse(invalid, 250)

    def test_multi_object_exact_frame_exports_and_preflight(self):
        self.cmd.read_pdbstr(
            "MODEL        1\n" + pdb_frame(20., 18.) + "ENDMDL\n"
            "MODEL        2\n" + pdb_frame(30., 28.) + "ENDMDL\nEND\n", "distance_other")
        try:
            self.window.refresh_objects()
            self.select_mode("both")
            self.cmd.frame(2)
            original_aggregate = self.window._aggregate_maps
            def no_other_aggregation(name, *args):
                self.assertEqual(name, "distance_test")
                return original_aggregate(name, *args)
            with tempfile.TemporaryDirectory() as destination, \
                    patch.object(self.window, "_aggregate_maps", side_effect=no_other_aggregation):
                results = self.window._export_objects(destination, ["distance_test", "distance_other"],
                                                      "2,1", True, True, image_format="jpg")
                self.assertEqual([result["completed"] for result in results], [2, 2])
                for result, expected in zip(results, (14., 30.)):
                    self.assertEqual(result["directory"].parent, Path(destination))
                    self.assertEqual(len(list(result["directory"].glob("*.jpg"))), 2)
                    self.assertEqual(len(list(result["directory"].glob("*.csv"))), 2)
                    image = next(result["directory"].glob("*.jpg"))
                    self.assertFalse(QtGui.QImage(str(image)).isNull())
                    with next(result["directory"].glob("*000002.csv")).open(newline="", encoding="utf-8") as handle:
                        self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), expected)
            self.assertEqual(self.window.object_combo.currentText(), "distance_test")
            self.assertEqual(self.cmd.get_frame(), 2)
            self.assertTrue(self.window.timer.isActive())
            self.cmd.delete_states("distance_other", "2")
            with tempfile.TemporaryDirectory() as destination:
                with self.assertRaisesRegex(ValueError, "distance_other"):
                    self.window._export_objects(destination, ["distance_test", "distance_other"], "2")
                self.assertEqual(list(Path(destination).iterdir()), [])
        finally:
            self.cmd.delete("distance_other")

    def test_export_dialog_frame_slider_and_executable_copy_command(self):
        self.cmd.read_pdbstr(
            "MODEL        1\n" + pdb_frame(20., 18.) + "ENDMDL\n"
            "MODEL        2\n" + pdb_frame(30., 28.) + "ENDMDL\nEND\n", "distance_other")
        commands, problems = [], []
        try:
            self.window.refresh_objects()
            self.select_mode("frame_centers")
            self.cmd.frame(2)
            with tempfile.TemporaryDirectory() as temporary:
                destination = Path(temporary) / "chosen folder, sample's" / "nested"
                def copy_options():
                    dialog = self.application.activeModalWidget()
                    try:
                        dialog.findChild(QtWidgets.QComboBox, "export_scope").setCurrentIndex(1)
                        objects = dialog.findChild(QtWidgets.QListWidget, "export_objects")
                        for index in range(objects.count()):
                            objects.item(index).setCheckState(QtCore.Qt.Checked)
                        dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(str(destination))
                        dialog.findChild(QtWidgets.QCheckBox, "export_png").setChecked(False)
                        slider = dialog.findChild(QtWidgets.QWidget, "export_frame_slider")
                        self.assertTrue(slider.isEnabled())
                        slider.active_handle = "minimum"
                        slider._drag(slider.width() - 1)
                        slider.active_handle = None
                        self.assertEqual(dialog.findChild(QtWidgets.QSpinBox, "export_frame_first").value(), 2)
                        dialog.findChild(QtWidgets.QRadioButton, "export_specific_option").click()
                        dialog.findChild(QtWidgets.QLineEdit, "export_exact_frames").setText("2,1")
                        command = dialog.findChild(QtWidgets.QLineEdit, "export_command")
                        self.assertTrue(command.isReadOnly())
                        commands.append(command.text())
                        copied = []
                        with patch.object(QtWidgets.QApplication, "clipboard", return_value=SimpleNamespace(setText=copied.append)):
                            next(button for button in dialog.findChildren(QtWidgets.QPushButton) if button.text() == "Copy").click()
                        self.assertEqual(copied, commands)
                    except Exception as exc:
                        problems.append(exc)
                    finally:
                        if dialog is not None:
                            dialog.reject()
                QtCore.QTimer.singleShot(0, copy_options)
                self.window.export_dialog()
                self.assertFalse(problems)
                self.assertEqual(len(commands), 1)
                self.assertTrue(commands[0].startswith("/"))
                compile(commands[0][1:], "<export command>", "exec")
                original_factory = self.plugin.ContactMapWindow
                def frame_export_window():
                    window = original_factory()
                    def unexpected_aggregate(*unused):
                        error = AssertionError("Frame-only command attempted temporal aggregation")
                        problems.append(error)
                        raise error
                    window._aggregate_maps = unexpected_aggregate
                    return window
                with patch.object(self.plugin, "ContactMapWindow", side_effect=frame_export_window):
                    self.cmd.do(commands[0], flush=1)
                self.assertFalse(problems)
                paths = list(destination.glob("*/*.csv"))
                self.assertEqual(len(paths), 4)
                self.assertEqual(self.cmd.get_frame(), 2)
                self.assertEqual(self.window.object_combo.currentText(), "distance_test")
        finally:
            self.cmd.delete("distance_other")

    def test_v4_current_export_retains_jpeg_support(self):
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        with tempfile.TemporaryDirectory() as destination:
            output = self.window._export_current_bundle(destination, True, ["std"], image_format="jpg")
            images = list(output.glob("*.jpg"))
            self.assertEqual(len(images), 1)
            self.assertFalse(QtGui.QImage(str(images[0])).isNull())
            self.assertEqual(len(list(output.glob("*.csv"))), 1)


    def test_v8_exclusive_frame_selection_modes(self):
        for selection, expected in (("range", [1, 2]), ("specific", [1, 2]), ("single", [2])):
            problems = []
            with self.subTest(selection=selection), tempfile.TemporaryDirectory() as destination:
                def choose_selection():
                    dialog = self.application.activeModalWidget()
                    try:
                        range_option = dialog.findChild(QtWidgets.QRadioButton, "export_range_option")
                        specific_option = dialog.findChild(QtWidgets.QRadioButton, "export_specific_option")
                        specific = dialog.findChild(QtWidgets.QLineEdit, "export_exact_frames")
                        slider = dialog.findChild(QtWidgets.QWidget, "export_frame_slider")
                        bounds = [dialog.findChild(QtWidgets.QSpinBox, "export_frame_" + name)
                                  for name in ("first", "last", "interval")]
                        command = dialog.findChild(QtWidgets.QLineEdit, "export_command")
                        self.assertTrue(range_option.isChecked())
                        self.assertFalse(specific.isEnabled())
                        self.assertTrue(slider.isEnabled())
                        specific_option.click()
                        self.assertTrue(specific.isEnabled())
                        self.assertFalse(slider.isEnabled())
                        self.assertTrue(all(not bound.isEnabled() for bound in bounds))
                        specific.setText("2" if selection == "single" else "2,1")
                        self.assertIn("frames=" + repr(specific.text()), command.text())
                        range_option.click()
                        self.assertFalse(specific.isEnabled())
                        self.assertTrue(slider.isEnabled())
                        self.assertTrue(all(bound.isEnabled() for bound in bounds))
                        self.assertIn("frames='1-2:1'", command.text())
                        if selection != "range":
                            specific_option.click()
                            self.assertEqual(specific.text(), "2" if selection == "single" else "2,1")
                        dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                        dialog.findChild(QtWidgets.QCheckBox, "export_png").setChecked(False)
                        dialog.accept()
                    except Exception as exc:
                        problems.append(exc)
                        if dialog is not None:
                            dialog.reject()
                QtCore.QTimer.singleShot(0, choose_selection)
                self.window.export_dialog()
                self.assertFalse(problems)
                frames = sorted(int(path.stem.rsplit("_", 1)[1]) for path in Path(destination).glob("*/*.csv"))
                self.assertEqual(frames, expected)
            self.assertEqual(self.plugin.__version__, "v9")


    def test_v8_playback_current_frame_views_and_looping(self):
        for mode, expected in (("contacts", 0.), ("frame_centers", 11.)):
            with self.subTest(mode=mode):
                self.select_mode(mode)
                self.window.current_frame_spin.setValue(1)
                self.window.cutoff_spin.setValue(11.)
                self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1],
                                       1. if mode == "contacts" else 7.)
                self.window.first_frame_spin.setValue(2)
                self.assertTrue(self.window.play_button.isEnabled())
                self.window.play_button.click()
                self.assertTrue(self.window.play_button.isChecked())
                self.assertEqual(self.window.play_button.text(), "Pause")
                self.window.timer.timeout.emit()
                self.assertEqual(self.cmd.get_frame(), 2)
                self.assertEqual(self.window.frame_slider.value(), 2)
                self.assertEqual(self.window.current_frame_spin.value(), 2)
                self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], expected)
                self.window.timer.timeout.emit()
                self.assertEqual(self.cmd.get_frame(), 1)
                self.assertEqual(self.window.first_frame_spin.value(), 2)
                self.window.play_button.click()
                self.assertFalse(self.window.play_button.isChecked())
                self.assertEqual(self.window.play_button.text(), "Play")
                self.window.timer.timeout.emit()
                self.assertEqual(self.cmd.get_frame(), 1)
                self.window.play_button.click()
                self.window.frame_slider.setValue(2)
                self.assertFalse(self.window.play_button.isChecked())

    def test_v8_playback_disabled_for_statistics_and_single_state(self):
        self.window.play_button.click()
        self.assertTrue(self.window.play_button.isChecked())
        for mode in ("both", "split"):
            self.select_mode(mode)
            self.assertFalse(self.window.play_button.isChecked())
            self.assertFalse(self.window.play_button.isEnabled())
        self.select_mode("contacts")
        self.cmd.create("single_state", "distance_test", 1, 1)
        try:
            self.window.refresh_objects()
            self.window.play_button.click()
            self.window.object_combo.setCurrentText("single_state")
            self.assertFalse(self.window.play_button.isChecked())
            self.assertFalse(self.window.play_button.isEnabled())
            self.window.object_combo.setCurrentText("distance_test")
            self.window.play_button.click()
            self.window.close()
            self.assertFalse(self.window.timer.isActive())
            self.assertFalse(self.window.play_button.isChecked())
        finally:
            self.cmd.delete("single_state")


    def test_v8_playback_pauses_while_export_dialog_is_open(self):
        self.window.play_button.click()
        original_frame = self.cmd.get_frame()
        problems = []
        def cancel_export():
            dialog = self.application.activeModalWidget()
            try:
                self.assertFalse(self.window.timer.isActive())
                self.assertTrue(self.window.play_button.isChecked())
                self.application.processEvents()
                self.assertEqual(self.cmd.get_frame(), original_frame)
            except Exception as exc:
                problems.append(exc)
            finally:
                if dialog is not None:
                    dialog.reject()
        QtCore.QTimer.singleShot(0, cancel_export)
        self.window.export_dialog()
        self.assertFalse(problems)
        self.assertTrue(self.window.timer.isActive())
        self.assertTrue(self.window.play_button.isChecked())
        self.window.timer.timeout.emit()
        self.assertEqual(self.cmd.get_frame(), original_frame % 2 + 1)
        self.window.play_button.click()
        self.assertEqual(self.plugin.__version__, "v9")


    def test_fixed_palette_range_does_not_follow_cutoffs(self):
        matrix = np.array([[0., 5., 10.], [5., 0., 15.], [10., 15., 0.]])
        for mode in ("linear", "log"):
            settings = {"mean": (mode, 0., 20.), "std": (mode, 0., 4.),
                        "mean_range": 20., "std_range": 4.}
            before = self.plugin._statistics_display(matrix, matrix / 5, settings)
            settings["mean"] = (mode, 3., 12.)
            settings["std"] = (mode, 1.2, 2.5)
            after = self.plugin._statistics_display(matrix, matrix / 5, settings)
            self.assertEqual((after["mean"]["range"], after["std"]["range"]), (20., 4.))
            np.testing.assert_array_equal(before["mean"]["pixels"][0, 1], after["mean"]["pixels"][0, 1])
            np.testing.assert_array_equal(before["std"]["pixels"][0, 2], after["std"]["pixels"][0, 2])
            self.assertEqual(tuple(after["mean"]["pixels"][0, 0]), (80, 80, 80))
            self.assertEqual(tuple(after["mean"]["pixels"][1, 2]), (210, 210, 210))


    def test_fixed_palette_survives_ui_changes_and_frame_exports(self):
        for mode in ("both", "frame_centers"):
            with self.subTest(mode=mode):
                self.select_mode(mode)
                controls = self.window.mean_split_controls if mode == "both" else self.window.distance_controls
                controls["range"].setValue(30.)
                controls["auto"].setChecked(False)
                controls["maximum"].setValue(25.)
                controls["minimum"].setValue(0.)
                original_color = self.window.canvas._panels[0]["image"].pixelColor(1, 0).getRgb()
                controls["maximum"].setValue(15.)
                controls["minimum"].setValue(5.)
                self.assertEqual(self.window.canvas._panels[0]["image"].pixelColor(1, 0).getRgb(), original_color)
                self.assertEqual(controls["range"].value(), 30.)
                controls["maximum"].setValue(25.)
                static_pair_color = self.window.canvas._panels[0]["image"].pixelColor(2, 0).getRgb()
                if mode == "both":
                    self.window.first_frame_spin.setValue(2)
                else:
                    self.window.current_frame_spin.setValue(2)
                self.assertEqual(controls["range"].value(), 30.)
                self.assertEqual(self.window.canvas._panels[0]["image"].pixelColor(2, 0).getRgb(), static_pair_color)
                saved_ranges = []
                save_image = self.window._save_canvas_image
                def capture_range(path):
                    saved_ranges.append(self.window.canvas._panels[0]["color_range"])
                    save_image(path)
                with tempfile.TemporaryDirectory() as destination, \
                        patch.object(self.window, "_save_canvas_image", side_effect=capture_range):
                    result = self.window._export_frame_series(destination, 1, 2, 1, True, False)
                    self.assertEqual(len(list(result["directory"].glob("*.png"))), 2)
                self.assertEqual(saved_ranges, [30., 30.])
                self.assertEqual(controls["range"].value(), 30.)
            self.assertEqual(self.plugin.__version__, "v9")


    def test_auto_range_never_changes_upper_cutoff(self):
        for mode, name in (("both", "mean"), ("both", "std"), ("frame_centers", "distance")):
            with self.subTest(mode=mode, statistic=name):
                self.select_mode(mode)
                controls = self.window.canvas.scale_controls[name]
                controls["auto"].setChecked(False)
                controls["range"].setValue(30.)
                controls["maximum"].setValue(1. if name == "std" else 15.)
                upper = controls["maximum"].value()
                controls["auto"].setChecked(True)
                self.assertEqual(controls["auto"].text(), "Auto range")
                self.assertEqual(controls["maximum"].value(), upper)
                self.assertTrue(controls["maximum"].isEnabled())
                self.assertFalse(controls["range"].isEnabled())
                original_range = controls["range"].value()
                controls["range_slider"].on_change("maximum", upper / 2)
                self.assertTrue(controls["auto"].isChecked())
                self.assertEqual(controls["maximum"].value(), upper / 2)
                self.assertEqual(controls["range"].value(), original_range)
                controls["auto"].setChecked(False)
                self.assertTrue(controls["range"].isEnabled())
                self.assertEqual(controls["maximum"].value(), upper / 2)


    def test_auto_range_tracks_data_but_preserves_cutoff_and_export_settings(self):
        self.select_mode("both")
        controls = self.window.std_split_controls
        controls["maximum"].setValue(0.5)
        self.assertTrue(controls["auto"].isChecked())
        self.assertEqual(controls["range"].value(), 2.)
        self.window.first_frame_spin.setValue(2)
        self.assertEqual(controls["range"].value(), 0.)
        self.assertEqual(controls["maximum"].value(), 0.)
        self.window.first_frame_spin.setValue(1)
        self.assertEqual(controls["range"].value(), 2.)
        self.assertEqual(controls["maximum"].value(), 0.)
        settings = self.window._export_settings()
        self.assertTrue(settings["statistics"]["std_auto_range"])
        self.assertEqual(settings["statistics"]["std"][2], 0.)
        other = self.plugin.ContactMapWindow()
        try:
            other._apply_export_settings("both", "ca", settings)
            self.assertTrue(other.std_split_controls["auto"].isChecked())
            self.assertEqual(other.std_split_controls["maximum"].value(), 0.)
            self.assertEqual(other.std_split_controls["range"].value(), 2.)
        finally:
            other.close()


    def test_v9_navigation_and_cutoff_bounds(self):
        self.window.show()
        for mode in ("contacts", "frame_centers", "both", "split"):
            self.select_mode(mode)
            self.application.processEvents()
            statistics = mode in {"both", "split"}
            self.assertEqual(self.window.range_controls.isVisible(), statistics)
            self.assertEqual(self.window.frame_controls.isVisible(), not statistics)
        for controls in (self.window.mean_split_controls, self.window.std_split_controls, self.window.distance_controls):
            controls["auto"].setChecked(False)
            controls["range"].setValue(30.)
            controls["maximum"].setValue(25.)
            controls["minimum"].setValue(20.)
            controls["range"].setValue(10.)
            self.assertEqual(controls["maximum"].value(), 10.)
            self.assertEqual(controls["minimum"].value(), 10.)
            controls["minimum"].setValue(15.)
            self.assertEqual(controls["minimum"].value(), 10.)
        self.assertEqual(self.plugin.__version__, "v9")


    def test_v9_statistics_selection_updates_all_datasets(self):
        self.select_mode("both")
        self.window._set_statistics_selection([2])
        self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 14.)
        self.assertAlmostEqual(self.window._csv_datasets["std"]["values"][0, 1], 0.)
        np.testing.assert_allclose(self.window._csv_datasets["mean_profile"]["values"], [17., 10., 13.])
        self.window._set_statistics_selection([2, 1])
        self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 12.)
        self.assertAlmostEqual(self.window._csv_datasets["std"]["values"][0, 1], 2.)


    def test_v9_independent_statistics_images_and_range_csvs(self):
        self.select_mode("both")
        with tempfile.TemporaryDirectory() as destination:
            results = self.window._export_objects(destination, ["distance_test"], "2", True, True,
                ["mean", "std", "mean_profile", "std_profile"], aggregate_frames=True,
                images=["mean", "std", "mean_profile", "std_profile"])
            output = results[0]["directory"]
            self.assertEqual(len(list(output.glob("*.png"))), 4)
            self.assertEqual(len(list(output.glob("*.csv"))), 4)
            for path in output.glob("*.png"):
                self.assertFalse(QtGui.QImage(str(path)).isNull())
            heatmap = QtGui.QImage(str(next(output.glob("*_mean.png"))))
            profile = QtGui.QImage(str(next(output.glob("*_mean_profile.png"))))
            self.assertEqual(heatmap.width(), 1920)
            self.assertEqual(profile.width(), 1920)
            self.assertEqual(heatmap.pixelColor(160, 180).getRgb()[:3], (68, 1, 84))
            self.assertNotEqual(profile.pixelColor(160, 180).getRgb()[:3], (68, 1, 84))
            with next(output.glob("*mean_distance.csv")).open(newline="", encoding="utf-8") as handle:
                self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), 14.)
            with next(output.glob("*distance_standard_deviation.csv")).open(newline="", encoding="utf-8") as handle:
                self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), 0.)


    def test_v9_export_dialog_previews_and_exports_statistics(self):
        self.window.show()
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        self.window.std_profile_check.setChecked(True)
        problems = []
        commands = []
        with tempfile.TemporaryDirectory() as destination:
            def choose_statistics():
                dialog = self.application.activeModalWidget()
                try:
                    scope = dialog.findChild(QtWidgets.QComboBox, "export_scope")
                    self.assertEqual(scope.currentData(), "frames")
                    dialog.findChild(QtWidgets.QSpinBox, "export_frame_first").setValue(2)
                    preview = dialog.findChild(QtCore.QTimer, "export_preview_timer")
                    preview.stop()
                    preview.timeout.emit()
                    self.assertAlmostEqual(self.window.canvas._panels[0]["matrix"][0, 1], 14.)
                    self.assertAlmostEqual(self.window.canvas._panels[1]["matrix"][0, 1], 0.)
                    np.testing.assert_allclose(self.window._csv_datasets["mean_profile"]["values"], [17., 10., 13.])
                    np.testing.assert_allclose(self.window.canvas._panels[0]["profile"], [17., 10., 13.])
                    np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [0., 0., 0.])
                    tabs = dialog.findChild(QtWidgets.QTabWidget, "export_preview_tabs")
                    self.assertIsNotNone(tabs)
                    mean_preview = dialog.findChild(QtWidgets.QWidget, "export_preview_mean")
                    self.assertAlmostEqual(mean_preview._panels[0]["matrix"][0, 1], 14.)
                    tabs.setCurrentIndex(1)
                    sd_preview = dialog.findChild(QtWidgets.QWidget, "export_preview_std")
                    self.assertAlmostEqual(sd_preview._panels[0]["matrix"][0, 1], 0.)
                    for checkbox in dialog.findChildren(QtWidgets.QCheckBox):
                        if checkbox.objectName().startswith("export_image_"):
                            checkbox.setChecked(checkbox.objectName() in {"export_image_std", "export_image_mean_profile"})
                        if checkbox.objectName().startswith("export_csv_"):
                            self.assertTrue(checkbox.isVisible())
                            checkbox.setChecked(True)
                    dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                    command = dialog.findChild(QtWidgets.QLineEdit, "export_command").text()
                    self.assertIn("aggregate_frames=True", command)
                    self.assertIn("images=['std', 'mean_profile']", command)
                    commands.append(command)
                    dialog.accept()
                except Exception as exc:
                    problems.append(exc)
                    if dialog is not None:
                        dialog.reject()
            QtCore.QTimer.singleShot(0, choose_statistics)
            with patch.object(QtWidgets.QMessageBox, "warning") as warning:
                self.window.export_dialog()
                warning.assert_not_called()
            self.assertFalse(problems)
            self.assertEqual(len(list(Path(destination).glob("*/*.png"))), 2)
            self.assertEqual(len(list(Path(destination).glob("*/*.csv"))), 4)
            with next(Path(destination).glob("*/*mean_distance.csv")).open(newline="", encoding="utf-8") as handle:
                self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), 14.)
            with next(Path(destination).glob("*/*mean_residue_profile.csv")).open(newline="", encoding="utf-8") as handle:
                np.testing.assert_allclose([float(row[1]) for row in list(csv.reader(handle))[1:]], [17., 10., 13.])
            self.cmd.do(commands[0], flush=1)
            self.assertEqual(len(list(Path(destination).glob("*/*.png"))), 4)
            self.assertEqual(len(list(Path(destination).glob("*/*.csv"))), 8)


    def test_v9_frame_edit_during_statistics_calculation_is_recomputed(self):
        self.select_mode("both")
        self.window._aggregate_cache.clear()
        get_model = self.cmd.get_model
        changed = []
        def edit_range(*args, **kwargs):
            if not changed:
                changed.append(True)
                self.window.first_frame_spin.setValue(2)
            return get_model(*args, **kwargs)
        with patch.object(self.cmd, "get_model", side_effect=edit_range):
            self.window.update_contact_map()
            self.application.processEvents()
        self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 14.)
        self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 0.)


    def test_v9_specific_statistics_frames_exclude_intervening_states(self):
        self.cmd.read_pdbstr(
            "MODEL        1\n" + pdb_frame(10., 6.) + "ENDMDL\n"
            "MODEL        2\n" + pdb_frame(100., 96.) + "ENDMDL\n"
            "MODEL        3\n" + pdb_frame(18., 14.) + "ENDMDL\nEND\n", "selection_test")
        try:
            self.window.refresh_objects()
            self.window.object_combo.setCurrentText("selection_test")
            self.select_mode("both")
            self.window._set_statistics_selection([1, 3])
            self.assertAlmostEqual(self.window._raw_matrices["mean"][0, 1], 14.)
            self.assertAlmostEqual(self.window._raw_matrices["std"][0, 1], 4.)
            with tempfile.TemporaryDirectory() as destination:
                result = self.window._export_objects(destination, ["selection_test"], "1,3", False, True,
                                                     ["mean", "std"], aggregate_frames=True)
                for suffix, expected in (("mean_distance", 14.), ("distance_standard_deviation", 4.)):
                    with next(result[0]["directory"].glob("*" + suffix + ".csv")).open(newline="", encoding="utf-8") as handle:
                        self.assertAlmostEqual(float(list(csv.reader(handle))[1][2]), expected)
        finally:
            self.window.object_combo.setCurrentText("distance_test")
            self.cmd.delete("selection_test")


    def test_v9_split_range_preview_and_selected_split_outputs(self):
        self.select_mode("split")
        problems = []
        with tempfile.TemporaryDirectory() as destination:
            def choose_split():
                dialog = self.application.activeModalWidget()
                try:
                    dialog.findChild(QtWidgets.QSpinBox, "export_frame_first").setValue(2)
                    timer = dialog.findChild(QtCore.QTimer, "export_preview_timer")
                    timer.stop()
                    timer.timeout.emit()
                    tabs = dialog.findChild(QtWidgets.QTabWidget, "export_preview_tabs")
                    tabs.setCurrentIndex(tabs.count() - 1)
                    preview = dialog.findChild(QtWidgets.QWidget, "export_preview_split")
                    self.assertGreaterEqual(preview.minimumWidth(), 960)
                    self.assertAlmostEqual(preview._panels[0]["mean_matrix"][0, 1], 14.)
                    self.assertAlmostEqual(preview._panels[0]["std_matrix"][0, 1], 0.)
                    for checkbox in dialog.findChildren(QtWidgets.QCheckBox):
                        if checkbox.objectName().startswith("export_image_"):
                            checkbox.setChecked(checkbox.objectName() == "export_image_split")
                        if checkbox.objectName().startswith("export_csv_"):
                            checkbox.setChecked(checkbox.objectName() in {"export_csv_split", "export_csv_std_profile"})
                    dialog.findChild(QtWidgets.QLineEdit, "export_destination").setText(destination)
                    dialog.accept()
                except Exception as exc:
                    problems.append(exc)
                    if dialog is not None:
                        dialog.reject()
            QtCore.QTimer.singleShot(0, choose_split)
            with patch.object(QtWidgets.QMessageBox, "warning") as warning:
                self.window.export_dialog()
                warning.assert_not_called()
            self.assertFalse(problems)
            self.assertEqual(len(list(Path(destination).glob("*/*.png"))), 1)
            self.assertEqual(len(list(Path(destination).glob("*/*.csv"))), 2)
            with next(Path(destination).glob("*/*split_mean_sd.csv")).open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
                self.assertAlmostEqual(float(rows[1][2]), 14.)
                self.assertAlmostEqual(float(rows[2][1]), 0.)


if __name__ == "__main__":
    unittest.main()