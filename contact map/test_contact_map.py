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
        self.assertEqual(self.plugin.__version__, "v5")
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
        np.testing.assert_allclose(panels[1]["profile"], [2., 2., np.nan], equal_nan=True)
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
        self.window.mean_split_controls["minimum"].setValue(100.)
        self.assertEqual(self.window._csv_datasets, {})
        self.assertFalse(self.window.export_csv_button.isEnabled())

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
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], [1., 2., 1.])
        np.testing.assert_allclose(self.window._csv_datasets["std_profile"]["values"], [1., 2., 1.])
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
            frame_paths = list(Path(destination).glob("*/*frame_000002.csv"))
            self.assertEqual(len(frame_paths), 1)

    def test_colorbar_handle_drags_preserve_sd_profile(self):
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
        self.assertFalse(mean["auto"].isChecked())
        self.assertAlmostEqual(mean["maximum"].value(), 15., delta=0.2)
        self.assertAlmostEqual(mean["minimum"].value(), 10., delta=0.2)
        original = self.window.canvas._panels[1]["profile"].copy()
        sd = self.window.std_split_controls
        drag(sd, "maximum", 0.5)
        drag(sd, "minimum", 0.25)
        np.testing.assert_allclose(self.window.canvas._panels[1]["profile"], original, equal_nan=True)
        np.testing.assert_allclose(self.window._csv_datasets["std_profile"]["values"], original, equal_nan=True)
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
        self.assertLess(chart.right(), plot.right() + 335)
        self.window.profile_check.setChecked(False)
        chart = self.window.canvas._profile_rect(self.window.canvas._panels[0], plot)
        self.assertLess(chart.right(), plot.right() + 14)

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

    def test_v4_current_export_retains_jpeg_support(self):
        self.select_mode("both")
        self.window.profile_check.setChecked(True)
        with tempfile.TemporaryDirectory() as destination:
            output = self.window._export_current_bundle(destination, True, ["std"], image_format="jpg")
            images = list(output.glob("*.jpg"))
            self.assertEqual(len(images), 1)
            self.assertFalse(QtGui.QImage(str(images[0])).isNull())
            self.assertEqual(len(list(output.glob("*.csv"))), 1)


if __name__ == "__main__":
    unittest.main()