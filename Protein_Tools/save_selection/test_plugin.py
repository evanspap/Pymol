import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch


spec = importlib.util.spec_from_file_location("save_selection_plugin", Path(__file__).with_name("__init__.py"))
plugin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plugin)


class SaveSelectionTests(unittest.TestCase):
    def setUp(self):
        self.original_cmd = plugin.cmd
        self.cmd = Mock()
        self.cmd.get_names.return_value = ["C1", "empty"]

        def iterate(selection, expression, space=None, quiet=1):
            if selection.startswith("(C1)"):
                space["rows"].append(("A", "42", "ALA"))
                space["rows"].append(("B", "7A", "GLY"))
            return 0

        self.cmd.iterate.side_effect = iterate
        plugin.cmd = self.cmd
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)

    def tearDown(self):
        plugin.cmd = self.original_cmd

    def test_writes_one_csv_for_each_selection(self):
        paths = plugin.export_selections(["C1", "empty"], self.temp_dir.name)

        self.assertEqual([Path(path).name for path in paths], [
            "C1_residues.csv",
            "empty_residues.csv",
        ])
        with open(paths[0], newline="", encoding="utf-8") as handle:
            self.assertEqual(list(csv.reader(handle)), [
                ["Chain", "Res_Num", "Res_Name"],
                ["A", "42", "ALA"],
                ["B", "7A", "GLY"],
            ])
        with open(paths[1], newline="", encoding="utf-8") as handle:
            self.assertEqual(list(csv.reader(handle)), [["Chain", "Res_Num", "Res_Name"]])
        self.cmd.get_names.assert_not_called()

    def test_dialog_receives_all_named_selections(self):
        with patch.object(plugin, "_open_selection_dialog") as open_dialog:
            plugin.save_selections(self.temp_dir.name)

        self.cmd.get_names.assert_called_once_with("selections", enabled_only=0)
        open_dialog.assert_called_once_with(["C1", "empty"], self.temp_dir.name)

    def test_can_export_only_a_chosen_selection(self):
        paths = plugin.export_selections(["C1"], self.temp_dir.name)

        self.assertEqual([Path(path).name for path in paths], ["C1_residues.csv"])

    def test_iterates_only_ca_atoms(self):
        plugin.export_selections(["C1", "empty"], self.temp_dir.name)

        self.assertEqual(self.cmd.iterate.call_args_list[0].args[0], "(C1) and name CA")
        self.assertEqual(self.cmd.iterate.call_args_list[1].args[0], "(empty) and name CA")

    def test_no_selections_returns_no_files(self):
        self.cmd.get_names.return_value = []

        self.assertEqual(plugin.save_selections(self.temp_dir.name), [])
        self.cmd.get_names.assert_called_once_with("selections", enabled_only=0)
        self.assertEqual(list(Path(self.temp_dir.name).iterdir()), [])


if __name__ == "__main__":
    unittest.main()