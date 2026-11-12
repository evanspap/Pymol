import importlib.util
from pathlib import Path
import unittest
from unittest.mock import Mock


spec = importlib.util.spec_from_file_location("show_around_plugin", Path(__file__).with_name("__init__.py"))
plugin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plugin)


class ShowAroundTests(unittest.TestCase):
    def setUp(self):
        self.original_cmd = plugin.cmd
        self.cmd = Mock()
        self.cmd.count_atoms.return_value = 3
        self.cmd.select.return_value = 12
        plugin.cmd = self.cmd

    def tearDown(self):
        plugin.cmd = self.original_cmd

    def test_defaults_to_sele_six_angstroms_and_whole_protein_residues(self):
        result = plugin.show_around()

        self.assertEqual(result, 12)
        self.cmd.count_atoms.assert_called_once_with("sele")
        self.cmd.select.assert_called_once_with(
            "nearby_residues",
            "byres ((polymer.protein) within 6 of (sele))",
            quiet=1,
        )

    def test_accepts_custom_selection_distance_and_name(self):
        plugin.show_around("ligand and name C1", 8.5, "nearby")

        self.cmd.select.assert_called_once_with(
            "nearby",
            "byres ((polymer.protein) within 8.5 of (ligand and name C1))",
            quiet=1,
        )

    def test_rejects_empty_source_and_invalid_distance(self):
        self.cmd.count_atoms.return_value = 0
        with self.assertRaisesRegex(ValueError, "empty or does not exist"):
            plugin.show_around()

        self.cmd.count_atoms.return_value = 1
        for distance in (0, -1, float("inf"), "invalid"):
            with self.subTest(distance=distance), self.assertRaises(ValueError):
                plugin.show_around(distance=distance)

    def test_rejects_invalid_output_name(self):
        with self.assertRaisesRegex(ValueError, "Selection name"):
            plugin.show_around(name="nearby selection")


if __name__ == "__main__":
    unittest.main()