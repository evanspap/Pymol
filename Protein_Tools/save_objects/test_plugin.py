from pathlib import Path
import sys
import unittest
from unittest.mock import Mock, patch


PACKAGE_ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(PACKAGE_ROOT.parent))
from Protein_Tools import save_objects as plugin


class SaveObjectsTests(unittest.TestCase):
    def setUp(self):
        self.original_cmd = plugin.cmd
        plugin.cmd = Mock()
        plugin.cmd.get_names.return_value = ["4TPW", "group1"]
        plugin.cmd.get_type.side_effect = {
            "4TPW": "object:molecule",
            "group1": "object:group",
        }.__getitem__

    def tearDown(self):
        plugin.cmd = self.original_cmd

    def test_exports_only_enabled_molecular_objects(self):
        with patch.object(plugin, "export_selections", return_value=["4TPW_residues.csv"]) as export:
            paths = plugin.save_objects("/tmp/output")

        self.assertEqual(paths, ["4TPW_residues.csv"])
        plugin.cmd.get_names.assert_called_once_with("objects", enabled_only=1)
        export.assert_called_once_with(["4TPW"], "/tmp/output")

    def test_no_enabled_molecules_does_not_export_groups(self):
        plugin.cmd.get_names.return_value = ["group1"]

        with patch.object(plugin, "export_selections") as export:
            self.assertEqual(plugin.save_objects(), [])
        export.assert_not_called()


if __name__ == "__main__":
    unittest.main()