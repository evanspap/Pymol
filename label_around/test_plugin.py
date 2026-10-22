import importlib.util
from pathlib import Path
import unittest
from unittest.mock import Mock


spec = importlib.util.spec_from_file_location("label_around_plugin", Path(__file__).with_name("__init__.py"))
plugin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plugin)


class LabelAroundTests(unittest.TestCase):
    def setUp(self):
        self.original_cmd = plugin.cmd
        self.cmd = Mock()
        self.cmd.count_atoms.return_value = 3
        plugin.cmd = self.cmd

    def tearDown(self):
        plugin.cmd = self.original_cmd

    def test_labels_alpha_carbons_without_a_position_offset(self):
        result = plugin.label_around("chain A")

        self.assertEqual(result, 3)
        ca_selection = "(chain A) and name CA"
        self.cmd.count_atoms.assert_called_once_with(ca_selection)
        label_call = self.cmd.label.call_args
        self.assertEqual(label_call.args[0], ca_selection)
        self.assertEqual(label_call.kwargs, {"quiet": 1})
        self.assertEqual(eval(label_call.args[1], {}, {"resn": "ALA", "resi": "42"}), "A42")
        self.assertEqual(
            self.cmd.unset.call_args_list,
            [
                unittest.mock.call("label_position", ca_selection, quiet=1),
                unittest.mock.call("label_placement_offset", ca_selection, quiet=1),
            ],
        )

    def test_defaults_to_current_selection(self):
        plugin.label_around()

        self.cmd.count_atoms.assert_called_once_with("(sele) and name CA")

    def test_unknown_residue_uses_x(self):
        expression = f"{plugin._AMINO_ACID_CODES!r}.get(resn, 'X') + str(resi)"
        self.assertEqual(eval(expression, {}, {"resn": "UNK", "resi": "7A"}), "X7A")

    def test_rejects_empty_ca_selection(self):
        self.cmd.count_atoms.return_value = 0
        with self.assertRaisesRegex(ValueError, "contains no CA atoms"):
            plugin.label_around()


if __name__ == "__main__":
    unittest.main()