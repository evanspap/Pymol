import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).parent
spec = importlib.util.spec_from_file_location(
    "protein_tools_test",
    ROOT / "__init__.py",
    submodule_search_locations=[str(ROOT)],
)
plugin = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = plugin
spec.loader.exec_module(plugin)


class PackageTests(unittest.TestCase):
    def test_initializes_all_child_plugins(self):
        app = object()
        with patch.object(plugin.label_around, "__init_plugin__") as label_init, \
                patch.object(plugin.ShowAround, "__init_plugin__") as around_init, \
                patch.object(plugin.save_selection, "__init_plugin__") as save_init, \
                patch.object(plugin.save_objects, "__init_plugin__") as objects_init:
            plugin.__init_plugin__(app)

        label_init.assert_called_once_with(app)
        around_init.assert_called_once_with(app)
        save_init.assert_called_once_with(app)
        objects_init.assert_called_once_with(app)


if __name__ == "__main__":
    unittest.main()