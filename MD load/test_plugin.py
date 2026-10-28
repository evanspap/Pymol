"""Integration checks for MD Analysis using a Python interpreter with PyMOL.
Author: OpenAI Codex. Summary: verify actual states, coordinates and safe loading.
Usage/example: python -B test_plugin.py
"""
import contextlib
import importlib.util
import io
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pymol
from pymol2 import SingletonPyMOL
_instance = SingletonPyMOL()
_instance.start()
from pymol import cmd

spec = importlib.util.spec_from_file_location("md_analysis", Path(__file__).with_name("__init__.py"))
plugin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plugin)


def pdb_frame(offset):
    return "".join(
        f"ATOM  {i:5d}  CA  ALA A{i:4d}    {x + offset:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
        for i, (x, y, z) in enumerate(((0., 0., 0.), (2., 1., 0.), (4., 0., 1.)), 1)
    ) + "TER\n"


def write_dcd(path):
    # Minimal little-endian CHARMM DCD, three frames, three atoms, no unit cell.
    def record(handle, data):
        handle.write(struct.pack("<i", len(data)) + data + struct.pack("<i", len(data)))
    with path.open("wb") as handle:
        controls = [3, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24]
        header = bytearray(b"CORD" + struct.pack("<20i", *controls))
        header[40:44] = struct.pack("<f", 1.0)
        record(handle, header)
        record(handle, struct.pack("<i", 1) + b"MD Analysis test".ljust(80))
        record(handle, struct.pack("<i", 3))
        for offset in (0., 5., 10.):
            for values in ((offset, offset + 2, offset + 4), (0., 1., 0.), (0., 0., 1.)):
                record(handle, struct.pack("<3f", *values))


class LoaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.topology = self.root / "protein with spaces.pdb"
        self.topology.write_text(pdb_frame(0) + "END\n")
        self.ensemble = self.root / "ensemble.pdb"
        self.ensemble.write_text("".join(
            f"MODEL     {i}\n" + pdb_frame(offset) + "ENDMDL\n"
            for i, offset in enumerate((0, 5, 10), 1)) + "END\n")
        self.dcd = self.root / "trajectory.dcd"
        write_dcd(self.dcd)
        self.names = set(cmd.get_names("all"))

    def tearDown(self):
        for name in set(cmd.get_names("all")) - self.names:
            cmd.delete(name)
        self.tmp.cleanup()

    def test_multimodel_and_unique_names(self):
        first = plugin.md_load(str(self.ensemble), name="md_test")
        second = plugin.md_load(str(self.topology), name="md_test")
        self.assertNotEqual(first, second)
        self.assertEqual(cmd.count_states(first), 3)
        self.assertEqual(cmd.count_states(second), 1)
        self.assertAlmostEqual(float(cmd.get_coords(first, 3)[0][0]), 10., places=3)

    def test_separate_trajectory(self):
        name = plugin.md_load(str(self.topology), str(self.dcd), name="md_test")
        self.assertEqual(cmd.count_states(name), 3)
        self.assertAlmostEqual(float(cmd.get_coords(name, 3)[0][0]), 10., places=3)

    def test_principal_axis_alignment_centers_and_orders_each_state(self):
        name = plugin.md_load(str(self.topology), str(self.dcd), name="md_test",
                              align_principal_axes=1)
        self.assertEqual(cmd.count_states(name), 3)
        for state in range(1, 4):
            coordinates = np.asarray(cmd.get_coords(name, state=state), dtype=float)
            self.assertTrue(np.allclose(coordinates.mean(axis=0), 0., atol=1e-5))
            covariance = coordinates.T @ coordinates
            self.assertTrue(np.allclose(covariance, np.diag(np.diag(covariance)), atol=1e-4))
            variances = np.diag(covariance)
            self.assertGreaterEqual(variances[0] + 1e-5, variances[1])
            self.assertGreaterEqual(variances[1] + 1e-5, variances[2])

    def test_trajectory_frame_selection(self):
        name = plugin.md_load(str(self.topology), str(self.dcd), name="md_test",
                              start=2, stop=3, interval=1)
        self.assertEqual(cmd.count_states(name), 2)
        self.assertAlmostEqual(float(cmd.get_coords(name, 1)[0][0]), 5., places=3)

    def test_failed_load_removes_only_new_object(self):
        original = plugin.md_load(str(self.topology), name="md_test")
        with self.assertRaises(ValueError):
            plugin.md_load(str(self.ensemble), str(self.dcd), name="md_test")
        self.assertEqual(set(cmd.get_names("all")), self.names | {original})

    def test_dry_run(self):
        before = set(cmd.get_names("all"))
        with contextlib.redirect_stdout(io.StringIO()) as output:
            plugin.md_load(str(self.topology), str(self.dcd), dry_run=1)
        self.assertEqual(set(cmd.get_names("all")), before)
        self.assertIn("Running:", output.getvalue())
        self.assertIn("Dry run:", output.getvalue())

    def test_bad_options(self):
        for options in ({"start": 0}, {"interval": 0}, {"stop": 0}, {"name": "bad name"}):
            with self.assertRaises(ValueError):
                plugin.md_load(str(self.topology), str(self.dcd), **options)

    def test_cli_header_and_dry_run(self):
        result = subprocess.run([sys.executable, "-B", plugin.__file__],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Author:", result.stdout)
        self.assertIn("Usage", result.stdout)
        result = subprocess.run([sys.executable, "-B", plugin.__file__, "--structure",
                                 str(self.topology), "--dry-run", "--align-principal-axes"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Dry run:", result.stdout)
        self.assertIn("align_principal_axes=1", result.stdout)


if __name__ == "__main__":
    unittest.main()

