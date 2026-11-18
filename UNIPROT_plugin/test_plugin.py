"""Regression tests; uses an isolated headless PyMOL instance when available."""
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("uniprot_test_plugin", Path(__file__).with_name("__init__.py"))
plugin = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = plugin
spec.loader.exec_module(plugin)


class Response(io.BytesIO):
    def __init__(self, body, size=None):
        super().__init__(body)
        self.headers = {"Content-Length": str(len(body) if size is None else size)}


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)

    def parse(self, contents):
        path = self.root / "list.txt"
        path.write_text(contents, encoding="utf-8-sig")
        return plugin.read_uniprot_file(str(path))

    def test_list_formats(self):
        self.assertEqual(self.parse("# comment\nP06730 IF4E_HUMAN\nP06730\n"), ["P06730", "IF4E_HUMAN"])
        self.assertEqual(self.parse("note,UniProt accession\nx,P06730\ny,P04637\n"), ["P06730", "P04637"])
        self.assertEqual(self.parse("Entry Name\nIF4E_HUMAN\n"), ["IF4E_HUMAN"])
        self.assertEqual(self.parse("accession\tnote\nP06730\tx\n"), ["P06730"])
        self.assertEqual(self.parse("P06730,x\nP04637,y\n"), ["P06730", "P04637"])
        with self.assertRaises(ValueError):
            self.parse("# empty\n")

    def test_pasted_batch_names(self):
        expected = ["BCL2_HUMAN", "CAB39_HUMAN", "PYR1_HUMAN", "CARM1_HUMAN", "PEBB_HUMAN"]
        for separator in (" ", "\n", "\r\n", "\t", ", ", ";"):
            self.assertEqual(plugin.parse_uniprot_queries(separator.join(expected)), expected)
        self.assertEqual(plugin.parse_uniprot_queries(" BCL2_HUMAN, CAB39_HUMAN; BCL2_HUMAN "), expected[:2])
        self.assertEqual(plugin.parse_uniprot_queries(" "), [])
        summary = {"errors": [], "downloaded": 0, "cached": 0}
        with patch.object(plugin, "_load_uniprot_structures", return_value=summary) as load:
            result = plugin.load_uniprot_batch(plugin.parse_uniprot_queries("\n".join(expected)))
        self.assertEqual([call.args[0] for call in load.call_args_list], expected)
        self.assertEqual(len(result["summaries"]), 5)

    def test_resolver_rejects_unsplit_batch(self):
        with patch.object(plugin, "_http_json") as network:
            with self.assertRaisesRegex(ValueError, "Expected one UniProt entry"):
                plugin.resolve_uniprot("BCL2_HUMAN\nCAB39_HUMAN")
            network.assert_not_called()

    def test_experimental_category_filter_before_limit(self):
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", [
            plugin.PDBReference("XR01", ["A"], "X-ray diffraction"),
            plugin.PDBReference("NM01", ["A"], "NMR"),
            plugin.PDBReference("EM01", ["A"], "Electron microscopy"),
            plugin.PDBReference("NM02", ["A"], "NMR")])
        for selected, expected in [("NMR", ["NM01"]), ("X-ray", ["XR01"]), ("EM", ["EM01"]),
                                   ("NMR+EM", ["NM01"])]:
            with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif", return_value="cached.pdb") as fetch:
                plugin.load_uniprot_structures("TFE2_HUMAN", include_alphafold=False,
                    download_only=True, pdb_methods=selected, limit=1, cache_dir=str(self.root))
            self.assertEqual([call.args[0] for call in fetch.call_args_list], expected)
        with self.assertRaises(ValueError):
            plugin._parse_pdb_methods("invalid")

    def test_absent_category_has_no_download_or_empty_group(self):
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", [plugin.PDBReference("XR01", ["A"], "X-ray")])
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif") as fetch:
            result = plugin.load_uniprot_structures("TFE2_HUMAN", include_alphafold=False,
                download_only=True, pdb_methods="NMR", cache_dir=str(self.root))
        fetch.assert_not_called()
        self.assertEqual(result["method_groups"], {})
        self.assertEqual(result["files"], [])

    def test_alphafold_only_skips_experimental_downloads(self):
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", [plugin.PDBReference("NM01", ["A"], "NMR")])
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "get_alphafold_pdb_url", return_value="https://example.org/af.pdb"), patch.object(plugin, "_download"), patch.object(plugin, "_download_pdb_or_cif") as fetch:
            result = plugin.load_uniprot_structures("TFE2_HUMAN", pdb_methods="none",
                download_only=True, cache_dir=str(self.root))
        fetch.assert_not_called()
        self.assertEqual(len(result["files"]), 1)

    def test_atomic_download_and_cache(self):
        path = self.root / "model.pdb"
        with patch.object(plugin.urllib.request, "urlopen", return_value=Response(b"coordinates")):
            plugin._download("https://example.org/model", str(path))
        self.assertEqual(path.read_bytes(), b"coordinates")
        self.assertEqual(list(self.root.glob("*.part")), [])
        with patch.object(plugin.urllib.request, "urlopen", return_value=Response(b"bad", size=20)):
            with self.assertRaises(IOError):
                plugin._download("https://example.org/model", str(path))
        self.assertEqual(path.read_bytes(), b"coordinates")
        self.assertEqual(list(self.root.glob("*.part")), [])
        with patch.object(plugin.urllib.request, "urlopen", side_effect=AssertionError("cache should avoid network")):
            (self.root / "TEST.pdb").write_bytes(b"cached")
            self.assertEqual(plugin._download_pdb_or_cif("TEST", str(self.root)), str(self.root / "TEST.pdb"))

    def test_cancelled_download_leaves_no_cache(self):
        path = self.root / "model.pdb"
        cancel = threading.Event()
        def report(event):
            cancel.set()
        with patch.object(plugin.urllib.request, "urlopen", return_value=Response(b"x" * 600000)):
            with self.assertRaises(plugin.Cancelled):
                plugin._download("https://example.org/model", str(path), progress=report, cancel=cancel)
        self.assertFalse(path.exists())
        self.assertEqual(list(self.root.iterdir()), [])

    def test_progress_distinguishes_cached_downloaded_failed(self):
        refs = [plugin.PDBReference(code, ["A"]) for code in ("OLD1", "NEW1", "BAD1")]
        entry = plugin.UniProtEntry("P06730", "IF4E_HUMAN", refs)
        (self.root / "OLD1.pdb").write_text("old")
        def fetch(code, folder, **kwargs):
            if code == "BAD1":
                raise IOError("missing")
            path = Path(folder) / (code + ".pdb")
            path.write_text("coordinates")
            return str(path)
        events = []
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif", side_effect=fetch):
            result = plugin.load_uniprot_structures("P06730", include_alphafold=False,
                download_only=True, cache_dir=str(self.root), progress=events.append)
        self.assertEqual((result["downloaded"], result["cached"], result["failed"]), (1, 1, 1))
        self.assertEqual((events[-1]["completed"], events[-1]["total"]), (3, 3))
        self.assertIn("0 remaining", events[-1]["message"])
        self.assertEqual(result["loaded"], [])

    def test_worker_reports_memory_error(self):
        with patch.object(plugin, "load_uniprot_batch", side_effect=MemoryError):
            thread, cancel, events = plugin._start_job(["P06730"], {})
            thread.join(timeout=2)
        self.assertFalse(thread.is_alive())
        kind, message = events.get_nowait()
        self.assertEqual(kind, "error")
        self.assertIn("Download only", message)

    def test_batch_recovers_and_releases_lock(self):
        summary = {"errors": [], "downloaded": 1, "cached": 0}
        with patch.object(plugin, "_load_uniprot_structures", side_effect=[ValueError("bad entry"), summary]) as load:
            result = plugin.load_uniprot_batch(["bad", "good", "good"])
        self.assertEqual(load.call_count, 2)
        self.assertEqual(len(result["summaries"]), 1)
        self.assertEqual(len(result["errors"]), 1)
        self.assertFalse(plugin._RUN_LOCK.locked())
        cancelled = threading.Event()
        cancelled.set()
        with self.assertRaises(plugin.Cancelled):
            plugin.load_uniprot_batch(["P06730"], cancel=cancelled)
        self.assertFalse(plugin._RUN_LOCK.locked())


try:
    import pymol2
except ImportError:
    pymol2 = None


@unittest.skipIf(pymol2 is None, "PyMOL not available")
class PyMOLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.instance = pymol2.PyMOL()
        cls.instance.start()
        cls.old_cmd = plugin.cmd
        plugin.cmd = cls.instance.cmd

    @classmethod
    def tearDownClass(cls):
        plugin.cmd = cls.old_cmd
        cls.instance.stop()

    def setUp(self):
        self.cmd = self.instance.cmd
        self.cmd.delete("all")
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.path = str(Path(self.folder.name) / "full.pdb")
        self.cmd.fab("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY", "mapped")
        self.cmd.alter("mapped", "chain='B'")
        self.cmd.create("unrelated", "mapped")
        self.cmd.alter("unrelated", "chain='A'")
        self.cmd.translate([30, 0, 0], "unrelated")
        self.cmd.create("full", "mapped or unrelated")
        self.cmd.save(self.path, "full")
        self.cmd.delete("all")

    def test_isolates_one_of_116_groups_with_disabled_children(self):
        self.cmd.delete("all")
        for index in range(116):
            group = f"PROT{index:03d}_HUMAN"
            obj = f"pdb_{index}"
            subgroup = f"method_{index}"
            self.cmd.pseudoatom(obj)
            self.cmd.group(subgroup, obj)
            self.cmd.group(group, subgroup)
        self.cmd.disable("all")
        objects = plugin.uniprot_show("PROT007_HUMAN")
        self.assertEqual(objects, ["pdb_7"])
        enabled = set(self.cmd.get_names("objects", enabled_only=1))
        self.assertEqual(enabled, {"PROT007_HUMAN", "method_7", "pdb_7"})
        self.assertEqual(len(plugin.list_uniprot_groups()), 116)
        plugin.uniprot_show("PROT008_HUMAN")
        self.assertEqual(set(self.cmd.get_names("objects", enabled_only=1)),
                         {"PROT008_HUMAN", "method_8", "pdb_8"})
        plugin.uniprot_show("PROT007_HUMAN", only="0")
        self.assertIn("PROT008_HUMAN", self.cmd.get_names("objects", enabled_only=1))
        plugin.uniprot_hide("PROT007_HUMAN")
        self.assertNotIn("PROT007_HUMAN", self.cmd.get_names("objects", enabled_only=1))

    def test_registered_visibility_command(self):
        self.cmd.pseudoatom("target")
        self.cmd.pseudoatom("other")
        self.cmd.group("X-ray", "target")
        self.cmd.group("BCL2_HUMAN", "X-ray")
        self.cmd.group("CAB39_HUMAN", "other")
        self.cmd.disable("all")
        plugin._register_commands()
        self.cmd.do("uniprot_show BCL2_HUMAN")
        enabled = set(self.cmd.get_names("objects", enabled_only=1))
        self.assertEqual(enabled, {"BCL2_HUMAN", "X-ray", "target"})

    def test_invalid_visibility_target_does_not_blank_scene(self):
        self.cmd.pseudoatom("existing")
        before = self.cmd.get_names("objects", enabled_only=1)
        with self.assertRaises(ValueError):
            plugin.uniprot_show("MISSING_HUMAN")
        self.assertEqual(before, self.cmd.get_names("objects", enabled_only=1))
        self.assertFalse(plugin._RUN_LOCK.locked())

    def test_visibility_preserves_hidden_representations(self):
        self.cmd.pseudoatom("ligand")
        self.cmd.group("BCL2_HUMAN", "ligand")
        self.cmd.hide("everything", "ligand")
        plugin.uniprot_show("BCL2_HUMAN")
        self.assertEqual(self.cmd.count_atoms("ligand and rep spheres"), 0)
        self.assertEqual(self.cmd.count_atoms("ligand and rep nonbonded"), 0)

    def test_nonempty_counted_method_groups_and_repeated_load(self):
        refs = [plugin.PDBReference("NM01", ["B"], "NMR")] + [
            plugin.PDBReference(f"XR0{i}", ["B"], "X-ray diffraction") for i in range(4)]
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", refs)
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif", return_value=self.path):
            result = plugin.load_uniprot_structures("TFE2_HUMAN", include_alphafold=False,
                                                  align=False, cache_dir=self.folder.name)
            repeated = plugin.load_uniprot_structures("TFE2_HUMAN", include_alphafold=False,
                                                    align=False, cache_dir=self.folder.name)
        self.assertEqual(result["method_counts"], {"NMR": 1, "X-ray": 4})
        self.assertEqual(result["method_groups"], {"NMR": "TFE2_NMR_1", "X-ray": "TFE2_X-ray_4"})
        for category, count in result["method_counts"].items():
            self.assertEqual(len(self.cmd.get_object_list(result["method_groups"][category])), count)
        self.assertEqual(len(plugin.list_uniprot_groups()), 2)
        self.assertNotEqual(repeated["method_groups"]["NMR"], result["method_groups"]["NMR"])
        self.assertFalse(any("_EM_" in name or "Other_methods" in name for name in self.cmd.get_names("objects")))

    def test_failed_load_creates_no_empty_method_group(self):
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", [plugin.PDBReference("NM01", ["Z"], "NMR")])
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif", return_value=self.path):
            result = plugin.load_uniprot_structures("TFE2_HUMAN", include_alphafold=False,
                                                  align=False, cache_dir=self.folder.name)
        self.assertEqual(result["method_groups"], {})
        self.assertEqual(self.cmd.get_names("objects"), ["TFE2_HUMAN"])

    def test_alphafold_confidence_selection_all_methods(self):
        ref = plugin._load_related_chains(self.path, "AF_REFERENCE", ["B"], "group")
        mobile = plugin._load_related_chains(self.path, "MOBILE", ["B"], "group")
        self.cmd.alter(ref, "b=90")
        self.cmd.alter(ref + " and resi 1-10", "b=10")
        self.cmd.alter(ref + " and resi 11-20", "b=50")
        self.cmd.alter(mobile, "b=0")
        original_atoms = self.cmd.count_atoms(ref)
        for method in ("align", "super", "cealign"):
            self.cmd.translate([8, -4, 3], mobile)
            events = []
            plugin._align_object(mobile, ref, method, events.append, reference_b_min=50)
            self.assertTrue(any("and b > 50" in event["message"] for event in events))
        self.assertEqual(self.cmd.count_atoms(ref), original_atoms)
        self.assertGreater(self.cmd.count_atoms(mobile + " and b = 0"), 0)
        self.cmd.alter(ref, "b=50")
        before = self.cmd.get_coords(mobile).copy()
        with self.assertRaisesRegex(ValueError, "fewer than 3"):
            plugin._align_object(mobile, ref, "align", None, reference_b_min=50)
        self.assertTrue((before == self.cmd.get_coords(mobile)).all())

    def test_confidence_filter_does_not_fall_back_to_experimental_b(self):
        entry = plugin.UniProtEntry("P06730", "TFE2_HUMAN", [plugin.PDBReference("NM01", ["B"], "NMR")])
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "get_alphafold_pdb_url", return_value=None), patch.object(plugin, "_download_pdb_or_cif", return_value=self.path), patch.object(plugin, "_align_object") as fit:
            result = plugin.load_uniprot_structures("TFE2_HUMAN", af_confidence_only=True, cache_dir=self.folder.name)
        fit.assert_not_called()
        self.assertEqual(result["visible"], [])
        self.assertIn("No AlphaFold reference", result["errors"][0])

    def test_keeps_mapped_B_instead_of_unrelated_A(self):
        obj = plugin._load_related_chains(self.path, "TEST", ["B"], "group")
        self.assertEqual(self.cmd.get_chains(obj), ["B"])
        self.assertFalse(any(name.endswith("_full") for name in self.cmd.get_names("objects")))

    def test_retains_displays_and_transforms_ligands_waters_and_ions(self):
        self.cmd.load(self.path, "full")
        # Ligand deliberately shares the discarded protein's chain A.
        self.cmd.pseudoatom("extras", name="C1", resn="LIG", chain="A", elem="C", pos=[1, 2, 3])
        self.cmd.pseudoatom("extras", name="C2", resn="LIG", chain="A", elem="C", pos=[2.4, 2, 3])
        self.cmd.bond("extras and name C1", "extras and name C2")
        self.cmd.pseudoatom("extras", name="O", resn="HOH", chain="W", resi="2", elem="O", pos=[5, 6, 7])
        self.cmd.pseudoatom("extras", name="ZN", resn="ZN", chain="", resi="3", elem="Zn", pos=[8, 9, 10])
        self.cmd.save(self.path, "full or extras")
        self.cmd.delete("all")
        ref = plugin._load_related_chains(self.path, "REF", ["B"], "group")
        mobile = plugin._load_related_chains(self.path, "MOB", ["B"], "group")
        self.assertEqual(self.cmd.get_chains(mobile + " and polymer.protein"), ["B"])
        self.assertEqual(self.cmd.count_atoms(mobile + " and resn LIG"), 2)
        self.assertEqual(self.cmd.count_atoms(mobile + " and solvent"), 1)
        self.assertEqual(self.cmd.count_atoms(mobile + " and elem Zn"), 1)
        plugin._show_structure(mobile)
        self.assertEqual(self.cmd.count_atoms(mobile + " and resn LIG and rep sticks"), 2)
        self.assertEqual(self.cmd.count_atoms(mobile + " and solvent and rep nonbonded"), 1)
        for method in ("align", "super", "cealign"):
            self.cmd.translate([8, -4, 3], mobile)
            plugin._align_object(mobile, ref, method, None)
            actual = self.cmd.get_atom_coords(mobile + " and resn LIG and name C1")
            expected = self.cmd.get_atom_coords(ref + " and resn LIG and name C1")
            for a, b in zip(actual, expected):
                self.assertAlmostEqual(a, b, delta=0.002)

    def test_prefers_A_only_when_mapped(self):
        obj = plugin._load_related_chains(self.path, "TEST", ["A", "B"], "group")
        self.assertEqual(self.cmd.get_chains(obj), ["A"])

    def test_no_matching_mapping_cleans_temporary_object(self):
        with self.assertRaises(ValueError):
            plugin._load_related_chains(self.path, "TEST", ["Z"], "group")
        self.assertEqual(self.cmd.get_object_list(), [])

    def test_preserves_ensemble_states(self):
        self.cmd.load(self.path, "ensemble")
        self.cmd.create("ensemble", "ensemble", source_state=1, target_state=2)
        self.cmd.save(self.path, "ensemble", state=0)
        self.cmd.delete("all")
        obj = plugin._load_related_chains(self.path, "TEST", ["B"], "group")
        self.assertEqual(self.cmd.count_states(obj), 2)
        self.assertEqual(self.cmd.get_chains(obj), ["B"])

    def test_real_alignment_to_reference_chain(self):
        ref = plugin._load_related_chains(self.path, "REF", ["B"], "group")
        mobile = plugin._load_related_chains(self.path, "MOB", ["B"], "group")
        self.cmd.translate([8, -4, 3], mobile)
        for method in ("align", "super", "cealign"):
            plugin._align_object(mobile, ref, method, None)
            self.assertLess(self.cmd.rms_cur(mobile + " and name CA", ref + " and name CA", matchmaker=-1), 0.1)

    def test_fallback_reference_and_failed_fit_hidden(self):
        entry = plugin.UniProtEntry("P06730", "IF4E_HUMAN", [plugin.PDBReference("ONE1", ["B"]), plugin.PDBReference("TWO2", ["B"])])
        with patch.object(plugin, "resolve_uniprot", return_value=entry), patch.object(plugin, "_download_pdb_or_cif", return_value=self.path), patch.object(plugin, "_align_object", side_effect=ValueError("no fit")):
            summary = plugin.load_uniprot_structures("P06730", include_alphafold=False, cache_dir=self.folder.name)
        self.assertEqual(summary["visible"], ["ONE1"])
        self.assertEqual(summary["loaded"], ["ONE1", "TWO2"])
        self.assertNotIn("TWO2", self.cmd.get_names("objects", enabled_only=1))
        self.assertEqual(len(summary["errors"]), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
