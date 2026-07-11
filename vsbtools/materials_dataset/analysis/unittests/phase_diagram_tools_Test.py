import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from ...io.preset_loaders import load_from_materials_project
from ..phase_diagram_tools import PhaseDiagramTools
from ..summary import collect_summary_df, print_pretty_df

PATH_WITH_TESTS = Path(__file__).parent

class symmetry_tools_Test(unittest.TestCase):
    def setUp(self):
        self.ds = load_from_materials_project({"Mo", "Si", "B", "P"}, message='')
        self.pd_tools = PhaseDiagramTools(self.ds)
        self.e_hull_fn = self.pd_tools.height_above_hull_pa

    def test_pretty_phasediag(self):
        callables = {"e_hull": self.e_hull_fn}
        summary_df = collect_summary_df(self.ds, callables=callables)
        self.assertIn("e_hull", summary_df.columns)
        self.assertEqual(len(summary_df), len(self.ds))
        with TemporaryDirectory() as tmpdir:
            table_path = Path(tmpdir) / "table.txt"
            print_pretty_df(summary_df, table_path, sort_by='e_hull')
            self.assertTrue(table_path.exists())
            self.assertIn("e_hull", table_path.read_text())
