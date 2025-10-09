import unittest
from pathlib import Path
from ase.io import read
from ..coordination_tools import compute_average_coordination_atoms, _to_torch_from_atoms, _soft_neighbor_counts_per_A, compute_target_share_atoms

PATH_WITH_TESTS = Path(__file__).parent

class symmetry_tools_Test(unittest.TestCase):
    def setUp(self):
        self.system = read(PATH_WITH_TESTS / "Na3B7O12_POSCAR")

    def test_pretty_phasediag(self):
        self.assertAlmostEqual(compute_average_coordination_atoms(self.system, type_A="B", type_B="O", kernel='sigmoid'),
                               3.43, delta=0.1)
        print(_soft_neighbor_counts_per_A(*_to_torch_from_atoms(self.system), type_A="B", type_B="O"))
        self.assertAlmostEqual(compute_target_share_atoms(self.system, type_A="B", type_B="O", target=3.),
                               0.57, delta=0.05)
