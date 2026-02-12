import unittest
from pathlib import Path

from pymatgen.core import Structure

from ..structure_checks import check_min_dist_pmg

PATH_WITH_TESTS = Path(__file__).parent

BAD_POSCAR_PATH = PATH_WITH_TESTS / "149POSCAR"

LICOO_POSCARS_DIR = PATH_WITH_TESTS / "../../materials_dataset/unittests_datasets/alex_LiCoO_POSCARS"


class StructureChecks_Test(unittest.TestCase):
    def test_check_min_dist_pmg_invalid_structure(self):
        structure = Structure.from_file(BAD_POSCAR_PATH)
        is_ok, info = check_min_dist_pmg(structure)
        self.assertFalse(is_ok, msg=info.get("reason"))

    def test_check_min_dist_pmg_alex_licoo_poscars(self):
        for poscar_path in sorted(LICOO_POSCARS_DIR.iterdir()):
            if not poscar_path.is_file():
                continue
            with self.subTest(poscar=poscar_path.name):
                structure = Structure.from_file(poscar_path)
                is_ok, info = check_min_dist_pmg(structure)
                self.assertTrue(is_ok, msg=f"{poscar_path.name}: {info.get('reason')}")
