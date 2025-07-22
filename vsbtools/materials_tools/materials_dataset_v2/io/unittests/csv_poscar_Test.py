import unittest
from pathlib import Path
from ..csv_poscars import read, write
from ..preset_loaders import load_mattersim_estimated_set

PATH_WITH_TESTS = Path(__file__).parent


class csv_poscars_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.csv = PATH_WITH_TESTS / 'mattersim_res_for_POSCARS.csv'
        self.poscars_folder = PATH_WITH_TESTS / "POSCARS"


    def test_read_write(self):
        ds = load_mattersim_estimated_set(self.csv, self.poscars_folder)
        write(ds,PATH_WITH_TESTS / 'description.yaml')