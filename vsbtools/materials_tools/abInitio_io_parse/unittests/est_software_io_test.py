import unittest
from pathlib import Path
from my_packages.materials_tools.abInitio_io_parse import ext_software_io


class CrystalDataset_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.multiimage_poscars = Path("multiimage_POSCARS")
        self.cifs_folder = Path("cifs")

    def test_readPoscars(self):
        poscars = ext_software_io.read_poscars(self.multiimage_poscars)
        self.assertEqual(len(poscars), 1051)
