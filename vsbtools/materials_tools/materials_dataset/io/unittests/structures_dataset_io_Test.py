import unittest
import tempfile
from pathlib import Path
from ..structures_dataset_io import StructureDatasetIO, exploded_zip_tree

PATH_WITH_TESTS = Path(__file__).parent
PATH_TEST_DATASET = PATH_WITH_TESTS / "../../unittests_datasets"


class yaml_csv_poscars_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.zipped_cifs = PATH_TEST_DATASET / 'zipped_cifs'

    def test_unzip(self):
        with exploded_zip_tree(self.zipped_cifs) as tmp_path:
            self.utils = StructureDatasetIO(tmp_path)
            ds = self.utils.load_from_directory()
            self.assertEqual(len(ds), 696)


class extxyz_Test(unittest.TestCase):

    def setUp(self) -> None:
        self.extxyz_file = (PATH_WITH_TESTS /
                            "../../unittests_datasets/Cu-Si-P-Ca_guided.extxyz").resolve()

    def test_extxyz_to_poscars(self):
        utils = StructureDatasetIO(PATH_WITH_TESTS)
        ds = utils.load_from_extxyz(self.extxyz_file)
        self.assertGreater(len(ds), 0)
        self.assertEqual(ds[0].natoms, 15)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir)
            utils.dump_poscar_files(ds, out_dir)
            poscars = list(out_dir.rglob("*POSCAR"))
            self.assertEqual(len(poscars), len(ds))
