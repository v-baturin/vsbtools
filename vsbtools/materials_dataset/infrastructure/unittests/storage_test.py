import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from pymatgen.core import Structure
from ...crystal_entry import CrystalEntry
from ...crystal_dataset import CrystalDataset
from ..storage import load, save

PATH_WITH_TESTS = Path(__file__).parent

class storage_Test(unittest.TestCase):
    def setUp(self):
        sample_poscars_file = PATH_WITH_TESTS /  "../../io/unittests/POSCARS"
        entries = [CrystalEntry(id = str(i), structure=Structure.from_file(f), energy = -i/10) \
                      for i, f in enumerate(sample_poscars_file.rglob('*POSCAR'))]
        self.ds = CrystalDataset(entries)

    def test_saveload(self):
        with TemporaryDirectory() as tmpdir:
            pickle_path = Path(tmpdir) / "dataset.pkl.gz"
            save(self.ds, pickle_path)
            ds = load(pickle_path)

        self.assertEqual(len(ds), len(self.ds))
        self.assertEqual([entry.id for entry in ds], [entry.id for entry in self.ds])
        self.assertEqual([entry.energy for entry in ds], [entry.energy for entry in self.ds])

