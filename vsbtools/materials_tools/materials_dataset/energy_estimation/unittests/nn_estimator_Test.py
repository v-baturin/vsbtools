from pathlib import Path
import unittest
from ..nn_estimator import NNEstimator
from .. import mattersim_bridge, grace_bridge
from ...io.structures_dataset_io import StructureDatasetIO

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"
AVAILABLE_MODELS = {'mattersim': mattersim_bridge, 'grace': grace_bridge}

class NNEstimator_Test(unittest.TestCase):

    def setUp(self):
        self.poscars = PATH_WITH_DATASETS / "POSCARS"
        self.dataset = StructureDatasetIO(self.poscars).load_from_directory()
        for name, bridge in AVAILABLE_MODELS.items():
            NNEstimator.register_model(name, bridge)
        self.estimator = NNEstimator()

    def test_single_entry_estimation(self):
        print(self.dataset[2])
        self.assertAlmostEqual(self.estimator.estimate_entry_energy(self.dataset[2]), -44.0794, places=3)
        print(self.estimator.estimate_entry_energy(self.dataset[2], model='grace'))

    def test_batch_estimation(self):
        ds2 = self.estimator.estimate_dataset_energies(self.dataset)
        print('\n'.join([str(e.energy) for e in ds2]))

    def test_batch_relaxation(self):
        ds2 = self.estimator.relax_dataset(self.dataset)
