from pathlib import Path
import unittest

import numpy as np

from ..nn_estimator import NNEstimator
from .. import mattersim_bridge, grace_bridge
from ...crystal_entry import CrystalEntry
from ...io.structures_dataset_io import StructureDatasetIO, _safe_structure_from_file

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"
AVAILABLE_ESTIMATORS = {"mattersim": mattersim_bridge, "grace": grace_bridge}

class NNEstimator_Test(unittest.TestCase):

    def setUp(self):
        # self.poscars = PATH_WITH_DATASETS / "POSCARS"
        self.poscars = PATH_WITH_DATASETS / "two_systems_MgAlO"
        self.chosen_fname = 'systemRDU2.POSCAR'
        self.test_entry = CrystalEntry(id = '0',
                                       structure=_safe_structure_from_file(self.poscars / self.chosen_fname))
        self.entry_energies = {"mattersim": -44.0794, "grace": -44.1207}
        self.dataset = StructureDatasetIO(self.poscars).load_from_directory()
        for name, bridge in AVAILABLE_ESTIMATORS.items():
            NNEstimator.register_model(name, bridge)
        self.estimator = NNEstimator()

    def test_single_entry_estimation(self):
        for model, energy in self.entry_energies.items():
            self.assertAlmostEqual(self.estimator.estimate_entry_energy(self.test_entry, model=model), energy, places=3)

    # def test_batch_estimation(self):
    #     for model in AVAILABLE_ESTIMATORS:
    #         ds = self.estimator.estimate_dataset_energies(self.dataset, model_name=model)
    #         energies = np.array([str(e.energy) for e in ds], dtype=float)
    #         print(f"{model}_mean = {energies.mean()}")

    # def test_single_relaxation(self):


    def test_batch_relaxation(self):
        ds2 = self.estimator.relax_dataset(self.dataset, model_name='grace')
