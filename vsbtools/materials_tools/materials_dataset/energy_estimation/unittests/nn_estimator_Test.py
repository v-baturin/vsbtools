import os
from pathlib import Path
import types
import unittest

import numpy as np

from ..nn_estimator import NNEstimator
from .. import mattersim_bridge, grace_bridge
from ...crystal_entry import CrystalEntry
from ...io.structures_dataset_io import StructureDatasetIO, _safe_structure_from_file
from ...io import write

PATH_WITH_TESTS = Path(__file__).parent
PATH_WITH_DATASETS = PATH_WITH_TESTS / "../../unittests_datasets"
AVAILABLE_ESTIMATORS = {"mattersim": mattersim_bridge, "grace": grace_bridge}

class NNEstimator_Test(unittest.TestCase):

    def setUp(self):
        # self.poscars = PATH_WITH_DATASETS / "POSCARS"
        self.poscars = PATH_WITH_DATASETS / "two_systems_MgAlO"
        self.chosen_fname = 'systemRDU2_POSCAR'
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
        ds2 = self.estimator.relax_dataset(self.dataset, model_name='grace', force_gpu=0)
        relaxed_path = PATH_WITH_TESTS / 'relaxed_grace'
        os.makedirs(relaxed_path, exist_ok=True)
        write(ds2,relaxed_path)

    def test_kwargs_forwarded_to_model_bridge(self):
        calls = {}

        def estimate_batch(dataset, **kwargs):
            calls["estimate_batch"] = kwargs
            return np.zeros(len(dataset))

        def estimate_entry_energy(entry, **kwargs):
            calls["estimate_entry_energy"] = kwargs
            return -1.0

        def relax_batch(dataset, **kwargs):
            calls["relax_batch"] = kwargs
            return [e.structure.to_ase_atoms() for e in dataset]

        def relax_entry(entry, **kwargs):
            calls["relax_entry"] = kwargs
            return entry.structure.to_ase_atoms()

        bridge = types.SimpleNamespace(
            estimate_batch=estimate_batch,
            estimate_entry_energy=estimate_entry_energy,
            relax_batch=relax_batch,
            relax_entry=relax_entry,
        )
        NNEstimator.register_model("dummy_forward", bridge)
        est = NNEstimator(default_model="dummy_forward")

        est.estimate_dataset_energies(self.dataset, force_gpu=1)
        est.estimate_entry_energy(self.test_entry, model="dummy_forward", force_gpu=1)
        est.relax_dataset(self.dataset, model_name="dummy_forward", force_gpu=1)
        est.relax_entry(self.test_entry, model="dummy_forward", force_gpu=1)

        self.assertEqual(calls["estimate_batch"]["force_gpu"], 1)
        self.assertEqual(calls["estimate_entry_energy"]["force_gpu"], 1)
        self.assertEqual(calls["relax_batch"]["force_gpu"], 1)
        self.assertEqual(calls["relax_entry"]["force_gpu"], 1)
