from dataclasses import dataclass
from typing import Callable, ClassVar, Dict
import types
from ..crystal_dataset import CrystalDataset, CrystalEntry
from ...geom_utils.structure_checks import check_structure_sanity_ase
from pymatgen.core import Structure


@dataclass(slots=True)
class NNEstimator:

    known_models: ClassVar[Dict[str, types.ModuleType]] = {}
    default_model: str = "mattersim"

    def estimate_dataset_energies(self, dataset: CrystalDataset, model_name: str | None = None):
        model_name = self.default_model if model_name is None else model_name
        if model_name not in self.known_models:
            raise ValueError(f"Unknown model: '{model_name}', available: {list(self.known_models.keys())}")
        estimator = self.known_models[model_name]
        new_energies = estimator.estimate_batch(dataset)
        new_entries = []
        for entry, estimation in zip(dataset, new_energies):
            new_entries.append(entry.copy_with(**{"energy": estimation}))
        msg = f"Energies estimated with {model_name}"
        return CrystalDataset.from_parents(new_entries, (dataset,), message=msg)

    def estimate_entry_energy(self, e: CrystalEntry, model: str | None = "mattersim"):
        return self.known_models[model].estimate_entry_energy(e)

    def relax_dataset(self, dataset: CrystalDataset, model_name: str | None = None):
        model_name = self.default_model if model_name is None else model_name
        if model_name not in self.known_models:
            raise ValueError(f"Unknown model: '{model_name}', available: {list(self.known_models.keys())}")
        estimator = self.known_models[model_name]
        new_structures = estimator.relax_batch(dataset)
        new_entries = []
        for entry, atoms in zip(dataset, new_structures):
            md = entry.metadata
            if atoms is None:
                is_OK = False
                info = "Relaxation_failed"
            else:
                is_OK, info = check_structure_sanity_ase(atoms)
            if is_OK:
                new_entries.append(entry.copy_with(**{"structure": Structure.from_ase_atoms(atoms),
                                                      "metadata": {**md, **{'relaxation_done': True}}}))
            else:
                print(info)
                new_entries.append(entry.copy_with(**{"metadata": {**md, **{'relaxation_done': False}}}))

        msg = f"Structures relaxed with {model_name}"
        return CrystalDataset.from_parents(new_entries, (dataset,), message=msg)

    def relax_entry(self, e: CrystalEntry, model: str | None = "mattersim"):
        return self.known_models[model].relax_entry(e)

    @classmethod
    def register_model(cls, name: str, module: types.ModuleType):
        cls.known_models[name] = module
