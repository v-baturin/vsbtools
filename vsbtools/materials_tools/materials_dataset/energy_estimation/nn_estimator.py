from dataclasses import dataclass
from typing import Callable, ClassVar, Dict
import types
from ..crystal_dataset import CrystalDataset, CrystalEntry


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

    @classmethod
    def register_model(cls, name: str, module: types.ModuleType):
        cls.known_models[name] = module
