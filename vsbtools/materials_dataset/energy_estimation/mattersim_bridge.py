import numpy as np
import time
from ..crystal_dataset import CrystalDataset, CrystalEntry
from ..NN_energy_estimators import mattersim_client

def estimate_batch(dataset: CrystalDataset, force_gpu: int | None = None, **kwargs):
    progress_interval_s = kwargs.get("progress_interval_s", 0.5)
    estimator_name = kwargs.get("estimator_name", "mattersim")
    with mattersim_client.EnergyStream(force_gpu=force_gpu) as es:
        n_entries = len(dataset)
        energies = np.zeros(n_entries)
        last_progress = 0.0
        status_len = 0
        for i, e in enumerate(dataset):
            try:
                source = e.metadata.get("file", e.metadata.get("source", "unknown source"))
                status = f"Now estimating energy with {estimator_name} {i + 1}/{n_entries} of entry with id={e.id} from {source}"
                now = time.monotonic()
                if now - last_progress >= progress_interval_s:
                    print(f"\r{status}{' ' * max(0, status_len - len(status))}", end="", flush=True)
                    last_progress, status_len = now, len(status)
                energies[i] = es.calc(e.structure.to_ase_atoms())
            except RuntimeError as err:
                print()
                print(f"Couldn't get energy for {e.id}, set to infinity")
                energies[i] = np.inf
        if n_entries:
            status = f"Now estimating energy with {estimator_name} {n_entries}/{n_entries}"
            print(f"\r{status}{' ' * max(0, status_len - len(status))}")
    return energies

def estimate_entry_energy(e: CrystalEntry, force_gpu: int | None = None, **kwargs):
    return mattersim_client.get_energy(e.structure.to_ase_atoms(), force_gpu=force_gpu)

def relax_batch(dataset: CrystalDataset, force_gpu: int | None = None, **kwargs):
    with mattersim_client.RelaxStream(force_gpu=force_gpu) as rs:
        new_structures = [None] * len(dataset)
        for i, e in enumerate(dataset):
            try:
                new_structures[i] = rs.relax(e.structure.to_ase_atoms())
            except RuntimeError as err:
                print (f"Couldn't relax structure with id = {e.id} ({e.poscarname}) from dataset id = {dataset.dataset_id}")
    return new_structures

def relax_entry(e: CrystalEntry, force_gpu: int | None = None, **kwargs):
    return mattersim_client.relax(e.structure.to_ase_atoms(), force_gpu=force_gpu)
