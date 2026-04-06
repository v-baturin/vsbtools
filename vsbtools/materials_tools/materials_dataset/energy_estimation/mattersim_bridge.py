import numpy as np
from ..crystal_dataset import CrystalDataset, CrystalEntry
from ...NN_energy_estimators import mattersim_client

def estimate_batch(dataset: CrystalDataset, force_gpu: int | None = None, **kwargs):
    with mattersim_client.EnergyStream(force_gpu=force_gpu) as es:
        energies = np.zeros(len(dataset))
        for i, e in enumerate(dataset):
            try:
                energies[i] = es.calc(e.structure.to_ase_atoms())
            except RuntimeError as err:
                print (f"Couldn't get energy for {e.id}, set to infinity")
                energies[i] = np.inf
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
