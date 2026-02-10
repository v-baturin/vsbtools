import numpy as np
from ..crystal_dataset import CrystalDataset, CrystalEntry
from ...NN_energy_estimators import grace_client

def estimate_batch(dataset: CrystalDataset):
    with grace_client.EnergyStream() as es:
        energies = np.zeros(len(dataset))
        for i, e in enumerate(dataset):
            try:
                energies[i] = es.calc(e.structure.to_ase_atoms())
            except RuntimeError as err:
                print (f"Couldn't get energy for {e.id}, set to infinity")
                energies[i] = np.inf
    return energies

def estimate_entry_energy(e: CrystalEntry):
    return grace_client.get_energy(e.structure.to_ase_atoms())

def relax_batch(dataset: CrystalDataset):
    with grace_client.RelaxStream() as rs:
        new_structures = [None] * len(dataset)
        for i, e in enumerate(dataset):
            try:
                print(f"Now relaxing structure with id = {e.id} ({e.poscarname}) from dataset id = {dataset.dataset_id}")
                new_structures[i] = rs.relax(e.structure.to_ase_atoms())
            except RuntimeError as err:
                print (f"Couldn't relax structure with id = {e.id} ({e.poscarname}) from dataset id = {dataset.dataset_id}")
    return new_structures

def relax_entry(e: CrystalEntry):
    return grace_client.relax(e.structure.to_ase_atoms())