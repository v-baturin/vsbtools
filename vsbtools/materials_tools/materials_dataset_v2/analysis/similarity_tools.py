from __future__ import annotations
from pathlib import Path
from typing import Callable, Any, Tuple
from ..crystal_entry import  CrystalEntry
from ..crystal_dataset import CrystalDataset
from my_packages.genutils.duplicate_analysis import remove_duplicates

class SimilarityTools:
    def __init__(self, dist_fn: Callable[[CrystalEntry, CrystalEntry], float], tol_fp: float = 0.08) -> None:
        self.tol_fp = tol_fp
        self.dist = dist_fn

    def is_duplicate(self, a: CrystalEntry, b: CrystalEntry) -> bool:
        """Fingerprint-distance based criterion"""
        return self.dist(a,b) < self.tol_fp

    def contains_structure(self, entry: CrystalEntry, ds: CrystalDataset) -> Tuple[list, list]:
        true_idcs = [i for i, ds_entry in enumerate(ds) if self.dist(entry, ds_entry) <= self.tol_fp]
        return true_idcs, [ds[i].id for i in true_idcs]

    def deduplicated(self, ds: CrystalDataset,
                     check_clusters_file=False, clusters_file: Path = None, check_dist_matrix_file=False,
                     dist_matrix_file=None, fitness_list=None,
                     enforce_compositions_separation=False, **kwargs) -> tuple[CrystalDataset, Any, Any]:

        """
        Remove duplicates from the dataset using USPEX's remove_duplicates function.
        :param check_clusters_file: If True, will write clusters to a file.
        :param check_dist_matrix_file: If True, will write distance matrix to a file.
        :param tol_FP: Tolerance for fingerprint distance.
        :param enforce_compositions_separation: If True, will enforce separation of compositions in clusters.
        :param fitness_list: List of fitness values (e.g., energies) for each entry.
        """
        try:
            fitness_list = fitness_list or [e.energy / e.natoms for e in ds]
        except TypeError:
            fitness_list = None

        if enforce_compositions_separation:
            reduced_compositions = [e.composition.reduced_formula for e in ds]
        else:
            reduced_compositions = None

        clusters_path = clusters_file or ds.base_path / f"{ds.dataset_id}_clusters.pkl"
        dist_matrix_path = dist_matrix_file or ds.base_path / f"{ds.dataset_id}_dist_matrix.pkl"
        best_representatives, clusters, best_idx = remove_duplicates(ds, dist_fn=self.dist,
                                                                     fitness_list=fitness_list,
                                                                     intercluster_mindistance=self.tol_fp,
                                                                     check_clusters_file=check_clusters_file,
                                                                     check_dist_matrix_file=check_dist_matrix_file,
                                                                     dist_matrix_file=dist_matrix_path,
                                                                     clusters_file=clusters_path,
                                                                     do_split_clusters_by_labels=enforce_compositions_separation,
                                                                     labels_list=reduced_compositions)
        filtered_list = [ds[i] for i in best_idx]
        message = f"Parent deduplicated with tol_FP={self.tol_fp} "
        return CrystalDataset.from_parents(filtered_list, parents=(ds.dataset_id,), message=message, **kwargs), clusters, best_idx


    def get_unseen_in_ref(self, ds: CrystalDataset, ref_ds: CrystalDataset, verbose=True):
        new_entries = []
        duplicates_counter = set()
        rho = self.dist
        for i, e in enumerate(ds):
            for j, e2 in enumerate(ref_ds):
                if rho(e, e2) <= self.tol_fp:
                    duplicates_counter.add(j)
                    break
            else:
                new_entries.append(ds[i])

        reproductibility = len(duplicates_counter) / len(ref_ds)
        if verbose:
            print(f"{reproductibility:.2%} of ref present in ds")
        message = f"Structures of parent1 unseen in parent2. parent1 has {reproductibility:.2%} of parent2"
        return CrystalDataset.from_parents(new_entries, parents=(ref_ds, ds), message=message)