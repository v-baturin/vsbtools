from __future__ import annotations
import time
from pathlib import Path
from typing import Callable, Any, Tuple
from ..crystal_entry import  CrystalEntry
from ..crystal_dataset import CrystalDataset
from ....genutils.duplicate_analysis import remove_duplicates

DEBUG = False

class SimilarityTools:
    def __init__(self, dist_fn: Callable[[CrystalEntry, CrystalEntry], float], tol_FP: float = 0.08) -> None:
        self.tol_FP = tol_FP
        self.dist = dist_fn

    def is_duplicate(self, a: CrystalEntry, b: CrystalEntry) -> bool:
        """Fingerprint-distance based criterion"""
        return self.dist(a,b) < self.tol_FP

    def contains_structure(self, entry: CrystalEntry, ds: CrystalDataset) -> Tuple[list, list]:
        true_idcs = [i for i, ds_entry in enumerate(ds) if self.dist(entry, ds_entry) <= self.tol_FP]
        return true_idcs, [ds[i].id for i in true_idcs]

    def deduplicate(self, ds: CrystalDataset,
                    check_clusters_file=False, clusters_file: Path = None, check_dist_matrix_file=False,
                    dist_matrix_file=None, fitness_list=None, tol_FP: float = None,
                    enforce_compositions_separation=True, **kwargs) -> tuple[CrystalDataset, Any, Any]:

        """
        Remove duplicates from the dataset using USPEX's remove_duplicates function.
        :param check_clusters_file: If True, will write clusters to a file.
        :param check_dist_matrix_file: If True, will write distance matrix to a file.
        :param tol_FP: Tolerance for fingerprint distance.
        :param enforce_compositions_separation: If True, will enforce separation of compositions in clusters.
        :param fitness_list: List of fitness values (e.g., energies) for each entry.
        """
        tol_FP = self.tol_FP if tol_FP is None else tol_FP
        try:
            fitness_list = fitness_list or [e.energy / e.natoms for e in ds]
        except TypeError:
            fitness_list = None

        if enforce_compositions_separation:
            reduced_compositions = [e.composition.reduced_formula for e in ds]
        else:
            reduced_compositions = None

        clusters_file = clusters_file or ds.base_path / f"{ds.dataset_id}_clusters.pkl"
        dist_matrix_file = dist_matrix_file or ds.base_path / f"{ds.dataset_id}_dist_matrix.pkl"
        best_representatives, clusters, best_idx = remove_duplicates(ds, dist_fn=self.dist,
                                                                     fitness_list=fitness_list,
                                                                     intercluster_mindistance=tol_FP,
                                                                     check_clusters_file=check_clusters_file,
                                                                     check_dist_matrix_file=check_dist_matrix_file,
                                                                     dist_matrix_file=dist_matrix_file,
                                                                     clusters_file=clusters_file,
                                                                     do_split_clusters_by_labels=enforce_compositions_separation,
                                                                     labels_list=reduced_compositions)
        for no, idx in enumerate(best_idx):
            if "duplicates" not in ds[idx].metadata or ds[idx].metadata["duplicates"] is None:
                ds[idx].metadata["duplicates"] = set([ds[j].id for j in clusters[no] if j != idx])
            else:
                ds[idx].metadata["duplicates"] |= set([ds[j].id for j in clusters[no] if j != idx])
        filtered_list = [ds[i] for i in best_idx]

        message = f"Parent deduplicated with tol_FP={tol_FP} "
        return CrystalDataset.from_parents(filtered_list, parents=(ds,), message=message, **kwargs), clusters, best_idx


    def get_unseen_in_ref(self, ds: CrystalDataset, ref_ds: CrystalDataset, tol_FP=None):
        tol_FP = self.tol_FP if tol_FP is None else tol_FP
        new_entries = []
        duplicates_counter = set()
        reproduced = set()
        rho = self.dist
        for i, examined in enumerate(ds):
            for j, reference in enumerate(ref_ds):
                tic = time.time()
                dist = rho(examined, reference)
                duration = time.time() - tic
                if duration > 10:
                    print(f"Fingerprint for i = {i} Took too long")
                    examined.structure.to_file(f"PROBLEM_POSCAR{i}", fmt='vasp')
                if dist <= tol_FP:
                    duplicates_counter.add(j)
                    reproduced.add(reference.id)
                    if not "duplicates" in reference.metadata:
                        reference.metadata["duplicates"] = {examined.id}
                    else:
                        reference.metadata["duplicates"].add(examined.id)
                    if DEBUG:
                        print(f"{examined.id} in {examined.metadata['source']} is a duplicate of {reference.id} in {reference.metadata['source']}")
                    break
            else:
                new_entries.append(ds[i])

        reproducibility = len(duplicates_counter) / len(ref_ds) if len(ref_ds) > 0 else 1.
        message = (f"{duplicates_counter} out of {len(ds)} in dataset {ds.dataset_id} are duplicates of {ref_ds.dataset_id}\n"
                   f"{ds.dataset_id} reproduces {reproducibility:.2%} of {ref_ds.dataset_id}")
        res = CrystalDataset.from_parents(new_entries, parents=(ref_ds, ds), message=message)
        res.metadata["reproducibility"] = reproducibility
        res.metadata["reproduced"] = list(reproduced)
        return res