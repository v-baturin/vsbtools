from __future__ import annotations
import time
from pathlib import Path
from typing import Callable, Any, Tuple, Iterable, TypeVar
from ..crystal_entry import  CrystalEntry
from ..crystal_dataset import CrystalDataset
from ..genutils.duplicate_analysis import remove_duplicates

DEBUG = False
T = TypeVar("T")


def describe_similarity_tool(similarity_tk: "SimilarityTools") -> str:
    dist_fn = getattr(similarity_tk, "dist", None)
    owner = getattr(dist_fn, "__self__", None)
    if owner is not None:
        return owner.__class__.__name__
    return getattr(dist_fn, "__qualname__", repr(dist_fn))


class SimilarityTools:
    def __init__(
            self,
            dist_fn: Callable[[CrystalEntry, CrystalEntry], float],
            tol_FP: float = 0.08,
            tol_fp: float | None = None,
    ) -> None:
        if tol_fp is not None:
            tol_FP = tol_fp
        self.tol_FP = tol_FP
        self.dist = dist_fn

    def is_duplicate(self, a: CrystalEntry, b: CrystalEntry) -> bool:
        """Fingerprint-distance based criterion"""
        return self.dist(a,b) < self.tol_FP

    def get_unseen_successively(
            self,
            items: Iterable[T],
            reference_items: Iterable[T] | None = None,
            *,
            entry_factory: Callable[[T], CrystalEntry] | None = None,
            tol_FP: float | None = None,
            match_formula: bool = True,
            progress: Callable[[str], None] | None = None,
            progress_prefix: str | None = None,
            progress_interval_s: float = 0.5,
            on_duplicate: Callable[[T, T], None] | None = None,
            on_slow_distance: Callable[[int, T], None] | None = None,
            skip_errors: bool = False,
    ) -> list[T]:
        """Return items that are unseen in the accepted reference set."""
        entry_factory = entry_factory or (lambda item: item)
        items = list(items)
        accepted = list(reference_items or [])
        unseen = []
        tol_FP = self.tol_FP if tol_FP is None else tol_FP
        last_progress = 0.0

        for item_no, item in enumerate(items, start=1):
            now = time.monotonic()
            if (
                    progress is not None
                    and progress_prefix is not None
                    and (
                            item_no == 1
                            or item_no == len(items)
                            or now - last_progress >= progress_interval_s
                    )
            ):
                progress(f"{progress_prefix}: merging entry {item_no}/{len(items)}")
                last_progress = now
            entry = entry_factory(item)
            duplicate = None
            for reference in accepted:
                reference_entry = entry_factory(reference)
                if match_formula and entry.formula and reference_entry.formula and entry.formula != reference_entry.formula:
                    continue
                try:
                    tic = time.time()
                    dist = self.dist(entry, reference_entry)
                except Exception:
                    if skip_errors:
                        continue
                    raise
                if time.time() - tic > 10 and on_slow_distance is not None:
                    on_slow_distance(item_no, item)
                if dist <= tol_FP:
                    duplicate = reference
                    break
            if duplicate is not None:
                if on_duplicate is not None:
                    on_duplicate(item, duplicate)
                continue
            unseen.append(item)
            accepted.append(item)

        return unseen

    def contains_structure(self, entry: CrystalEntry, ds: CrystalDataset) -> Tuple[list, list]:
        true_idcs = [i for i, ds_entry in enumerate(ds) if self.dist(entry, ds_entry) <= self.tol_FP]
        return true_idcs, [ds[i].id for i in true_idcs]

    def deduplicate(self, ds: CrystalDataset,
                    check_clusters_file=False, clusters_file: Path = None, check_dist_matrix_file=False,
                    dist_matrix_file=None, fitness_list=None, tol_FP: float = None,
                    enforce_compositions_separation=True,
                    save_clusters_file: bool = False,
                    save_dist_matrix_file: bool = False,
                    **kwargs) -> tuple[CrystalDataset, Any, Any]:

        """
        Remove duplicates from the dataset using USPEX's remove_duplicates function.
        :param check_clusters_file: If True, try to read clusters from a file.
        :param check_dist_matrix_file: If True, try to read distance matrix from a file.
        :param save_clusters_file: If True, write clusters to a file.
        :param save_dist_matrix_file: If True, write distance matrix to a file.
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

        if (check_clusters_file or save_clusters_file) and clusters_file is None:
            clusters_file = ds.base_path / f"{ds.dataset_id}_clusters.pkl"
        if (check_dist_matrix_file or save_dist_matrix_file) and dist_matrix_file is None:
            dist_matrix_file = ds.base_path / f"{ds.dataset_id}_dist_matrix.pkl"
        best_representatives, clusters, best_idx = remove_duplicates(ds, dist_fn=self.dist,
                                                                     fitness_list=fitness_list,
                                                                     intercluster_mindistance=tol_FP,
                                                                     check_clusters_file=check_clusters_file,
                                                                     check_dist_matrix_file=check_dist_matrix_file,
                                                                     dist_matrix_file=dist_matrix_file,
                                                                     clusters_file=clusters_file,
                                                                     do_split_clusters_by_labels=enforce_compositions_separation,
                                                                     labels_list=reduced_compositions,
                                                                     save_clusters_file=save_clusters_file,
                                                                     save_dist_matrix_file=save_dist_matrix_file)
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
        duplicates_counter = set()
        reproduced = set()

        def on_duplicate(examined: CrystalEntry, reference: CrystalEntry):
            reference_idx = ref_index_by_identity.get(id(reference))
            if reference_idx is not None:
                duplicates_counter.add(reference_idx)
                reproduced.add(reference.id)
                if not "duplicates" in reference.metadata:
                    reference.metadata["duplicates"] = {examined.id}
                else:
                    reference.metadata["duplicates"].add(examined.id)
            examined.metadata["known"] = True
            if DEBUG:
                print(f"{examined.id} in {examined.metadata['source']} is a duplicate of {reference.id} in {reference.metadata['source']}")

        def on_slow_distance(i: int, examined: CrystalEntry):
            print(f"Fingerprint for i = {i - 1} Took too long")
            examined.structure.to_file(f"PROBLEM_POSCAR{i - 1}", fmt='poscar')

        ref_index_by_identity = {id(entry): i for i, entry in enumerate(ref_ds)}
        for examined in ds:
            examined.metadata["known"] = False
        new_entries = self.get_unseen_successively(
            ds,
            ref_ds,
            tol_FP=tol_FP,
            on_duplicate=on_duplicate,
            on_slow_distance=on_slow_distance,
        )

        reproducibility = len(duplicates_counter) / len(ref_ds) if len(ref_ds) > 0 else 1.
        message = (f"{duplicates_counter} out of {len(ds)} in dataset {ds.dataset_id} are duplicates of {ref_ds.dataset_id}\n"
                   f"{ds.dataset_id} reproduces {reproducibility:.2%} of {ref_ds.dataset_id}")
        res = CrystalDataset.from_parents(new_entries, parents=(ref_ds, ds), message=message)
        res.metadata["reproducibility"] = reproducibility
        res.metadata["reproduced"] = list(reproduced)
        return res
