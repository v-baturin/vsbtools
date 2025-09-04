import warnings
import os
from pathlib import Path
from typing import Callable, Any
import numpy as np
import pickle as pkl
from .clustering_chunk import make_dist_matrix, clusterize_dist_matrix, separate_by_labels, select_best_representatives
from datetime import datetime

tdy = datetime.today().strftime('%Y%m%d')

def remove_duplicates(entries, dist_fn: Callable[[Any, Any], float],
                      intercluster_mindistance=0.08,
                      fitness_list=None, threshold=np.inf,
                      data_path: Path | None = None, clusters_file=None, dist_matrix_file=None,
                      check_clusters_file=False, check_dist_matrix_file=False, do_split_clusters_by_labels=False,
                      labels_list=None):
    """
    Abstract deduplication of entries
    """

    data_path = data_path or Path(os.getcwd())
    clusters_file = clusters_file or data_path / f"clusters_{tdy}.pkl"
    dist_matrix_file = dist_matrix_file or data_path / f"dist_mat_{tdy}.pkl"

    clusters = None
    if check_clusters_file:
        if Path(clusters_file).is_file():
            with open(clusters_file, 'rb') as f:
                clusters = pkl.load(f)
            warnings.warn(f"Clusters are loaded from {clusters_file}. tolFP is ignored")
        else:
            warnings.warn("Clusters are not loaded, will be created from the distance matrix.")

    # If clusters weren’t loaded, ensure we have a distance matrix
    if clusters is None:
        # Try loading existing distance matrix if requested
        if check_dist_matrix_file:
            if Path(dist_matrix_file).is_file():
                with open(dist_matrix_file, 'rb') as f:
                    dist_matrix = pkl.load(f)
                warnings.warn(f"Distance matrix is loaded from {dist_matrix_file}.")
            else:
                warnings.warn("Distance matrix is not loaded, will be created from the entries.")
                dist_matrix = make_dist_matrix(entries, dist_fn, dist_matrix_file)
        else:
            # Always create if not checking for a pre-existing file
            dist_matrix = make_dist_matrix(entries, dist_fn, dist_matrix_file)

        # Finally clusterize
        clusters = clusterize_dist_matrix(
            dist_matrix,
            clusters_out_file=clusters_file,
            tolFP=intercluster_mindistance,
            separation_labels=(labels_list if labels_list and do_split_clusters_by_labels else None)
        )

    best_representatives, best_idc = select_best_representatives(clusters, entries, fitness_list, threshold)
    return best_representatives, clusters, best_idc