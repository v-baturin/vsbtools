import os
import numpy as np
from collections import deque
import pickle as pkl
from datetime import datetime
tdy = datetime.today().strftime('%Y%m%d')

def make_dist_matrix(points, dist_fun, dist_matrix_out_file=None):
    matrix_dim = len(points)
    dist_mat = np.zeros((matrix_dim, matrix_dim), dtype=float)
    for i in range(matrix_dim):
        for j in range(i, matrix_dim):
            print(f"i = {i}, j = {j}")
            dist_mat[i, j] = dist_mat[j, i] = dist_fun(points[i], points[j])
    if dist_matrix_out_file is not None:
        with open(dist_matrix_out_file, 'wb') as dm_file:
            pkl.dump(dist_mat, dm_file)
            print('dist matrix saved')
    return dist_mat

def adj_matrix_from_dist(dist_matrix, delta):
    return dist_matrix <= delta

def get_clusters_from_adj_mat(adj_mat):
    clusters = []
    rest = list(range(len(adj_mat)))
    while rest:
        to_visit = deque([rest[0]])
        curr_cluster = []
        while to_visit:
            this_node = to_visit.pop()
            if this_node not in curr_cluster:
                adj_mat[:, this_node] = False
                curr_cluster.append(this_node)
                rest.remove(this_node)
                next_to_visit = np.where(adj_mat[this_node])[0]
                to_visit.extendleft(next_to_visit)
        clusters.append(curr_cluster)
    return clusters


def clusterize_dist_matrix(dist_matrix=None, clusters_out_file=None, separation_labels=None, tolFP=0.16, **kwargs):
    adj_mat = dist_matrix <= tolFP
    clusters = get_clusters_from_adj_mat(adj_mat)
    if separation_labels:
        clusters = separate_by_labels(clusters, separation_labels)
    if clusters_out_file is None:
        clusters_out_file = f'{os.getcwd()}/clusters_{tdy}.pkl'
    with open(clusters_out_file, 'wb') as cl_file:
        pkl.dump(clusters, cl_file)
        print(f'clustering saved in {clusters_out_file}')
    return clusters

def separate_by_labels(clusters, labels_list):
    new_clusters  = []
    print(f"Separating {len(clusters)} clusters by labels list")
    for c in clusters:
        class_comps = [labels_list[i] for i in c]
        eq_classes = [[c[i] for i, v in enumerate(class_comps) if v == key]
                      for key in dict.fromkeys(class_comps)]
        new_clusters.extend(eq_classes)
    print(f"Separated into {len(new_clusters)} clusters to have same labels per clusters")
    return new_clusters

def select_best_representatives(clusters, entries, fitness_list = None, max_fitness_delta = np.inf, **kwargs):
    print(f"Processing {len(clusters)} clusters")
    if fitness_list is not None:
        fitness_list = np.array(fitness_list)
        ref_fitness = np.min(fitness_list)
    else:
        fitness_list = np.zeros(len(entries), dtype=float)
        ref_fitness = 0.
    good_fitnesses = []
    i = 0
    best_of_each = []
    global_indices_of_best = []
    for cluster in clusters:
        cl_en = fitness_list[cluster]
        min_en = np.min(cl_en)
        idx_min = cluster[np.where(cl_en == np.min(cl_en))[0][0]]
        struct = entries[idx_min]
        if min_en is not None and ref_fitness is not None and (min_en - ref_fitness > max_fitness_delta):
            continue
        if hasattr(struct, 'setProperty'):
            struct.setProperty("ID", i)
        good_fitnesses.append(min_en)
        best_of_each.append(struct)
        i += 1
        global_indices_of_best.append(idx_min)
    print(f"{len(best_of_each)} good ones. Writing")

    sorting_idx = np.argsort(good_fitnesses)
    sorted_best = []
    for i in sorting_idx:
        sorted_best.append(best_of_each[i])
    return sorted_best, global_indices_of_best


def clustering_by_dist(fingerprints, rho, delta):
    """
    Clustering
    splits set of nodes 0...N into clusters based on the distance mtrx and the threshold distance,
    below which the nodes are considered as belonging to one
    @param fingerprints: list of fingerprint objects
    @param rho: a function object calculating metrics in FP space
    @param delta: fingerprint tolerance
    @return: list of lists of classes
    :Authors:
    Vladimir Baturin <v.baturin@skoltech.ru>
    """

    return get_clusters_from_adj_mat(adj_matrix_from_dist(make_dist_matrix(fingerprints, rho), delta))
