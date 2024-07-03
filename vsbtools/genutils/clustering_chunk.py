import numpy as np
from collections import deque
import pickle as pkl


def make_dist_matrix(points, dist_fun, dist_matrix_out_file=None):
    matrix_dim = len(points)
    dist_mat = np.zeros((matrix_dim, matrix_dim), dtype=float)
    adj_mat = np.zeros((matrix_dim, matrix_dim), dtype=bool)  # Boolean adjacency mtrx in FP-space
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
