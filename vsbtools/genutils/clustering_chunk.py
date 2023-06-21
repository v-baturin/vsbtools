import numpy as np
from collections import deque


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
    matrix_dim = len(fingerprints)

    adj_mat = np.zeros((matrix_dim, matrix_dim), dtype=bool)  # Boolean adjacency mtrx in FP-space
    for i in range(matrix_dim):
        for j in range(matrix_dim):
            adj_mat[i, j] = rho(fingerprints[i], fingerprints[j]) < delta

    clusters = []
    rest = list(range(matrix_dim))
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
