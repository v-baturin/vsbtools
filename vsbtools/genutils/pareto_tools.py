import numpy as np


# Faster than is_pareto_efficient_simple, but less readable.
def pareto_gen(costs, return_mask=False):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    init_cost = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    residue = costs
    while len(residue) > 0:
        residue = []
        non_efficient = []
        is_efficient = init_cost
        while next_point_index < len(costs):
            nondominated_point_mask = np.any(costs < costs[next_point_index], axis=1)
            nondominated_point_mask[next_point_index] = True
            non_efficient.append(is_efficient[np.logical_not(nondominated_point_mask)])
            residue.append(costs[np.logical_not(nondominated_point_mask)])
            is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
            costs = costs[nondominated_point_mask]
            next_point_index = np.sum(nondominated_point_mask[:next_point_index]) + 1
        if return_mask:
            is_efficient_mask = np.zeros(n_points, dtype=bool)
            is_efficient_mask[is_efficient] = True
            yield is_efficient_mask
        else:
            yield is_efficient
        if residue:
            costs = np.concatenate(residue)
            init_cost = np.hstack(non_efficient)


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    # pts = np.array([[1, 4], [4, 4], [3, 3], [4, 2], [5, 1], [2, 1], [5, 3]])
    pts = 10 * np.random.rand(120, 2)
    pfronts = pareto_gen(pts)
    plt.plot(pts[:, 0], pts[:, 1], 'o')
    for pf in pfronts:
        current_front = pts[pf]
        sorted_front = current_front[current_front[:, 0].argsort()]
        plt.plot(sorted_front[:, 0], sorted_front[:, 1], '-')
    plt.show()
