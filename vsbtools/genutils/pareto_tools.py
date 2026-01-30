import numpy as np


# Faster than is_pareto_efficient_simple, but less readable.
def pareto_gen(costs, return_mask=False, max_front = None):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    i_front= 0
    init_cost = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    if max_front is None:
        max_front = np.inf
     # Next index in the is_efficient array to search for
    residue = costs
    while len(residue) > 0 and i_front < max_front:
        next_point_index = 0
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
            if i_front == 57:
                pass
        else:
            yield is_efficient
            if i_front == 57:
                pass
        if residue:
            costs = np.concatenate(residue)
            init_cost = np.hstack(non_efficient)
            i_front += 1

def pareto_subdataframe_indices(df, cols, max_front=None):
    """
    Non-dominated sorting (all objectives in `cols` are minimized).

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.
    cols : list of str
        Columns to use as objectives, e.g. ['a0', 'a1', ..., 'aN'].
    max_front : int
        Highest pareto front to output

    Returns
    -------
    fronts_idx : list of np.ndarray
        fronts_idx[0] is the index positions (iloc) of the first Pareto front,
        fronts_idx[1] is the second front, etc.
    rank : np.ndarray of int, shape (len(df),)
        rank[i] is the Pareto front number (1, 2, 3, ...) of row i.
    """
    data = df[cols].to_numpy()
    n = data.shape[0]

    remaining = np.arange(n)       # positions 0..n-1
    fronts_idx = []
    rank = np.empty(n, dtype=int)  # will hold front number for each point

    current_front = 1
    while remaining.size > 0 and (max_front is None or current_front <= max_front):
        M = data[remaining]  # currently remaining objectives

        # i dominates j if i <= j in all dims and < in at least one (minimization)
        all_le = np.all(M[:, None, :] <= M[None, :, :], axis=2)
        any_lt = np.any(M[:, None, :] < M[None, :, :], axis=2)
        dominates = all_le & any_lt

        # j is dominated if exists i that dominates j
        is_dominated = dominates.any(axis=0)

        # non-dominated subset among 'remaining'
        front = remaining[~is_dominated]

        fronts_idx.append(front)
        rank[front] = current_front

        # keep only dominated points for next iteration
        remaining = remaining[is_dominated]
        current_front += 1

    return fronts_idx, rank

def add_pf_idx(df, col1, col2, pf_col="pf_idx"):
    """
    Add Pareto front index for 2D minimization (col1, col2).

    Pareto dominance definition (strong dominance):
      Point j is dominated if there exists i such that:
          df[col1][i] < df[col1][j]  AND  df[col2][i] < df[col2][j]

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe (modified in place).
    col1, col2 : str
        Names of the two columns to minimize.
    pf_col : str, default "pf_idx"
        Name of the output column with Pareto front numbers.

    Returns
    -------
    df : pandas.DataFrame
        Same dataframe, with an extra integer column `pf_col`
        where 1 = first (best) front, 2 = second, etc.
    """
    # Extract the 2D data as a NumPy array (in current row order)
    vals = df[[col1, col2]].to_numpy()
    n = vals.shape[0]

    remaining = np.arange(n)             # positions 0..n-1 (iloc-style)
    pf_idx = np.empty(n, dtype=int)      # Pareto front index for each row
    current_front = 1

    while remaining.size > 0:
        M = vals[remaining]              # shape (k, 2) for currently remaining points

        # i dominates j if both coordinates are strictly smaller (strong dominance)
        better_1 = M[:, None, 0] < M[None, :, 0]
        better_2 = M[:, None, 1] < M[None, :, 1]
        dominates = better_1 & better_2  # (i, j) True when i dominates j

        # j is dominated if any i dominates j
        is_dominated = dominates.any(axis=0)

        # non-dominated among remaining form the current front
        front = remaining[~is_dominated]
        pf_idx[front] = current_front

        # prepare for next front
        remaining = remaining[is_dominated]
        current_front += 1

    # Add / overwrite the Pareto-front column (position-wise assignment)
    df[pf_col] = pf_idx
    return df


if __name__ == '__main__':
    from matplotlib import pyplot as plt

    # pts = np.array([[1, 4], [4, 4], [3, 3], [4, 2], [5, 1], [2, 1], [5, 3]])
    pts = 10 * np.random.rand(120, 2)
    pfronts = pareto_gen(pts)
    plt.plot(pts[:, 0], pts[:, 1], 'o')
    for pf in pfronts:
        curr_front = pts[pf]
        sorted_front = curr_front[curr_front[:, 0].argsort()]
        plt.plot(sorted_front[:, 0], sorted_front[:, 1], '-')
    plt.show()
