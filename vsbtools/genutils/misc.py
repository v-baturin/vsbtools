import numpy as np
import functools

def rgetattr(obj, attr, default=None):
    try: left, right = attr.split('.', 1)
    except: return getattr(obj, attr, default)
    return rgetattr(getattr(obj, left), right, default)

def rsetattr(obj, attr, val):
    try: left, right = attr.split('.', 1)
    except: return setattr(obj, attr, val)
    return rsetattr(getattr(obj, left), right, val)

def rhasattr(obj, attr):
    try: left, right = attr.split('.', 1)
    except: return hasattr(obj, attr)
    return rhasattr(getattr(obj, left), right)

def odometer(maxcounts, mincounts=None):
    maxcounts = np.array(maxcounts)
    if not mincounts:
        mincounts = 0 * maxcounts
    else:
        mincounts = np.array(mincounts)
    max_tmp = maxcounts - mincounts
    moduli = [np.prod(max_tmp[k:]) for k in range(1, len(max_tmp))]
    moduli = moduli + [1]

    total = np.prod(max_tmp)

    for k in range(0, total):
        q, p = divmod(k, moduli[0])
        r = [q]
        for i in range(1, len(max_tmp)):
            q, p = divmod(p, moduli[i])
            r.append(q)
        yield np.array(r) + mincounts

def get_sorted_compositions(master_dict):
    """
    Utility for getting sorted list of compositions, where the latest index is sorted first
    Example: dictionary consists of compositions [1,3], [2,4], [1,2], [2,1]
    Output: [[1,2], [1,3], [2,1], [2,4]]
    @param master_dict: master-dictionary extracted from database(s)
    @return: list of lists
    """
    composition_array = np.array([list(key) for key in master_dict.keys()])
    composition_array = composition_array[composition_array[:, -1].argsort()]
    for ind in range(-2, -composition_array.shape[1] - 1, -1):
        composition_array = composition_array[composition_array[:, ind].argsort(kind='mergesort')]
    return composition_array

def merge(a: dict, b: dict, path=[]):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge(a[key], b[key], path + [str(key)])
            elif a[key] != b[key]:
                raise Exception('Conflict at ' + '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

if __name__ == '__main__':
    for o in odometer((2,2,2), (-1,-1,-1)):
        print(o)