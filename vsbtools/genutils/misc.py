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