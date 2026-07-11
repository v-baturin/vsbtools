from typing import Any
import re
import numbers
import numpy as np
import textwrap, inspect, ast

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

def describe_predicate(fn):
    """Return a clean string like  'lambda x: x % 2 == 0'."""
    try:
        src = textwrap.dedent(inspect.getsource(fn))
        tree = ast.parse(src)

        # find the first (and usually only) Lambda node
        for node in ast.walk(tree):
            if isinstance(node, ast.Lambda):
                # Python ≥3.8: returns the exact slice of source for *node*
                expr = ast.get_source_segment(src, node)
                return expr.strip()
        # fallback: maybe it was a normal def
        first_line = src.strip().splitlines()[0]
        return first_line
    except (OSError, TypeError, SyntaxError):
        # interactive / C-level / no source
        qname = getattr(fn, "__qualname__", None)
        return qname if qname and "<lambda>" not in qname else repr(fn)

_MULTI_US = re.compile(r"_+")

def _fmt_num(x: Any) -> str:
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)

def _sanitize_token(s: str) -> str:
    # Turn any run of braces/brackets/parentheses/quotes/spaces into "_"
    s = re.sub(r"[{}\[\]\(\)\"'\s]+", "_", s.strip())
    s = _MULTI_US.sub("_", s)
    return s

def serialize_structure(obj: Any) -> str:
    """
    Dicts: key and value concatenated with '_'
    Lists/Tuples: elements joined with '-'
    Scalars: sanitized; hyphens and dots are preserved
    Leading/trailing '_' are stripped at the end.
    """
    def ser(x: Any) -> str:
        if isinstance(x, dict):
            parts: list[str] = []
            for k, v in x.items():  # preserves insertion order
                k_str = _sanitize_token(str(k))
                v_str = ser(v)
                parts.append(k_str if not v_str else f"{k_str}_{v_str}")
            return "_".join(p for p in parts if p)
        if isinstance(x, (list, tuple)):
            items = [ser(e) for e in x]
            items = [i for i in items if i]  # drop empties
            return "-".join(items)
        if x is None:
            return "None"
        if isinstance(x, bool):
            return "True" if x else "False"
        if isinstance(x, (int, float)):
            return _fmt_num(x)
        return _sanitize_token(str(x))

    s = ser(obj)
    s = _MULTI_US.sub("_", s).strip("_")
    return s

def is_substructure(tree, pattern, *, epsilon=1e-8):
    """
    Return True if `pattern` is a substructure of `tree`.

    - Dict pattern:
        matches some dict node in `tree` that has at least those keys.
        For each key k, pattern[k] is a substructure of node[k].
    - List pattern:
        matches a contiguous sublist of some list node in `tree`.
    - Leaf pattern:
        if both leaves are numeric (non-bool reals), they match if
        abs(a - b) < epsilon; otherwise they match iff a == b.
    """

    def is_real_number(x):
        # Real number, but not bool
        return isinstance(x, numbers.Real) and not isinstance(x, bool)

    def leaves_equal(a, b):
        # Bool vs non-bool: require exact equality
        if isinstance(a, bool) or isinstance(b, bool):
            return a is b

        # Numeric comparison with tolerance
        if is_real_number(a) and is_real_number(b):
            return abs(a - b) < epsilon

        # Fallback: exact equality
        return a == b

    def match_here(node, pat):
        """Check if `pat` matches the subtree rooted exactly at `node`."""

        # dict pattern: subset-of-keys constraint
        if isinstance(pat, dict):
            if not isinstance(node, dict):
                return False
            for k, v in pat.items():
                if k not in node:
                    return False
                if not is_substructure(node[k], v, epsilon=epsilon):
                    return False
            return True

        # list pattern: must match a contiguous sublist
        if isinstance(pat, list):
            if not isinstance(node, list):
                return False
            if not pat:  # empty list pattern matches trivially
                return True
            n, m = len(node), len(pat)
            if m > n:
                return False
            for start in range(n - m + 1):
                if all(is_substructure(node[start + i], pat[i], epsilon=epsilon)
                       for i in range(m)):
                    return True
            return False

        # leaf pattern: tolerant numeric comparison, otherwise ==
        return leaves_equal(node, pat)

    # First, try to match with the root aligned at `tree`
    if match_here(tree, pattern):
        return True

    # Otherwise, search deeper in children
    if isinstance(tree, dict):
        return any(is_substructure(v, pattern, epsilon=epsilon)
                   for v in tree.values())
    if isinstance(tree, list):
        return any(is_substructure(x, pattern, epsilon=epsilon)
                   for x in tree)
    return False  # leaf: nowhere deeper to look

# If you still want the original name:
is_subtree = is_substructure



if __name__ == '__main__':
    for o in odometer((2,2,2), (-1,-1,-1)):
        print(o)