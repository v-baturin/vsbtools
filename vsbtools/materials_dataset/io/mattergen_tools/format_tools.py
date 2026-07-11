import re

def build_pairs_from_guidance(s: str) -> dict:
    """
    Parse the 'guidance' block of a config-like string and return:
    {"A-B neigh": {"type_A": "A", "type_B": "B"}, ...}
    Returns {} if no A-B keys are found.
    """
    # Locate start of 'guidance' and the brace-enclosed block
    m = re.search(r'guidance\s*:\s*', s)
    if not m:
        return {}
    i = s.find('{', m.end())
    if i == -1:
        return {}
    depth, end = 0, -1
    for j, ch in enumerate(s[i:], start=i):
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    if end == -1:
        return {}

    block = s[i:end]

    # Find all "A-B" keys inside the guidance block (quotes optional, single/double/doubled)
    pairs = re.findall(r"""['"]{0,2}([A-Za-z][A-Za-z0-9]*-[A-Za-z][A-Za-z0-9]*)['"]{0,2}\s*:""", block)

    # Build result (deduplicate while preserving order)
    out, seen = {}, set()
    for cs in pairs:
        if cs in seen:
            continue
        seen.add(cs)
        a, b = cs.split('-', 1)
        out[f"{cs} neigh"] = {'type_A': a, 'type_B': b}

    return out