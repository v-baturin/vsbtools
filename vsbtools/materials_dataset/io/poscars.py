from __future__ import annotations

from ase.io import read as ase_read


def read_poscars(poscars_fname):
    """Read a concatenated POSCARS file into ASE Atoms objects."""
    poscars_list = []
    with open(poscars_fname) as poscars_fid:
        while True:
            try:
                poscars_list.append(ase_read(poscars_fid, format="vasp"))
            except (IndexError, StopIteration):
                break
    return poscars_list
