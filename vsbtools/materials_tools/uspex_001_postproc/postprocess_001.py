import os
import numpy as np
from .gather_best import parse_001_results
from ..tools_stability.aux_routines import list_fmt2table
from genutils.misc import get_sorted_compositions
from genutils.filesystem_tools import add_index
from ase.io import write

def write_best_n_isom(all_gathered_data, format='xyz', n_best=1, out_dir='gathered_isoms'):
    os.makedirs(out_dir, exist_ok=True)
    if isinstance(all_gathered_data, str):
        all_gathered_data = parse_001_results(all_gathered_data)
    if format == 'xyz':
        for k, v in all_gathered_data.items():
            for i in range(min(len(v['atoms']), n_best)):
                xyz_fname = out_dir + v['atoms'][i].get_chemical_formula(mode='reduce') + f'_{i}.xyz'
                write(xyz_fname, v['atoms'][i], format='xyz')
    elif format.casefold() in ['poscar', 'vasp']:
        sorted_cmp_array = get_sorted_compositions(all_gathered_data)
        poscar_fname = add_index(out_dir + '/POSCARS')
        for cmp in sorted_cmp_array:
            cmp = tuple(cmp)
            v = all_gathered_data[cmp]
            for i in range(min(len(v['atoms']), n_best)):
                write(poscar_fname, v['atoms'][i], format='vasp', vasp5=True, append=True)


def get_energy_matrix(all_gathered_data, outfile='en_matrix.txt', **kwargs):
    if isinstance(all_gathered_data, str):
        all_gathered_data = parse_001_results(all_gathered_data)
    energytable_formatted = []
    for k, v in all_gathered_data.items():
        energytable_formatted.append(list(k) + [v['Enthalpies'][0]])
    en_matrix_nopure = list_fmt2table(energytable_formatted, outfile=outfile, **kwargs)
    return en_matrix_nopure

def get_energy_table(all_gathered_data):
    if isinstance(all_gathered_data, str):
        all_gathered_data = parse_001_results(all_gathered_data)
    energytable_sorted = []
    sorted_cmp_array = get_sorted_compositions(all_gathered_data)
    for cmp in sorted_cmp_array:
        cmp = tuple(cmp)
        energytable_sorted.append(list(cmp) + [all_gathered_data[cmp]['Enthalpies'][0]])
    return np.array(energytable_sorted)
