import os
from gather_best import parse_001_results
from tools_stability.aux_routines import list_fmt2table
from genutils.misc import get_sorted_compositions
from genutils.filesystem_tools import add_index
from ase.io import write

def write_best_n_isom(all_gathered_data, format='xyz', n_best=1, out_dir='gathered_isoms'):
    os.makedirs(out_dir, exist_ok=True)
    if isinstance(all_gathered_data, str):
        all_gathered_data = parse_001_results(all_gathered_data)
    if format == 'xyz':
        for k, v in all_gathered_data.items():
            for k in range(max(len(v['atoms']), n_best)):
                xyz_fname = out_dir + v['atoms'][k].get_chemical_formula(mode='reduce') + f'_{k}.xyz'
                write(xyz_fname, v['atoms'][k], format='xyz')
    elif format.casefold() in ['poscar', 'vasp']:
        sorted_cmp_array = get_sorted_compositions(all_gathered_data)
        poscar_fname = add_index(out_dir + '/POSCARS')
        for cmp in sorted_cmp_array:
            cmp = tuple(cmp)
            for k in range(max(len(v['atoms']), n_best)):
                write(poscar_fname, v['atoms'][k], format='vasp', vasp5=True, append=True)

def get_energy_table(all_gathered_data, outfile='en_table.txt'):
    if isinstance(all_gathered_data, str):
        all_gathered_data = parse_001_results(all_gathered_data)
    energytable_formatted = []
    for k, v in all_struct_ordered_dict.items():
        energytable_formatted.append(list(k) + [v['Enthalpies'][0]])
    en_matrix_nopure = list_fmt2table(energytable_formatted, outfile=out_dir + 'en_table.txt', placeholder=0.0)
    return en_matrix_nopure
