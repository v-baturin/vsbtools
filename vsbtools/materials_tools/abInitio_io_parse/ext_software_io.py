import re
import os
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read as ase_read
from ase.io import write as ase_write
from cclib.io import ccread
from genutils.filesystem_tools import add_index
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar
# from .common_tools import recursive_map_to_keys, recursively_map_to_vals, try_numerize_string


def read_poscars(poscars_fname):
    """
    Reads POSCARS multiple images file into a list of Atoms objects
    @param poscars_fname: POSCARS filename
    @return: list of Atoms objects
    """
    poscars_list = list()
    with open(poscars_fname) as poscars_fid:
        while True:
            try:
                poscars_list.append(ase_read(poscars_fid, format='vasp'))
            except IndexError:
                break
    return poscars_list

def read_coords_g(poscars_fname, format='vasp', slc=None):
    """
    Generator,
    Reads POSCARS multiple images file into a list of Atoms objects
    @param poscars_fname: POSCARS filename
    @yield: Atoms object
    """
    poscars_list = list()
    with open(poscars_fname) as poscars_fid:
        while True:
            try:
                yield ase_read(poscars_fid, format=format, index=slc)
            except IndexError:
                break

def poscars2xyz_files(poscars_fname, folder_name=None):
    if folder_name is None:
        folder_name = poscars_fname + '_XYZ'
    os.mkdir(add_index(folder_name))
    for k, ats in enumerate(read_coords_g(poscars_fname)):
        chem_formula = ats.get_chemical_formula()
        filename = add_index(folder_name + '/' + chem_formula + '.xyz')
        ase_write(filename, ats, format='xyz')

def xyzfile2poscars(xyz_fname, box=None, poscars_fname=None):
    if poscars_fname is None:
        poscars_fname = f'{xyz_fname.split(".")[0]}_POSCARS'
    if box is None:
        box = np.array([20., 20., 20])
    ats_all = ase_read(xyz_fname, index=':', format='xyz')
    for k, ats in enumerate(ats_all):
        ats.positions = ats.positions - sum(ats.positions)/len(ats) + box/2
        ats.cell = box
        ats.pbc = [True, True, True]
        ase_write(poscars_fname, ats, format='vasp', append=True, direct=True, vasp5=True)

def dict2gauformat(dct_param: dict, gau_route=False):
    txt = ''
    dct = dct_param.copy()
    if gau_route:
        txt += '#p ' + dct['approach_basis']
        del dct['approach_basis']
    for upper_key, upper_val in dct.items():
        txt += (' ' + upper_key)
        if isinstance(upper_val, dict):
            if len(upper_val) == 1 and next(iter(upper_val.values())) is None:
                txt += ('=' + next(iter(upper_val.keys())))
            elif upper_key == 'iop':
                txt += ('(' + recursive_print(upper_val) + ')')
            else:
                txt += ('=(' + recursive_print(upper_val) + ')')
        elif bool(upper_val):
            txt += ('=' + upper_val)
    return txt.replace(',)', ')')

def recursive_print(dct: dict):
    txt1 = ''
    for k, v in dct.items():
        if isinstance(v, dict):
            if len(v) == 1 and next(iter(v.values())) is None:
                txt1 += (k + '=' + next(iter(v.keys())))
            else:
                txt1 += (k + '=(' + recursive_print(v) + '),')
        elif bool(v):
            txt1 += (k + '=' + str(v) + ',')
        else:
            txt1 += (k + ',')
    return txt1.replace(',)', ')')


def xyz_2atoms(xyz):
    symbols = []
    positions = []
    closefile = False
    if isinstance(xyz, str):
        if '\n' in xyz:
            xyz = xyz.split('\n')
        else:
            closefile = True
            xyz = open(xyz)

    for item in xyz:
        if re.match(r'(^\s*[A-Z][a-z]?\s+(?:[-0-9.]+.*){3})', item):
            tokens = item.split()
            symbol = tokens[0]
            pos = list(map(float, tokens[1:4]))
            symbols.append(symbol)
            positions.append(pos)
    atoms = Atoms(symbols, positions)
    if closefile:
        xyz.close()
    return atoms


def gaussian_in_2_Atoms(fd):
    """
    Borrowed from ASE package
    @param fd:
    @return:
    """
    re_chgmult = re.compile(r'^\s*[+-]?\d+(?:,\s*|\s+)[+-]?\d+\s*$')

    symbols = []
    positions = []
    pbc = np.zeros(3, dtype=bool)
    cell = np.zeros((3, 3))
    npbc = 0
    # We're looking for charge and multiplicity
    for line in fd:
        if re_chgmult.match(line) is not None:
            tokens = fd.readline().split()
            while tokens:
                symbol = tokens[0]
                pos = list(map(float, tokens[1:4]))
                if symbol.upper() == 'TV':
                    pbc[npbc] = True
                    cell[npbc] = pos
                    npbc += 1
                else:
                    symbols.append(symbol)
                    positions.append(pos)
                tokens = fd.readline().split()
            atoms = Atoms(symbols, positions, pbc=pbc, cell=cell)
            return atoms

def atoms_2_str(atms):
    output = []
    for at in atms:
        output.append('{:<10s}{:20.10f}{:20.10f}{:20.10f}'.format(at.symbol, *at.position))
    return '\n'.join(output)

def parse_gout(logfile):
    ccdata = ccread(logfile)
    try:
        last_positions = ccdata.atomcoords[-1]
    except:
        return None, None
    numbers = ccdata.atomnos
    return ccdata, Atoms(positions=last_positions, numbers=numbers)

def cif2poscar(in_cif_fname, out_poscar_fname):
    structure = CifParser(in_cif_fname).get_structures()[0]
    p = Poscar(structure)
    p.write_file(out_poscar_fname)

def cif_dir2poscars_dir(cif_dir, poscars_dir=None, sort_by_comp=False):
    if poscars_dir is None:
        if "cif" in cif_dir.casefold():
            poscars_dir = re.sub('cif', 'POSCAR', cif_dir, flags=re.I)
        else:
            poscars_dir = cif_dir + "_poscars"
    cif_dir = Path(cif_dir)
    for cif_path in cif_dir.glob("*.cif"):
        print(f"Trying to convert {cif_path.name}.")
        try:
            struct = Structure.from_file(cif_path)
            if sort_by_comp:
                comp = struct.composition.get_integer_formula_and_factor()[0]
                target_dir = Path(poscars_dir) / comp
            else:
                target_dir = Path(poscars_dir)
            target_dir.mkdir(parents=True, exist_ok=True)
            poscar_path = target_dir / f'{cif_path.name.replace(".cif","_POSCAR")}'
            struct.to(fmt="poscar", filename=str(poscar_path))
        except ValueError:
            print(f"Bad structure detected. Check {cif_path.name}!")

# def get_options(string, option):
#     """
#     :param string: 'option=val option1=(val1, val2) option3(val3, val4)
#     :param option: one of the options
#     :return: list of values for given option
#     """
#     # Preconditioning the string:
#     string = re.sub(' +', ' ', string)
#     string = re.sub('\s*=\s*', '=', string)
#     string = re.sub('\s*\)', ')', string)
#     string = re.sub('\(\s*', '(', string)
#     string = re.sub('\s*,\s*', ',', string)
#
#     # Analyzing the string
#     if option.casefold() not in string.casefold():
#         return None
#     option_opts_findres = re.findall('(?:' + option + '[=\(]+)([^ \)\n]+)', string, re.IGNORECASE)
#
#     if len(option_opts_findres) == 0:
#         return []
#     else:
#         option_opts = ','.join(option_opts_findres).split(',')
#     return [j.strip() for j in option_opts]

# if __name__ == '__main__':
#     upper = lambda x: x.casefold()
#     testdct = {'Mem': 10, 'HaHa': {'Lalal': 'AAA', 'uUuU': 777}}
#     resdct = recursive_map_to_keys(upper, testdct)
#     print(resdct)
#
#
#     def try_numerize_string(string):
#         if not isinstance(string, str):
#             return string
#         else:
#             string = string.strip()
#         if re.match('-?[0-9]+$', string):
#             return int(string)
#         else:
#             try:
#                 return float(string)
#             except:
#                 return string
#
#
#     testdct2 = {'a': 'Hello', 'b': '150.00', 'c': 10, 'd': '10', 'e': '1e5'}
#     resdct = recursively_map_to_vals(try_numerize_string, testdct2)
#     print(resdct)
#
#     g_in_fname = '/home/vsbat/SYNC/1_TALKS/20200604_Sk_lab_TALK/TXM-CAT_CHCl3_Step1.gjf'
#     with open(g_in_fname) as f:
#         ats = gaussian_in_2_Atoms(f)
#     xyz_at = xyz_2atoms(g_in_fname)
#
#     bad_g_out = '/home/vsbat/SYNC/00_Current_PyWork/refine_with_gaussian/testfolder/log737771.8835'
#     # output_g = '/home/vsbat/SYNC/00_Current_PyWork/refine_with_gaussian/testfolder/log737771.7237'
#     ccdata = ccread(bad_g_out)
#
#
#     pass
#     # poscar_fname = '/home/vsbat/SYNC/00_Current_PyWork/20200319_CuAu_project/CuAu/CuAu_varcomp_6_10_6_10/results2/POSCAR'
#     # poscar_atoms = ase_read(poscar_fname)

if __name__ == "__main__":
    # xyz_batch_fname = '/home/vsbat/SYNC/00__WORK/20201013_PdBi/restartPd1015/xyz_all.xyz'
    # poscars_fname = '/home/vsbat/SYNC/00__WORK/20201013_PdBi/restartPd1015/POSCARS_1'
    # xyzfile2poscars(xyz_batch_fname, np.array([20, 20, 20]), poscars_fname)
    cif_dir_mopb = "/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/STRUCTURE_GENERATION/gen_results/20250403_MoPB_SG140/CIF_Mo2B1P1_140"
    cif_dir2poscars_dir(cif_dir_mopb)

