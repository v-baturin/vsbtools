import os
import re
import shutil
from pathlib import Path
from genutils.filesystem_tools import sh_execute
import numpy as np
from USPEX.Atomistic.RadialDistributionUtility import RadialDistributionUtility
from USPEX.components import AtomisticRepresentation, AtomisticPoolEntry, Cell


data_vs_struct_files = {'Individuals': 'gatheredPOSCARS.uspex',
                        'BESTIndividuals': 'BESTgatheredPOSCARS.uspex',
                        'goodStructures': 'goodStructures_POSCARS.uspex'}
def read_Individuals_uspexPY(fname, out_dict=None):
    with open(fname, 'r') as ind_fid:
        ind_fid.readline()
        headers = [h.strip() for h in ind_fid.readline().split('|')[1:-1]]
        if out_dict is None:
            out_dict = {h: [] for h in headers}
            out_dict['structures'] = []
        for ln in ind_fid:
            if ln[0] == '+':
                continue
            for i, val in enumerate([h.strip() for h in ln.split('|')[1:-1]]):
                v = float(val) if headers[i] == 'Enthalpy (eV)' else val
                out_dict[headers[i]].append(v)
    return out_dict

def readResFolders_uspexPY(path, individuals_kind: str = 'Individuals'):
    out_dict = None
    for i, ind_file in enumerate(path.rglob(individuals_kind)):
            out_dict = read_Individuals_uspexPY(ind_file, out_dict)
            out_dict['structures'] += [AtomisticPoolEntry(**s) for s in AtomisticRepresentation.readAtomicStructures(ind_file.parent/data_vs_struct_files[individuals_kind])]
    return out_dict

def read_Individuals_uspexML(fname, out_dict=None, max_entries=200):
    def indiv_line_parser(input_string):
        elements = re.findall(r'\[[^\]]+\]|\S+', input_string)
        elements = [element.strip('[]').strip() if element.startswith('[') else element for element in elements]
        return elements
    with open(fname, 'r') as ind_fid:
        headers_line = ind_fid.readline()
        headers = [h.strip() for h in headers_line.split()]
        start_end_idcs = []
        for h in headers:
            for m in re.finditer(h, headers_line):
                span = [m.start(), m.end()]
                if span not in start_end_idcs:
                    start_end_idcs.append(span)
        start_end_idcs.append([len(headers_line), len(headers_line)])
        start_end_idcs.sort()
        second_line = ind_fid.readline() + ' ' * len(headers_line)
        for i in range(len(headers)):
            headers[i] += second_line[(start_end_idcs[i][0] - 1):(start_end_idcs[i+1][0] - 1)].strip()
        if out_dict is None:
            out_dict = {h: [] for h in headers}
            out_dict['structures'] = []
            out_dict['k_infile'] = []
        for k_ind, ln in enumerate(ind_fid):
            if k_ind < max_entries:
                out_dict['k_infile'].append(k_ind)
                for i, val in enumerate(indiv_line_parser(ln)):
                    v = float(val) if headers[i] == 'Enthalpies(eV)' else val
                    out_dict[headers[i]].append(v)

    return out_dict

def readResFolders_uspexML(path, individuals_fname_pattern: str = 'composition_*', max_entries=200, iscluster=False):
    out_dict = None
    for i, ind_file in enumerate(path.rglob(individuals_fname_pattern)):
        if 'POSCARS' not in ind_file.name:
            out_dict = read_Individuals_uspexML(ind_file, out_dict, max_entries=max_entries)
            for k_struct, s in enumerate(AtomisticRepresentation.readAtomicStructures(ind_file.parent/(ind_file.name + '_POSCARS'))):
                if k_struct < max_entries:
                    if iscluster:
                        s['cell'] = Cell.initFromCellVectors((0,0,0))
                    out_dict['structures'].append(AtomisticPoolEntry(**s))
    return out_dict

def create_uspexfold(dest_folder,
                     uspex_common_source,
                     template_path=None,
                     seeds_pool_entries=None,
                     varied_calc_components=None,
                     open_input=True,
                     **kwargs):
    if template_path is None:
        os.makedirs(dest_folder, exist_ok=True)
    else:
        shutil.copytree(template_path, dest_folder, dirs_exist_ok=True)
    sh_execute(f'ln -s {uspex_common_source} {dest_folder}/Common ')

    if varied_calc_components is not None:
        for component, path in varied_calc_components.items():
            if os.path.isfile(path):
                shutil.copy(path, dest_folder / component)
            elif os.path.isdir(path):
                shutil.copytree(path, dest_folder / component, dirs_exist_ok=True)

    input_path = dest_folder / 'input.uspex'

    if seeds_pool_entries is not None:
        seeds_path = dest_folder / 'Seeds/1'
        os.makedirs(seeds_path, exist_ok=True)
        AtomisticRepresentation.writeAtomicStructures(seeds_path / 'POSCARSeeds', seeds_pool_entries)

    if open_input:
        sh_execute(f'gedit {input_path}')

def calc_time_from_pool_entries(pool_entries, std_time, time_limit):
    for seeds_pool_entry in pool_entries:
        natoms = len(seeds_pool_entry.getAtomicStructure())
        newtime = std_time + 0.01 * natoms ** 3
        maxtime = max(newtime, maxtime)
        maxtime = min(maxtime, time_limit)
    return maxtime

if __name__ == '__main__':
    # file = '/20230414_Borohydrures/case_studies/fingerprint_choice/resCa/results2/goodStructures'
    # output = read_Individuals_uspexPY(file)
    # print('.')
    # file = '/home/vsbat/Downloads/Telegram_downloads/composition_5_5'
    # outp = read_Individuals_uspexML(file)
    path = Path('/home/vsbat/SYNC/00__WORK/20220324_LiP_PROJECT/02_LiP_clusters/filter_by_FP/5_5')
    outp = readResFolders_uspexML(path)
    print('x')
