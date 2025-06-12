import os
import re
import shutil
from pathlib import Path
from typing import Dict, Union, Callable
from my_packages.genutils.filesystem_tools import sh_execute
import numpy as np
from USPEX.Atomistic.RadialDistributionUtility import RadialDistributionUtility
from USPEX.components import AtomisticRepresentation, Atomistic
from USPEX.DataModel.Engine import Engine
from USPEX.DataModel.Flavour import Flavour
from USPEX.DataModel.Entry import Entry

atomistic = Atomistic()

data_vs_struct_files = {'Individuals': 'gatheredPOSCARS.uspex',
                        'BESTIndividuals': 'BESTgatheredPOSCARS.uspex',
                        'goodStructures': 'goodStructures_POSCARS.uspex'}


def read_Individuals_uspexPY(fname, out_dict=None):
    with open(fname, 'r') as ind_fid:
        s = '+'
        while s[0] != '|':
            s = ind_fid.readline()
        headers = [h.strip() for h in s.split('|')[1:-1]]
        if out_dict is None:
            out_dict = {h: [] for h in headers}
            out_dict['structures'] = []
        for ln in ind_fid:
            if ln[0] != '|':
                continue
            for i, val in enumerate([h.strip() for h in ln.split('|')[1:-1]]):
                if val.strip() != headers[i]:
                    if headers[i] == 'ID':
                        v = int(val)
                    else:
                        try:
                            v = float(val)
                        except ValueError:
                            v = val
                    out_dict[headers[i]].append(v)
    return out_dict


def readAtomicStructuresToPoolEntries(struc_file_path):
    Engine.createEngine(':memory:')
    extensions = dict(atomistic = (atomistic, atomistic.propertyExtension.propertyTable))  #.propertyExtension())
    systems = []
    for i, s in enumerate(Atomistic.readAtomicStructures(struc_file_path)):
        print(f"reading {i}")
        systems.append(Entry.newEntry(Flavour(extensions=extensions,
                                                       **{'.howCome': 'Seeds', '.parent': None, '.label': f"EA{i}"},
                                                       **s)))
    return systems

def atomsListToPoolEntries(atoms_list: list) -> list:
    Engine.createEngine(':memory:')
    extensions = dict(atomistic=(atomistic, atomistic.propertyExtension.propertyTable))  # .propertyExtension())
    systems = []
    for i, s in enumerate(Atomistic.atomsListToSystems(atoms_list)):
        systems.append(Entry.newEntry(Flavour(extensions=extensions,
                                              **{'.howCome': 'Seeds', '.parent': None, '.label': f"EA{i}"},
                                              **s)))
    return systems

def readResFolders_uspexPY(path, individuals_kind: str = 'Individuals'):
    out_dict = None
    for i, ind_file in enumerate(path.rglob(individuals_kind)):
            out_dict = read_Individuals_uspexPY(ind_file, out_dict)
            out_dict['structures'] += readAtomicStructuresToPoolEntries(ind_file.parent/data_vs_struct_files[individuals_kind])
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
            headers[i] += second_line[(start_end_idcs[i][0] - 1):(start_end_idcs[i + 1][0] - 1)].strip()
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


# def readResFolders_uspexML(path, individuals_fname_pattern: str = 'composition_*', max_entries=200, iscluster=False):
#     out_dict = None
#     for i, ind_file in enumerate(path.rglob(individuals_fname_pattern)):
#         if 'POSCARS' not in ind_file.name:
#             out_dict = read_Individuals_uspexML(ind_file, out_dict, max_entries=max_entries)
#             for k_struct, s in enumerate(
#                     AtomisticRepresentation.readAtomicStructures(ind_file.parent / (ind_file.name + '_POSCARS'))):
#                 if k_struct < max_entries:
#                     if iscluster:
#                         s['cell'] = Cell.initFromCellVectors((0, 0, 0))
#                     out_dict['structures'].append(AtomisticPoolEntry(**s))
#     return out_dict


def create_uspexfold(dest_folder,
                     uspex_common_source=None,
                     template_path=None,
                     seeds_pool_entries=None,
                     components_to_copy=None,
                     callbacks: Dict[Callable, Dict] = None,  # {f1 : {args: [], kwargs = {kw1:v1,...}}, f2 : ...}
                     open_input=True,
                     **kwargs):
    if template_path is None:
        os.makedirs(dest_folder, exist_ok=True)
    else:
        shutil.copytree(template_path, dest_folder, dirs_exist_ok=True)

    if uspex_common_source is not None:
        sh_execute(f'ln -s {uspex_common_source} {dest_folder}/Common ')

    if components_to_copy is not None:
        for component, path in components_to_copy.items():
            if os.path.isfile(path):
                shutil.copy(path, dest_folder / component)
            elif os.path.isdir(path):
                shutil.copytree(path, dest_folder / component, dirs_exist_ok=True)

    if callbacks:
        for f, argdict in callbacks.items():
            f(*argdict['args'], **argdict['kwargs'])

    if seeds_pool_entries is not None:
        seeds_path = dest_folder / 'Seeds/1'
        os.makedirs(seeds_path, exist_ok=True)
        seeds_flavours = [pe.getFlavour('origin') for pe in seeds_pool_entries]
        for i, ef in enumerate(seeds_flavours):
            ef.setProperty('label', f"EA{i}")
        atomistic.writeAtomicStructures(seeds_path / 'POSCARSeeds', seeds_flavours)

    if open_input:
        sh_execute(f'gedit {dest_folder}/input.uspex')


def calc_time_from_pool_entries(pool_entries, std_time, time_limit):
    maxtime = std_time
    for seeds_pool_entry in pool_entries:
        natoms = len(seeds_pool_entry.getAtomicStructure())
        newtime = std_time + 0.01 * natoms ** 3
        maxtime = max(newtime, maxtime)
        maxtime = min(maxtime, time_limit)
    return maxtime


def prepare_potcars(specific_path, element_symbols, potcars_source, potcars_spec: Union[None, Dict[str, str]] = None):
    if potcars_spec is None:
        potcars_spec = {k: k for k in element_symbols}
    for symbol in element_symbols:
        shutil.copy(potcars_source / potcars_spec[symbol] / 'POTCAR', specific_path / f"POTCAR_{symbol}")


def input_corrector(input_source_path, corrector_dict, input_result_path=None, custom_patterns=None):
    """
    Modifyer of input.uspex
    @param input_source_path: original input path
    @param corrector_dict: dict of the following format: {keyword: newvalue}
                           Example: {'blocks': [2, 1, 1]}
    @param input_result_path: modified input path
    @param custom_patterns: regex pattern for replacement input parameters
                            Example: r'heredity: \( *((?:[\d\.]+ *)+) *\)'
                            MIND that only the captured group is replaced
    @return: None
    """

    if input_result_path is None:
        input_result_path = input_source_path

    infile_patterns = {'symbols': r'symbols:\s+\[ *((?:\w+ *)+) *\]',
                       'blocks': r'blocks:\s+\[\[((?:\d+\s*)+)]\]',
                       'range': r'range: \[\(((?:\d+\s*)+)\)\]',
                       'ionDistances': r'ionDistances *: *\{((?:\'\D+\' *: *[\d\.]+ *)+)\}'}

    if custom_patterns is not None:
        infile_patterns.update(custom_patterns)

    with open(input_source_path) as in_fid:
        replaced_text = in_fid.read()

    for k, v in corrector_dict.items():
        m = re.search(infile_patterns[k], replaced_text)
        if m:
            newdatastr = ' '.join([str(i) for i in v])
            replacement = re.sub(m.group(1), newdatastr, m.group(0))
            replaced_text = replaced_text.replace(m.group(0), replacement)
        else:
            print(f'Pattern for {k} not found')

    with open(input_result_path, 'w') as in_fid:
        in_fid.write(replaced_text)


if __name__ == '__main__':
    # file = '/20230414_Borohydrures/case_studies/fingerprint_choice/resCa/results2/goodStructures'
    # output = read_Individuals_uspexPY(file)
    # print('.')
    # file = '/home/vsbat/Downloads/Telegram_downloads/composition_5_5'
    # outp = read_Individuals_uspexML(file)
    # path = Path('/home/vsbat/SYNC/00__WORK/20220324_LiP_PROJECT/02_LiP_clusters/filter_by_FP/5_5')
    # outp = readResFolders_uspexML(path)
    # print('x')
    input_example = Path(
        '/home/vsbat/SYNC/00__WORK/20231018_ALANATES/03_USPEX_calc_management/01_USPEX_debug/templates/template_1gen_seeds_local/input.uspex')
    input_res = Path(
        '/home/vsbat/SYNC/00__WORK/20231018_ALANATES/03_USPEX_calc_management/01_USPEX_debug/templates/template_1gen_seeds_local/input.uspex_m')
    input_corrector(input_example, {'symbols': 'Li molmol'}, input_res)
