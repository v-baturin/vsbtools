from pathlib import Path
from ase.io import read, write
import re
import numpy as np
from ab_initio_postprocessing.abInitio_io_parse.ext_software_io import read_coords_g

# res_folder = Path("/home/vsbat/SYNC/00__WORK/PtAu_PROJECT/gathered_results")
#
def parse_001_results(res_folder):
    all_data_dict = {}
    for fld in Path(res_folder).rglob("__goodStructures"):
        for indivs in fld.rglob('composition_*'):
            if 'POSCARS' in str(indivs):
                continue
            poscar_file = str(indivs) + '_POSCARS'
            all_data_dict = {}
            all_data_dict = add_results_data(all_data_dict, indivs, poscar_file)
    for k, v in all_data_dict.items():
        indices = np.array(v['Enthalpies']).argsort()
        for element, lst in v.items():
            all_data_dict[k][element] = [lst[i] for i in indices]

    return all_data_dict

# composition_fpath = "/home/vsbat/SYNC/00__WORK/PtAu_PROJECT/gathered_results/Pt1-4_Au1-4/results2/__goodStructures/composition_1_1"

def add_results_data(all_data_dict, indivs_fpath, poscars_fpath):
    if all_data_dict is None:
        all_data_dict = {}  # all_data_dict[(3, 8)] - data on X3Y8 compound
    with open(indivs_fpath, 'r') as comp_fid:
        header = comp_fid.readline().split()
        header = [x for x in header if x != 'Compositions']
        comp_fid.readline()
        poscars_gen = read_coords_g(poscars_fpath)
        for line in comp_fid:
            individual = {}
            detect_comp = re.findall('(.*)\[(.*)\](.*)', line)[0]
            assert len(detect_comp) == 3
            composition = tuple([int (i) for i in detect_comp[1].split()])
            formatted_info = detect_comp[0].split() + detect_comp[2].split()
            if composition not in all_data_dict:
                all_data_dict[composition] = {x:[] for x in header}
                all_data_dict[composition]['atoms'] = []
            enthalpies = np.inf
            for i,element in enumerate(header):
                if element in ['ID', 'SYMM']:
                    all_data_dict[composition][element].append(int(formatted_info[i]))
                elif element == 'Enthalpies':
                    enthalpies = np.min((enthalpies, float(formatted_info[i])))
                else:
                    all_data_dict[composition][element].append(float(formatted_info[i]))
            all_data_dict[composition]['Enthalpies'].append(enthalpies)
            all_data_dict[composition]['atoms'].append(next(poscars_gen))
    return all_data_dict


if __name__ == '__main__':
    #  test1
    # composition_fpath = "/home/vsbat/SYNC/00__WORK/PtAu_PROJECT/gathered_results/Pt1-4_Au1-4/results2/__goodStructures/composition_1_1"
    # all_data_dict = {}
    # all_data_dict = add_results_data(all_data_dict, composition_fpath, composition_fpath + '_POSCARS')
    # print(all_data_dict)
    #  test2

    res_dir_test = "/home/vsbat/SYNC/00__WORK/PtAu_PROJECT/gathered_results"
    all_data_dict = parse_001_results(res_dir_test)
    # print(all_data_dict)










