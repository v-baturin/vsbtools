from USPEX.Atomistic.RadialDistributionUtility import RadialDistributionUtility
from USPEX.components import AtomisticRepresentation, AtomisticPoolEntry
from pathlib import Path
import numpy as np

data_vs_struct_files = {'Individuals': 'gatheredPOSCARS.uspex',
                        'BESTIndividuals': 'BESTgatheredPOSCARS.uspex',
                        'goodStructures': 'goodStructures_POSCARS.uspex'}
def read_Individuals(fname, out_dict=None):
    with open(fname, 'r') as ind_fid:
        ind_fid.readline()
        headers = [h.strip() for h in ind_fid.readline().split('|')[1:-1]]
        if out_dict is None:
            out_dct = {h: [] for h in headers}
            out_dct['structures'] = []
        for ln in ind_fid:
            if ln[0] == '+':
                continue
            for i, val in enumerate([h.strip() for h in ln.split('|')[1:-1]]):
                v = float(val) if headers[i] == 'Enthalpy (eV)' else val
                out_dct[headers[i]].append(v)
    return out_dct

def readResFolders(path, individuals_kind: str = 'Individuals'):
    outdict = None
    for i, ind_file in enumerate(path.rglob(individuals_kind)):
            out_dict = read_Individuals(ind_file, outdict)
            out_dict['structures'] += AtomisticRepresentation.readAtomicStructures(ind_file.parent/data_vs_struct_files[individuals_kind])
    return out_dict




if __name__ == '__main__':
    file = '/20230414_Borohydrures/case_studies/fingerprint_choice/resCa/results2/goodStructures'
    output = read_Individuals(file)
    print('.')

