import sys
import numpy as np
import re
from .. import my_io
from cclib.method import MPA
from cclib.io import ccread
from cclib.parser.utils import PeriodicTable as pt

ptable = pt()
element_labels = np.array(ptable.element[:])


def get_kohn_sham(clusterdata):
    if isinstance(clusterdata, str):
        gaudata = ccread(clusterdata)
    else:
        gaudata = clusterdata

    homo_idx = gaudata.homos[0]
    up_states = gaudata.moenergies[0]
    dn_states = gaudata.moenergies[-1]
    fermi = 0.9 * up_states[homo_idx] + 0.1 * up_states[homo_idx + 1]

    return up_states, dn_states, fermi


def get_formula(clusterdata, extended_out=False):
    """

    @param clusterdata: name of out file or a ccData object containing cclib data attributes
    @param extended_out: if True returns not only formula but lists of atomic symbols and number of atoms in a dict
    @return:
    """
    if isinstance(clusterdata, str):
        gaudata = ccread(clusterdata)
    else:
        gaudata = clusterdata

    symbols_list = list(element_labels[gaudata.atomnos])

    # Now get different symbols and numbers of atoms with preserving order
    diff_symbols = []
    num_atoms = dict()
    for symbol in symbols_list:
        if symbol not in diff_symbols:
            diff_symbols.extend([symbol])
            num_atoms[symbol] = 1
        else:
            num_atoms[symbol] += 1

    formula = ''
    for symbol in diff_symbols:
        formula += symbol + str(num_atoms[symbol])
    if extended_out == True:
        return {'formula': formula, 'symbols': diff_symbols, 'fdict': num_atoms}
    else:
        return formula


def get_orbital_coefficients(clusterdata):
    if isinstance(clusterdata, str):
        gaudata = ccread(clusterdata)
    else:
        gaudata = clusterdata
    coeffs = gaudata.mocoeffs[0]
    # # Returns the mtrx coeffs(2X2) of expansion coefficients over atomic orbitals
    # # coeffs[2] is the array of expansion coefficients for the 3-rd MO
    # ext_formula = get_formula(log_filename, extended_out=True)
    # symbols = ext_formula['symbols']
    # symbols_set = set(symbols)
    # coeff_block = my_parse.getblock(log_filename, 'Molecular Orbital Coefficients:', 'Density Matrix:')
    # try:
    #     n_orbitals = int(coeff_block[-1].split()[0])
    # except:
    #     print('i am here')
    # del coeff_block[0::n_orbitals + 3]  # Deletes all strings with numbers of molecular orbitals like '1   2   3   4   5' and so on
    # del coeff_block[0::n_orbitals + 2]  # Deletes all occupancy flags like 'O   O   O   V   V' and so on
    # del coeff_block[0::n_orbitals + 1]  # Deletes all 'Eigenvalues --' strings
    #
    # i_atom = 0
    # atom_0_label = set(coeff_block[0].split()) & symbols_set
    # if len(atom_0_label) == 1:
    #     raw_coeffs = ' '.join(coeff_block[0::n_orbitals])
    #     raw_coeffs_split = raw_coeffs.split()
    #     for x in range(9, 5, -1):
    #         del raw_coeffs_split[::x]
    #     coeffs = [[float(x) for x in raw_coeffs_split]]
    # else:
    #     raise IOError("Unexpected file format")
    #
    # ranges = [{'label': atom_0_label, 'range': [0]}, ]
    #
    # for line_no in range(1, n_orbitals):
    #     curr_line_split = coeff_block[line_no].split()
    #     curr_line_split_set = set(curr_line_split)
    #     new_atom_label = curr_line_split_set & symbols_set
    #     # print(line_no)
    #     if len(new_atom_label) == 1:
    #         i_atom += 1
    #         ranges[-1]['range'].extend([line_no])
    #         ranges.extend([{'label': new_atom_label, 'range': [line_no]}])
    #         raw_coeffs = ' '.join(coeff_block[line_no::n_orbitals])
    #         raw_coeffs_split = raw_coeffs.split()
    #         for x in range(9, 5, -1):
    #             del raw_coeffs_split[::x]
    #         coeffs.append([float(x) for x in raw_coeffs_split])
    #     else:
    #         raw_coeffs = ' '.join(coeff_block[line_no::n_orbitals])
    #         raw_coeffs_split = raw_coeffs.split()
    #         # print(raw_coeffs_split[0:20])
    #         for x in range(len(curr_line_split), 5, -1):
    #             del raw_coeffs_split[::x]
    #         coeffs.append([float(x) for x in raw_coeffs_split])
    # ranges[-1]['range'].extend([line_no + 1])
    #
    # coeffs = np.array(coeffs).transpose()

    return coeffs


def get_atombasis(clusterdata):
    # [{'label': {'Cd'}, 'range': [0, 24]},
    #  {'label': {'Cd'}, 'range': [24, 48]},...]
    rgs = []
    if isinstance(clusterdata, str):
        gaudata = ccread(clusterdata)
    else:
        gaudata = clusterdata
    # basisbyatom = np.array(gaudata.atombasis)
    basisbyatom = gaudata.atombasis
    symbols_list = list(element_labels[gaudata.atomnos])
    for idx, rangeslist in enumerate(basisbyatom):
        rgs.append({'label': {symbols_list[idx]}, 'range': [rangeslist[0], rangeslist[-1] + 1]})
    return rgs


def get_mullik_contributions(clusterdata):
    if isinstance(clusterdata, str):
        gaudata = ccread(clusterdata)
    else:
        gaudata = clusterdata
    m = MPA(gaudata)
    m.calculate()
    return np.array(m.aoresults)


def get_pdos_on_atoms(log_filename, atomic_numbers, method='coeff'):
    atomic_numbers = np.array(atomic_numbers)
    coef, rgs = get_orbital_coefficients(log_filename)  # Here rgs (ranges) are important for both methods
    print("Attention! Atomic numbers will be interpreted as 1-based")
    atomic_numbers -= 1
    pdos_coeff = 0 * coef[0]
    # if method == 'coeff': # TO CHECK VERY CAREFULLY!!!
    #     coef /= np.sqrt(sum(coef ** 2, 0))  # Normalizing coefficients!
    #     for at_no in atomic_numbers:
    #         pdos_coeff += np.sum(coef[rgs[at_no]['range'][0]: rgs[at_no]['range'][1]] ** 2, 0)
    if method == 'mulliken':
        mulcontrib = get_mullik_contributions(log_filename)[0]
        for at_no in atomic_numbers:
            pdos_coeff += np.sum(mulcontrib[:, rgs[at_no]['range'][0]: rgs[at_no]['range'][1]], 1)

    return pdos_coeff


def getcontribs(clusterdata):
    # Returns atomic contributions to orbitals (Mulliken charges)
    if isinstance(clusterdata, str):
        parseddata = ccread(clusterdata)
    else:
        parseddata = clusterdata

    mulliken = get_mullik_contributions(clusterdata)
    rgs = get_atombasis(clusterdata)
    norb = parseddata.nmo
    natom = parseddata.natom
    contribs = np.zeros((norb, natom))
    for k_atom in range(natom):
        orbs_of_k_atom = range(rgs[k_atom]['range'][0], rgs[k_atom]['range'][1])
        for i_orbital in range(norb):
            contribs[i_orbital, k_atom] = np.sum(mulliken[0][i_orbital][orbs_of_k_atom])
    ipr = np.sum(contribs ** 2, 1)  # $\sum_i{c_i^4}$
    return ipr, contribs, parseddata.moenergies, parseddata.homos


def getpdos_general(gaudata, flags_list, method='mulliken'):
    # flaglist is a string of type '1Cd,S;Se,P'
    # Meaning - contribution of S-orbitals of atom Cd1 AND all P-orbitals of all Se-atoms
    # Nomeration is 1-based
    # 1Cd is equivalent to Cd1

    orbs = list(enumerate(gaudata.aonames))  # List of tuples (orbital_no, orbital_type_str)

    # Parsing indices of orbitals corresponding to flags
    orgroups = flags_list.split(';')
    orbrange = []
    for orgroup in orgroups:
        if len(orgroup) == 0:
            continue
        andgroup = orgroup.split(',')
        target_element = re.findall('[A-Z][a-z]*', andgroup[0])
        target_at_no = re.findall('\d+', andgroup[0])
        if len(target_element) == 0:
            target_element = ['']
        if len(target_at_no) == 0:
            search_pattern = target_element[0]
        else:
            search_pattern = target_element[0] + target_at_no[0] + '_'
        for orb in orbs:
            orbnamesplit = orb[1].split('_')
            if len(andgroup) == 2 and (search_pattern in (orbnamesplit[0] + '_')) and (andgroup[1] in orbnamesplit[1]):
                orbrange.append(orb[0])
            elif len(andgroup) == 1 and (search_pattern in (orbnamesplit[0] + '_')):
                orbrange.append(orb[0])

    # Constructing partial DOS coefficients
    # if method == 'coeff':
    #     coef /= np.sqrt(sum(coef ** 2, 0))  # Normalizing coefficients!
    #     for at_no in atomic_numbers:
    #         pdos_coeff += np.sum(coef[rgs[at_no]['range'][0]: rgs[at_no]['range'][1]] ** 2, 0)
    if method == 'mulliken':
        mc = MPA(gaudata)
        mc.calculate()
        mull_contrib = np.array(mc.aoresults)
        up_mull_contrib = mull_contrib[0]
        pdos_coeff = np.sum(up_mull_contrib[:, orbrange], 1)

    return pdos_coeff, gaudata.moenergies  # , orbrange, gaudata

def get_total_energy(clusterdata):
    if isinstance(clusterdata, str):
        clname = clusterdata.split('/')[0]
        parseddata = ccread(clusterdata)
    else:
        parseddata = clusterdata
        clname = get_formula(clusterdata)
    try:
        return parseddata.scfenergies[-1]
    except Exception as e:
        print('No scf done: ' + e.__class__ + ' ' + clname)
        return None
    

# def deprecated_get_kohn_sham(log_filename):
#     # The function works for spin-nonpolarized cases only
#     spec_raw = my_parse.getblock(log_filename, 'The electronic state', 'Condensed to atoms')
#     spectrum = {'up_occ': [], 'up_free': [], 'dn_occ': [], 'dn_free': []}
#     for spec_string in spec_raw:
#         splt_spec_string = spec_string.split()
#         evals = splt_spec_string[4:]
#         if len(evals) > 0 and len(evals[0]) > 10: # Sometimes Gaussian glues together the evals: -343.44321-342.44321...
#             evals_dummy = evals[0].strip().replace('-',' -').split()
#             if len(evals) == 1:
#                 evals = evals_dummy
#             else:
#                 evals = evals_dummy + evals[1:]
#         if 'Alpha  occ.' in spec_string:
#             spectrum['up_occ'].extend(evals)
#         elif ' Alpha virt.' in spec_string:
#             spectrum['up_free'].extend(evals)
#         elif 'Beta  occ.' in spec_string:
#             spectrum['dn_occ'].extend(evals)
#         elif 'Beta virt.' in spec_string:
#             spectrum['dn_free'].extend(evals)
#
#     for key in spectrum.keys():
#         spectrum[key] = np.array([float(i) for i in spectrum[key]])
#         # print(key, spectrum[key])
#
#     occ_states = np.append(spectrum['up_occ'], spectrum['dn_occ'])
#     free_states = np.append(spectrum['up_free'], spectrum['dn_free'])
#     fermi = 0.9 * np.max(occ_states) + 0.1 * np.min(free_states)
#
#     up_states = np.append(spectrum['up_occ'], spectrum['up_free'])
#     dn_states = np.append(spectrum['dn_occ'], spectrum['dn_free'])
#
#     return up_states, dn_states, fermi
