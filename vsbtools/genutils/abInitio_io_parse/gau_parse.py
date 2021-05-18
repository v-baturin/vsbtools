import sys
import numpy as np
import re
from .. import my_io
from cclib.method import MPA
from cclib.io import ccread
from cclib.parser.utils import PeriodicTable as pt

ptable = pt()
element_labels = np.array(ptable.element[:])

def path_or_ccobject(f):
    def wrapped(*args, **kwargs):
        if isinstance(args[0], str):
            args = list(args)
            args[0] = ccread(args[0])
        return f(*args, **kwargs)
    return wrapped

@path_or_ccobject
def get_kohn_sham(gaudata):

    homo_idx = gaudata.homos[0]
    up_states = gaudata.moenergies[0]
    dn_states = gaudata.moenergies[-1]
    fermi = 0.9 * up_states[homo_idx] + 0.1 * up_states[homo_idx + 1]

    return up_states, dn_states, fermi

@path_or_ccobject
def get_formula(gaudata, extended_out=False):
    """

    @param clusterdata: name of out file or a ccData object containing cclib data attributes
    @param extended_out: if True returns not only formula but lists of atomic symbols and number of atoms in a dict
    @return:
    """
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

@path_or_ccobject
def get_orbital_coefficients(gaudata):
    coeffs = gaudata.mocoeffs[0]
    return coeffs

@path_or_ccobject
def get_atombasis(gaudata):
    # [{'label': {'Cd'}, 'range': [0, 24]},
    #  {'label': {'Cd'}, 'range': [24, 48]},...]
    rgs = []
    # basisbyatom = np.array(gaudata.atombasis)
    basisbyatom = gaudata.atombasis
    symbols_list = list(element_labels[gaudata.atomnos])
    for idx, rangeslist in enumerate(basisbyatom):
        rgs.append({'label': {symbols_list[idx]}, 'range': [rangeslist[0], rangeslist[-1] + 1]})
    return rgs

@path_or_ccobject
def get_mullik_contributions(gaudata):
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

@path_or_ccobject
def getcontribs(parseddata):
    # Returns atomic contributions to orbitals (Mulliken charges)

    mulliken = get_mullik_contributions(parseddata)
    rgs = get_atombasis(parseddata)
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

@path_or_ccobject
def get_total_energy(parseddata):
    # if isinstance(clusterdata, str):
    #     clname = clusterdata.split('/')[0]
    #     parseddata = ccread(clusterdata)
    # else:
    #     parseddata = clusterdata
    clname = get_formula(parseddata)
    try:
        return parseddata.scfenergies[-1]
    except Exception as e:
        print('No scf done: ' + e.__class__ + ' ' + clname)
        return None

@path_or_ccobject
def get_gap(clusterdata):
    try:
        homo_indices = clusterdata.homos
    except AttributeError:
        return np.nan
    if len(homo_indices) == 2:
        gap = np.min([clusterdata.moenergies[0][homo_indices[0] + 1],
                      clusterdata.moenergies[1][homo_indices[1] + 1]]) - \
              np.max([clusterdata.moenergies[0][homo_indices[0]],
                      clusterdata.moenergies[1][homo_indices[1]]])
    else:
        gap = clusterdata.moenergies[0][homo_indices[0] + 1] - \
              clusterdata.moenergies[0][homo_indices[0]]
    return gap
    

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
