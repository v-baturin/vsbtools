import os
import sys
import numpy as np
import re
from genutils import my_io
from cclib.method import MPA
from cclib.io import ccread
from cclib.parser.utils import PeriodicTable as pt
from genutils.filesystem_tools import add_index
from functools import wraps

ptable = pt()
element_labels = np.array(ptable.element[:])


def path_or_ccobject(f):
    @wraps(f)
    def wrapped(*args, **kwargs):
        if isinstance(args[0], str):
            args = list(args)
            args[0] = ccread(args[0])
        return f(*args, **kwargs)
    return wrapped


@path_or_ccobject
def trace_relaxation(gaudata, outfolder, fname_root=None):
    os.makedirs(outfolder, exist_ok=True)
    gaudata.metadata['comments'] = []
    for k in range(len(gaudata.atomcoords)):
        gaudata.metadata['comments'].append('E_tot = {:6.5f}'.format(gaudata.scfenergies[k]))
        gaudata.writexyz(add_index(outfolder + '/' + fname_root.split('/')[-1].split('.')[0] + '.xyz', zerobased=True,
                                   respect_file_extension=True), indices=k)

@path_or_ccobject
def get_atom_symbols(gaudata):
    return element_labels[gaudata.atomnos]


@path_or_ccobject
def get_kohn_sham(gaudata, x_homo_fermi = 0.5):
    homo_idx = gaudata.homos[0]
    up_states = gaudata.moenergies[0]
    dn_states = gaudata.moenergies[-1]
    fermi = x_homo_fermi * up_states[homo_idx] + (1 - x_homo_fermi) * up_states[homo_idx + 1]

    return up_states, dn_states, fermi


@path_or_ccobject
def get_formula(gaudata, extended_out=False):
    """

    @param clusterdata: name of out file or a ccData object containing cclib data attributes
    @param extended_out: if True returns not only formula but lists of atomic symbols and number of atoms in a dict
    @return:
    """
    symbols_list = list(get_atom_symbols(gaudata))

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
    basis_ranges = get_atombasis(gaudata)
    return coeffs, basis_ranges


@path_or_ccobject
def get_atombasis(gaudata):
    # [{'label': {'Cd'}, 'range': [0, 24]},
    #  {'label': {'Cd'}, 'range': [24, 48]},...]
    ranges = []
    # basisbyatom = np.array(gaudata.atombasis)
    basisbyatom = gaudata.atombasis
    symbols_list = list(get_atom_symbols(gaudata))
    for idx, rangeslist in enumerate(basisbyatom):
        ranges.append({'label': {symbols_list[idx]}, 'range': [rangeslist[0], rangeslist[-1] + 1]})
    return ranges


@path_or_ccobject
def get_mullik_contributions(gaudata):
    """
    aoresults[s][MO][AO] = c(MO, AO) * \sum_{AO'}{c(MO, AO') * <AO|AO'>} -- Mullilken contribution of atomic orbital AO
    to molecular orbital MO. So that sum_{AO} {aoresults[s][MO][AO]} = <MO|MO> = 1.
    Here c(MO, AO) is axpansion coefficient: |MO>= \sum_{AO}{c(MO, AO) * |AO>}
    :param gaudata: result of parsing of gaussian file
    :return: aouresults[spin_component][MO][AO]
    """
    m = MPA(gaudata)
    m.calculate()
    return np.array(m.aoresults)

@path_or_ccobject
def get_atomic_contribs(gaudata):
    # Returns atomic contributions to orbitals (Mulliken charges)
    mulliken = get_mullik_contributions(gaudata)
    rgs = get_atombasis(gaudata)
    n_MO = gaudata.nmo
    n_atom = gaudata.natom
    atomic_contribs = np.zeros((n_MO, n_atom))
    for k_atom in range(n_atom):
        orbs_of_k_atom = range(rgs[k_atom]['range'][0], rgs[k_atom]['range'][1])
        for i_orbital in range(n_MO):
            atomic_contribs[i_orbital, k_atom] = np.sum(mulliken[0][i_orbital][orbs_of_k_atom])
    return atomic_contribs


def get_pdos_on_atoms(log_filename, atomic_numbers, method='mulliken'):
    atomic_numbers = np.array(atomic_numbers)
    coef, rgs = get_orbital_coefficients(log_filename)  # Here rgs (ranges) are important for both methods
    print("Attention! Atomic numbers will be interpreted as 1-based")
    atomic_numbers -= 1
    pdos_coeff = 0 * coef[0]
    mulcontrib = get_mullik_contributions(log_filename)[0]
    for at_no in atomic_numbers:
        pdos_coeff += np.sum(mulcontrib[:, rgs[at_no]['range'][0]: rgs[at_no]['range'][1]], 1)

    return pdos_coeff


@path_or_ccobject
def get_ipr(parseddata):
    # Returns inverse population ratio (IPR)
    contribs = get_atomic_contribs(parseddata)
    ipr = np.sum(contribs ** 2, 1)  # $\sum_i{c_i^4}$
    return ipr


def getpdos_general(gaudata, flags_list, method='mulliken'):
    # flaglist is a string of type '1Cd,S;Se,P'
    # Meaning - contribution of S-orbitals of atom Cd1 AND all P-orbitals of all Se-atoms
    # Numeration is 1-based
    # 1Cd is equivalent to Cd1

    orbs = list(enumerate(gaudata.aonames))  # List of tuples (orbital_no, orbital_type_str)

    # Parsing indices of orbitals corresponding to flags
    orgroups = flags_list.split(';')
    orbrange = []
    for orgroup in orgroups:
        if len(orgroup) == 0:
            continue
        andgroup = [x.strip() for x in orgroup.split(',')]
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
    mull_contrib = get_mullik_contributions(gaudata)
    up_mull_contrib = mull_contrib[0]
    if len(mull_contrib) == 2:
        dn_mull_contrib = mull_contrib[1]
        pdos_coeff = (np.sum(up_mull_contrib[:, orbrange], 1), np.sum(dn_mull_contrib[:, orbrange], 1))
    else:
        dn_mull_contrib = None
        pdos_coeff = (np.sum(up_mull_contrib[:, orbrange], 1), None)

    print(np.sum(up_mull_contrib))
    return pdos_coeff, gaudata.moenergies  # , orbrange, gaudata


@path_or_ccobject
def get_total_energy(parseddata):
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
