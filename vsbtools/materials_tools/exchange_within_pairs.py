import numpy as np
import sys

sys.path.append("..")
from genutils.format_tools import stoich2formula
from .aux_routines import bin_paths, list_fmt2table


def min_d2e_2d(energies_table, diags=False, shift=False):
    # shift variable corresponds to the exchange of two molecules of second type
    if energies_table.ndim != 2:
        raise ValueError('Energy array must be two-dimensional numpy.ndarray of floats')
    core_energies = energies_table[1:-1, 1:-1]
    v_d2E = energies_table[0:-2, 1:-1] + energies_table[2:, 1:-1] - 2 * core_energies
    h_d2E = energies_table[1:-1, 0:-2] + energies_table[1:-1, 2:] - 2 * core_energies
    min_d2E = np.minimum(h_d2E, v_d2E)
    if shift:
        h_offset2_d2E = 100 * np.ones(h_d2E.shape)
        h_offset2_d2E[:, 1:-1] = energies_table[1:-1, 0:-4] + energies_table[1:-1, 4:] - 2 * core_energies[:, 1:-1]
        min_d2E = np.minimum(min_d2E, h_offset2_d2E)
    if diags is True:
        d1_d2E = energies_table[0:-2, 0:-2] + energies_table[2:, 2:] - 2 * core_energies
        d2_d2E = energies_table[0:-2, 2:] + energies_table[2:, 0:-2] - 2 * core_energies
        min_d2E = np.minimum(np.minimum(min_d2E, d1_d2E), d2_d2E)
    return min_d2E


def min_d2e_nd(energies_fmt, ref_energies=None, for_CNH=False, shifts=None, astable=False, return_mask=None):
    """
    Second derivative for systems of any dimensionality
    Turns formatted list with elements of format [i_0, i_1, i_2, ..., i_n, E(i_0, i_1, i_2, ..., i_n)]
    into a list of [j_0, j_1, j_2, ..., j_n, d2E(j_0, j_1, j_2, ..., j_n)] for all j_k which allow calculation of
    second derivatives
    @param astable:
    @param ref_energies:
    @param energies_fmt: List
    @param diags:
    @param shift:
    @return:
    """
    if ref_energies is None:
        ref_energies = energies_fmt
    if return_mask:
        has_d2e_mask = np.zeros(len(energies_fmt), dtype=bool)

    all_comps = ref_energies[:, :-1]
    ndims = energies_fmt.shape[1] - 1

    d2e = []

    if shifts is None:
        shifts = np.eye(ndims)

    n_shifts = shifts.shape[0]

    for el_i, element in enumerate(energies_fmt):
        d2e_components = np.empty(n_shifts)
        d2e_components[:] = np.nan
        composition = element[:-1]
        for i in range(n_shifts):
            comp_plus = composition + shifts[i]
            comp_minus = composition - shifts[i]
            where_plus = np.where((all_comps == comp_plus).all(axis=1))[0]
            where_minus = np.where((all_comps == comp_minus).all(axis=1))[0]
            if where_plus and where_minus:
                d2e_components[i] = ref_energies[where_plus, -1] + ref_energies[where_minus, -1] - 2 * element[-1]

        if not np.any(np.isnan(d2e_components)):
            if for_CNH:
                where_compos = np.where((all_comps == composition).all(axis=1))[0]
                ref_energy_compos = ref_energies[where_compos, -1]
                below_ground = ref_energy_compos - element[-1]
                d2e.append(list(composition) + [np.min(np.hstack((0.5 * d2e_components, below_ground)))])
            else:
                d2e.append(list(composition) + [np.min(d2e_components)])
            if return_mask:
                has_d2e_mask[el_i] = True

    d2e_np = np.array(d2e)

    if astable:
        assert ndims == 2, "For table representation a system must be binary"
        outp = list_fmt2table(d2e_np)[0]
    else:
        outp = d2e_np

    if return_mask:
        outp = (outp, has_d2e_mask)

    return outp


def min_exch_en(fmt_data, atomlabels, outfile):
    avail_stoich = fmt_data[:, :-1].tolist()
    all_en = fmt_data[:, -1]
    fmt_frag = []
    with open(outfile, 'w') as out_fid:
        for entry in fmt_data:
            current_stoich = entry[:-1]
            current_formula = stoich2formula(current_stoich, atomlabels)
            cl_en = 2 * entry[-1]
            l_frags, r_frags = bin_paths(2 * current_stoich)
            reactives = []
            react_ens = []
            for l_frag, r_frag in zip(l_frags, r_frags):
                if (l_frag in avail_stoich) and (r_frag in avail_stoich) and (l_frag != list(current_stoich)):
                    reactives.append({'l': l_frag, 'r': r_frag})
                    l_idx = avail_stoich.index(l_frag)
                    r_idx = avail_stoich.index(r_frag)
                    e_react = all_en[l_idx] + all_en[r_idx] - cl_en  # E_parts - E_whole
                    react_ens.append([e_react])
            if len(react_ens) > 0:
                e_frag_min = np.min(react_ens)
                min_idx = react_ens.index(e_frag_min)
                l_formula = stoich2formula(reactives[min_idx]['l'], atomlabels)
                r_formula = stoich2formula(reactives[min_idx]['r'], atomlabels)
                out_fid.write(
                    '2 * %s -> %s + %s, E_exch = %6.3f\n' % (current_formula, l_formula, r_formula, e_frag_min))
                fmt_frag.append(list(current_stoich) + [e_frag_min])
    return fmt_frag
