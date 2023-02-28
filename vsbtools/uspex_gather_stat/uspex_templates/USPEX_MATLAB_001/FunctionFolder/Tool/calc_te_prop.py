# calc_te_prop.py ---
#
# Filename: calc_te_prop.py
# Author: Fred Qi
# Created: 2015-11-11 20:58:55(+0300)
#
# Last-Updated: 2015-11-20 12:57:08(+0300) [by Fred Qi]
#     Update #: 863
#

# Commentary:
# Calculate thermoelectronic properties using output files of BoltzTraP.
# The interface function is calculate_thermoelectric_properties()
# The required arguments are the file names of case.trace and case.condtens,
# and also the expected temperature.

# Change Log:
#
#
#
from __future__ import print_function
from os import path
import numpy as np
from scipy import signal
from scipy.signal import argrelextrema


def load_BoltzTraP_data(datafile, cols=None,
                        temperature=300.0, temperature_col=2,
                        eps=1e-6):
    """Load BoltzTraP data for given columns and given temperature.

    :param datafile: name of a BoltzTraP generated file.
    :param cols: columns to be retrieved, numbers are 1-based.
    :param temperature: float, expected temperature (K).
    :param temperature_col: the index of the column containing temperatures.
    :param eps: accuracy when comparing temperatures.
    :returns: data selected by the given cols and temperature.
    :rtype: numpy.ndarray

    """
    data = np.loadtxt(datafile)
    # assert max(cols) <= data.shape[1]
    if cols is None:
        cols_sel = np.ones((data.shape[1],), dtype=np.bool)
    else:
        cols_sel = [idx+1 in cols for idx in range(data.shape[1])]
        cols_sel = np.array(cols_sel, dtype=np.bool)
    tk_sel = np.abs(data[:, temperature_col-1] - temperature) <= eps
    return data[tk_sel][:, cols_sel]


def load_BoltzTraP_trace(tracefile, temperature, eps=1e-6):
    """Load data from a BoltzTraP generated trace file.

    :param tracefile: name of a BoltzTraP generated case.trace file.
    :param temperature: float, expect temperature (K).
    :param eps: accuracy when comparing temperatures.
    :returns: data selected by the given temperature.
    :rtype: numpy.ndarray

    """
    # Columns for calculating ZT
    # (1, mu), (3, N), (4, DOS)
    # (5, Seebeck), (6, sigma/tau), (8, kappa0)
    cols = (1, 4, 5, 6, 8)
    data = load_BoltzTraP_data(tracefile, cols, temperature, eps=eps)
    return data


def load_BoltzTraP_condtens(condtensfile, temperature, eps=1e-6):
    """Load data from a BoltzTraP generated condtens file.

    :param condtensfile: name of a BoltzTraP generated case.condtens file.
    :param temperature: float, expect temperature (K).
    :param eps: accuracy when comparing temperatures.
    :returns: data selected by the given temperature.
    :rtype: numpy.ndarray

    """
    cols = (1, 4, 12, 13, 21, 22, 30)  # Columns for calculating ZT
    data = load_BoltzTraP_data(condtensfile, cols, temperature, eps=eps)
    return data


def calculate_zt(data, c_sigma, c_seebeck, c_kappa, temperature):
    """To calculate the ZT property using given data.

    :param data: data loaded from a BoltzTraP generated file.
    :param c_sigma: index of the sigma column.
    :param c_seebeck: index of the seebeck coefficient column.
    :param c_kappa: index of the kappa0 column.
    :param temperature: the expected temperature.
    :returns: a tuple of sigma*S*S and ZT
    :rtype: (float, float)

    """
    cols = (c_sigma, c_seebeck, c_kappa)
    sigma, seebeck, kappa = (data[:, c-1] for c in cols)
    sigma_SS = sigma*seebeck*seebeck*1.0e-10
    # this is a magic
    zeros = np.abs(kappa) < 1e-6
    non_zeros = np.logical_not(zeros)
    # print(np.mean(kappa[non_zeros]), zeros.shape, np.sum(zeros), np.sum(non_zeros))
    kappa[zeros] = np.mean(kappa[non_zeros])

    ZT = sigma*seebeck*seebeck*temperature/kappa
    return np.vstack((sigma, seebeck, kappa, sigma_SS, ZT)).T


def calculate_te_prop_trace(tracefile, temperature):
    """Calculate thermoelectric properties using BoltzTraP output trace files.

    :param tracefile: name of casename.trace file with full path,
                      or a corresponding file object.
    :param temperature: expected temperature.
    :returns: data with columns: (mu, dos),
              (tr(sigma), tr(seebeck), tr(kappa0), tr(sigma_SS), tr(ZT))
    :rtype: tuple(numpy.ndarray)

    """
    data = load_BoltzTraP_trace(tracefile, temperature)
    prop_tr = calculate_zt(data, 4, 3, 5, temperature)
    return data[:, :2], prop_tr


def calculate_te_prop_condtens(condtensfile, temperature):
    """Calculate thermoelectric properties using BoltzTraP output files.

    :param tracefile: name of casename.trace file with full path,
                      or a corresponding file object.
    :param condtensfile: name of casename.condtens file with full path,
                         or a corresponding file object.
    :param temperature: expected temperature.
    :returns: data with columns: mu,
              (sigma_xx, S_xx, kappa0_xx, sigma_SS_xx, ZT_xx),
              (sigma_zz, S_zz, kappa0_zz, sigma_SS_zz, ZT_zz)
    :rtype: tuple(numpy.ndarray)

    """
    data = load_BoltzTraP_condtens(condtensfile, temperature)
    prop_xx = calculate_zt(data, 2, 4, 6, temperature)
    prop_zz = calculate_zt(data, 3, 5, 7, temperature)
    return data[:, 0], prop_xx, prop_zz


def calculate_te_prop(tracefile, condtensfile, temperature):
    """Calculate thermoelectric properties using BoltzTraP output files.

    :param tracefile: name of casename.trace file with full path,
                      or a corresponding file object.
    :param condtensfile: name of casename.condtens file with full path,
                         or a corresponding file object.
    :param temperature: expected temperature.
    :returns: (mu, dos),
              dict of results of (sigma, S, kappa0, sigma_SS, ZT)
              dict keys are:
                - tr: for trace
                - xx: for xx direction
                - ...
                - yy: for yy direction
                - ...
                - zz: for zz direction
    :rtype: numpy.ndarray, dict(numpy.ndarray)

    """
    prop_all = dict()
    data = load_BoltzTraP_data(tracefile, temperature=temperature)
    common = data[:, (0, 3)]
    prop_all['ZT'] = calculate_zt(data, 6, 5, 8, temperature)

    data = load_BoltzTraP_data(condtensfile, temperature=temperature)
    cols_xx = (4, 13, 22)
    offset = dict(ZT_xx=0, ZT_xy=1, ZT_xz=2, ZT_yx=3,
                  ZT_yy=4, ZT_yz=5, ZT_zx=6, ZT_zy=7, ZT_zz=8)
    for key in ['ZT_xx', 'ZT_yy', 'ZT_zz']:
        c_sigma, c_seebeck, c_kappa = (c + offset[key] for c in cols_xx)
        prop_all[key] = calculate_zt(data, c_sigma, c_seebeck, c_kappa,
                                     temperature)
    return common, prop_all


class column_mapping:
    """Mapping a column from the name to its corresponding data column."""

    def __init__(self):
        """Initialize a class instance."""
        properties = ('sigma', 'seebeck', 'kappa0', 'sigma_SS', 'ZT')
        self.index = dict()
        for idx, key in enumerate(properties):
            self.index[key] = idx

    def __getitem__(self, column):
        return self.index[column]


def zt_camel_hump_slice(zt):
    """Analyze the value of ZT and extract a slice of ZT.
    Typical ZT curve is with a camel hump likes shape, extract the slice
    containing the hump is very useful.

    :param zt: thermoelectric property figure of merit.
    :returns: positions of the two maxima and range of the hump slice
    :rtype: (int, int, int, int)

    """
    # Get all positions of local maxima
    maxpos = argrelextrema(zt, np.greater)[0]
    # minpos = argrelextrema(zt, np.less)[0]
    # Sort and checkout the two largest maxima
    maxima_sort = np.argsort(zt[maxpos])
    # Determine and order the positions of hump maxima
    pos_l, pos_r = maxpos[maxima_sort[-2:]]
    pos_l, pos_r = min(pos_l, pos_r), max(pos_l, pos_r)

    # Smoothing the zt curve with a Hann window in order to
    # remove small ripples of local minima.
    # Theoretically, this smoothing will not change positions of extrema.
    win = signal.hann(32)
    ztv = signal.convolve(zt, win, mode='same')/sum(win)
    # The first minimum left to the left hump
    l_sel = argrelextrema(ztv[:pos_l], np.less)[0]
    edge_l = l_sel[-1] if len(l_sel) > 0 else 0
    # The first minimum right to the right hump
    r_sel = argrelextrema(ztv[pos_r:], np.less)[0]
    edge_r = r_sel[0] + pos_r if len(r_sel) > 0 else len(zt)
    # print(l_sel, r_sel)
    # print(type(l_sel), l_sel, type(r_sel), r_sel, len(r_sel))
    # return the detected positions
    return pos_l, pos_r, edge_l, edge_r


def extract_log_info(common, prop_all, optgoal='ZT'):
    """Extrac useful information for printing out to a text log file.

    :param common: common information for (mu, dos)
    :param prop_all: dictionary of properties, keys are: tr, xx, yy, zz
                     each item has columns in the order:
                     (sigma, seebeck, kappa0, sigma_SS, ZT)
    :returns: (mu_max_1, ZT_max_1, mu_max_2, ZT_max_2) for trace, xx, yy, zz
              slices for saving to a file
              BUT currently, only the trace is supported.
    :rtype: numpy.ndarray, numpy.ndarray

    """
    ZT = prop_all[optgoal][:, -1]
    pl, pr, l_sel, r_sel = zt_camel_hump_slice(ZT)
    p1, p2 = (pl, pr) if ZT[pl] >= ZT[pr] else (pr, pl)
    data_log = np.array([common[p1, 0], ZT[p1],
                         common[p2, 0], ZT[p2]])
    # Range selected for storage
    sel = np.zeros_like(ZT, dtype=np.bool)
    sel[l_sel:r_sel+1] = True
    # Try to save all directions
    data_save = dict()
    for key in prop_all.keys():
        data_save[key] = np.hstack((common[sel, :], prop_all[key][sel, :]))
    # print(data_log, data_save)
    return data_log, data_save


def write_summary_with_order(filename, data, index=0):
    """Write summary data to the given log file.
       The data should be a list fit for writing in one line.
    """
    if isinstance(filename, str):
        if path.exists(filename):
            summary = open(filename, 'a')
        else:
            hdr_line1 = '#' + "trace".center(39) + '\n'
            hdr_items = ["mu_max_1", "ZT_max_1",
                         "mu_max_2", "ZT_max_2"]
            hdr_items = [item.rjust(9) for item in hdr_items]
            hdr_line2 = "# ID  " + ' '.join(hdr_items) + '\n'
            summary = open(filename, 'w')
            summary.write(hdr_line1 + hdr_line2)
    elif hasattr(filename, 'wirte'):
        # filename is already an opened file/stream.
        summary = filename

    
    if index is None:
        index_str = ' '*6
    else:
        index_str = '{index:5d} '.format(index=index)
    line = ' '.join(['{v:9f}'.format(v=v) for v in data])
    summary.write(index_str + line + '\n')
    summary.close()


def write_summary_header(filename):
    """Write the header line to a file."""
    if isinstance(filename, str):
        hdr_line1 = '#' + "trace".center(39) + '\n'
        hdr_items = ["mu_max_1", "ZT_max_1",
                     "mu_max_2", "ZT_max_2"]
        hdr_items = [item.rjust(9) for item in hdr_items]
        hdr_line2 = "# ID  " + ' '.join(hdr_items) + '\n'
        summary = open(filename, 'w')
        summary.write(hdr_line1 + hdr_line2)
        summary.close()


def write_summary(filename, data):
    """Write summary data to the given log file."""
    if isinstance(filename, str):
        summary = open(filename, 'w')
    elif hasattr(filename, 'write'):
        summary = filename

    line = ' '.join(['{v:9f}'.format(v=v) for v in data])
    summary.write(line + '\n')
    summary.close()


def make_storage_header(direction, field_width=16):
    """Get the header string for the storage of given direction."""
    key = 'trace'
    if direction.find('_') > 0:
        key = direction.split('_')[-1]
    header = list()
    header.append('mu'.center(field_width-1))
    header.append('DOS'.center(field_width))
    name_properties = ('sigma', 'seebeck', 'kappa0', 'sigma_SS', 'ZT')
    for name in name_properties:
        name_dir = name + '_' + key
        header.append(name_dir.center(field_width))
    return ' '.join(header)


def find_ticker_base(xmin, xmax, interval):
    from matplotlib.ticker import MultipleLocator

    delta = abs(xmax - xmin)/(interval + 1)
    base = round(delta/10.0)*10.0
    digits = 0
    while base < 1e-6:
        base = round(delta, digits)
        digits += 1

    return MultipleLocator(base)


def plot_te_figures(data, show_slice=True):
    """Plot curves of thermoelectronic properties."""
    import matplotlib.pyplot as plt

    index = column_mapping()
    if show_slice:
        ZT_xx, ZT_zz = data[:, index['ZT_xx']], data[:, index['ZT_zz']]
        _, _, pos_lxx, pos_rxx = zt_camel_hump_slice(ZT_xx)
        _, _, pos_lzz, pos_rzz = zt_camel_hump_slice(ZT_zz)
        ll, rr = min(pos_lxx, pos_lzz), max(pos_rxx, pos_rzz)
    else:
        ll, rr = 0, len(data[:, 0])

    mu = data[ll:rr, index['mu']]
    dos = data[ll:rr, index['dos']]
    S_xx, S_zz = data[ll:rr, index['S_xx']]*1e6, data[ll:rr, index['S_zz']]*1e6
    sigma_SS_xx = data[ll:rr, index['sigma_SS_xx']]
    sigma_SS_zz = data[ll:rr, index['sigma_SS_zz']]
    ZT_xx, ZT_zz = data[ll:rr, index['ZT_xx']], data[ll:rr, index['ZT_zz']]
    # mu = data[:, index['mu']]
    # dos = data[:, index['dos']]
    # S_xx, S_zz = data[:, index['S_xx']]*1e6, data[:, index['S_zz']]*1e6
    # sigma_SS_xx = data[:, index['sigma_SS_xx']]
    # sigma_SS_zz = data[:, index['sigma_SS_zz']]
    # ZT_xx, ZT_zz = data[:, index['ZT_xx']], data[:, index['ZT_zz']]

    ticker_locator = find_ticker_base(mu[0], mu[-1], 20)

    plt.figure()
    ax = plt.subplot(4, 1, 1)
    plt.plot(mu, dos)
    plt.xlim(mu[0], mu[-1])
    ax.xaxis.set_major_locator(ticker_locator)

    ax = plt.subplot(4, 1, 2)
    plt.plot(mu, S_xx)
    plt.plot(mu, S_zz)
    plt.xlim(mu[0], mu[-1])
    ax.xaxis.set_major_locator(ticker_locator)

    ax = plt.subplot(4, 1, 3)
    plt.plot(mu, sigma_SS_xx, 'm', label='$\sigma S^2_{xx}$')
    plt.plot(mu, sigma_SS_zz, 'k', label='$\sigma S^2_{zz}$')
    plt.xlim(mu[0], mu[-1])
    ax.xaxis.set_major_locator(ticker_locator)
    plt.legend()

    ax = plt.subplot(4, 1, 4)
    plt.plot(mu, ZT_xx, 'r', label='ZTxx')
    plt.plot(mu, ZT_zz, 'b', label='ZTzz')
    win = signal.hann(32)
    ZT_xx_p = signal.convolve(ZT_xx, win, mode='same')/np.sum(win)
    ZT_zz_p = signal.convolve(ZT_zz, win, mode='same')/np.sum(win)
    plt.plot(mu, ZT_xx_p, 'r--', lw=1, label='ZT_xx^p')
    plt.plot(mu, ZT_zz_p, 'b--', lw=1, label='ZT_zz^p')
    plt.xlim(mu[0], mu[-1])
    ax.xaxis.set_major_locator(ticker_locator)
    # print(ll, rr)
    plt.axvline(x=mu[pos_lxx], color='r', lw=1)
    plt.axvline(x=mu[pos_rxx], color='r', lw=1)
    plt.axvline(x=mu[pos_lzz], color='b', lw=1)
    plt.axvline(x=mu[pos_rzz], color='b', lw=1)
    plt.legend()
    plt.show()


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(prog="calc_te_prop")

    parser.add_argument('--casename', '-n', dest='casename',
                        type=str, required=True,
                        help="The casename of BoltzTraP files.")
    parser.add_argument('--temperature', '-t',
                        dest='T', type=float, default=300.0,
                        help="A temperature of interesets.")
    parser.add_argument('--threshold', '-r',
                        dest='threshold', type=float,
                        help="A threshold set to drop bad structures.")
    parser.add_argument('--optimization-goal', '-g',
                        dest='optgoal', default='ZT',
                        choices=['ZT', 'ZT_xx', 'ZT_yy', 'ZT_zz'],
                        help="The property to be optimized.")
    parser.add_argument('--output-directory', '-o',
                        dest='outfolder', type=str,
                        help="Directory to write BoltzTraP results.")
    parser.add_argument('--output-prefix', '-p',
                        dest='prefix', type=str,
                        help="Directory to write BoltzTraP results.")
    # parser.add_argument('--index', '-i',
    #                     dest='index', type=int,
    #                     help="The index of the structure in USPEX.")
    parser.add_argument('--summary', '-s',
                        dest='summary', default='summary.txt',
                        help="Summary of thermoelectric properties.")
    parser.add_argument('input_directory', type=str,
                        help="The directory to read BoltzTraP output files from.")

    args = parser.parse_args()
    # print(args)                # For debug purpose only

    fn_trace = path.join(args.input_directory, args.casename + '.trace')
    fn_condtens = path.join(args.input_directory, args.casename + '.condtens')
    dt, prop_all = calculate_te_prop(fn_trace, fn_condtens, args.T)

    prop = prop_all[args.optgoal]
    row = prop[:, -1].argmax()
    values = prop[row, :].tolist()
    # MATLAB interface, added '<CALLRESULT>' for python_uspex compatibility.
    ret = values[-1] if values[-1] < 10.0 else -1.0*values[-1]
    if np.isnan(ret) or np.isinf(ret):
        ret = -88888.0

    dt_log, dt_log_str = None, ''
    if args.outfolder and args.prefix:
        dt_log, dt_save = extract_log_info(dt, prop_all, args.optgoal)
        fn_summary = path.join(args.outfolder, args.summary)
        if not path.exists(fn_summary):
            write_summary_header(fn_summary)
        fn_log = path.join(args.outfolder, args.prefix + '-summary.txt')
        write_summary(fn_log, dt_log.tolist())
        if args.threshold is None or values[-1] > args.threshold:
            for direction in ('ZT', 'ZT_xx', 'ZT_yy', 'ZT_zz'):
                fn_output = path.join(args.outfolder, 
                                      args.prefix + direction + '.txt')
                np.savetxt(fn_output, dt_save[direction],
                           fmt='%.10e',
                           header=make_storage_header(direction, 16))
        # print(dt_log)

    if dt_log is not None:
        values = dt_log.tolist()
        values = [str(val) for val in values]
        dt_log_str = ', '.join(values) + ', '
    print('<CALLRESULT>' + dt_log_str + '{v}'.format(v=ret))
#
# calc_te_prop.py ends here
