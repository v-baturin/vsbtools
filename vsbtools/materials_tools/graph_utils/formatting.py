import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from collections.abc import Iterable


def cm2inch(*tupl):
    """
    Usage: plt.figure(figsize=cm2inch(12.8, 9.6))
    @param tupl: tuple of dimensions in centimeters
    @return: tuple of corresponding sizes in inches to be used with matplotlib
    """
    inch = 2.54
    if isinstance(tupl[0], Iterable):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def set_subplots_pos_cm(axdims_cm: list, figsize_cm=None, fig=None, subpl_shape=None, wspace_cm=None, hspace_cm=None):
    """
    Sets position and size of all subplots inside a figure (in centimeters)
    Usage: set_subplots_cm([leftpad, downpad, width, height])
    @param axdims_cm: [x0, y0, w, h] coordinates of left-down corner, followed by dimensions of subplots box
    @param figsize_cm: tuple
    @param fig: specifies the figure to which modifications are applied. if None, gcf is used
    @param subpl_shape: tuple of (nrows, ncols) of subfigures array
    @param hspace_cm: vertical padding (in cm): works only if all subfigures are the same
    @param wspace_cm: horizontal padding (in cm): works only if all subfigures are the same
    @return:
    """
    if fig is None:
        fig = plt.gcf()

    if figsize_cm is None:
        wfig, hfig = fig.get_size_inches()
    else:
        wfig, hfig = cm2inch(tuple(figsize_cm))
        fig.set_size_inches(wfig, hfig)

    x0, y0, w, h = cm2inch(axdims_cm)

    lpad = x0/wfig
    bpad = y0/hfig
    rpad = (x0 + w) / wfig
    tpad = (y0 + h) / hfig

    nx, ny = subpl_shape

    if wspace_cm is not None and nx > 1:
        dx_rel = cm2inch(wspace_cm)[0] / (wfig / (nx - 1) + cm2inch(wspace_cm)[0])
    else:
        dx_rel = 0.2

    if hspace_cm is not None and ny > 1:
        dy_rel = cm2inch(hspace_cm)[0] / (hfig / (ny - 1) + cm2inch(hspace_cm)[0])
    else:
        dy_rel = 0.2

    fig.subplots_adjust(left=lpad, right=rpad, bottom=bpad, top=tpad, wspace=dx_rel, hspace=dy_rel)


def set_ax_position_cm(ax, dims_cm):
    """

    :param ax: axes
    :param dims_cm: (left offset, low_offset, width, height)
    :return:
    """
    fig_obj = ax.figure
    figsize_in = fig_obj.get_size_inches()
    dims_in = cm2inch(dims_cm)
    fractional_left = dims_in[0] / figsize_in[0]
    fractional_btm  = dims_in[1] / figsize_in[1]
    fractional_w = dims_in[2] / figsize_in[0]
    fractional_h = dims_in[3] / figsize_in[1]
    ax.set_position([fractional_left, fractional_btm, fractional_w, fractional_h])
    ax.autoscale(enable=True)

