import copy
import numpy as np
from scipy import interpolate
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from .aux_routines import smeared_spectrum, heatmap, annotate_heatmap
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Spectra and spectral characteristics plots:

def smeared_dos(limits,
                spec_up,
                smear_sigma=None, smear_type='Gaussian',
                spec_dn=None,
                eaxis=None,
                weights_up=None, weights_dn=None,
                pdos_up_coef=None, pdos_dn_coef=None, normalization=True):
    """
    
    @param limits: Span of energies [emin, emax]
    @param spec_up: ndarray of energy values of energies for spin up states
    @param spec_dn: ndarray of energy values of energies for spin down states
    @param smear_sigma: smearing of dos graph
    @param eaxis: interpolation grid of energies
    @param weights_up: weights (degeneracies) of spin up states
    @param weights_dn: weights (degeneracies) of spin down states
    @param pdos_up_coef: partial dos coefficients for spin up states
    @param pdos_dn_coef: partial dos coefficients for spin down states
    @param normalization:
    @param smear_type:
    @return: eaxis, {'updos': updos, 'dndos': dndos}, {'pdos_up': pdos_up, 'pdos_dn': pdos_dn}
    """
    if smear_sigma is None:
        smear_sigma = np.min(np.sort(spec_up)[1:] - np.sort(spec_up)[:-1]) / 2
    if spec_dn is None:
        spec_dn = []
    if eaxis is None:
        eaxis = np.arange(limits[0] - smear_sigma, limits[1] + smear_sigma, (limits[1] - limits[0]) / 9999)
    if weights_up is None:
        weights_up = np.ones(spec_up.shape)
    updos, _ = smeared_spectrum(spec_up, weights_up, smear_sigma, x_grid=eaxis, sm_type=smear_type)
    if len(spec_dn) > 0:
        if weights_dn is None:
            weights_dn = np.ones(spec_dn.shape)
        dndos, _ = smeared_spectrum(spec_dn, weights_dn, smear_sigma, x_grid=eaxis, sm_type=smear_type)
    else:
        dndos = np.array([])

    if pdos_up_coef is not None:
        pdos_up, _ = smeared_spectrum(spec_up, pdos_up_coef, smear_sigma, x_grid=eaxis, sm_type=smear_type)
    else:
        pdos_up = np.array([])

    if pdos_dn_coef is not None:
        pdos_dn, _ = smeared_spectrum(spec_dn, pdos_dn_coef, smear_sigma, x_grid=eaxis, sm_type=smear_type)
    else:
        pdos_dn = np.array([])

    if normalization:

        if isinstance(normalization, bool):
            max_dos = np.max(updos)
            if len(spec_dn) > 0:
                max_dos = np.max((max_dos, np.max(dndos)))
        else:
            max_dos = normalization
        updos /= max_dos
        if len(spec_dn) > 0:
            dndos /= max_dos
        if pdos_up_coef is not None:
            pdos_up /= max_dos
        if pdos_dn_coef is not None:
            pdos_dn /= max_dos

    return eaxis, {'updos': updos, 'dndos': dndos}, {'pdos_up': pdos_up, 'pdos_dn': pdos_dn}


def draw_spectrum(span,
                  specup, specdn=None,
                  weights_up=None, weights_dn=None,
                  e_fermi=0, units='eV',
                  sigma=0.5, smear_type='Gaussian',
                  partial_up=None, partial_dn=None,
                  shareax=None,  # total_axes=1, curr_ax_no=1, curr_fig_no=None,
                  normalization=True,
                  dos_line_kwargs=None,
                  dos_area_fill=False, dos_area_kwargs=None,
                  pdosup_area_kwargs=None, pdosup_clr='r', pdosdn_area_kwargs=None, pdosdn_clr='b',
                  orientation='h',
                  label='', label_kwargs=None):
    """

    @param span: range of plot (e_min, e_max)
    @param specup: ndarray of spin-up levels
    @param specdn: ndarray of spin-down levels
    @param weights_up: weights for degeneracies of Kpoints or for intensities (in case of Raman/IR spectra), spin-up
    @param weights_dn: weights, spin-down
    @param e_fermi: Fermi energy, None for any non-electron spectra
    @param units: units of energies in specup and specdn
    @param smear_type: type of smearing: "Gaussian" or "Lorentzian"
    @param sigma: dispersion parameter for chosen type
    @param partial_up: weights for spin-up PDOS
    @param partial_dn: weights for spin-down PDOS
    @param shareax: base axes/figure/plt object where the graph will be plotted
    @param dos_line_kwargs: kwargs for DOS line
    @param dos_area_fill: Bool: filling
    @param dos_area_kwargs:
    @param pdosup_area_kwargs:
    @param pdosup_clr:
    @param pdosdn_area_kwargs:
    @param pdosdn_clr:
    @param label:
    @param label_kwargs:
    @param normalization:
    @param orientation:
    @return:
    """

    def o_plot(orient, base, *args, **kwargs):
        """
        Orientation-dependent plot
        @param orient: orientation, 'h' for horizontal, 'v' for vertical
        @param base: axes or figure
        @param args:
        @param kwargs:
        @return:
        """
        if orient == 'h':
            return base.plot(*args, **kwargs)
        elif orient == 'v':
            rev = list(args[::-1])
            x = rev.pop()
            y = rev.pop()
            newargs = [y, x] + rev[::-1]
            return base.plot(*newargs, **kwargs)

    def o_fill(orient, base, *args, **kwargs):
        if orient == 'h':
            return base.fill_between(*args, **kwargs)
        elif orient == 'v':
            return base.fill_betweenx(*args, **kwargs)

    if shareax is None:
        shareax = plt

    if e_fermi is None:
        e_fermi = 0
        plot_fermi = False
    else:
        plot_fermi = True

    if specdn is None:
        specdn = np.array([])

    if units not in ['ev', 'eV', 'EV']:
        if units in ['ha', 'Ha', 'HA']:
            specup_local = specup * 27.2114
            specdn_local = specdn * 27.2114
            e_fermi *= 27.2114
        elif units in ['Ry', 'ry', 'RY']:
            specup_local = specup * 13.6057
            specdn_local = specdn * 13.6057
            e_fermi *= 13.6057
        else:
            raise ValueError("units must be 'eV', 'Ha' or 'Ry' (case insensitive)")
    else:
        specup_local = specup.copy()
        specdn_local = specdn.copy()

    specup_local -= e_fermi
    specdn_local -= e_fermi

    eaxis, dos, pdos = smeared_dos(np.array(span), specup_local, spec_dn=specdn_local, smear_sigma=sigma,
                                   weights_up=weights_up, weights_dn=weights_dn,
                                   pdos_up_coef=partial_up, pdos_dn_coef=partial_dn, normalization=normalization,
                                   smear_type=smear_type)

    dos_handle = {'up': [], 'dn': []}
    pdos_handle = {'up': [], 'dn': []}

    # Spin up
    if dos_line_kwargs is None:
        dos_line_kwargs = {'color': 'k', 'linewidth': 0.2}
    dos_handle['up'] = o_plot(orientation, shareax, eaxis, dos['updos'], **dos_line_kwargs)
    if dos_area_fill:
        if dos_area_kwargs is None:
            dos_area_kwargs = {}
        dos_handle['up'] = o_fill(orientation, shareax, eaxis, 0, dos['updos'], **dos_area_kwargs)
    if partial_up is not None:
        if pdosup_area_kwargs is None:
            pdosup_area_kwargs = {'facecolor': pdosup_clr, 'alpha': 0.5}
        pdos_handle['up'] = o_fill(orientation, shareax, eaxis, 0, pdos['pdos_up'], **pdosup_area_kwargs)

    fermiline_args = (np.array([0, 0]), np.array([0, np.max(dos['updos'])]))

    # Spin down
    if len(dos['dndos']) > 0:
        fermiline_args = (np.array([0, 0]), np.array([-1, 1]) * np.max([dos['updos'].max(), dos['dndos'].max()]))
        if dos_line_kwargs is None:
            dos_line_kwargs = {'color': 'k', 'linewidth': 0.2}
        dos_handle['dn'] = o_plot(orientation, shareax, eaxis, -dos['dndos'], **dos_line_kwargs)
        if dos_area_fill:
            if dos_area_kwargs is None:
                dos_area_kwargs = {}
            dos_handle['dn'] = o_fill(orientation, shareax, eaxis, 0, -dos['dndos'], **dos_area_kwargs)
        if partial_dn is not None:
            if pdosdn_area_kwargs is None:
                pdosdn_area_kwargs = {'color': pdosdn_clr, 'alpha': 0.5}
            pdos_handle['dn'] = o_fill(orientation, shareax, eaxis, 0, -pdos['pdos_dn'], **pdosdn_area_kwargs)

    if plot_fermi:
        o_plot(orientation, shareax, *fermiline_args, '--', color='gray')

    if label_kwargs is None:
        label_kwargs = {'x': 0.1, 'y': 0.7, 'fontsize': 9}
    shareax.text(s=label, **label_kwargs)
    return shareax, dos_handle, pdos_handle


def draw_ipr(span, eaxup, iprup, eaxdn=None, eFermi=0, units='eV',
             currFigNo=None, totalAxs=1, currAxNo=1, shareax=None, fmtarg='r'):
    if eaxdn is None:
        eaxdn = np.array([])

    if units not in ['ev', 'eV', 'EV']:
        if units in ['ha', 'Ha', 'HA']:
            eaxup *= 27.2114
            eaxdn *= 27.2114
            eFermi *= 27.2114
        elif units in ['Ry', 'ry', 'RY']:
            eaxup *= 13.6057
            eaxdn *= 13.6057
            eFermi *= 13.6057
        else:
            raise ValueError("units must be 'eV', 'Ha' or 'Ry' (case insensitive)")
    else:
        eaxup = eaxup.copy()
        eaxdn = eaxdn.copy()
    eaxup -= eFermi
    eaxdn -= eFermi

    if not (plt.fignum_exists(currFigNo)):
        currFig = plt.figure(num=currFigNo, figsize=(3.375, totalAxs * 3.375 / 3))
    else:
        currFig = plt.figure(currFigNo)

    if shareax is None:
        shareax = currFig.add_subplot(totalAxs, 1, currAxNo)
        ipr_handle = plt.stem(eaxup[(eaxup > span[0]) & (eaxup < span[1])],
                              iprup[(eaxup > span[0]) & (eaxup < span[1])],
                              linefmt=fmtarg, basefmt=" ")
    else:
        ipr_handle = shareax.stem(eaxup[(eaxup > span[0]) & (eaxup < span[1])],
                                  iprup[(eaxup > span[0]) & (eaxup < span[1])],
                                  linefmt=fmtarg, markerfmt=" ", basefmt=" ", use_line_collection=True)
        shareax.plot(np.array([0, 0]), np.array([0, 1]), '--', color='gray')
    return shareax, ipr_handle


def interp_matrix(v_array, h_array, values_array, grid_res=0.02, interp_type='interp2d', interp_kwargs=None):
    v_array = np.array(v_array)
    h_array = np.array(h_array)
    values_array = np.array(values_array)
    fine_v_array = np.arange(v_array.min(), v_array.max() + grid_res, grid_res)
    fine_h_array = np.arange(h_array.min(), h_array.max() + grid_res, grid_res)
    if interp_type == 'RectBivariateSpline':
        if interp_kwargs is None:
            interp_kwargs = {'s': 0}
        interp_fun = interpolate.RectBivariateSpline(h_array, v_array, values_array.T, **interp_kwargs)
        return fine_v_array, fine_h_array, interp_fun(fine_h_array, fine_v_array).T
    else:
        if interp_kwargs is None:
            interp_kwargs = {'kind': 'linear'}
        interp_fun = interpolate.interp2d(h_array, v_array, values_array, **interp_kwargs)
        return fine_v_array, fine_h_array, interp_fun(fine_h_array, fine_v_array)


# 2D data plts
def prop_map(v_array, h_array, values_array, contours=None, cmap='jet',
             vcenter=None, vmin=None, vmax=None, low_cutoff=None, high_cutoff=None, bad_color_low='k', titlestr=None,
             x_ticks=None, y_ticks=None, plot_cbar=True, cbar_label=None, cbar_ticklabels=None,
             pcolormesh_kwargs=None, cbar_kwargs=None, contours_kwargs=None, shallowlevel=None, grid_res=0.02, **kwargs):
    """
    Draws 2D map of given data VALUES ARRAY(V_ARRAY, H_ARRAY)

    @type high_cutoff: object
    @param v_array:[1xN ndarray] y-axis array (vertical)
    @param h_array:[1xM ndarray] x-axis array (horizontal)
    @param values_array:[NxM ndarray] mtrx of values values_array[i,j], i in v_array, j in h_array
    @param grid_res: resolution of interpolation grid. If None is passed, no interpolation is performed
    @param contours: array of values for plotting contours, or just True for default values based on array values
    @param cmap: colormap
    @param vcenter: setting value to corresponding the "center" of colormap
    @param vmin: setting value for "lowest" color. Default is values_array.min()
    @param vmax: setting value for "highest" color. Default is values_array.max()
    @param low_cutoff: setting value, below which the values are considered "bad" (masked)
    @param bad_color_low: color for masked elements - the values below max_badvalue
    @param titlestr: title of figure
    @param cbar_label: title of colorbar
    @param cbar_ticklabels: tick labels of colorbar
    @param x_ticks: set of x-ticks
    @param y_ticks: set of y-ticks
    @param plot_cbar: colorbar flag, default is True
    @param pcolormesh_kwargs: keyword args for plt.pcolormesh
    @param cbar_kwargs: kwargs dictionary for plt.colorbar
    @param contours_kwargs: kwargs dictionary of plt.contour()
    @param shallowlevel: if not None and in the range of [low_cutoff...vcenter], then all values in this range
                         become equal to shallowlevel
    @param grid_res: resolution of interpolation grid. If None is passed, no interpolation is performed
    """

    # Processing kwargs:
    if pcolormesh_kwargs is None:
        pcolormesh_kwargs = {}
    if cbar_kwargs is None:
        cbar_kwargs = {}

    minval = np.nanmin(values_array)
    maxval = np.nanmax(values_array)

    if grid_res is None:
        fine_h_array, fine_v_array, interp_values = h_array, v_array, values_array
    else:
        assert 0 < grid_res <= 1.0
        fine_v_array, fine_h_array, interp_values = interp_matrix(v_array, h_array, values_array, grid_res=grid_res)

    if vmax is None:
        vmax = maxval
    if vmin is None:
        vmin = minval

    my_cmap = copy.copy(get_cmap(cmap))
    if (vcenter is not None) and (vcenter >= vmin):
        if low_cutoff is None and high_cutoff is None:
            interp_masked = interp_values
        elif high_cutoff is None:
            interp_masked = np.ma.array(interp_values, mask=interp_values <= low_cutoff)
            my_cmap.set_bad(bad_color_low)
        elif low_cutoff is None:
            interp_masked = np.ma.array(interp_values, mask=interp_values >= high_cutoff)
            my_cmap.set_bad(bad_color_low)
        else:
            interp_masked = np.ma.array(interp_values,
                                        mask=(np.logical_or(interp_values >= high_cutoff, interp_values <= low_cutoff)))
            my_cmap.set_bad(bad_color_low)

        norm = TwoSlopeNorm(vcenter=vcenter, vmax=vmax, vmin=vmin)
        norm.clip = True
        if shallowlevel is not None and (low_cutoff <= shallowlevel <= vcenter):
            interp_values[(interp_values < vcenter) & (interp_values > low_cutoff)] = shallowlevel
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
        interp_masked = interp_values

    mesh_obj = plt.pcolormesh(fine_h_array, fine_v_array, interp_masked, cmap=my_cmap, norm=norm, shading='auto',
                              **pcolormesh_kwargs)

    ax = plt.gca()

    if x_ticks is None:
        ax.set_xticks(h_array)
    else:
        ax.set_xticks(x_ticks)
    if y_ticks is None:
        ax.set_yticks(v_array)
    else:
        ax.set_yticks(y_ticks)

    ax.set_aspect(1)

    if titlestr is not None:
        plt.title(titlestr)

    plt.grid(linestyle='--', linewidth=0.5, zorder=10)

    if contours is not None:
        if contours_kwargs is None:
            contours_kwargs = {'colors': 'k', 'linewidths': 0.5}
        CS = plt.contour(fine_h_array, fine_v_array, interp_values, contours, zorder=1.2, **contours_kwargs)
        if 'clabel' in kwargs and kwargs['clabel'] is not None:
            if kwargs['clabel'] is True:
                kwargs['clabel'] = {'inline': 1, 'fontsize': 7, 'fmt': '%1.2g'}
            plt.clabel(CS, **kwargs['clabel'])

    if plot_cbar:
        if cbar_kwargs is None:
            cbar_kwargs = {}
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = plt.colorbar(mesh_obj, cax=cax, **cbar_kwargs)
        if cbar_label is not None:
            cbar.set_label(cbar_label)

    plt.sca(ax)
    return mesh_obj


def my_heatmap(data, row_labels, col_labels, cbar_ttl='', title='', badvalue=None, clr_ranges=None, **kwargs):
    if clr_ranges is None:
        norm = None
    else:
        norm = TwoSlopeNorm(**clr_ranges)
    fig, ax = plt.figure(), plt.gca()

    im, _ = heatmap(data[::-1], row_labels[::-1], col_labels, ax=ax, cmap='jet', cbarlabel=cbar_ttl, aspect=0.4,
                       cbar_kw={'orientation': 'horizontal', 'pad': 0.08, 'aspect': 40}, norm=norm)
    annotate_heatmap(im, badvalue=badvalue, valfmt="{x:6.2f}", size=7)
    plt.title(title)
    plt.margins(0.0001)
