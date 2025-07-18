import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


# Routines for spectra plots
def smeared_spectrum(x_array, y_array=None, sigma=None, x_grid=None, normalized_to_n=False, sm_type='Gaussian'):
    """
    Given Xi, Yi, sigma and x_grid = {x} returns array
    z = sum_i{Yi*exp(-(x-Xi)^2/(2*sigma^2))}

    @param x_array: Xi
    @param y_array: Yi
    @param sigma: 
    @param x_grid: {x}
    @param normalized_to_n:
    @param sm_type:
    @return: z 
    """

    def smear_fun(x, x_i, weight_i, dispersion, f):
        x = x.reshape(-1)
        x_i = x_i.reshape(-1, 1)
        weight_i = weight_i.reshape(-1, 1)
        return np.sum(weight_i * f(x - x_i, dispersion), 0)

    assert sm_type[0].casefold() in ('l', 'g'), " Wrong smearing type"

    if sm_type[0].casefold() == 'l':
        def f_smear(x, s):
            return (1 / np.pi) * s / (x ** 2 + s ** 2)
    else:
        def f_smear(x, s):
            return np.exp(- x ** 2 / (2 * s ** 2))

    if y_array is None:
        y_array = np.ones(x_array.shape())
    if sigma is None:
        dx = abs(x_array.reshape(x_array.size, 1) - x_array)
        sigma = max(min(dx[dx > 0]) / 2, x_grid[1] - x_grid[0])
    if x_grid is None:
        x_grid = np.arange(np.min(x_array[:]) - sigma, np.max(x_array[:]) + sigma,
                           (np.max(x_array[:]) - np.min(x_array[:])) / 9999)
    smeared = smear_fun(x_grid, x_array, y_array, sigma, f_smear)
    if normalized_to_n:
        smeared /= sigma * np.sqrt(2 * np.pi)
    return smeared, x_grid


# Routines for heatmap tables
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_xlabel(cbarlabel, rotation=0, va="bottom", labelpad=10)  # before: rotation=-90

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on bottom.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), ha="center")  # before: rotation=-30, rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     badvalue=None,
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = 0.7

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    #       1. Preparing an array of colors:
    mindata = im.norm.vmin
    maxdata = im.norm.vmax

    color_scale = (256 * ((data - mindata) / (maxdata - mindata)))
    color_scale[color_scale < 0] = 0
    color_scale[color_scale > 255] = 256
    color_scale = color_scale.astype(int)

    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j] == badvalue:
                continue
            brightness = np.dot(im.cmap(color_scale[i, j])[0:3], np.array([0.299, 0.587, 0.114]))
            kw.update(color=textcolors[brightness < threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        v_ext = np.max([np.abs(self.vmin), np.abs(self.vmax)])
        x, y = [-v_ext, self.midpoint, v_ext], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
