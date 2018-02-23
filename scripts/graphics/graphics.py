"""General plotting functions"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np

from ..config import EXT
from ..fileio import read_calcsfh_param
from ..utils import center_grid

try:
    import seaborn
    seaborn.set()
except ImportError:
    pass

import matplotlib as mpl


def corner_setup(ndim):
    """
    Setup a corner plot
    Only keep ticklabels on outter left and bottom plots, only set visible
    lower triangle of the plot grid.

    ndim: int
        number of rows (= number of columns)
    Returns
    fig, axs : output of plt.subplots
    """
    fig, axs = plt.subplots(nrows=ndim, ncols=ndim,
                            figsize=(ndim * 2, ndim * 2))
    # Turn off ticks on diagonal
    offs = ['top', 'right', 'left', 'labelleft']
    lnp_ticks = dict(zip(offs, ['off'] * len(offs)))
    [axs[k, k].tick_params(**lnp_ticks) for k in range(ndim)]

    # Turn off bottom tick labels on all but bottom axes
    [ax.tick_params(labelbottom='off') for ax in axs[:ndim-1, :].ravel()]

    # Turn off left tick labels on all but left axes
    [ax.tick_params(labelleft='off') for ax in axs[:, 1:].ravel()]

    # Turn off upper triangle axes
    [[ax.set_visible(False) for ax in axs[k, k+1:]] for k in range(ndim)]

    if ndim > 2:
        fig.subplots_adjust(left=0.1, bottom=0.1, top=0.95, right=0.95)
    return fig, axs


def fix_diagonal_axes(raxs, ndim):
    """
    set diagonal xaxis limits to that of the off diagonal xaxis limits.

    With corner_meshgrid, the axes limits for the 2d plots will be slightly
    different than the 1d plots.
    """
    nplots = len(raxs)
    idiag = [i * (ndim + 1) for i in range(nplots // (ndim + 1))]
    [raxs[i].set_xlim(raxs[i+1].get_xlim()) for i in idiag]
    [raxs[nplots - 1].set_xlim(raxs[ndim - 1].get_ylim()) for i in idiag]
    return


def key2label(string, gyr=False):
    """latex labels for different strings"""
    def_fmt = r'$\rm{{{}}}$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu_0$',
               'lage': r'$\log\ \rm{Age\ (yr)}$',
               'logZ': r'$\log\ \rm{Z}$',
               'fit': r'$-2 \ln\ \rm{P}$',
               'ov': r'$\Lambda_c$',
               'IMF': r'$\Gamma$',
               'chi2': r'$\chi^2$',
               'bf': r'$\rm{Binary\ Fraction}$',
               'dav': r'$dA_V$',
               'trueov': r'$\Lambda_c\ \rm{In}$'}
    if string not in convert.keys():
        convstr = def_fmt.format(string)
    else:
        convstr = convert[string]

    if gyr and 'age' in string:
        convstr = convstr.replace('yr', 'Gyr').replace(r'\log\ ', '')
    return convstr


def put_a_line_on_it(ax, val, consty=False, color='black',
                     ls='--', lw=2, yfilt='I', steps=20):
    """
    if consty is True: plots a constant y value across ax.xlims().
    if consty is False: plots a constant x on a plot of y vs x-y
    """
    (xmin, xmax) = ax.get_xlim()
    (ymin, ymax) = ax.get_ylim()

    xarr = np.linspace(xmin, xmax, steps)
    yarr = np.linspace(ymin, ymax, steps)
    if consty is True:
        # just a contsant y value over the plot range of x.
        # e.g., a LF
        ax.axhline(val, color=color, lw=lw)
    if consty is False:
        # a plot of y vs x-y and we want to mark
        # where a constant value of x is
        # e.g, f814w vs f555-f814; val is f555
        new_xarr = val - yarr
        # e.g, f555w vs f555-f814; val is f814
        if yfilt.upper() != 'I':
            yarr = xarr + val
            new_xarr = xarr
        ax.plot(new_xarr, yarr, ls, color=color, lw=lw)
    return new_xarr, yarr


def square_aspect(ax):
    """
    Set aspect ratio of a plot to have equal sized x and y axes length.
    (I.e, a square figure)
    """
    ax.set_aspect(np.abs((np.diff(ax.get_xlim()) / (np.diff(ax.get_ylim())))))
    # ax.set_aspect(np.abs((extent[1] - extent[0])/(extent[3] - extent[2])))
    return ax


def stitch_cmap(cmap1, cmap2, stitch_frac=0.5, dfrac=0.001, transparent=False):
    '''
    Code adapted from Dr. Adrienne Stilp
    Stitch two color maps together:
        cmap1 from 0 and stitch_frac
        and
        cmap2 from stitch_frac to 1
        with dfrac spacing inbetween

    ex: stitch black to white to white to red:
    stitched = stitch_cmap(cm.Greys_r, cm.Reds, stitch_frac=0.525, dfrac=0.05)
    '''
    from matplotlib.colors import LinearSegmentedColormap

    def left(seg):
        """left color segment"""
        return [(i * (stitch_frac - dfrac), j, k) for i, j, k in seg]

    def right(seg):
        """right color segment"""
        frac = stitch_frac + dfrac
        return [(i * (1 - frac) + frac, j, k) for i, j, k in seg]

    def new_seg(color):
        """combine left and right segments"""
        seg = left(cmap1._segmentdata[color]) + \
            right(cmap2._segmentdata[color])
        return seg

    rgb = ['blue', 'red', 'green']
    cname = '_'.join((cmap1.name, cmap2.name))
    cdict = dict([(key, new_seg(key)) for key in rgb])
    ncmap = LinearSegmentedColormap(cname, cdict, 1024)

    if transparent:
        # set the middle value to zero transparency.
        # it's probably better if you set alpha on the call using the
        # color map rather than change a single value.
        ncmap._init()
        ind = np.max([np.argmax(ncmap._lut.T[i])
                      for i in range(len(ncmap._lut.T)-1)])
        ncmap._lut[ind][-1] = 0
    return ncmap


def add_inner_title(ax, title, loc, size=None):
    '''
    add a label to an axes as if it were a legend (loc must be 1-11)
        'upper right'     1
        'upper left'      2
        'lower left'      3
        'lower right'     4
        'right'           5
        'center left'     6
        'center right'    7
        'lower center'    8
        'upper center'    9
        'center'          10
    '''
    from matplotlib.patheffects import withStroke
    from matplotlib.offsetbox import AnchoredText
    mpl.rc('text', usetex=True)

    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    anct = AnchoredText(title, loc=loc, prop=size, pad=0., borderpad=0.5,
                        frameon=False)
    ax.add_artist(anct)
    anct.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    anct.patch.set_alpha(0.5)
    mpl.rc('text', usetex=False)
    return anct


def zeroed_cmap(hess, cmap1=plt.cm.Greys_r, cmap2=plt.cm.Blues, dfrac=0.001,
                transparent=False):
    """make a diverging color map with white set to 0.0"""
    fhess = hess[np.isfinite(hess)]
    minfhess = np.abs(np.min(fhess))
    # stitch to make a diverging color map with white set to 0.0
    frac = minfhess / (minfhess + np.abs(np.max(fhess)))
    return stitch_cmap(cmap1, cmap2, stitch_frac=frac, dfrac=dfrac,
                       transparent=transparent)
