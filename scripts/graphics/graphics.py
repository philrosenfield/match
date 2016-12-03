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
               'dmod': r'$\mu$',
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
        seg = left(cmap1._segmentdata[color]) + right(cmap2._segmentdata[color])
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


def match_diagnostic(param, phot, fake=None, save=True):
    """
    make two panel cmd figure (ymag = mag1, mag2)
    from match photometery and cmd parameter space from the match param drawn
    """
    moff = 0.1
    off = 0.1
    mag1, mag2 = np.loadtxt(phot, unpack=True)
    color = mag1 - mag2

    params = read_calcsfh_param(param)

    cmin, cmax = params['vimin'], params['vimax']
    m1min, m1max = params['vmin'], params['vmax']
    m2min, m2max = params['imin'], params['imax']

    filters = [params['v'], params['i']]

    verts = [np.array([[cmin, m1min], [cmin, m1max], [cmax, m1max],
                       [cmax, m1min], [cmin, m1min]]),
             np.array([[cmin, m2min], [cmin, m2max], [cmax, m2max],
                       [cmax, m2min], [cmin, m2min]])]

    magcuts, = np.nonzero((mag1 < m1max) & (mag1 > m1min) &
                          (mag2 < m2max) & (mag2 > m2min))

    _, axs = plt.subplots(ncols=2, sharex=True, figsize=(12, 6))

    for i, ymag in enumerate([mag1, mag2]):
        axs[i].plot(color, ymag, '.', label='all phot')
        axs[i].plot(color[magcuts], ymag[magcuts], '.', label='mag limits')
        axs[i].plot(verts[i][:, 0], verts[i][:, 1], label='param limits')
        axs[i].set_ylabel(r'${}$'.format(filters[i]))
        axs[i].set_xlabel(r'${}-{}$'.format(*filters))
        axs[i].invert_yaxis()

    if fake is not None:
        mag1in, mag2in, _, _ = np.loadtxt(fake, unpack=True)
        colin = mag1in - mag2in
        for i, ymag in enumerate([mag1in, mag2in]):
            fverts = np.array([[colin.min(), ymag.min()],
                               [colin.min(), ymag.max()],
                               [colin.max(), ymag.max()],
                               [colin.max(), ymag.min()],
                               [colin.min(), ymag.min()]])
            axs[i].plot(fverts[:, 0], fverts[:, 1], label='fake limits')
            plt.legend(loc='best')

    if save:
        plt.savefig(param + EXT)
        print('wrote', param + EXT)
        plt.close()
    return axs


def add_inner_title(ax, title, loc, size=None):
    '''add a label to an axes as if it were a legend (loc must be 1-11)'''
    from matplotlib.patheffects import withStroke
    from matplotlib.offsetbox import AnchoredText

    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    anct = AnchoredText(title, loc=loc, prop=size, pad=0., borderpad=0.5,
                        frameon=False)
    ax.add_artist(anct)
    anct.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    anct.patch.set_alpha(0.5)
    return anct


def zeroed_cmap(hess, cmap1=plt.cm.Reds_r, cmap2=plt.cm.Blues, dfrac=0.05,
                transparent=False):
    """make a diverging color map with white set to 0.0"""
    fhess = hess[np.isfinite(hess)]
    minfhess = np.abs(np.min(fhess))
    # stitch to make a diverging color map with white set to 0.0
    frac = minfhess / (minfhess + np.abs(np.max(fhess)))
    return stitch_cmap(cmap1, cmap2, stitch_frac=frac, dfrac=dfrac,
                       transparent=transparent)


def setup_imgrid(figsize=[12, 3], nrows=1, ncols=4):
    from mpl_toolkits.axes_grid1 import ImageGrid
    igkw = {'nrows_ncols': (nrows, ncols),
            'axes_pad': .7,
            'label_mode': "all",
            'share_all': True,
            'cbar_location': "top",
            'cbar_mode': "each",
            'cbar_size': "7%",
            'cbar_pad': "2%"}

    fig = plt.figure(figsize=figsize)
    grid = ImageGrid(fig, 111, **igkw)
    return grid


def match_plot(hesslist, extent, labels=None, twobytwo=True,
               xlabel=None, ylabel=None, cmap=None, logcounts=False):
    '''plot four hess diagrams with indivdual color bars using ImageGrid'''
    from mpl_toolkits.axes_grid1 import ImageGrid

    if twobytwo:
        figsize=[9, 9]
        nrows = 2
        ncols = 2
    else:
        figsize=[12, 3]
        nrows = 1
        ncols = 4

    grid = setup_imgrid(figsize=figsize, nrows=nrows, ncols=ncols)

    for i, (ax, hess) in enumerate(zip(grid, hesslist)):
        if cmap is None:
            if i > 1:
                # bottom row: diff, sig
                colors = zeroed_cmap(hess)
            else:
                # first row: data, model. White will be on the left of color bar
                if i == 0:
                    colors = plt.cm.Blues
                if i == 1:
                    colors = plt.cm.Reds
                # colors = plt.cm.get_cmap('binary', 11)
        else:
            colors = cmap
        if logcounts:
            hess = np.log10(hess)
        img = ax.imshow(hess, origin='upper', extent=extent,
                        interpolation="nearest", cmap=colors)
        ax.cax.colorbar(img)
        if labels is not None:
            _ = add_inner_title(ax, labels[i], loc=1)

        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax = square_aspect(ax)

    if xlabel is not None:
        ind = 0
        if twobytwo:
            ind = 1
        _ = [ax.set_xlabel(xlabel) for ax in grid.axes_row[ind]]
        grid.axes_all[0].xaxis.label.set_visible(True)
    if ylabel is not None:
        _ = [ax.set_ylabel(ylabel) for ax in grid.axes_column[0]]

    return grid
