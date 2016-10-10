"""General plotting functions"""
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np

from .config import EXT
from .fileio import read_calcsfh_param
from .utils import center_grid

try:
    import seaborn
    seaborn.set()
except ImportError:
    pass


def pdf_plot(SSP, xattr, yattr=None, ax=None, sub=None, save=False,
             truth=None, cmap=None, plt_kw=None):
    """Plot prob vs marginalized attributes

    SSP : SSP class instance

    xattr, yattr : string, string
        column names to marginalize and plot
    ax : plt.Axes
    sub : string
        if save, add this string to the filename

    save : bool
        save the plot to file with axes lables.
    truth : dict(attr : val)
        truth dictionary to overplot as hline and vline on axes.
    cmap : if yattr isn't None, call pcolor with this cmap.
    plt_kw : if yattr is None, pass these kwargs to plt.plot

    Returns
    ax
    """
    def validate_arr(atr):
        v = SSP.data[atr]
        if len(np.unique(v)) == 1:
            print('{} not varied.'.format(atr))
            plt.close()
            v = None
        return v

    x = validate_arr(xattr)
    if x is None:
        return
    if yattr is not None:
        y = validate_arr(yattr)
        if y is None:
            return
    z = SSP.data['fit']

    plt_kw = plt_kw or {'lw': 4, 'color': 'k'}
    cmap = cmap or plt.cm.viridis_r
    sub = sub or ''
    truth = truth or {}
    do_cbar = False

    if ax is None:
        _, ax = plt.subplots()
        if yattr is not None:
            do_cbar = True

    if yattr is None:
        # plot type is marginal probability. Attribute vs -2 ln P
        # array bins
        l = ax.plot(*SSP.marginalize(xattr), **plt_kw)
        if save:
            ptype = 'marginal'
            ax.set_ylabel(key2label('fit'))
    else:
        # plot type is joint probability.
        # Attribute1 vs Attribute2 colored by fit
        [X, Y], prob = SSP.marginalize(xattr, yattr=yattr)
        l = ax.pcolor(X, Y, prob, cmap=cmap)
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())

        if do_cbar:
            cbar = plt.colorbar(l)
            cbar.set_label(key2label('fit'))

        if save:
            ptype = '{}_joint'.format(yattr)
            ax.set_ylabel(key2label(yattr, gyr=SSP.gyr))

        if yattr in truth:
            ax.axhline(truth[yattr], color='darkred', lw=3)

    if xattr in truth:
        ax.axvline(truth[attr], color='darkred', lw=3)

    if save:
        ax.set_xlabel(key2label(xattr, gyr=SSP.gyr))
        # add subdirectory to filename
        if len(sub) > 0:
            sub = '_' + sub
        outfmt = '{}_{}{}_{}{}'
        outname = outfmt.format(SSP.name.replace('.csv', ''),
                                xattr, sub, ptype, EXT)
        plt.savefig(outname, bbox_inches='tight')
        print('wrote {}'.format(outname))
        plt.close()
    return ax


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


def pdf_plots(SSP, marginals=None, sub=None, twod=False, truth=None,
              text=None):
    """Call pdf_plot for a list of xattr and yattr"""
    text = text or ''
    sub = sub or ''
    truth = truth or {}
    marginals = marginals or SSP._getmarginals()

    ndim = len(marginals)
    if twod:
        fig, axs = corner_setup(ndim)
        raxs = []
        for c, mx in enumerate(marginals):
            for r, my in enumerate(marginals):
                ax = axs[r, c]
                if r == c:
                    # diagonal
                    my = 'fit'  # my is reset for ylabel call
                    raxs.append(SSP.pdf_plot(mx, ax=ax, truth=truth))
                else:
                    # off-diagonal
                    raxs.append(SSP.pdf_plot(mx, yattr=my, ax=ax, truth=truth))

                if c == 0:
                    # left most column
                    ax.set_ylabel(key2label(my))

                if r == ndim - 1:
                    # bottom row
                    ax.set_xlabel(key2label(mx))
        [ax.locator_params(axis='y', nbins=6) for ax in axs.ravel()]
    else:
        fig, axs = plt.subplots(ncols=ndim, figsize=(15, 3))
        [ax.tick_params(left='off', labelleft='off') for ax in axs[1:]]
        [ax.tick_params(right='off', labelright='off') for ax in axs[:-1]]
        [ax.tick_params(top='off') for ax in axs]
        axs[-1].tick_params(labelright='on')
        raxs = [SSP.pdf_plot(i, truth=truth, ax=axs[marginals.index(i)])
                for i in marginals]
        # need to add xlabels!
        if text:
            axs[-1].text(0.10, 0.90, '${}$'.format(text),
                         transform=axs[-1].transAxes)
        fig.subplots_adjust(bottom=0.2, left=0.05)

    [ax.locator_params(axis='x', nbins=6) for ax in axs.ravel()]

    return fig, raxs


def key2label(string, gyr=False):
    """latex labels for different strings"""
    def_fmt = r'${}$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu$',
               'lage': r'$\log\ \rm{Age\ (yr)}$',
               'logZ': r'$\log\ \rm{Z}$',
               'fit': r'$-2 \ln\ \rm{P}$',
               'ov': r'$\Lambda_c$',
               'chi2': r'$\chi^2$',
               'bf': r'$\rm{Binary\ Fraction}$',
               'dav': r'$dA_V$',
               'trueov': r'$\Lambda_c\ \rm{In}$'}
    if string not in convert.keys():
        convstr = def_fmt.format(string)
    else:
        convstr = convert[string]

    if gyr:
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


def stitch_cmap(cmap1, cmap2, stitch_frac=0.5, dfrac=0.001):
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
    return LinearSegmentedColormap('_'.join((cmap1.name, cmap2.name)),
                                   dict([(key, new_seg(key)) for key in rgb]),
                                   1024)


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
    '''
    add a title to an ax inside to a location loc, which follows plt.legends
    locations.
    '''
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


def match_plot(hesslist, extent, labels=None, imagegrid_kw=None,
               xlabel=None, ylabel=None):
    '''plot four hess diagrams with indivdual color bars using ImageGrid'''
    from mpl_toolkits.axes_grid1 import ImageGrid

    imagegrid_kw = imagegrid_kw or {}

    defaults = {'nrows_ncols': (2, 2),
                'axes_pad': .7,
                'label_mode': "all",
                'share_all': True,
                'cbar_location': "top",
                'cbar_mode': "each",
                'cbar_size': "7%",
                'cbar_pad': "2%"}

    defaults.update(imagegrid_kw)

    fig = plt.figure(figsize=(9, 9))
    grid = ImageGrid(fig, 111, **defaults)

    for i, (ax, hess) in enumerate(zip(grid, hesslist)):
        if i > 1:
            # bottom row: diff, sig
            fhess = hess[np.isfinite(hess)]
            minfhess = np.abs(np.min(fhess))
            # stitch to make a diverging color map with white set to 0.0
            frac = minfhess / (minfhess + np.abs(np.max(fhess)))
            colors = stitch_cmap(plt.cm.Reds_r, plt.cm.Blues, stitch_frac=frac,
                                 dfrac=0.05)
        else:
            # first row: data, model. White will be on the left of color bar
            if i == 0:
                colors = plt.cm.Reds
            if i == 1:
                colors = plt.cm.Blues
            # colors = plt.cm.get_cmap('binary', 11)

        aspect = abs((extent[1] - extent[0]) / (extent[3] - extent[2]))
        img = ax.imshow(hess, origin='upper', extent=extent, aspect=aspect,
                        interpolation="nearest", cmap=colors)
        ax.cax.colorbar(img)
        if labels is not None:
            _ = add_inner_title(ax, labels[i], loc=1)

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    if xlabel is not None:
        _ = [ax.set_xlabel(xlabel) for ax in grid.axes_row[1]]
        grid.axes_all[0].xaxis.label.set_visible(True)
    if ylabel is not None:
        _ = [ax.set_ylabel(ylabel) for ax in grid.axes_column[0]]

    return grid
