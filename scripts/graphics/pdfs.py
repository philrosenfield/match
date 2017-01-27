from __future__ import print_function, absolute_import
import matplotlib.pyplot as plt
import numpy as np
from .graphics import corner_setup, fix_diagonal_axes, key2label

def add_quantiles(SSP, ax, attrs, uvalss=None, probs=None,
                  twod=False):
    """
    attrs must go [xattr, yattr]
    """
    attrs = np.atleast_1d(attrs)
    if uvalss is None:
        uvalss = [None] * len(attrs)
    else:
        uvalss = np.atleast_1d(uvalss)

    if probs is None:
        probs = [None] * len(attrs)
    else:
        probs = np.atleast_1d(probs)

    pltkw = {'color': 'k'}
    linefuncs = [ax.axvline, ax.axhline]
    pts = []
    for i, (attr, uvals, prob) in enumerate(zip(attrs, uvalss, probs)):
        if attr is None:
            continue
        gatr = '{:s}g'.format(attr)
        if not hasattr(SSP, gatr):
            g = SSP.fitgauss1D(attr, uvals, prob)
        else:
            g = SSP.__getattribute__(gatr)

        lines = [g.mean, g.mean + g.stddev / 2, g.mean - g.stddev / 2]
        lstys = ['-', '--', '--']

        if twod:
            lines = [g.mean + g.stddev / 2, g.mean - g.stddev / 2]
            lstys = ['--', '--']
            pts.append(g.mean)

        [linefuncs[i](l, ls=ls, **pltkw) for (l, ls) in zip(lines, lstys)]

    if twod:
        ax.plot(*pts, 'o', color='white', mec='k', mew=1)
    else:
        X = uvalss[0]
        prob = probs[0]
        # plot the Gaussian 10 steps beyond the calculated limits.
        dx = 10 * np.diff(X)[0]
        xx = np.linspace(X.min() - dx, X.max() + dx, 100)
        ax.plot(xx, g(xx), color='darkred')
    return ax


def pdf_plot(SSP, xattr, yattr=None, ax=None, sub=None, save=False,
             truth=None, cmap=None, plt_kw=None, X=None, prob=None,
             logp=True, gauss1D=False):
    """Plot -2 ln P vs marginalized attributes

    SSP : SSP class instance

    xattr, yattr : str
        column names to marginalize and plot

    ax : plt.Axes
        add plot to this axis

    sub : str
        if save, add this string to the filename

    save : bool
        save the plot to file with axes lables.

    truth : dict
        truth dictionary with attributes as keys and truth values as values
        overplot as hline and vline on axes.

    cmap : cm.cmap instance
        if yattr isn't None, call plt.pcolor with this cmap.

    plt_kw : dict
        if yattr is None, pass these kwargs to plt.plot

    Returns
    ax : plt.Axes
        new or updated axes instance
    """
    plt_kw = plt_kw or {'lw': 4, 'color': 'k'}
    cmap = cmap or plt.cm.viridis_r
    sub = sub or ''
    truth = truth or {}
    do_cbar = False
    pstr = ''
    if logp:
        pstr = '\ln\ '

    if ax is None:
        _, ax = plt.subplots()
        if yattr is not None:
            do_cbar = True

    if yattr is None:
        # plot type is marginal probability. Attribute vs -2 ln P
        if X is None and prob is None:
            X, prob = SSP.marginalize(xattr, log=logp)
            if not SSP.vdict[xattr]:
                return ax
            SSP.build_posterior(xattr, X, prob)
        l = ax.plot(X, prob, **plt_kw)
        if gauss1D:
            ax = add_quantiles(SSP, ax, xattr, uvalss=[X], probs=[prob])
        ax.set_xlim(X.min(), X.max())
        # yaxis max is the larger of 10% higher than the max val or current ylim.
        ymax = np.max([prob[prob>0].max() + (prob[prob>0].max() * 0.1),
                                             ax.get_ylim()[1]])
        ax.set_ylim(prob[prob>0].min(), ymax)
        if save:
            ptype = 'marginal'
            # ax.set_ylabel(key2label('fit'))
            ax.set_ylabel(key2label(pstr+'Probability'))
    else:
        # plot type is joint probability.
        # Attribute1 vs Attribute2 colored by fit
        [X, Y], prob = SSP.marginalize(xattr, yattr=yattr, log=logp)
        if not SSP.vdict[xattr] or not SSP.vdict[yattr]:
            return ax
        l = ax.pcolor(X, Y, prob, cmap=cmap)
        if gauss1D:
            add_quantiles(SSP, ax, [xattr, yattr], twod=True)
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())

        if do_cbar:
            cbar = plt.colorbar(l)
            # cbar.set_label(key2label('fit'))
            cbar.set_label(key2label(pstr+'Probability'))

        if save:
            ptype = '{}_joint'.format(yattr)
            ax.set_ylabel(key2label(yattr, gyr=SSP.gyr))

        if yattr in truth:
            ax.axhline(truth[yattr], color='k', lw=3)

    if xattr in truth:
        ax.axvline(truth[xattr], color='k', lw=3)

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


def pdf_plots(SSP, marginals=None, sub=None, twod=False, truth=None,
              text=None, cmap=None, fig=None, axs=None, frompost=False,
              logp=True, gauss1D=False):
    """Call pdf_plot for a list of xattr and yattr"""
    text = text or ''
    sub = sub or ''
    truth = truth or {}
    marginals = marginals or SSP._getmarginals()
    pstr = ''
    if logp:
        pstr = '\ln\ '
    if not hasattr(SSP, 'vdict'):
        SSP.check_grid()
    valid_margs = [k for (k, v) in list(SSP.vdict.items()) if v]
    ndim = len(marginals)
    if ndim != len(valid_margs):
        bad_margs = [m for m in marginals if m not in valid_margs]
        marginals = [m for m in marginals if m in valid_margs]
        print('Warning: {} does not vary and will be skipped.'.format(bad_margs))
        ndim = len(marginals)

    if twod:
        fig, axs = corner_setup(ndim)
        raxs = []
        for c, mx in enumerate(marginals):
            for r, my in enumerate(marginals):
                ax = axs[r, c]
                if r == c:
                    # diagonal
                    # my = 'fit'  # my is reset for ylabel call
                    my = pstr+'Probability'
                    raxs.append(SSP.pdf_plot(mx, ax=ax, truth=truth, logp=logp,
                                             gauss1D=gauss1D))
                else:
                    # off-diagonal
                    raxs.append(SSP.pdf_plot(mx, yattr=my, ax=ax, truth=truth,
                                             cmap=cmap, logp=logp,
                                             gauss1D=gauss1D))

                if c == 0:
                    # left most column
                    ax.set_ylabel(key2label(my))

                if r == ndim - 1:
                    # bottom row
                    ax.set_xlabel(key2label(mx))
        [ax.locator_params(axis='y', nbins=6) for ax in axs.ravel()]
        fix_diagonal_axes(raxs, ndim)
    else:
        raxs = []
        if fig is None and axs is None:
            fig, axs = plt.subplots(ncols=ndim, figsize=(15, 3))
        [ax.tick_params(left='off', labelleft='off', right='off', top='off')
         for ax in axs]
        X = None
        prob = None
        for i in marginals:
            if frompost:
                pattr = '{:s}prob'.format(i)
                X = SSP.data[i][np.isfinite(SSP.data[i])]
                prob = SSP.data[pattr][np.isfinite(SSP.data[pattr])]
            ax = axs[marginals.index(i)]
            ax = SSP.pdf_plot(i, truth=truth, ax=ax, X=X, prob=prob, logp=logp,
                              gauss1D=gauss1D)
            ax.set_xlabel(key2label(i, gyr=SSP.gyr))
            raxs.append(ax)

        if text:
            axs[-1].text(0.90, 0.10, '${}$'.format(text), ha='right',
                         transform=axs[-1].transAxes)
        fig.subplots_adjust(bottom=0.22, left=0.05)

        raxs[0].set_ylabel(key2label('Probability'))
    [ax.locator_params(axis='x', nbins=6) for ax in axs.ravel()]
    return fig, raxs
