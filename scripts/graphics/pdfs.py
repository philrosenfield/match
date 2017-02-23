from __future__ import print_function, absolute_import
import matplotlib.pyplot as plt
import numpy as np
from .graphics import (corner_setup, fix_diagonal_axes, key2label,
                       add_inner_title, square_aspect)


def add_quantiles(SSP, ax, attrs, uvalss=None, probs=None,
                  twod=False, gauss=False, interpolate=True):
    """
    Add some lines!

    Parameters
    ----------
    SSP : ssp.SSP instance
    ax : plt.axes instance
        add lines/points to this plot
    attrs: 1 or 2D str array
        (marginalized) SSP.data columns to plot
        NB: attrs must go [xattr, yattr] for 2d plotting
    uvalls: 1 or 2D float array
        unique values of attrs
    probs: 1 or 2D float array
        marginalized probabilities corresponding to attrs
    twod:
        add lines
    gauss: bool (False)
        fit or access gaussian fit of the attribute(s) and add lines
        for the mean and +/- stdev / 2
        or (default)
        add 16 and 84 percentile as well as maximum posterior probability.
        if twod, the mean or max post. prob will be plotted as point.

    Returns
    -------
    ax : plt.axes instance
    (may add attributes to SSP if quantiles or fittgauss1D had not been called)
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
    # why order matteres: xattr yattr will plot vline and hline
    linefuncs = [ax.axvline, ax.axhline]
    # mean or max post. prob to plot as points
    pts = []
    for i, (attr, uvals, prob) in enumerate(zip(attrs, uvalss, probs)):
        if attr is None:
            # for second loop 1D
            continue
        # look up value
        qatr = '{:s}q'.format(attr)
        if hasattr(SSP, qatr):
            q = SSP.__getattribute__(qatr)
        else:
            if gauss:
                # fit 1D Gaussian
                SSP.fitgauss1D(attr, uvals, prob)
                # go with quantiles (default 0.16, 0.84)
                # g = SSP.quantiles(attr, uvals, prob, maxp=True, k=1, ax=ax)
            q = SSP.quantiles(attr, uvals, prob, maxp=True, k=1,
                              interpolate=interpolate)

        try:
            lines = [g.mean, g.mean + g.stddev / 2, g.mean - g.stddev / 2]
        except:
            # if maxp=False when SSP.quantiles called
            # this will raise a value error because g will be length 2.
            lines = [q[2], q[0], q[1]]
        lstys = ['-', '--', '--']

        if twod:
            if gauss:
                lines = [g.mean + g.stddev / 2, g.mean - g.stddev / 2]
                pts.append(g.mean)
            else:
                lines = [q[0], q[1]]
                pts.append(q[2])
            lstys = ['--', '--']

        # plot.
        [linefuncs[i](l, ls=ls, **pltkw) for (l, ls) in zip(lines, lstys)]

    if twod:
        # plot mean or max post prob
        ax.plot(pts[0], pts[1], 'o', color='white', mec='k', mew=1)

    else:
        if gauss:
            # over plot Gaussian fit
            X = uvalss[0]
            prob = probs[0]
            # plot the Gaussian 10 steps beyond the calculated limits.
            dx = 10 * np.diff(X)[0]
            xx = np.linspace(X.min() - dx, X.max() + dx, 100)
            ax.plot(xx, g(xx), color='darkred')

    return ax, g


def pdf_plot(SSP, xattr, yattr=None, ax=None, sub=None, save=False,
             truth=None, cmap=None, plt_kw=None, X=None, prob=None,
             logp=True, quantile=False, gauss1D=False, plotfit=False,
             interpolateq=True):
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
    plt_kw = plt_kw or {}
    def_kw = {'lw': 4, 'color': 'k'}
    def_kw.update(plt_kw)
    plt_kw = def_kw

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
            if not SSP.vdict[xattr]:
                return ax
            X, prob = SSP.marginalize(xattr, log=logp)
            SSP.build_posterior(xattr, X, prob)

        if gauss1D:
            gx = SSP.fitgauss1D(xattr, X, prob)

        l = ax.plot(X, prob, **plt_kw)

        if quantile:
            ax = add_quantiles(SSP, ax, xattr, uvalss=[X], probs=[prob],
                               gauss=gauss1D, interpolate=interpolateq)

        ax.set_xlim(X.min(), X.max())
        # yaxis max is the larger of 10% higher than the max val or current ylim.
        ymax = np.max([prob.max() + (prob.max() * 0.1), ax.get_ylim()[1]])
        ax.set_ylim(prob.min(), ymax)

        if gauss1D and plotfit:
            # over plot Gaussian fit
            # plot the Gaussian 10 steps beyond the calculated limits.
            dx = 10 * np.diff(X)[0]
            xx = np.linspace(X.min() - dx, X.max() + dx, 100)
            ax.plot(xx, gx(xx), color='darkred')

        if save:
            ptype = 'marginal'
            # ax.set_ylabel(key2label('fit'))
            ax.set_ylabel(key2label(pstr+'Probability'))
    else:

        # plot type is joint probability.
        # Attribute1 vs Attribute2 colored by fit
        if not SSP.vdict[xattr] or not SSP.vdict[yattr]:
            return ax

        [X, Y], prob = SSP.marginalize(xattr, yattr=yattr, log=logp)

        if gauss1D:
            gy = SSP.fitgauss1D(yattr, Y, prob)

        l = ax.pcolor(X, Y, prob, cmap=cmap)
        # use imshow instead of pcolor, has strange aspect ratio...
        # ux = SSP.__getattribute__('u{0:s}'.format(xattr))
        # uy = SSP.__getattribute__('u{0:s}'.format(yattr))
        # l = ax.imshow(prob.T, extent=[ux[0], ux[-1], uy[0], uy[-1]], cmap=cmap,
        #               origin='lower')
        # ax = square_aspect(ax)

        if quantile:
            add_quantiles(SSP, ax, [xattr, yattr], twod=True, gauss=gauss1D,
                          interpolate=interpolateq)

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
            ax.axhline(truth[yattr], color='darkred', lw=3, zorder=0)

    if xattr in truth:
        ax.axvline(truth[xattr], color='darkred', lw=3, zorder=0)

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
              logp=True, gauss1D=False, quantile=True, interpolateq=True):
    """Call pdf_plot for a list of xattr and yattr"""
    text = text or ''
    sub = sub or ''
    truth = truth or {}
    marginals = marginals or SSP._getmarginals()
    pstr = ''
    plkw = {'logp': logp, 'gauss1D': gauss1D, 'quantile': quantile,
            'interpolateq': interpolateq}

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

    raxs = []
    if twod:
        fig, axs = corner_setup(ndim)
        for c, mx in enumerate(marginals):
            for r, my in enumerate(marginals):
                ax = axs[r, c]
                if r == c:
                    # diagonal
                    # my = 'fit'  # my is reset for ylabel call
                    my = pstr+'Probability'
                    raxs.append(SSP.pdf_plot(mx, ax=ax, truth=truth, **plkw))
                else:
                    # off-diagonal
                    raxs.append(SSP.pdf_plot(mx, yattr=my, ax=ax, truth=truth,
                                             cmap=cmap, **plkw))

                if c == 0:
                    # left most column
                    ax.set_ylabel(key2label(my, gyr=SSP.gyr))

                if r == ndim - 1:
                    # bottom row
                    ax.set_xlabel(key2label(mx, gyr=SSP.gyr))
        [ax.locator_params(axis='y', nbins=6) for ax in axs.ravel()]
        fix_diagonal_axes(raxs, ndim)
    else:
        if fig is None and axs is None:
            fig, axs = plt.subplots(ncols=ndim, figsize=(ndim * 3., ndim * 0.6))
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
            ax = SSP.pdf_plot(i, truth=truth, ax=ax, X=X, prob=prob, **plkw)
            ax.set_xlabel(key2label(i, gyr=SSP.gyr))
            raxs.append(ax)

        if text:
            add_inner_title(raxs[-1], '${}$'.format(text), 3, size=None)
        fig.subplots_adjust(bottom=0.22, left=0.05)
        raxs[0].set_ylabel(key2label('Probability'))
    [ax.locator_params(axis='x', nbins=5) for ax in axs.ravel()]
    return fig, raxs
