"""Stats and visualization of calcsfh -ssp runs"""
from __future__ import print_function
import argparse
import itertools
import os
import sys

try:
    import seaborn
    seaborn.set()
except ImportError:
    pass

import pandas as pd
import matplotlib.pylab as plt
import numpy as np

from .config import EXT, match_base
from .fileio import filename_data, add_filename_info_to_file
from .utils import strip_header
from .cmd import call_pgcmd_byfit

__all__ = ['SSP']


def combine_files(fnames, outfile='combined_files.csv', best=False,
                  stats=False):
    """add files together including columns based on params in filename"""
    all_data = pd.DataFrame()
    for fname in fnames:
        dframe = add_filename_info_to_file(fname, best=best, stats=stats)
        all_data = all_data.append(dframe, ignore_index=True)

    all_data.to_csv(outfile, index=False)
    return outfile


def sspcombine(fname, dry_run=True, outfile=None):
    """
    call bin/sspcombine
    fname : string
        input filename
    outfile : None default: [fname].stats
        output filename
    dry_run : bool
        do not actually run sspcombine
    return string command to run sspcombine
    """
    sspname = strip_header(fname)
    if outfile is None:
        outfile = '> {}.stats'.format(fname)
    else:
        outfile = '>> {}'.format(outfile)
    cmd = '{} {} {}'.format(os.path.join(match_base, 'bin/sspcombine'),
                            sspname, outfile)
    if not dry_run:
        print('excecuting: {}'.format(cmd))
        os.system(cmd)
    return cmd


class SSP(object):
    """
    Class for calcsfh -ssp outputs
    """
    def __init__(self, filename, data=None, filterby=None):
        """
        filenames are the calcsfh -ssp terminal or console output.
        They do not need to be stripped of their header or footer or
        be concatenated as is typical in MATCH useage.
        """
        self.base, self.name = os.path.split(filename)
        if data is None:
            if filename.endswith('.csv'):
                data = pd.read_csv(filename)
            else:
                dat = np.genfromtxt(filename, skip_header=10, skip_footer=2,
                                    names=['Av', 'IMF', 'dmod', 'lage', 'logZ',
                                           'fit', 'sfr'])
                data = pd.DataFrame(dat)

        if filterby is not None:
            for key, val in filterby.items():
                data = data.loc[data[key] == val].copy(deep=True)

        self.data = data
        self.ibest = np.argmin(self.data['fit'])

        self.absprob = \
            np.exp(0.5 * (self.data['fit'].min() - self.data['fit']))

    def _getmarginals(self):
        """get the values to marginalize over that exist in the data"""
        # marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav',
        # 'ov', 'bf'])
        marg = np.array([k for k in self.data.keys() if k != 'fit'])
        inds = [i for i, m in enumerate(marg) if self._haskey(m)]
        return marg[inds]

    def _haskey(self, key):
        """test if the key requested is available"""
        ecode = True
        if key not in self.data.keys():
            ecode = False
        return ecode

    def marginalize(self, attr, attr2=None):
        """Find the best fit for each unique value of attr"""
        assert self._haskey(attr), '{} not found'.format(attr)
        x = self.data[attr]
        unq_x = np.unique(x)

        if attr2 is not None:
            assert self._haskey(attr2), '{} not found'.format(attr)
            y = self.data[attr2]
            unq_y = np.unique(y)

            size = len(unq_x) * len(unq_y)
            prob = np.zeros(size)
            vals_x = np.zeros(size)
            vals_y = np.zeros(size)

            k = 0
            for ix in unq_x:
                for iy in unq_y:
                    inds, = np.nonzero((x == ix) & (y == iy))
                    prob[k] = np.sum(self.absprob.iloc[inds])
                    vals_x[k] = ix
                    vals_y[k] = iy
                    k += 1
            vals = [vals_x, vals_y]
        else:
            # compute linear probabilites
            # sum over probabilites for each unique grid value
            prob = np.zeros(len(unq_x))
            for i, ix in enumerate(unq_x):
                prob[i] = np.sum(self.absprob[x == ix])
            vals = unq_x

        prob /= prob.sum()
        return vals, prob, np.log(prob)

    def pdf_plot(self, *args, **kwargs):
        return pdf_plot(self, *args, **kwargs)

    def pdf_plots(self, *args, **kwargs):
        return pdf_plots(self, *args, **kwargs)


def pdf_plot(SSP, attr, attr2=None, ax=None, sub='', save=False,
             truth=None):
    """Plot prob vs marginalized attr"""
    def center_grid(a):
        """
        uniquify and shift a uniform array half a bin maintaining its size
        """
        x = np.unique(a)
        dx = np.diff(x)[0]
        x = np.append(x, x[-1] + dx)
        x -= dx / 2
        return x

    truth = truth or {}

    if len(sub) > 0:
        sub = '_' + sub

    vals, prob, _ = SSP.marginalize(attr, attr2=attr2)

    if len(np.unique(vals)) == 1:
        print('{} not varied.'.format(attr))
        return

    if attr2 is None:
        if ax is None:
            _, ax = plt.subplots()
        ax.hist(vals, weights=prob, bins=center_grid(vals), histtype='step',
                lw=4, color='k')
        ax.set_ylabel(r'$\rm{Probability}$')
        ptype = 'marginal'
        if attr in truth:
            ax.axvline(truth[attr], color='darkred')
    else:
        [vals, vals2] = vals
        if len(np.unique(vals2)) == 1 or len(np.unique(vals)) == 1:
            print('{} not varied.'.format(attr))
            return

        if ax is None:
            _, ax = plt.subplots()

        ybins = center_grid(vals2)
        hist, xedge, yedge = np.histogram2d(vals, vals2,
                                            bins=[center_grid(vals),
                                                  center_grid(vals2)],
                                            weights=prob)
        img = ax.imshow(hist.T, origin='low', interpolation='nearest',
                        extent=[xedge[0], xedge[-1], yedge[0], yedge[-1]],
                        cmap=plt.cm.Blues, aspect='auto')

        # cbar = plt.colorbar(img)
        # cbar.set_label(r'$\rm{Probability}$')
        # cbar.set_clim(0, 1)
        ax.set_ylabel(key2label(attr2))
        ptype = '{}_joint'.format(attr2)
        if attr in truth:
            ax.axvline(truth[attr], color='darkred', lw=3)
        if attr2 in truth:
            ax.axhline(truth[attr2], color='darkred', lw=3)
            ax.set_xlim(xedge[0], xedge[-1])
            ax.set_ylim(yedge[0], yedge[-1])

    ax.set_xlabel(key2label(attr))
    if save:
        outfmt = '{}_{}{}_{}_gamma{}'
        outname = outfmt.format(SSP.name.replace('.csv', ''),
                                attr, sub, ptype, EXT)
        plt.savefig(outname, bbox_inches='tight')
        print('wrote {}'.format(outname))
        plt.close()
    return ax


def pdf_plots(SSP, marginals='default', sub='', twod=False, truth=None):
    """Call pdf_plot2d for a list of attr and attr2"""
    truth = truth or {}
    if marginals == 'default':
        marg = SSP._getmarginals()
    else:
        marg = marginals

    ndim = len(marg)
    if twod:
        raxs = []

        fig, axs = plt.subplots(nrows=ndim, ncols=ndim)
        [[ax.set_visible(False) for ax in axs[i, i:]] for i in range(ndim)]
        [[ax.tick_params(left='off', labelleft='off') for ax in axs.T[i]]
         for i in np.arange(ndim-1)+1]
        [ax.tick_params(bottom='off', labelbottom='off')
         for ax in axs[:-1].ravel()]
        [[ax.tick_params(right='off', top='off') for ax in axs[i+1, :i]]
         for i in range(ndim-1)]

        for j in marg:
            for i in marg:
                # e.g., skip Av vs Av and Av vs IMF
                # if already plotted IMF vs Av
                iy, ix = np.sort([marg.index(j), marg.index(i)])
                if i >= j:
                    continue

                    ax.tick_params(labelleft='off', labelbottom='off')
                raxs.append(SSP.pdf_plot(i, attr2=j, sub=sub, truth=truth,
                                         ax=axs[ix, iy]))
    else:
        fig, axs = plt.subplots(ncols=ndim, figsize=(15, 3))
        [ax.tick_params(left='off', labelleft='off') for ax in axs[1:]]
        [ax.tick_params(right='off', labelright='off') for ax in axs[:-1]]
        [ax.tick_params(top='off') for ax in axs]
        axs[-1].tick_params(labelright='on')
        raxs = [SSP.pdf_plot(i, sub=sub, truth=truth, ax=axs[marg.index(i)])
                for i in marg]
        [ax.set_ylabel('') for ax in axs[1:]]
    return fig, raxs


def key2label(string):
    """latex labels for different strings"""
    def_fmt = r'${}$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu$',
               'lage': r'$\log\ \rm{Age\ (yr)}$',
               'logZ': r'$\log\ \rm{Z}$',
               'fit': r'$\rm{Fit\ Parameter}$',
               'ov': r'$\Lambda_c$',
               'chi2': r'$\chi^2$',
               'bf': r'$\rm{Binary\ Fraction}$',
               'dav': r'$dA_V$'}
    if string not in convert.keys():
        return def_fmt.format(string)
    return convert[string]


def main(argv):
    """
    Main function for ssp.py plot or reformat ssp output.

    e.g., Reformat and then plot a OV=0.30 run:
    python -m match.scripts.ssp -fot --sub=ov0.30 *scrn
    """
    parser = argparse.ArgumentParser(description="stats for calcsfh -ssp")

    parser.add_argument('-f', '--format', action='store_true',
                        help='combine the files including data from filename')

    parser.add_argument('-d', '--oned', action='store_true',
                        help='make val vs prob plots')

    parser.add_argument('-t', '--twod', action='store_true',
                        help='make val vs val vs prob plots')

    parser.add_argument('-b', '--best', action='store_true',
                        help='include best fits only')

    parser.add_argument('--outfile', type=str,
                        default='combined_files.csv',
                        help='if -f file name to write to')

    parser.add_argument('--sub', type=str, default='',
                        help='add substring to figure names')

    parser.add_argument('-c', '--sspcombine', action='store_true',
                        help='run sspcombine on the file(s) (and exit)')

    parser.add_argument('-s', '--stats', action='store_true',
                        help='with -f, use stats (see .cmd.main -c)')

    parser.add_argument('-p', '--plotcmd', action='store_true',
                        help='run pgcmd (need .out.cmd file) and exit')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='invoke pdb')

    parser.add_argument('fnames', nargs='*', type=str,
                        help='ssp output(s) or formated output(s)')

    args = parser.parse_args(argv)

    if args.verbose:
        import pdb
        pdb.set_trace()

    if args.plotcmd:
        call_pgcmd_byfit(args.fnames, nmax=16)
        sys.exit(0)

    filtdict = {}
    if args.sub is not '':
        filtdict = filename_data(args.sub, skip=0)

    if args.format:
        fname = combine_files(args.fnames, outfile=args.outfile,
                              best=args.best)
    elif args.sspcombine:
        _ = [sspcombine(f, dry_run=False) for f in args.fnames]
        sys.exit(0)
    else:
        fname = args.fnames[0]

    ssp = SSP(filename=fname, filterby=filtdict)

    if args.oned:
        ssp.pdf_plots(sub=args.sub)

    if args.twod:
        ssp.pdf_plots(sub=args.sub, twod=args.twod)


if __name__ == "__main__":
    main(sys.argv[1:])
