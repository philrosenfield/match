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
from .utils import strip_header, center_grid
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
    def __init__(self, filename, data=None, filterby=None, gyr=False):
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

        self.gyr = gyr
        if gyr:
            data['lage'] = (10 ** (data['lage'] - 9))

        if filterby is not None:
            for key, val in filterby.items():
                data = data.loc[data[key] == val].copy(deep=True)

        self.data = data
        self.absprob = np.zeros(len(data))
        relprob = None
        if relprob is not None:
            x = data[relprob]
            un_x = np.unique(x)
            for ix in un_x:
                inds, = np.nonzero(x == ix)
                self.absprob[inds] = \
                    np.sum(np.exp(0.5 * self.data['fit'].iloc[inds].min() -
                                  self.data['fit'].iloc[inds]))
        else:
            self.absprob = \
                np.exp(0.5 * (data['fit'].min() - data['fit']))
        self.ibest = np.argmax(self.absprob)
        self.data = data

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

    def marginalize(self, attr, attr2=None, absprob=True):
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
                    if absprob:
                        prob[k] = np.sum(self.absprob.iloc[inds])
                    else:
                        prob[k] = \
                            np.sum(np.exp(0.5 *
                                          (self.data['fit'].iloc[inds].min() -
                                           self.data['fit'].iloc[inds])))
                    vals_x[k] = ix
                    vals_y[k] = iy
                    k += 1
            vals = [vals_x, vals_y]
        else:
            # compute linear probabilites
            # sum over probabilites for each unique grid value
            prob = np.zeros(len(unq_x))
            for i, ix in enumerate(unq_x):
                inds, = np.nonzero(x == ix)
                if absprob:
                    prob[i] = np.sum(self.absprob[inds])
                else:
                    prob[i] = \
                        np.sum(np.exp(0.5 * self.data['fit'].iloc[inds].min() -
                               self.data['fit'].iloc[inds]))
            vals = unq_x

        prob /= prob.sum()
        return vals, prob

    def pdf_plot(self, *args, **kwargs):
        return pdf_plot(self, *args, **kwargs)

    def pdf_plots(self, *args, **kwargs):
        return pdf_plots(self, *args, **kwargs)


def pdf_plot(SSP, attr, attr2=None, ax=None, sub=None, save=False,
             truth=None, absprob=True, useweights=False, plt_kw=None):
    """Plot prob vs marginalized attr"""
    plt_kw = plt_kw or {}
    do_cbar = False
    sub = sub or ''
    truth = truth or {}

    vals, prob = SSP.marginalize(attr, attr2=attr2, absprob=absprob)

    weights = None
    if useweights:
        weights = prob

    if len(np.unique(vals)) == 1:
        print('{} not varied.'.format(attr))
        return

    if attr2 is None:
        # plot type is marginal probability.
        ptype = 'marginal'
        if ax is None:
            _, ax = plt.subplots()
            do_cbar = True

        # grid edges
        xedge = center_grid(vals)

        # 1D Historgram
        ax.hist(vals, weights=weights, bins=xedge, histtype='step',
                lw=4, color='k')
        #ax.plot(vals, prob, lw=4, color='k')
        ax.set_ylabel(r'$\rm{Probability}$')
    else:
        [vals, vals2] = vals
        if len(np.unique(vals2)) == 1 or len(np.unique(vals)) == 1:
            print('{} not varied.'.format(attr))
            return

        # plot type is joint probability.
        ptype = '{}_joint'.format(attr2)

        # grid edges
        xedge = center_grid(vals)
        yedge = center_grid(vals2)

        if ax is None:
            _, ax = plt.subplots()
            docbar = True

        # 2D histogram weighted by probabibily
        hist, xedge, yedge = np.histogram2d(vals, vals2, bins=[xedge, yedge],
                                            weights=weights)
        img = ax.imshow(hist.T, origin='low', interpolation='nearest',
                        extent=[xedge[0], xedge[-1], yedge[0], yedge[-1]],
                        cmap=plt.cm.Blues, aspect='auto')

        if do_cbar:
            cbar = plt.colorbar(img)
            cbar.set_label(r'$\rm{Probability}$')

        ax.set_ylabel(key2label(attr2, gyr=SSP.gyr))
        if attr2 in truth:
            ax.axhline(truth[attr2], color='darkred', lw=3)
        ax.set_ylim(yedge[0], yedge[-1])

    if attr in truth:
        ax.axvline(truth[attr], color='darkred', lw=3)

    ax.set_xlim(xedge[0], xedge[-1])
    ax.set_xlabel(key2label(attr, gyr=SSP.gyr))
    if save:
        # add subdirectory to filename
        if len(sub) > 0:
            sub = '_' + sub
        outfmt = '{}_{}{}_{}_gamma{}'
        outname = outfmt.format(SSP.name.replace('.csv', ''),
                                attr, sub, ptype, EXT)
        plt.savefig(outname, bbox_inches='tight')
        print('wrote {}'.format(outname))
        plt.close()
    return ax


def pdf_plots(SSP, marginals=None, sub=None, twod=False, truth=None,
              text=None, absprob=True):
    """Call pdf_plot2d for a list of attr and attr2"""
    text = text or ''
    sub = sub or ''
    truth = truth or {}
    marg = marginals or SSP._getmarginals()

    ndim = len(marg)
    if twod:
        fig, axs = plt.subplots(nrows=ndim, ncols=ndim)
        [[ax.tick_params(left='off', labelleft='off') for ax in axs.T[i]]
         for i in np.arange(ndim-1)+1]
        [ax.tick_params(bottom='off', labelbottom='off')
         for ax in axs[:-1].ravel()]
        [[ax.tick_params(right='off', top='off') for ax in axs[i+1, :i]]
         for i in range(ndim-1)]
        raxs = []
        for j, mj in enumerate(marg):
            for i, mi in enumerate(marg):
                # e.g., skip Av vs Av and Av vs IMF
                # if already plotted IMF vs Av
                if i == j:
                    raxs.append(ssp.pdf_plot(mj, ax=axs[i, j],
                                             sub=sub, truth=truth,
                                             absprob=False))
                    # ax.tick_params(labelleft='off', labelbottom='off')
                elif i > j:
                    raxs.append(ssp.pdf_plot(mj, attr2=mi, ax=axs[i, j],
                                             sub=sub, truth=truth,
                                             absprob=True))
                else:
                    axs[i, j].set_visible(False)
    else:
        fig, axs = plt.subplots(ncols=ndim, figsize=(15, 3))
        [ax.tick_params(left='off', labelleft='off') for ax in axs[1:]]
        [ax.tick_params(right='off', labelright='off') for ax in axs[:-1]]
        [ax.tick_params(top='off') for ax in axs]
        axs[-1].tick_params(labelright='on')
        raxs = [SSP.pdf_plot(i, sub=sub, truth=truth, ax=axs[marg.index(i)],
                             absprob=absprob)
                for i in marg]
        [ax.set_ylabel('') for ax in axs[1:]]
        [ax.locator_params(axis='x', nbins=6) for ax in axs]
        if text:
            axs[-1].text(0.10, 0.90, '${}$'.format(text),
                         transform=axs[-1].transAxes)
        fig.subplots_adjust(bottom=0.2, left=0.05)
    return fig, raxs


def key2label(string, gyr=False):
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
               'dav': r'$dA_V$',
               'trueov': r'$\Lambda_c\ \rm{In}$'}
    if string not in convert.keys():
        convstr = def_fmt.format(string)
    else:
        convstr = convert[string]

    if gyr:
        convstr = convstr.replace('yr', 'Gyr').replace(r'\log\ ', '')
    return convstr


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

    parser.add_argument('--list', action='store_true',
                        help='fnames is one file with a list of filenames')

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

    # load filenames is it is a list (e.g., if $(ls *) is too long)
    if args.list:
        with open(args.fnames[0]) as inp:
            fnames = [l.strip() for l in inp.readlines()]
        args.fnames = fnames

    # plot cmd files by best fit nmax is the max number of images to make
    if args.plotcmd:
        call_pgcmd_byfit(args.fnames, nmax=16)
        sys.exit(0)

    # filter the ssp to only include sub
    filtdict = {}
    if args.sub is not '':
        filtdict = filename_data(args.sub, skip=0)

    # combine the match ssp output into one file
    if args.format:
        fname = combine_files(args.fnames, outfile=args.outfile,
                              best=args.best)
    # call match/sspcombine
    elif args.sspcombine:
        _ = [sspcombine(f, dry_run=False) for f in args.fnames]
        sys.exit(0)
    else:
        # Load one.
        fname = args.fnames[0]

    # load the combined ssp file or one match ssp output file
    ssp = SSP(filename=fname, filterby=filtdict)

    # call 1-plotting
    if args.oned:
        ssp.pdf_plots(sub=args.sub)

    # call 2-plotting
    if args.twod:
        ssp.pdf_plots(sub=args.sub, twod=args.twod)

if __name__ == "__main__":
    main(sys.argv[1:])
