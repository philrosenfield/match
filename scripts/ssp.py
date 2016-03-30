"""Stats and visualization of calcsfh -ssp runs"""
import argparse
import itertools
import os
import sys

import pandas as pd
import matplotlib.pylab as plt
import numpy as np

from .config import EXT, match_base
from .fileio import filename_data, add_filename_info_to_file
from .utils import strip_header

__all__ = ['SSP']

try:
    plt.style.use('presentation')
except:
    pass


def combine_files(fnames, outfile='combined_files.csv'):
    all_data = pd.DataFrame()
    for fname in fnames:
        df = add_filename_info_to_file(fname)
        all_data = all_data.append(df, ignore_index=True)

    all_data.to_csv(outfile)
    return


def sspcombine(fname, dry_run=True, outfile=None):
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
    def __init__(self, filename, data=None, filterby=None):
        """
        filenames are the calcsfh -ssp terminal or console output.
        They do not need to be stripped of their header or footer or
        be concatenated as is typical in MATCH useage.
        """
        self.base, self.name = os.path.split(filename)
        if data is None:
            data = pd.read_csv(filename)

        if filterby is not None:
            for k, v in filterby.items():
                data = data[data[k] == v].copy(deep=True)
        
        self.data = data

        self.ibest = np.argmin(self.data['fit'])

        self.absprob = \
            np.exp(0.5 * (self.data['fit'].min() - self.data['fit']))

    def _getmarginals(self):
        """get the values to marginalize over that exist in the data"""
        marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav', 'ov',
                         'bf'])
        inds = [i for i, m in enumerate(marg) if self._haskey(m)]
        return marg[inds]


    def _haskey(self, key):
        """test if the key requested is available"""
        ecode = True
        if not key in self.data.keys():
            ecode = False
        return ecode

    def marginalize(self, attr, attr2=None):
        """Find the best fit for each unique value of attr"""
        assert self._haskey(attr), '{} not found'.format(attr)
        x = self.data[attr]
        unq_x = np.unique(x)
        
        size = len(unq_x)
        prob = np.zeros(size)
        vals_x = np.zeros(size)
        
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
            for i, ix in enumerate(vals_x):
                prob[i] = np.sum(self.absprob[x == ix])
            vals = vals_x

        prob /= prob.sum()
        return vals, prob, np.log(prob)


    def pdf_plot(self, attr, attr2=None, ax=None, sub=''):
        """Plot prob vs marginalized attr"""
        save = False
        if len(sub) > 0:
            sub = '_' + sub

        if ax is None:
            fig, ax = plt.subplots()
            save = True

        vals, prob, _ = self.marginalize(attr, attr2=attr2)

        if attr2 is None:
            ax.hist(vals, weights=prob, bins=31, histtype='step',
                    lw=4, color='k')
            ax.set_ylabel(r'$\rm{Probability}$')
            ptype = 'marginal'
        else:
            [vals, vals2] = vals

            h, xe, ye = np.histogram2d(vals, vals2, weights=prob)
            c = plt.imshow(h.T, origin='low', interpolation='nearest',
                           extent=[xe[0], xe[-1], ye[0], ye[-1]],
                           cmap=plt.cm.Blues, aspect='auto')

            cb = plt.colorbar(c)
            cb.set_label(r'$\rm{Probability}$')
            ax.set_ylabel(key2label(attr2))
            ptype = 'joint'
        ax.set_xlabel(key2label(attr))
        
        if save:
            outname ='{}_{}{}_{}_gamma{}'.format(self.name.split('_')[0],
                                                 attr, sub, ptye, EXT)
            plt.savefig(outname, bbox_inches='tight')
            plt.close()
        return ax
    

    def pdf_plots(self, marginals='default', sub='', twod=False):
        """Call pdf_plot2d for a list of attr and attr2"""
        if marginals=='default':
            marg = self._getmarginals()
        else:
            marg = marginals

        if twod:
            for i, j in itertools.product(marg, marg):
                # e.g., skip Av vs Av and Av vs IMF if already plotted IMF vs Av
                if i >= j:
                    continue
                             
                self.pdf_plot(i, attr2=j, sub=sub)
        else:
            [self.pdf_plot(i, sub=sub) for i in marg]
        
        return

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
    if not string in convert.keys():
        return def_fmt.format(string)
    return convert[string]

def main(argv):
    """
    Main function for ssp.py plot or reformat ssp output.

    e.g., Reformat and then plot a OV=0.30 run:
    python -m match.scripts.ssp -fot --sub=ov0.30 *scrn
    """
    parser = argparse.ArgumentParser(description="Plot or reformat calcsfh -ssp output")

    parser.add_argument('-f', '--format', action='store_true',
                        help='combine the files including data in the filename')

    parser.add_argument('-o', '--oned', action='store_true',
                        help='make val vs prob plots')

    parser.add_argument('-t', '--twod', action='store_true',
                        help='make val vs val vs prob plots')

    parser.add_argument('-s', '--sub', type=str, default='',
                        help='add substring to figure names')

    parser.add_argument('-l', '--list', action='store_true',
                        help='fnames is a file with a list of files to read')

    parser.add_argument('-c', '--sspcombine', action='store_true',
                        help='run sspcombine on the file(s) (and exit)')

    parser.add_argument('fnames', nargs='*', type=str,
                        help='ssp output(s) or formated output(s)')

    args = parser.parse_args(argv)

    if args.list:
        args.fnames = map(str.strip, open(args.fnames[0], 'r').readlines())

    filtdict = {}
    if args.sub is not '':
        filtdict = filename_data(args.sub, skip=0)

    if args.format:
        fname = combine_files(args.fnames)        
    elif args.sspcombine:
        [sspcombine(f, dry_run=False) for f in args.fnames]
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
