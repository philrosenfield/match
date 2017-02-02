"""Stats and visualization of calcsfh -ssp runs"""
from __future__ import print_function, absolute_import
import argparse
import itertools
import os
import sys

import matplotlib.pylab as plt
import numpy as np

from . import utils
from .fileio import filename_data, read_ssp_output, combine_files
from .cmd import call_pgcmd_byfit
from .graphics.pdfs import pdf_plot, pdf_plots
from .wrappers.sspcombine import sspcombine

try:
    import seaborn
    seaborn.set()
except ImportError:
    pass

yeahpd = True
try:
    import pandas as pd
except ImportError:
    print('some SSP functions (-f, pdf plots) will be missing without pandas')
    yeahpd = False

__all__ = ['SSP']


def get_absprob(data):
    """absprob is the posterior since the fit parameter is -2 ln (posterior)"""
    data['absprob'] = np.exp(0.5 * (data['fit'].min() - data['fit']))
    return data


class SSP(object):
    """
    Class for calcsfh -ssp outputs
    """
    def __init__(self, filename=None, data=None, filterby=None, gyr=False):
        """
        filenames are the calcsfh -ssp terminal or console output.
        They do not need to be stripped of their header or footer or
        be concatenated as is typical in MATCH useage.
        """
        self.gyr = gyr
        self.frompost = False
        if filename is not None:
            data = self.load_ssp(filename)
            if 'post' in filename:
                self.frompost = True

        if data is not None:
            if gyr:
                data['lage'] = (10 ** (np.array(data['lage'],
                                                dtype=float) - 9))

            if filterby is not None:
                # Perhaps this should split into a dict instead of culling...
                for key, val in filterby.items():
                    if len(np.nonzero(data[key] == val)[0]) == 0:
                        print('can not filter by {0:s}={1:g}: no matching values'.format(key, val))
                        print('available values:', np.unique(data[key]))
                        sys.exit(1)
                    data = data.loc[data[key] == val].copy(deep=True)

            self.data = data

    def check_grid(self, skip_cols=None):
        """
        utils.marg and utils.marg2d assume a uniform grid to marginalize
        over. If, for example, the ssp file is a combination of calcsfh
        calls that overlap parts of parameter space or are missing parts of
        parameter space, the marginalized probabilities will not be correct,
        even in a relative way.

        This is a call to self.unique_ to check if the number of unique
        elements in each attribute is the same.

        Since the unique array is saved as an attribute, this will not add
        much time to a pdf_plots call.
        """
        skip_cols = skip_cols or []
        cols = [c for c in self.data.columns if c not in skip_cols or 'prob' in c]
        [self.unique_(c, check=True) for c in cols]

    def _getmarginals(self, avoid_list=['fit']):
        """get the values to marginalize over that exist in the data"""
        # marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav',
        # 'ov', 'bf'])
        marg_ = np.array([k for k in self.data.columns if k not in avoid_list])
        inds = [i for i, m in enumerate(marg_) if self._haskey(m)]
        return marg_[inds]

    def load_ssp(self, filename):
        """call fileio.read_ssp_output add file base, name to self"""
        self.base, self.name = os.path.split(filename)
        return read_ssp_output(filename)

    def _haskey(self, key):
        """test if the key requested is in self.data.columns"""
        ecode = True
        if key not in self.data.columns:
            ecode = False
        return ecode

    def build_posterior(self, xattr, arr, prob):
        """save the posterior in a DataFrame"""
        if not yeahpd:
            print('posterior functions need pandas')
            return
        if not hasattr(self, 'posterior'):
            self.posterior = pd.DataFrame()
        df = pd.DataFrame()
        df[xattr] = arr
        df['{:s}prob'.format(xattr)] = prob
        self.posterior = self.posterior.append(df, ignore_index=True)
        return

    def fitgauss1D(self, xattr, ux, prob):
        """Fit a 1D Gaussian to a marginalized probability
        see .utils.fitgauss1D
        sets attribute 'xattr'g
        """
        g = utils.fitgauss1D(ux, prob)
        self.__setattr__('{0:s}g'.format(xattr), g)
        return g

    def quantiles(self, xattr, ux, prob, qs=[0.16, 0.84], res=200, maxp=False,
                  ax=None, k=3):
        g = utils.quantiles(ux, prob, qs=qs, res=res, maxp=maxp, ax=ax,
                            k=k)
        self.__setattr__('{0:s}g'.format(xattr), g)
        return g

    def marg_table(filename=None, gauss1D=True):
        if filename is not None:
            self.load_posterior(filename)


    def write_posterior(self, filename='post.dat'):
        """write the posterior to a csv"""
        if not yeahpd:
            print('posterior functions need pandas')
        else:
            self.posterior.to_csv(filename, index=False)

    def load_posterior(self, filename):
        """read the posterior csv, see bulid_posterior and wrote_posterior"""
        self.base, self.name = os.path.split(filename)
        if not yeahpd:
            print('posterior functions need pandas, but a csv reader could go here.')
            self.data = None
            return
        self.data = pd.read_csv(filename)

    def unique_(self, attr, uniq_attr='u{:s}', check=False):
        """
        call np.unique on self.data[attr] if it has not already been called.

        Will store the unique array as an attribute passed with uniq_attr.

        check: print warning if there are number of unique array values
        are not equal.

        vdict is also added to self. It is a dictionary of attr: bool
        where the bool is True if the attr has more than one value.
        """
        if not hasattr(self, 'vdict'):
            self.vdict = {}
        if self.frompost:
            self.vdict[attr] = True
            return self.data[attr]
        self.vdict[attr] = False
        uatr = uniq_attr.format(attr)
        if not hasattr(self, uatr):
            self.data[attr] = np.array(self.data[attr], dtype=float)
            uns, idx = np.unique(self.data[attr], return_counts=True)
            if check:
                unc, cidx = np.unique(idx, return_index=True)
                if len(unc) > 1:
                    print('{} grid is uneven: {}'.format(attr, uns[cidx[1:]]))
                    print(unc)
            self.__setattr__(uatr, uns)
        u = self.__getattribute__(uatr)
        if len(u) > 1:
            self.vdict[attr] = True
        return u

    def marginalize(self, xattr, yattr=None, **kwargs):
        """
        Marginalize over one or two quanitities
        xattr, yattr : string, string
            data column to marginalize over

        Returns
        vals : array or list of arrays
            if yattr is passed, vals is output of utils.cantered_meshgrid
            otherwise it's the unique values of data[xattr]
        prob : return from utils.marg or utils.marg2d

        NOTE:
        marg/marg2d only work for values calculated with an equal spaced grid.
        see self.check_grid
        """
        assert self._haskey(xattr), '{} not found'.format(xattr)
        if not self._haskey('absprob'):
            self.data = get_absprob(self.data)

        z = self.data.absprob
        x = self.data[xattr]
        ux = self.unique_(xattr)

        if yattr is not None:
            assert self._haskey(yattr), '{} not found'.format(yattr)
            y = self.data[yattr]
            uy = self.unique_(yattr)
            prob, ux, uy = utils.marg2d(x, y, z, unx=ux, uny=uy, **kwargs)
            vals = utils.centered_meshgrid(x, y, unx=ux, uny=uy)
        else:
            prob, ux = utils.marg(x, z, unx=ux, **kwargs)
            vals = ux

        return vals, prob

    def pdf_plot(self, *args, **kwargs):
        return pdf_plot(self, *args, **kwargs)

    def pdf_plots(self, *args, **kwargs):
        return pdf_plots(self, *args, **kwargs)


def parse_args(argv=None):
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
                        help='ssp output(s) or combined output')

    return parser.parse_args(argv)


def main(argv=None):
    """
    Main function for ssp.py plot or reformat ssp output.

    e.g., Reformat and then plot a OV=0.30 run:
    python -m match.scripts.ssp -fot --sub=ov0.30 *scrn
    """

    args = parse_args(argv)

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
        return

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
        return
    else:
        # Load one file
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
    sys.exit(main())
