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

from .config import match_base
from .fileio import filename_data, add_filename_info_to_file
from .utils import strip_header, centered_meshgrid, marg, marg2d
from .cmd import call_pgcmd_byfit
from .graphics import pdf_plot, pdf_plots

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
                # combined file of many calcsfh outputs
                data = pd.read_csv(filename)
            else:
                # one calcsfh output
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
        self.absprob = np.exp(0.5 * (data['fit'].min() - data['fit']))
        self.ibest = np.argmax(self.absprob)

    def _getmarginals(self, avoid_list=['fit']):
        """get the values to marginalize over that exist in the data"""
        # marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav',
        # 'ov', 'bf'])
        marg_ = np.array([k for k in self.data.columns if k not in avoid_list])
        inds = [i for i, m in enumerate(marg_) if self._haskey(m)]
        return marg_[inds]

    def _haskey(self, key):
        """test if the key requested is available"""
        ecode = True
        if key not in self.data.columns:
            ecode = False
        return ecode

    def marginalize(self, xattr, yattr=None):
        """
        Marginalize over one or two quanitities
        xattr : string
            data column to marginalize over
        yattr : string
            data column to marginalize over

        Returns
        vals : array or list of arrays
            if yattr is passed, vals is output of utils.cantered_meshgrid
            otherwise it's the unique values of data[xattr]
        prob : the minimum -2 ln P in each bin
        """
        assert self._haskey(xattr), '{} not found'.format(xattr)
        z = self.data.fit
        x = self.data[xattr]

        if yattr is not None:
            assert self._haskey(yattr), '{} not found'.format(yattr)
            y = self.data[yattr]
            vals = centered_meshgrid(x, y)
            prob = marg2d(x, y, z)
        else:
            vals = np.unique(x)
            prob = marg(x, z)
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
                        help='ssp output(s) or formated output(s)')

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
