"""Class for reading the output.cmd file from calcsfh"""
from __future__ import print_function, absolute_import
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from match.scripts.config import EXT
from match.scripts.fileio import read_match_cmd
from match.scripts.graphics.graphics import match_plot
from match.scripts.utils import parse_pipeline
from match.scripts.wrappers.stats import call_stats

__all__ = ['CMD']


class CMD(object):
    """
    A quikly made object to read the MATCH CMD file and hold paramters to
    automatically make plots with the same color scale as other MATCH CMD
    files.
    """
    def __init__(self, filename=None):
        if filename is not None:
            self.base, self.name = os.path.split(os.path.abspath(filename))
            self.cmd, self.fit, self.colors, self.yfilter = \
                read_match_cmd(filename)
            self.load_match_cmd()

    def set_labels(self):
        """Set up list of labels for pgpro"""
        strfmt = r'${{\rm {:s}}}$'
        labels = [strfmt.format(i) for i in ['data', 'model', 'd-m']]
        try:
            target, _ = parse_pipeline(self.name)
            labels[0] = strfmt.format(target)
        except:
            pass
        labels.append(r'$-2\ln P = {:g}$'.format(self.fit))
        return labels

    def set_axis_labels(self, ax=None):
        xlabel = r'$\rm{{{}}}$'.format(self.colors.replace('WFC', 'F').replace('UVIS', 'F'))
        ylabel = r'$\rm{{{}}}$'.format(self.yfilter.replace('WFC', 'F').replace('UVIS', 'F'))
        if ax is not None:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        return xlabel, ylabel

    def load_cmd_fromfake(self, filename, dcol=0.05, dmag=0.1, ymag='V',
                          yfilter=None, colors=None, xlim=None, ylim=None):
        self.data = None
        self.diff = None
        self.sig = None
        self.colors = colors
        self.yfilter = yfilter

        mag1, mag2 = np.loadtxt(filename, unpack=True)
        inds, = np.nonzero((mag1 < 30) & (mag2 < 30))
        mag1 = mag1[inds]
        mag2 = mag2[inds]
        color = mag1 - mag2
        if ymag.upper() == 'V':
            mag = mag1
        elif ymag.upper() == 'I':
            mag = mag2
        else:
            print('Error: ymag needs to be V or I')
            return

        if xlim is None:
            cmin = np.min(color)
            cmax = np.max(color)
        else:
            cmin, cmax = xlim
        if ylim is None:
            mmax = np.max(mag)
            mmin = np.min(mag)
        else:
            mmin, mmax = np.sort(ylim)

        cbins = np.arange(cmin, cmax, dcol)
        mbins = np.arange(mmin, mmax, dmag)
        H, ce, me = np.histogram2d(color, mag, bins=[cbins, mbins])

        self.model = H.T
        self.extent = [cmin, cmax, mmax, mmin]
        self.nmagbin = len(mbins)
        self.ncolbin = len(cbins)

    def load_match_cmd(self):
        """
        pgcmd needs hesses and extent. Optional are max_* which set the vmins
        and vmaxs.
        """
        self.nmagbin = len(np.unique(self.cmd['mag']))
        self.ncolbin = len(np.unique(self.cmd['color']))
        self.data = self.cmd['Nobs'].reshape(self.nmagbin, self.ncolbin)
        self.model = self.cmd['Nsim'].reshape(self.nmagbin, self.ncolbin)
        self.diff = self.cmd['diff'].reshape(self.nmagbin, self.ncolbin)
        self.sig = self.cmd['sig'].reshape(self.nmagbin, self.ncolbin)
        self.hesses = [self.data, self.model, self.diff, self.sig]

        self.extent = [np.min(self.cmd['color']), np.max(self.cmd['color']),
                       np.max(self.cmd['mag']), np.min(self.cmd['mag'])]
        # self.max_counts = np.nanmax(np.concatenate([self.data, self.model]))
        # self.max_diff = np.nanmax(np.abs(self.diff))
        # self.max_sig = np.nanmax(np.abs(self.sig))

    def pgcmd(self, labels=None, outdir=None, logcounts=False, figname=None,
              twobytwo=True):
        '''produce the image that pgcmd.pro makes'''
        labels = labels or self.set_labels()

        if figname is None:
            base = outdir or self.base
            assert os.path.isdir(base), \
                '{} directory not found'.format(base)
            figname = os.path.join(base, os.path.split(self.name)[1] + EXT)

        hesses = self.hesses
        if logcounts:
            hesses[0] = np.log10(hesses[0])
            hesses[1] = np.log10(hesses[1])

        xlabel, ylabel = self.set_axis_labels()

        grid = match_plot(hesses, self.extent, labels=labels, ylabel=ylabel,
                          xlabel=xlabel, twobytwo=twobytwo)

        gates = self.cmd['gate']
        ugates = np.unique(gates)
        if len(ugates) > 1:
            dinds = np.digitize(gates, bins=np.unique(gates, right=True))
            _ = [grid.axes_all[0].plot(self.cmd['color'][dinds == i],
                                       self.cmd['mag'][dinds == i],
                                       '.', alpha=0.3)
                 for i in range(len(dinds)) if i > 0]

        for ax in grid.axes_all:
            ax.locator_params(axis='x', nbins=6)

        plt.savefig(figname)
        plt.close()
        print('wrote {}'.format(figname))
        return grid


def call_pgcmd_byfit(cmdfns, nmax=5, outdir=None, logcounts=False):
    """
    Call pgcmd with filename order them by increasing best fit value.
    cmdfns : string or list
        .cmd filename or list of filenames.

    nmax : int
        make best nmax plots.
    """
    if not isinstance(cmdfns, list):
        cmd = CMD(cmdfns)
        cmd.pgcmd(outdir=outdir, logcounts=logcounts)
        return

    cmds = np.array([CMD(cmdfn) for cmdfn in cmdfns])
    try:
        icmd = np.argsort(np.concatenate([cmd.fit for cmd in cmds]))
    except ValueError:
        icmd = np.argsort([cmd.fit for cmd in cmds])
    for j, cmd in enumerate(cmds[icmd]):
        if j > nmax:
            break
        jstr = ('{}'.format(j)).zfill(4)
        cmd.pgcmd(figname='{}{}{}'.format(cmd.name, jstr, EXT),
                  outdir=outdir, logcounts=logcounts)
    return


def main(argv):
    """main function for cmd"""
    parser = argparse.ArgumentParser(description="plot cmd file")

    parser.add_argument('-o', '--outdir', type=str,
                        help='directory to place images')

    parser.add_argument('-f', '--byfit', action='store_true',
                        help='number filenames by best fit')

    parser.add_argument('-c', '--stats', action='store_true',
                        help='call match/bin/stats and exit')

    parser.add_argument('--logcounts', action='store_true',
                        help='use log binning for data and model')

    parser.add_argument('--nmax', type=int,
                        help='max best cmds to plot with --byfit')

    parser.add_argument('cmdfiles', type=str, nargs='*',
                        help='.cmd files to plot')

    args = parser.parse_args(argv)

    if args.outdir is not None:
        assert os.path.isdir(args.outdir), \
            'directory {} not found'.format(args.outdir)

    if args.stats:
        call_stats(args.cmdfiles, outdir=args.outdir)
        sys.exit()

    if args.byfit:
        if args.nmax is None:
            args.nmax = len(args.cmdfiles)
        call_pgcmd_byfit(args.cmdfiles, nmax=args.nmax,
                         outdir=args.outdir, logcounts=args.logcounts)
    else:
        for cmdfile in args.cmdfiles:
            cmd = CMD(cmdfile)
            cmd.pgcmd(outdir=args.outdir, logcounts=args.logcounts)

if __name__ == "__main__":
    main(sys.argv[1:])
