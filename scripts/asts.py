""" All about Artificial star tests """
from __future__ import print_function
import argparse
import logging
import os
from astropy.io import fits
import re
import sys

import matplotlib.pylab as plt
import numpy as np
from .config import EXT
from scipy.interpolate import interp1d
from .fileio import parse_pipeline

logger = logging.getLogger(__name__)

__all__ = ['ast_correct_starpop', 'ASTs', 'parse_pipeline']

plt.style.use('ggplot')

def hess(color, mag, binsize, **kw):
    """Compute a hess diagram (surface-density CMD) on photometry data."""
    defaults = dict(mbin=None, cbin=None, verbose=False)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    if kw['mbin'] is None:
        mbin = np.arange(mag.min(), mag.max(), binsize)
    else:
        mbin = np.array(kw['mbin']).copy()
    if kw['cbin'] is None:
        cbinsize = kw.get('cbinsize')
        if cbinsize is None:
            cbinsize = binsize
        cbin = np.arange(color.min(), color.max(), cbinsize)
    else:
        cbin = np.array(kw['cbin']).copy()

    hesst, cbin, mbin = np.histogram2d(color, mag, bins=[cbin, mbin])
    hess = hesst.T
    return (cbin, mbin, hess)


def ast_correct_starpop(sgal, fake_file=None, outfile=None, overwrite=False,
                        asts_obj=None, correct_kw={}, diag_plot=False,
                        plt_kw={}, hdf5=True):
    '''
    correct mags with artificial star tests, finds filters by fake_file name

    Parameters
    ----------
    sgal : galaxies.SimGalaxy or StarPop instance
        must have apparent mags (corrected for dmod and Av)

    fake_file : string
         matchfake file

    outfile : string
        if sgal, a place to write the table with ast_corrections

    overwrite : bool
        if sgal and outfile, overwite if outfile exists

    asts_obj : AST instance
        if not loading from fake_file

    correct_kw : dict
        passed to ASTs.correct important to consider, dxy, xrange, yrange
        see AST.correct.__doc__

    diag_plot : bool
        make a mag vs mag diff plot

    plt_kw :
        kwargs to pass to pylab.plot

    Returns
    -------
    adds corrected mag1 and mag2

    If sgal, adds columns to sgal.data
    '''
    fmt = '{}_cor'
    if asts_obj is None:
        sgal.fake_file = fake_file
        _, filter1, filter2 = parse_pipeline(fake_file)
        if fmt.format(filter1) in sgal.data.keys() or fmt.format(filter2) in sgal.data.keys():
            errfmt = '{}, {} ast corrections already in file.'
            logger.warning(errfmt.format(filter1, filter2))
            return sgal.data[fmt.format(filter1)], sgal.data[fmt.format(filter2)]
        ast = ASTs(fake_file)
    else:
        ast = asts_obj

    mag1 = sgal.data[ast.filter1]
    mag2 = sgal.data[ast.filter2]

    correct_kw = dict({'dxy': (0.2, 0.15)}.items() + correct_kw.items())
    cor_mag1, cor_mag2 = ast.correct(mag1, mag2, **correct_kw)
    for name, data in zip([fmt.format(ast.filter1), fmt.format(ast.filter2)],
                          [cor_mag1, cor_mag2]):
        sgal.add_data(name, data)

        if outfile is not None:
            sgal.write_data(outfile, overwrite=overwrite, hdf5=hdf5)

    if diag_plot:
        from ..fileio.fileIO import replace_ext
        plt_kw = dict({'color': 'navy', 'alpha': 0.3, 'label': 'sim'}.items() \
                      + plt_kw.items())
        axs = ast.magdiff_plot()
        mag1diff = cor_mag1 - mag1
        mag2diff = cor_mag2 - mag2
        rec, = np.nonzero((np.abs(mag1diff) < 10) & (np.abs(mag2diff) < 10))
        axs[0].plot(mag1[rec], mag1diff[rec], '.', **plt_kw)
        axs[1].plot(mag2[rec], mag2diff[rec], '.', **plt_kw)
        if 'label' in plt_kw.keys():
            [ax.legend(loc=0, frameon=False) for ax in axs]
        plt.savefig(replace_ext(outfile, '_ast_correction{}'.format(EXT)))
    return cor_mag1, cor_mag2


class ASTs(object):
    '''class for reading and using artificial stars'''
    def __init__(self, filename, filter1=None, filter2=None, filt_extra=''):
        '''
        if filename has 'match' in it will assume this is a matchfake file.
        if filename has .fits extention will assume it's a binary fits table.
        '''
        self.base, self.name = os.path.split(filename)
        self.filter1 = filter1
        self.filter2 = filter2
        self.filt_extra = filt_extra

        self.target, filters = parse_pipeline(filename)

        try:
            self.filter1, self.filter2 = filters
        except:
            self.filter1, self.filter2, self.filter3 = filters
        self.read_file(filename)

    def recovered(self, threshold=9.99):
        '''
        find indicies of stars with magdiff < threshold

        Parameters
        ----------
        threshold: float
            [9.99] magin - magout threshold for recovery

        Returns
        -------
        self.rec: list
            recovered stars in both filters
        rec1, rec2: list, list
            recovered stars in filter1, filter2
        '''
        rec1, = np.nonzero(np.abs(self.mag1diff) < threshold)
        rec2, = np.nonzero(np.abs(self.mag2diff) < threshold)
        self.rec = list(set(rec1) & set(rec2))
        if len(self.rec) == len(self.mag1diff):
            logger.warning('all stars recovered')
        return rec1, rec2

    def make_hess(self, binsize=0.1, yattr='mag2diff', hess_kw={}):
        '''make hess grid'''
        self.colordiff = self.mag1diff - self.mag2diff
        mag = self.__getattribute__(yattr)
        self.hess = hess(self.colordiff, mag, binsize, **hess_kw)

    def read_file(self, filename):
        '''
        read MATCH fake file into attributes
        format is mag1in mag1diff mag2in mag2diff
        mag1 is assumed to be mag1in
        mag2 is assumed to be mag2in
        mag1diff is assumed to be mag1in-mag1out
        mag2diff is assumed to be mag2in-mag2out
        '''
        if not filename.endswith('.fits'):
            names = ['mag1', 'mag2', 'mag1diff', 'mag2diff']
            self.data = np.genfromtxt(filename, names=names)
            # unpack into attribues
            for name in names:
                self.__setattr__(name, self.data[name])
        else:
            assert not None in [self.filter1, self.filter2], \
                'Must specify filter strings'
            self.data = fits.getdata(filename)
            self.mag1 = self.data['{}_IN'.format(self.filter1)]
            self.mag2 = self.data['{}_IN'.format(self.filter2)]
            mag1out = self.data['{}{}'.format(self.filter1, self.filt_extra)]
            mag2out = self.data['{}{}'.format(self.filter2, self.filt_extra)]
            self.mag1diff = self.mag1 - mag1out
            self.mag2diff = self.mag2 - mag2out

    def write_matchfake(self, newfile):
        '''write matchfake file'''
        dat = np.array([self.mag1, self.mag2, self.mag1diff, self.mag2diff]).T
        np.savetxt(newfile, dat, fmt='%.3f')

    def bin_asts(self, binsize=0.2, bins=None):
        '''
        bin the artificial star tests

        Parameters
        ----------
        bins: bins for the asts
        binsize: width of bins for the asts

        Returns
        -------
        self.am1_inds, self.am2_inds: the indices of the bins to
            which each value in mag1 and mag2 belong (see np.digitize).
        self.ast_bins: bins used for the asts.
        '''
        if bins is None:
            ast_max = np.max(np.concatenate((self.mag1, self.mag2)))
            ast_min = np.min(np.concatenate((self.mag1, self.mag2)))
            self.ast_bins = np.arange(ast_min, ast_max, binsize)
        else:
            self.ast_bins = bins

        self.am1_inds = np.digitize(self.mag1, self.ast_bins)
        self.am2_inds = np.digitize(self.mag2, self.ast_bins)

    def _random_select(self, arr, nselections):
        '''
        randomly sample arr nselections times

        Parameters
        ----------
        arr : array or list
            input to sample
        nselections : int
            number of times to sample

        Returns
        -------
        rands : array
            len(nselections) of randomly selected from arr (duplicates included)
        '''
        rands = np.array([np.random.choice(arr) for i in range(nselections)])
        return rands

    def ast_correction(self, obs_mag1, obs_mag2, binsize=0.2, bins=None,
                       not_rec_val=np.nan, missing_data1=0., missing_data2=0.):
        '''
        Apply ast correction to input mags.

        Corrections are made by going through obs_mag1 in bins of
        bin_asts and randomly selecting magdiff values in that ast_bin.
        obs_mag2 simply follows along since it is tied to obs_mag1.

        Random selection was chosen because of the spatial nature of
        artificial star tests. If there are 400 asts in one mag bin,
        and 30 are not recovered, random selection should match the
        distribution (if there are many obs stars).

        If there are obs stars in a mag bin where there are no asts,
        will throw the star out unless the completeness in that mag bin
        is more than 50%.
        Parameters
        ----------
        obs_mag1, obs_mag2 : N, 1 arrays
            input observerd mags

        binsize, bins : sent to bin_asts

        not_rec_val : float
            value for not recovered ast
        missing_data1, missing_data2 : float, float
            value for data outside ast limits per filter (include=0)

        Returns
        -------
        cor_mag1, cor_mag2: array, array
            ast corrected magnitudes

        Raises:
            returns -1 if obs_mag1 and obs_mag2 are different sizes

        To do:
        possibly return magXdiff rather than magX + magXdiff?
        reason not to: using AST results from one filter to another isn't
        kosher. At least not glatt kosher.
        '''
        self.completeness(combined_filters=True, interpolate=True)

        nstars = obs_mag1.size
        if obs_mag1.size != obs_mag2.size:
            logger.error('mag arrays of different lengths')
            return -1

        # corrected mags are filled with nan.
        cor_mag1 = np.empty(nstars)
        cor_mag1.fill(not_rec_val)
        cor_mag2 = np.empty(nstars)
        cor_mag2.fill(not_rec_val)

        # need asts to be binned for this method.
        if not hasattr(self, 'ast_bins'):
            self.bin_asts(binsize=binsize, bins=bins)
        om1_inds = np.digitize(obs_mag1, self.ast_bins)

        for i in range(len(self.ast_bins)):
            # the obs and artificial stars in each bin
            obsbin, = np.nonzero(om1_inds == i)
            astbin, = np.nonzero(self.am1_inds == i)

            nobs = len(obsbin)
            nast = len(astbin)
            if nobs == 0:
                # no stars in this mag bin to correct
                continue
            if nast == 0:
                # no asts in this bin, probably means the simulation
                # is too deep
                if self.fcomp2(self.ast_bins[i]) < 0.5:
                    continue
                else:
                    # model is producing stars where there was no data.
                    # assign correction for missing data
                    cor1 = missing_data1
                    cor2 = missing_data2
            else:
                # randomly select the appropriate ast correction for obs stars
                # in this bin
                cor1 = self._random_select(self.mag1diff[astbin], nobs)
                cor2 = self._random_select(self.mag2diff[astbin], nobs)

            # apply corrections
            cor_mag1[obsbin] = obs_mag1[obsbin] + cor1
            cor_mag2[obsbin] = obs_mag2[obsbin] + cor2
            # finite values only: not implemented because trilegal array should
            # maintain the same size.
            #fin1, = np.nonzero(np.isfinite(cor_mag1))
            #fin2, = np.nonzero(np.isfinite(cor_mag2))
            #fin = list(set(fin1) & set(fin2))
        return cor_mag1, cor_mag2

    def correct(self, obs_mag1, obs_mag2, bins=[100,200], xrange=[-0.5, 5.],
                yrange=[15., 27.], not_rec_val=0., dxy=None):
        """
        apply AST correction to obs_mag1 and obs_mag2

        Parameters
        ----------
        obs_mag1, obs_mag2 : arrays
            input mags to correct

        bins : [int, int]
            bins to pass to graphics.plotting.crazy_histogram2d

        xrange, yrange : shape 2, arrays
            limits of cmd space send to graphics.plotting.crazy_histogram2d
            since graphics.plotting.crazy_histogram2d is called twice it is
            important to have same bin sizes

        not_rec_val : float or nan
            value to fill output arrays where obs cmd does not overlap with
            ast cmd.

        dxy : array shape 2,
            color and mag step size to make graphics.plotting.crazy_histogram2d

        Returns
        -------
        cor_mag1, cor_mag2 : arrays len obs_mag1, obs_mag2
            corrections to obs_mag1 and obs_mag2
        """
        from ..graphics.plotting import crazy_histogram2d as chist

        nstars = obs_mag1.size
        if obs_mag1.size != obs_mag2.size:
            logger.error('mag arrays of different lengths')
            return -1, -1

        # corrected mags are filled with nan.
        cor_mag1 = np.empty(nstars)
        cor_mag1.fill(not_rec_val)
        cor_mag2 = np.empty(nstars)
        cor_mag2.fill(not_rec_val)

        obs_color = obs_mag1 - obs_mag2
        ast_color = self.mag1 - self.mag2

        if dxy is not None:
            # approx number of bins.
            bins[0] = len(np.arange(*xrange, step=dxy[0]))
            bins[1] = len(np.arange(*yrange, step=dxy[1]))

        ckw = {'bins': bins, 'reverse_indices': True, 'xrange': xrange,
                    'yrange': yrange}
        SH, _, _, sixy, sinds = chist(ast_color, self.mag2, **ckw)
        H, _, _, ixy, inds = chist(obs_color, obs_mag2, **ckw)

        x, y = np.nonzero(SH * H > 0)
        # there is a way to do this with masking ...
        for i, j in zip(x, y):
            sind, = np.nonzero((sixy[:, 0] == i) & (sixy[:, 1] == j))
            hind, = np.nonzero((ixy[:, 0] == i) & (ixy[:, 1] == j))
            nobs = int(H[i, j])
            xinds = self._random_select(sinds[sind], nobs)
            cor_mag1[inds[hind]] = self.mag1diff[xinds]
            cor_mag2[inds[hind]] = self.mag2diff[xinds]

        return obs_mag1 + cor_mag1, obs_mag2 + cor_mag2

    def completeness(self, combined_filters=False, interpolate=False,
                     binsize=0.2):
        '''
        calculate the completeness of the data in each filter

        Parameters
        ----------
        combined_filters : bool
            Use individual or combined ast recovery

        interpolate : bool
            add a 1d spline the completeness function to self

        Returns
        -------
        self.comp1, self.comp2 : array, array
            the completeness per filter binned with self.ast_bins
        '''
        # calculate stars recovered, could pass theshold here.
        rec1, rec2 = self.recovered()

        # make sure ast_bins are good to go
        if not hasattr(self, 'ast_bins'):
            self.bin_asts(binsize=binsize)

        # gst uses both filters for recovery.
        if combined_filters is True:
            rec1 = rec2 = self.rec

        # historgram of all artificial stars
        qhist1 = np.array(np.histogram(self.mag1, bins=self.ast_bins)[0],
                          dtype=float)

        # histogram of recovered artificial stars
        rhist1 = np.array(np.histogram(self.mag1[rec1], bins=self.ast_bins)[0],
                          dtype=float)

        # completeness histogram
        self.comp1 = rhist1 / qhist1

        qhist2 = np.array(np.histogram(self.mag2, bins=self.ast_bins)[0],
                          dtype=float)
        rhist2 = np.array(np.histogram(self.mag2[rec2], bins=self.ast_bins)[0],
                          dtype=float)
        self.comp2 = rhist2 / qhist2

        if interpolate is True:
            # sometimes the histogram isn't as useful as the a spline
            # function... add the interp1d function to self.
            self.fcomp1 = interp1d(self.ast_bins[1:], self.comp1,
                                   bounds_error=False)
            self.fcomp2 = interp1d(self.ast_bins[1:], self.comp2,
                                   bounds_error=False)
        return

    def get_completeness_fraction(self, frac, dmag=0.001, bright_lim=18):
        """Find the completeness magnitude at a given fraction"""
        assert hasattr(self, 'fcomp1'), \
            'need to run completeness with interpolate=True'

        # set up array to evaluate interpolation
        # sometimes with few asts at bright mags the curve starts with low
        # completeness, reaches toward 1, and then declines as expected.
        # To get around taking a value too bright, I search for values beginning
        # at the faint end
        search_arr = np.arange(bright_lim, 31, dmag)[::-1]

        # completeness in each filter, and the finite vals
        # (frac - nan = frac)
        cfrac1 = self.fcomp1(search_arr)
        ifin1 = np.isfinite(cfrac1)

        cfrac2 = self.fcomp2(search_arr)
        ifin2 = np.isfinite(cfrac2)

        # closest completeness fraction to passed fraction
        icomp1 = np.argmin(np.abs(frac - cfrac1[ifin1]))
        icomp2 = np.argmin(np.abs(frac - cfrac2[ifin2]))

        # mag associated with completeness
        comp1 = search_arr[ifin1][icomp1]
        comp2 = search_arr[ifin2][icomp2]

        if comp1 == bright_lim or comp2 == bright_lim:
            logger.warning('Completeness fraction is at mag search limit and probably wrong. '
                           'Try adjusting bright_lim')
        return comp1, comp2

    def magdiff_plot(self, axs=None):
        """Make a plot of input mag - output mag vs input mag"""
        if not hasattr(self, 'rec'):
            self.completeness(combined_filters=True)
        if axs is None:
            fig, axs = plt.subplots(ncols=2, figsize=(12, 6))

        axs[0].plot(self.mag1[self.rec], self.mag1diff[self.rec], '.',
                    color='k', alpha=0.5)
        axs[1].plot(self.mag2[self.rec], self.mag2diff[self.rec], '.',
                    color='k', alpha=0.5)

        xlab = r'${{\rm Input}}\ {}$'

        axs[0].set_xlabel(xlab.format(self.filter1), fontsize=20)
        axs[1].set_xlabel(xlab.format(self.filter2), fontsize=20)

        axs[0].set_ylabel(r'${{\rm Input}} - {{\rm Ouput}}$', fontsize=20)
        return axs

    def completeness_plot(self, ax=None, comp_fracs=None):
        """Make a plot of completeness vs mag"""
        assert hasattr(self, 'fcomp1'), \
            'need to run completeness with interpolate=True'

        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(self.ast_bins, self.fcomp1(self.ast_bins),
                label=r'${}$'.format(self.filter1))
        ax.plot(self.ast_bins, self.fcomp2(self.ast_bins),
                label=r'${}$'.format(self.filter2))

        if comp_fracs is not None:
            self.add_complines()
        ax.set_xlabel(r'${{\rm mag}}$', fontsize=20)
        ax.set_ylabel(r'${{\rm Completeness\ Fraction}}$', fontsize=20)
        plt.legend(loc='lower left', frameon=False)
        return ax

    def add_complines(self, ax, *fracs, **get_comp_frac_kw):
        """add verticle lines to a plot at given completeness fractions"""
        lblfmt = r'${frac}\ {filt}:\ {comp: .2f}$'
        for frac in fracs:
            ax.hlines(frac, *ax.get_xlim(), alpha=0.5)
            comp1, comp2 = self.get_completeness_fraction(frac,
                                                          **get_comp_frac_kw)
            for comp, filt in zip((comp1, comp2), (self.filter1, self.filter2)):
                lab = lblfmt.format(frac=frac, filt=filt, comp=comp)
                ax.vlines(comp, 0, 1, label=lab,
                          color=next(ax._get_lines.color_cycle))
        plt.legend(loc='lower left', frameon=False)
        return ax


def main(argv):
    parser = argparse.ArgumentParser(description="Calculate completeness fraction, make AST plots")

    parser.add_argument('-c', '--comp_frac', type=float, default=0.9,
                        help='completeness fraction to calculate')

    parser.add_argument('-p', '--makeplots', action='store_true',
                        help='make AST plots')

    parser.add_argument('-m', '--bright_mag', type=float, default=20.,
                        help='brighest mag to consider for completeness frac')

    parser.add_argument('-f', '--plot_fracs', type=str, default=None,
                        help='comma separated completeness fractions to overplot')

    parser.add_argument('fake', type=str, nargs='*', help='match AST file(s)')

    args = parser.parse_args(argv)
    for fake in args.fake:
        ast = ASTs(fake)
        ast.completeness(combined_filters=True, interpolate=True,
                         binsize=0.15)
        comp1, comp2 = ast.get_completeness_fraction(args.comp_frac,
                                                     bright_lim=args.bright_mag)
        print('{} {} completeness fraction:'.format(fake, args.comp_frac))
        print('{0:20s} {1:.4f} {2:.4f}'.format(ast.target, comp1, comp2))

        if args.makeplots:
            comp_name = os.path.join(ast.base, ast.name + '_comp{}'.format(EXT))
            ast_name = os.path.join(ast.base, ast.name + '_ast{}'.format(EXT))

            ax = ast.completeness_plot()
            if args.plot_fracs is not None:
                fracs = map(float, args.plot_fracs.split(','))
                ast.add_complines(ax, *fracs, **{'bright_lim': args.bright_mag})
            plt.savefig(comp_name)
            plt.close()

            ast.magdiff_plot()
            plt.savefig(ast_name)
            plt.close()


if __name__ == "__main__":
    main(sys.argv[1:])
