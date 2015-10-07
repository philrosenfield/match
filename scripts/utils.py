from __future__ import print_function
import logging
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger()

from astropy.table import Table
from .fileio import read_match_cmd, read_binned_sfh
#from .. import graphics


__all__ = ['check_boundaries', 'grab_val',
           'strip_header', 'MatchCMD', 'MatchSFH',
           'convertz']

def convertz(z=None, oh=None, mh=None, feh=None, oh_sun=8.76, z_sun=0.01524,
             y0=.2485, dy_dz=1.80):
    '''
    input:
    metallicity as z
    [O/H] as oh
    [M/H] as mh
    [Fe/H] as feh

    initial args can be oh_sun, z_sun, y0, and dy_dz

    returns oh, z, y, x, feh, mh where y = He and X = H mass fractions
    '''

    if oh is not None:
        feh = oh - oh_sun
        z = z_sun * 10 **(feh)

    if mh is not None:
        z = (1 - y0) / ((10**(-1. * mh) / 0.0207) + (1. + dy_dz))

    if z is not None:
        feh = np.log10(z / z_sun)

    if feh is not None:
        z = z_sun * 10**feh

    oh = feh + oh_sun
    y = y0 + dy_dz * z
    x = 1. - z - y
    if mh is None:
        mh = np.log10((z / x) / 0.0207)

    if __name__ == "__main__":
        print('''
                 [O/H] = %2f
                 z = %.4f
                 y = %.4f
                 x = %.4f
                 [Fe/H] = %.4f
                 [M/H] = %.4f''' % (oh, z, y, x, feh, mh))
    return oh, z, y, x, feh, mh


def check_boundaries(param, scrn):
    """
    check if the best fit file hits the Av or dmod edge of parameter search
    space.

    print information to terminal, nothing is printed if dmod and Av are
    within bounds.

    Parameters
    ----------
    param : match parameter file
    scrn : match console output (saved as a file)
    """
    def betweenie(val, upper, lower, retval=0, msg=''):
        if val >= upper:
            msg += 'extend boundary higher, %f >= %f\n' % (val, upper)
            retval += 1
        if val <= lower:
            msg += 'extend boundary lower, %f <= %f\n' % (val, lower)
            retval += 1
        return retval, msg

    msg = '{} / {}\n'.format(os.path.split(param)[1], os.path.split(scrn)[1])
    # parse scrn
    bfit = open(scrn).readlines()[-1]
    if not 'Best' in bfit:
        msg += 'error calcsfh not finished'
        retval = 1
    else:
        av, dmod, _ = bfit.split(',')
        dmod = float(dmod.replace(' dmod=', ''))
        av = float(av.replace('Best fit: Av=', ''))

        # parse param
        pars = open(param).readline()
        dmod0, dmod1, ddmod, av0, av1, dav = np.array(pars.split(), dtype=float)

        retval, msg = betweenie(dmod, dmod1, dmod0, msg=msg)
        retval, msg = betweenie(av, av1, av0, retval=retval, msg=msg)

    if retval > 0:
        print(msg)
    return


def make_matchfake(fname):
    """
    make four-column Vin, Iin, Vdiff, Idiff artificial stars

    made to work with pipeline fake.fits files with 3 filters, should work
    for two but not tested
    assumes _F*W_F*W_ etc in the file name label the filters.
    """
    try:
        tbl = Table.read(fname, format='fits')
    except:
        logger.error('problem with {}'.format(fname))
        return
    filters = [f for f in fname.split('_') if f.startswith('F')]
    pref = fname.split('F')[0]
    sufx = fname.split('W')[-1].replace('fits', 'matchfake')
    for i in range(len(filters)-1):
        mag1in_col = 'MAG{}IN'.format(i+1)
        mag2in_col = 'MAG{}IN'.format(len(filters))

        if mag1in_col == mag2in_col:
            continue

        mag1out_col = 'MAG{}OUT'.format(i+1)
        mag2out_col = 'MAG{}OUT'.format(len(filters))

        try:
            mag1in = tbl[mag1in_col]
            mag2in = tbl[mag2in_col]
            mag1diff = tbl[mag1in_col] - tbl[mag1out_col]
            mag2diff = tbl[mag2in_col] - tbl[mag2out_col]
        except:
            logger.error('problem with column formats in {}'.format(fname))
            return

        fout = pref + '{}_{}'.format(filters[i],filters[-1]) + sufx
        np.savetxt(fout, np.column_stack((mag1in, mag2in, mag1diff, mag2diff)),
                   fmt='%.4f')
        logger.info('wrote {}'.format(fout))
    return


def float2sci(num):
    return r'$%s}$' % ('%.0E' % num).replace('E', '0').replace('-0', '^{-').replace('+0', '^{').replace('O', '0')


def strip_header(ssp_file, skip_header=10):
    outfile = ssp_file + '.dat'
    with open(ssp_file, 'r') as infile:
        lines = [l.strip() for l in infile.readlines()]
    np.savetxt(outfile, np.array(lines[skip_header:], dtype=str), fmt='%s')


def cheat_fake(infakefile, outfakefile):
    """
    Increase the resolution of a match fake file by repeating the
    fake star entry with slight shifts within the same hess diagram
    cell of mag2.

    Parameters
    ----------
    infakefile, outfakefile : str, str
        input and output fake file names
    """
    # infake format is mag1in, mag2in, mag1idff, mag2diff
    infake = np.loadtxt(infakefile)

    offsets = [0.06, 0.03, -0.03, -0.06]

    outfake = np.copy(infake)
    for offset in offsets:
        tmp = np.copy(infake)
        tmp.T[1] += offset
        outfake = np.concatenate([outfake, tmp])
    np.savetxt(outfakefile, outfake, '%.3f')
    return


class MatchCMD(object):
    """
    A quikly made object to read the MATCH CMD file and hold paramters to
    automatically make plots with the same color scale as other MATCH CMD files.
    """
    def __init__(self, filename):
        self.cmd = read_match_cmd(filename)
        self.figname = os.path.split(filename)[1] + '.png'
        labels = ['${\\rm %s}$' % i for i in ('data', 'model', 'diff', 'sig')]
        labels[1] = '${\\rm %s}$' % self.figname.split('.')[0].replace('_', '\ ')
        self.labels = labels
        self.load_match_cmd(filename)

    def load_match_cmd(self, filename):
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
        self.extent = [self.cmd['color'][0], self.cmd['color'][-1],
                       self.cmd['mag'][-1], self.cmd['mag'][0]]
        self.max_counts = np.nanmax(np.concatenate([self.data, self.model]))
        self.max_diff = np.nanmax(np.abs(self.diff))
        self.max_sig = np.nanmax(np.abs(self.sig))


class MatchSFH(object):
    '''
    load the match sfh solution as a class with attributes set by the
    best fits from the sfh file.
    '''
    def __init__(self, filename):
        self.base, self.name = os.path.split(filename)
        self.data = read_binned_sfh(filename)
        self.load_match_header(filename)

    def load_match_header(self, filename):
        '''
        assumes header is from line 0 to 6 and sets footer to be the final
        line of the file

        header formatting is important:
        Line # format requirement
        first  Ends with "= %f (%s)"
        N      is the string "Best fit:\n"
        N+1    has ',' separated strings of "%s=%f+%f-%f"
        last   is formatted "%s %f %f %f"
        '''
        def set_value_err_attr(key, attr, pattr, mattr):
            '''
            set attributes [key], [key]_perr, [key]_merr
            to attr, pattr, mattr (must be floats)
            '''
            self.__setattr__(key, float(attr))
            self.__setattr__(key + '_perr', float(pattr))
            self.__setattr__(key + '_merr', float(mattr))

        with open(filename, 'r') as infile:
            lines = infile.readlines()

        if len(lines) == 0:
            print('empty file: %s' % filename)
            self.header = []
            self.footer = []
            self.bestfit = np.nan
            self.match_out = ''
            self.data = []
            return

        self.header = lines[0:6]
        self.footer = lines[-1]
        bestfit, fout = self.header[0].replace(' ', '').split('=')[1].split('(')
        self.bestfit = float(bestfit)
        self.match_out = fout.split(')')[0]

        try:
            iline = self.header.index('Best fit:\n') + 1
        except ValueError:
            print('Need Best fit line to assign attributes')
            raise ValueError

        line = self.header[iline].strip().replace(' ', '').split(',')
        for l in line:
            key, attrs = l.split('=')
            attr, pmattr = attrs.split('+')
            pattr, mattr = pmattr.split('-')
            set_value_err_attr(key, attr, pattr, mattr)
        # the final line has totalSF
        key, attr, pattr, mattr = self.header[-1].strip().split()
        set_value_err_attr(key, attr, pattr, mattr)

        self.flag = None
        if np.sum(np.diff(self.data.mh)) == 0:
            self.flag = 'setz'
        if len(np.nonzero(np.diff(self.data.mh) >= 0)[0]) == len(self.data.mh):
            self.flag = 'zinc'
        return

    def mh2z(self, num):
        return 0.02 * 10 ** num

    def plot_bins(self, val='sfr', err=False, convertz=False, offset=1.):
        '''make SFH bins for plotting'''
        if type(val) == str:
            if err:
                #import pdb; pdb.set_trace()
                valm = self.data['%s_errm' % val] * offset
                valp = self.data['%s_errp' % val] * offset
            val = self.data[val] * offset
            if convertz:
                val = self.mh2z(val)
                if err:
                    valm = self.mh2z(valm)
                    valp = self.mh2z(valp)
        lagei = self.data.lagei
        lagef = self.data.lagef

        # double up value
        # lagei_i, lagef_i, lagei_i+1, lagef_i+1 ...
        lages = np.ravel([(lagei[i], lagef[i]) for i in range(len(lagei))])
        vals = np.ravel([(val[i], val[i]) for i in range(len(val))])
        if err:
            valm = np.ravel([(valm[i], valm[i]) for i in range(len(val))])
            valp = np.ravel([(valp[i], valp[i]) for i in range(len(val))])
            data = (vals, valm, valp)
        else:
            data = vals
        return lages, data

    def age_plot(self, val='sfr', ax=None, plt_kw={}, errors=True,
                 convertz=False, xlabel=None, ylabel=None,
                 sfr_offset=1e3):
        plt_kw = dict({'lw': 3, 'color': 'black'}.items() + plt_kw.items())
        eplt_kw = plt_kw.copy()
        eplt_kw.update({'linestyle': 'None'})

        lages, sfrs = self.plot_bins(offset=sfr_offset)
        rlages, (rsfrs, sfr_merrs, sfr_perrs) = self.plot_bins(err=True,
                                                               offset=sfr_offset)

        if val != 'sfr':
            lages, vals = self.plot_bins(val=val, convertz=convertz)
            # mask values with no SF
            isfr, = np.nonzero(sfrs==0)
            vals[isfr] = np.nan
            if self.flag != 'setz':
                rlages, (rvals, val_merrs, val_perrs) = self.plot_bins(val=val,
                                                                       err=True)
                # mask values with no SF
                irsfr, = np.nonzero(rsfrs==0)
                val_merrs[irsfr] = 0.
                val_perrs[irsfr] = 0.
                if np.sum(val_merrs) == 0 or np.sum(val_perrs) == 0:
                    errors = False
            else:
                errors = False
            if 'mh' in val:
                if ylabel is not None:
                    ylabel = r'$\rm{[M/H]}$'
                if convertz:
                    ylabel = r'$Z$'
        else:
            ylabel = r'$SFR\ %s (\rm{M_\odot/yr})$' % \
                     float2sci(1./sfr_offset).replace('$','')
            vals = sfrs
            rvals = rsfrs
            val_merrs = sfr_merrs
            val_perrs = sfr_perrs
        if ax is None:
            fig, ax = plt.subplots()
            xlabel = r'$\log Age\ \rm{(yr)}$'

        ax.plot(lages, vals, **plt_kw)
        if errors:
            ax.errorbar(rlages, rvals, yerr=[val_merrs, val_perrs], **eplt_kw)

        if xlabel is not None:
            ax.set_xlabel(xlabel, fontsize=20)
        if ylabel is not None:
            ax.set_ylabel(ylabel, fontsize=20)
        return ax

    def plot_csfr(self, ax=None, errors=True, plt_kw={}, fill_between_kw={},
                  xlim=(13.4, -0.01), ylim=(-0.01, 1.01)):
        '''cumulative sfr plot from match'''
        one_off = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
            plt.subplots_adjust(right=0.95, left=0.1, bottom=0.11, top=0.95)
            one_off = True

        fill_between_kw = dict({'alpha': 0.2, 'color': 'gray'}.items() \
                               + fill_between_kw.items())

        plt_kw = dict({'lw': 3}.items() + plt_kw.items())

        lages, (csfh, csfh_errm, csfh_errp) = self.plot_bins(val='csfr',
                                                             err=True)
        age = 10 ** (lages - 9.)

        age = np.append(age, age[-1])
        csfh = np.append(csfh, 0)
        csfh_errm = np.append(csfh_errm, 0)
        csfh_errp = np.append(csfh_errp, 0)

        if errors:
            ax.fill_between(age, csfh - csfh_errm, csfh + csfh_errp,
                            **fill_between_kw)

        ax.plot(age, csfh, **plt_kw)

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if one_off:
            ax.set_xlabel('$\\rm{Time\ (Gyr)}$', fontsize=20)
            ax.set_ylabel('$\\rm{Culmulative\ SF}$', fontsize=20)
            plt.legend(loc=0, frameon=False)
            if 'label' in plt_kw.keys():
                outfile = '{}_csfr.png'.format(plt_kw['label'].replace('$', '').lower())
            else:
                outfile = '{}_csfr.png'.format(os.path.join(self.base, self.name))
            plt.savefig(outfile)
            print('wrote {}'.format(outfile))
        return ax

    def param_table(self, angst=True, agesplit=[1., 3.], target='',
                    filters=['','']):
        d = {'bestfit': self.bestfit, 'Av': self.Av, 'dmod': self.dmod}

        d['header'] = \
            (r'Galaxy & Optical Filters & A$_V$ & $(m\!-\!M)_0$ &'
             r'$\% \frac{{\rm{{SF}}}}{{\rm{{SF_{{TOT}}}}}}$ &'
             r'$\langle \mbox{{[Fe/H]}} \rangle$ &'
             r'$\% \frac{{\rm{{SF}}}}{{\rm{{SF_{{TOT}}}}}}$ &'
             r'$\langle \mbox{{[Fe/H]}} \rangle$ & $bestfit$ \\ & & & & '
             r'\multicolumn{{2}}{{c}}{{$<{0}\rm{{Gyr}}$}} & '
             r'\multicolumn{{2}}{{c}}{{${0}-{1}\rm{{Gyr}}$}} & \\ \hline'
             '\n'.format(*agesplit))

        if angst:
            d['target'], filters = parse_pipeline(self.name)
        else:
            d['target'] = target

        d['filters'] = ','.join(filters)

        iyoung = np.argmin(abs(agesplit[0] - 10 **(self.data.lagef - 8)))
        iinter = np.argmin(abs(agesplit[1] - 10 **(self.data.lagef - 8)))

        sf = self.data['sfr'] * \
            (10 ** self.data['lagef'] - 10 ** self.data['lagei'])
        fcsf = np.cumsum(sf)/np.sum(sf)

        d['fyoung'] = 100 * fcsf[iyoung]
        d['finter'] = 100 * fcsf[iinter] - fcsf[iyoung]

        # logZ = 0 if there is no SF, that will add error to mean Fe/H
        iyoungs, = np.nonzero(self.data.mh[:iyoung + 1] != 0)
        iinters, = np.nonzero(self.data.mh[:iinter + 1] != 0)
        iinters = list(set(iinters) - set(iyoungs))

        d['feh_young'] = convertz(z=0.02 * 10 ** np.mean(self.data.mh[iyoungs]))[-2]
        d['feh_inter'] = convertz(z=0.02 * 10 ** np.mean(self.data.mh[iinters]))[-2]

        line = ['{target}', '{filters}', '{Av: .2f}', '{dmod: .2f}',
                '{fyoung: .2f}', '{feh_young: .2f}', '{finter: .2f}',
                '{feh_inter: .2f}', '{bestfit: .1f}']

        d['fmt'] = '%s \\\\ \n' % (' & '.join(line))
        return d

    def mass_fraction(self, lagei, lagef):
        """
        Return the fraction of total mass formed between lagei and lagef.
        lage[] units can be log yr or yr.
        Multiply by self.totalSF to obtain the mass formed.
        """
        if lagei > 1e6:
            lagei = np.log10(lagei)
            print('Warning: converting input age to log age')
        if lagef > 1e6:
            lagef = np.log10(lagef)
            print('Warning: converting input age to log age')

        # min age bin size, will trigger warning if ages requested are
        # higher than the min binsize.
        tol = np.min(np.diff(self.data.lagei))

        agebins = (10 ** self.data.lagef - 10 ** self.data.lagei)

        # higher precision than self.totalSF
        totalSF = np.sum(self.data.sfr * agebins)

        #  find closest age bin to lagei
        idxi = np.argmin(np.abs(self.data.lagei - lagei))
        difi = np.abs(self.data.lagei[idxi] - lagei)
        if difi > tol:
            print('Warning: input lagei={} not found. Using {}'.format(lagei, self.data.lagei[idxi]))

        #  find closest age bin to lagef
        idxf = np.argmin(np.abs(self.data.lagef - lagef))
        dif = np.abs(self.data.lagef[idxf] - lagef)
        if dif > tol:
            print('Warning: input lagef={} not found using {}'.format(lagef, self.data.lagef[idxf]))

        fracsfr = np.sum(self.data.sfr[idxi:idxf + 1]  * agebins[idxi:idxf + 1])# +1 to include final bin
        return fracsfr / totalSF

