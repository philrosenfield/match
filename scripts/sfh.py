from __future__ import print_function
import logging
import os

import numpy as np
import matplotlib.pyplot as plt

from .config import EXT
from .fileio import read_binned_sfh, parse_pipeline
from .utils import convertz

logger = logging.getLogger()


class SFH(object):
    '''
    load the match sfh solution as a class with attributes set by the
    best fits from the sfh file.
    '''
    def __init__(self, filename, hmc_file=None):
        self.base, self.name = os.path.split(filename)
        self.data = read_binned_sfh(filename, hmc_file)
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
        def float2sci(num):
            """mpl has a better way of doing this now..."""
            return r'$%s}$' % ('%.0E' % num).replace('E', '0').replace('-0', '^{-').replace('+0', '^{').replace('O', '0')

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
            ylabel = r'$SFR\ %s\ (\rm{M_\odot/yr})$' % \
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
                  xlim=(10 ** (10.15-9), 10 **(6.5-9)), ylim=(-0.01, 1.01)):
        '''cumulative sfr plot from match'''
        one_off = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))
            plt.subplots_adjust(right=0.95, left=0.1, bottom=0.1, top=0.95)
            ax.tick_params(direction='in')
            one_off = True

        fill_between_kw = dict({'alpha': 1, 'color': 'gray'}.items() \
                               + fill_between_kw.items())

        plt_kw = dict({'lw': 3}.items() + plt_kw.items())

        lages, (csfh, csfh_errm, csfh_errp) = self.plot_bins(val='csfr',
                                                             err=True)
        age = 10 ** (lages - 9.)
        #age = lages
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
        #ax.set_xscale('log')
        #ax.xaxis.set_major_locator(LogNLocator)
        if one_off:
            ax.set_xlabel('$\\rm{Star\ Formation\ Time\ (Gyr)}$', fontsize=20)
            ax.set_ylabel('$\\rm{Culmulative\ Star\ Formation}$', fontsize=20)
            plt.legend(loc=0, frameon=False)
            if 'label' in plt_kw.keys():
                outfile = '{}_csfr'.format(plt_kw['label'].replace('$', '').lower(), EXT)
            else:
                outfile = '{}_csfr{}'.format(os.path.join(self.base, self.name, EXT))
            plt.savefig(outfile)
            print('wrote {}'.format(outfile))
        return ax

    def param_table(self, angst=True, agesplit=[1e9, 3e9], target='',
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

        fyoung, fyoung_errp, fyoung_errm = self.mass_fraction(0, agesplit[0])
        finter, finter_errp, finter_errm = self.mass_fraction(agesplit[0], agesplit[1])

        # logZ = 0 if there is no SF, that will add error to mean Fe/H
        iyoung = self.nearest_age(agesplit[0], i=False)
        iinter = self.nearest_age(agesplit[1], i=False)

        iyoungs, = np.nonzero(self.data.mh[:iyoung + 1] != 0)
        iinters, = np.nonzero(self.data.mh[:iinter + 1] != 0)
        iinters = list(set(iinters) - set(iyoungs))

        feh_young = convertz(z=0.02 * 10 ** np.mean(self.data.mh[iyoungs]))[-2]
        feh_inter = convertz(z=0.02 * 10 ** np.mean(self.data.mh[iinters]))[-2]
        feh_young_errp = convertz(z=0.02 * 10 ** quadriture(self.data.mh_errp[iyoungs]))[-2]
        feh_young_errm = convertz(z=0.02 * 10 ** quadriture(self.data.mh_errm[iyoungs]))[-2]
        feh_inter_errp = convertz(z=0.02 * 10 ** quadriture(self.data.mh_errp[iinters]))[-2]
        feh_inter_errm = convertz(z=0.02 * 10 ** quadriture(self.data.mh_errm[iinters]))[-2]

        maf = '${0: .2f}^{{+{1: .2f}}}_{{-{2: .2f}}}$'

        d['fyoung'], d['finter'] = [maf.format(v, p, m)
                                    for v,p,m in zip([fyoung, finter],
                                                     [fyoung_errp, finter_errp],
                                                     [fyoung_errm, finter_errm])]
        d['feh_young'], d['feh_inter'] = [maf.format(v, p, m)
                                          for v,p,m in zip([feh_young, feh_inter],
                                                           [feh_young_errp, feh_inter_errp],
                                                           [feh_young_errm, feh_inter_errm])]

        line = ['{target}', '{filters}', '{Av: .2f}', '{dmod: .2f}',
                '{fyoung}', '{feh_young}','{finter}',
                '{feh_inter}', '{bestfit: .1f}']

        d['fmt'] = '%s \\\\ \n' % (' & '.join(line))
        return d

    def nearest_age(self, lage, i=True):
        if lage > 10.15:
            lage = np.log10(lage)
            logger.warning('converting input age to log age')

        age_arr = self.data.lagef
        msg = 'lagef'
        if i:
            age_arr = self.data.lagei
            msg = 'lagei'
        # min age bin size, will trigger warning if ages requested are
        # higher than the min binsize.
        tol = np.min(np.diff(age_arr))

        #  find closest age bin to lage
        idx = np.argmin(np.abs(age_arr - lage))
        difi = np.abs(age_arr[idx] - lage)
        if difi > tol:
            logger.warning('input {}={} not found. Using {}'.format(msg, lage, age_arr[idx]))
        return idx

    def mass_fraction(self, lagei, lagef):
        """
        Return the fraction of total mass formed between lagei and lagef.
        lage[] units can be log yr or yr.
        Multiply by self.totalSF to obtain the mass formed.
        """
        agebins = (10 ** self.data.lagef - 10 ** self.data.lagei)
        if lagef-lagei < np.min(np.diff(self.data.lagei)):
            logger.error('Age difference smaller than bin sizes (or negative)')
            return 0, 0, 0

        # higher precision than self.totalSF
        totalSF = np.sum(self.data.sfr * agebins)
        idxi = self.nearest_age(lagei)
        idxf = self.nearest_age(lagef, i=False) + 1 # +1 is to include final bin

        fracsfr = np.sum(self.data.sfr[idxi:idxf] * agebins[idxi:idxf]) / totalSF
        fracsfr_errp = quadriture(self.data.sfr_errp[idxi:idxf] * agebins[idxi:idxf]) / totalSF
        fracsfr_errm = quadriture(self.data.sfr_errm[idxi:idxf] * agebins[idxi:idxf]) / totalSF

        return fracsfr, fracsfr_errp, fracsfr_errm

def quadriture(x):
    return np.sqrt(np.sum(x * x))
