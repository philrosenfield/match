from __future__ import print_function
import logging
import os

import numpy as np
import matplotlib.pyplot as plt

from .fileio import read_binned_sfh

logger = logging.getLogger()



class SFH(object):
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
            logger.warning('converting input age to log age')
        if lagef > 1e6:
            lagef = np.log10(lagef)
            logger.warning('converting input age to log age')

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
            logger.warning('input lagei={} not found. Using {}'.format(lagei, self.data.lagei[idxi]))

        #  find closest age bin to lagef
        idxf = np.argmin(np.abs(self.data.lagef - lagef))
        dif = np.abs(self.data.lagef[idxf] - lagef)
        if dif > tol:
            logger.warning('input lagef={} not found using {}'.format(lagef, self.data.lagef[idxf]))

        fracsfr = np.sum(self.data.sfr[idxi:idxf + 1]  * agebins[idxi:idxf + 1])# +1 to include final bin
        return fracsfr / totalSF

