from __future__ import print_function
import logging
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger()

from .fileio import InputParameters
from astropy.table import Table
#from .. import graphics


__all__ = ['check_boundaries', 'calcsfh_dict', 'call_match', 'grab_val',
           'check_for_bg_file', 'make_calcsfh_param_file', 'strip_header',
           'match_param_default_dict', 'match_param_fmt', 'process_match_sfh',
           'read_binned_sfh', 'read_match_cmd', 'write_match_bg', 'cheat_fake',
           'parse_pipeline', 'convertz']

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


def parse_pipeline(filename):
    '''find target and filters from the filename'''
    name = os.path.split(filename)[1].upper()

    # filters are assumed to be F???W
    starts = np.array([m.start() for m in re.finditer('_F', name)])
    starts += 1
    if len(starts) == 1:
        starts = np.append(starts, starts+6)
    filters = [name[s: s+5] for s in starts]

    # the target name is assumed to be before the filters in the filename
    pref = name[:starts[0]-1]
    for t in pref.split('_'):
        if t == 'IR':
            continue
        try:
            # this could be the proposal ID
            int(t)
        except:
            # a mix of str and int should be the target
            target = t
    return target, filters

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


def grab_val(s, val, v2=None, v3=None):
    def split_str(s, val):
        return float('.'.join(s.split('.%s' % val)[1].split('.')[:2]))
    d = np.nan
    try:
        d = split_str(s, val)
    except:
        if v2 is not None:
            try:
                d = split_str(s, v2)
            except:
                if v3 is not None:
                    d = split_str(s, v3)
    return d


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


def read_binned_sfh(filename):
    '''
    reads the file created using zcombine or HybridMC from match
    into a np.recarray.

    NOTE
    calls genfromtext up to 3 times. There may be a better way to figure out
    how many background lines/what if there is a header... (it's a small file)
    '''
    dtype = [('lagei', '<f8'),
             ('lagef', '<f8'),
             ('dmod', '<f8'),
             ('sfr', '<f8'),
             ('sfr_errp', '<f8'),
             ('sfr_errm', '<f8'),
             ('mh', '<f8'),
             ('mh_errp', '<f8'),
             ('mh_errm', '<f8'),
             ('mh_disp', '<f8'),
             ('mh_disp_errp', '<f8'),
             ('mh_disp_errm', '<f8'),
             ('csfr', '<f8'),
             ('csfr_errp', '<f8'),
             ('csfr_errm', '<f8')]
    try:
        data = np.genfromtxt(filename, dtype=dtype)
    except ValueError:
        try:
            data = np.genfromtxt(filename, dtype=dtype, skip_header=6,
                                 skip_footer=1)
        except ValueError:
            data = np.genfromtxt(filename, dtype=dtype, skip_header=6,
                                 skip_footer=2)
    return data.view(np.recarray)


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

def match_param_default_dict():
    ''' default params for match param file'''

    dd = {'ddmod': 0.05,
          'dav': 0.05,
          'logzmin': -2.3,
          'logzmax': 0.1,
          'dlogz': 0.1,
          'logzmin0': -2.3,
          'logzmax0': -1.0,
          'logzmin1': -1.3,
          'logzmax1': -0.1,
          'BF': 0.35,
          'bad0': 1e-6,
          'bad1': 1e-6,
          'ncmds': 1,
          'Vstep': 0.1,
          'V-Istep': 0.05,
          'fake_sm': 5,
          'nexclude_gates': 0,
          'exclude_gates': '',
          'ninclude_gates': 0,
          'include_gates': ''}

    therest = ['imf', 'dmod1', 'dmod2', 'av1', 'av2', 'V-Imin', 'V-Imax', 'V',
               'I', 'Vmin', 'Vmax', 'Imin', 'Imax']
    for key in therest:
        dd[key] = np.nan
    return dd



def match_param_fmt(set_z=False, zinc=True):
    '''
    calcsfh parameter format, set up for dan's runs and parsec M<12.
    NOTE exclude and include gates are strings and must have a space at
    their beginning.
    '''

    return '''%(imf)s %(dmod1).3f %(dmod2).3f %(ddmod).3f %(av1).3f %(av2).3f %(dav).3f
%(logzmin).2f %(logzmax).2f %(dlogz).2f %(logzmin0).2f %(logzmax0).2f %(logzmin1).2f %(logzmax1).2f
%(BF).2f %(bad0).6f %(bad1).6f
%(ncmds)i
%(Vstep).2f %(V-Istep).2f %(fake_sm)i %(V-Imin).2f %(V-Imax).2f %(V)s,%(I)s
%(Vmin).2f %(Vmax).2f %(V)s
%(Imin).2f %(Imax).2f %(I)s
%(nexclude_gates)i%(exclude_gates)s %(ninclude_gates)i%(include_gates)s
50
6.60 6.70
6.70 6.80
6.80 6.90
6.90 7.00
7.00 7.10
7.10 7.20
7.20 7.30
7.30 7.40
7.40 7.50
7.50 7.60
7.60 7.70
7.70 7.80
7.80 7.90
7.90 8.00
8.00 8.10
8.10 8.20
8.20 8.30
8.30 8.40
8.40 8.50
8.50 8.60
8.60 8.70
8.70 8.75
8.75 8.80
8.80 8.85
8.85 8.90
8.90 8.95
8.95 9.00
9.00 9.05
9.05 9.10
9.10 9.15
9.15 9.20
9.20 9.25
9.25 9.30
9.30 9.35
9.35 9.40
9.40 9.45
9.45 9.50
9.50 9.55
9.55 9.60
9.60 9.65
9.65 9.70
9.70 9.75
9.75 9.80
9.80 9.85
9.85 9.90
9.90 9.95
9.95 10.00
10.00 10.05
10.05 10.10
10.10 10.15
-1 5 -1bg.dat
-1  1 -1
'''


def process_match_sfh(sfhfile, outfile='processed_sfh.out', sarah_sim=False,
                      zdisp=0.):
    '''
    turn a match sfh output file into a sfr-z table for trilegal.

    todo: add possibility for z-dispersion.
    '''

    fmt = '%.6g %.6g %.4g %s \n'

    data = read_binned_sfh(sfhfile)
    sfr = data['sfr']
    # Trilegal only needs populated time bins, not fixed age array
    inds, = np.nonzero(sfr > 0)
    sfr = sfr[inds]
    to = data['lagei'][inds]
    tf = data['lagef'][inds]
    dlogz = data['mh'][inds]
    half_bin = np.diff(dlogz[0: 2])[0] / 2.
    if zdisp > 0:
        zdisp = '%.4g' % (0.02 * 10 ** zdisp)
    else:
        zdisp = ''

    # correct age for trilegal isochrones.
    # with PARSEC V1.1 and V1.2 no need!
    #tf[tf == 10.15] = 10.13

    with open(outfile, 'w') as out:
        for i in range(len(to)):
            if sarah_sim is True:
                z1 = dlogz[i] - half_bin
                z2 = dlogz[i] + half_bin
                sfr[i] /= 2.
            else:
                sfr[i] *= 1e3  # sfr is normalized in trilegal
                # MATCH conversion:
                z1 = 0.02 * 10 ** (dlogz[i] - half_bin)
                z2 = 0.02 * 10 ** (dlogz[i] + half_bin)
            age1a = 1.0 * 10 ** to[i]
            age1p = 1.0 * 10 ** (to[i] + 0.0001)
            age2a = 1.0 * 10 ** tf[i]
            age2p = 1.0 * 10 ** (tf[i] + 0.0001)

            out.write(fmt % (age1a, 0.0, z1, zdisp))
            out.write(fmt % (age1p, sfr[i], z1, zdisp))
            out.write(fmt % (age2a, sfr[i], z2, zdisp))
            out.write(fmt % (age2p, 0.0, z2, zdisp))
            out.write(fmt % (age1a, 0.0, z2, zdisp))
            out.write(fmt % (age1p, sfr[i], z2, zdisp))
            out.write(fmt % (age2a, sfr[i], z1, zdisp))
            out.write(fmt % (age2p, 0.0, z1, zdisp))

    print('wrote', outfile)
    return outfile


def read_match_cmd(filename):
    '''
    reads MATCH .cmd file
    '''
    # mc = open(filename, 'r').readlines()
    # I don't know what the 7th column is, so I call it lixo.
    names = ['mag', 'color', 'Nobs', 'Nsim', 'diff', 'sig', 'lixo']
    cmd = np.genfromtxt(filename, skip_header=4, names=names, invalid_raise=False)
    return cmd


def calcsfh_dict():
    '''
    default dictionary for calcsfh.
    '''
    return {'dmod': 10.,
            'Av': 0.,
            'filter1': None,
            'filter2': None,
            'bright1': None,
            'faint1': None,
            'bright2': None,
            'faint2': None,
            'color': None,
            'mag': None,
            'dmod2': None,
            'colmin': None,
            'colmax': None,
            'Av2': None,
            'imf': 1.30,
            'ddmod': 0.050,
            'dAv': 0.050,
            'logzmin': -2.3,
            'logzmax': 0.1,
            'dlogz': 0.1,
            'zinc': True,
            'bf': 0.35,
            'bad0': 1e-6,
            'bad1': 1e-6,
            'Ncmds': 1,
            'dmag': 0.1,
            'dcol': 0.05,
            'fake_sm': 5,
            'nexclude_gates': 0,
            'exclude_poly': None,
            'ncombine_gates': 0,
            'combine_poly': None,
            'ntbins': 0,
            'dobg': -1,
            'bg_hess': .0,   # neg if it's a .CMD, else it's same fmt as match_phot
            'smooth': 1,
            'ilogzmin': -2.3,
            'ilogzmax': -1.3,
            'flogzmin': -1.9,
            'flogzmax': -1.1,
            'match_bg': ''}


# moved from starpop
def make_match_param(gal, more_gal_kw=None):
    '''
    Make param.sfh input file for match
    see rsp.match_utils.match_param_fmt()

    takes calcsfh search limits to be the photometric limits of the stars in
    the cmd.
    gal is assumed to be angst galaxy, so make sure attr dmod, Av, comp50mag1,
    comp50mag2 are there.

    only set up for acs and wfpc, if other photsystems need to check syntax
    with match filters.

    All values passed to more_gal_kw overwrite defaults.
    '''

    more_gal_kw = more_gal_kw or {}

    # load parameters
    inp = input_parameters(default_dict=match_param_default_dict())

    # add parameteres
    cmin = gal.color.min()
    cmax = gal.color.max()
    vmin = gal.mag1.min()
    imin = gal.mag2.min()

    if 'acs' in gal.photsys:
        V = gal.filter1.replace('F', 'WFC')
        I = gal.filter2.replace('F', 'WFC')
    elif 'wfpc' in gal.photsys:
        V = gal.filter1.lower()
        I = gal.filter2.lower()
    else:
        print(gal.photsys, gal.name, gal.filter1, gal.filter2)

    # default doesn't move dmod or av.
    gal_kw = {'dmod1': gal.dmod, 'dmod2': gal.dmod, 'av1': gal.Av,
              'av2': gal.Av, 'V': V, 'I': I, 'Vmax': gal.comp50mag1,
              'Imax': gal.comp50mag2, 'V-Imin': cmin, 'V-Imax': cmax,
              'Vmin': vmin, 'Imin': imin}

    # combine sources of params
    phot_kw = dict(match_param_default_dict().items() \
                   + gal_kw.items() + more_gal_kw.items())

    inp.add_params(phot_kw)

    # write out
    inp.write_params('param.sfh', match_param_fmt())
    return inp

def is_numeric(lit):
    """
    value of numeric: literal, string, int, float, hex, binary
    From http://rosettacode.org/wiki/Determine_if_a_string_is_numeric#Python
    """
    # Empty String
    if len(lit) <= 0:
        return lit
    # Handle '0'
    if lit == '0':
        return 0
    # Hex/Binary
    if len(lit) > 1:  # sometimes just '-' means no data...
        litneg = lit[1:] if lit[0] == '-' else lit
        if litneg[0] == '0':
            if litneg[1] in 'xX':
                return int(lit, 16)
            elif litneg[1] in 'bB':
                return int(lit, 2)
            else:
                try:
                    return int(lit, 8)
                except ValueError:
                    pass
    # Int/Float/Complex
    try:
        return int(lit)
    except ValueError:
        pass
    try:
        return float(lit)
    except ValueError:
        pass
    try:
        return complex(lit)
    except ValueError:
        pass
    return lit
