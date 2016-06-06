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


__all__ = ['check_boundaries', 'strip_header', 'convertz']


def parse_argrange(strarr, arg):
    """Reade a comma separated string into np.arange.
    strarr : string
    e.g.,
    "0.,1,0.1"
    array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9])

    arg:
    if no comma in strarr, return an array of this value.
    """
    if strarr is None:
        return [strarr]

    if ',' in strarr:
        try:
            arr = np.arange(*map(float, strarr.split(',')))
        except ValueError:
            #not floats.
            arr = strarr.split(',')
    else:
        arr = np.array([arg])
    return arr


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
        if upper == lower:
            msg += 'value unchanging'
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
        try:
            dmod0, dmod1, ddmod, av0, av1, dav = np.array(pars.split(), dtype=float)
        except:
            imf, dmod0, dmod1, ddmod, av0, av1, dav = np.array(pars.split(), dtype=float)
            #print(sys.exc_info()[1], param)
            #raise
        retval, msg = betweenie(dmod, dmod1, dmod0, msg=msg)
        retval, msg = betweenie(av, av1, av0, retval=retval, msg=msg)

    if retval > 0:
        print(msg)
    return


def float2sci(num):
    """mpl has a better way of doing this now..."""
    return r'$%s}$' % ('%.0E' % num).replace('E', '0').replace('-0', '^{-').replace('+0', '^{').replace('O', '0')


def strip_header(ssp_file, skip_header=10):
    outfile = ssp_file + '.dat'
    with open(ssp_file, 'r') as infile:
        lines = [l.strip() for l in infile.readlines()]

    try:
        footer, = [i for i, l in enumerate(lines) if 'Best' in l]
    except:
        footer = None
    np.savetxt(outfile, np.array(lines[skip_header:footer], dtype=str), fmt='%s')
    return outfile


def writeorappend(outfile, line):
    """If outfile exists, append line to it, else write line to outfile"""
    wstr = 'w'
    wrote = 'wrote'
    if os.path.isfile(outfile):
        wstr = 'a'
        wrote = 'appended'
    with open(outfile, wstr) as outp:
        outp.write(line)
    print('{} {}'.format(wrote, outfile))
    return


def replaceext(filename, newext):
    '''
    Replace the extension of a filename

    Parameters:
    filename : string
    new_ext : string
        string to replace current extension

    e.g,:
        $ replaceext('data.02.SSS.v4.dat', '.log')
        data.02.SSS.v4.log
    '''
    return splitext(filename)[0] + ext


def splitext(filename):
    '''split the filename from its extension'''
    return '.'.join(filename.split('.')[:-1]), filename.split('.')[-1]

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
