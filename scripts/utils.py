"""Utility functions for match.scripts"""
from __future__ import print_function
import logging
import os
import re
import sys

import numpy as np

try:
    from .fileio import read_binned_sfh
except SystemError:
    from fileio import read_binned_sfh

logger = logging.getLogger()


__all__ = ['strip_header', 'convertz', 'center_grid', 'quantiles',
           'fitgauss1D']


def quantiles(ux, prob, qs=[0.16, 0.84], res=200, maxp=False,
              ax=None, k=3, interpolate=True):
    """
    Calculate quantiles, or fraction of total area under prob curve.

    Parameters
    ----------
    ux : 1D array
        vector assciated with prob
    prob : 1D array
        probability of ux
    qs : list or array
        quantiles
    res :
        number of resolution points for interpolation
    maxp : bool (False)
        append maximum posterior probability to returned quantile array
    ax : plt.Axes instance
        plot the interpolation
    Returns
    -------
    interpolated ux evaluated at qs
    """
    if interpolate:
        from scipy.interpolate import splprep, splev
        ((tckp, u), fp, ier, msg) = splprep([ux, prob], k=k, full_output=1)

        if ier != 0:
            k = 1

            ((tckp, u), fp, ier, msg) = splprep([ux, prob], k=k, full_output=1)
        iux, iprob = splev(np.linspace(0, 1, res), tckp)
    else:
        iux, iprob = ux, prob

    fac = np.cumsum(iprob) / np.sum(iprob)
    ipts = [np.argmin(np.abs(fac - q)) for q in qs]
    if maxp:
        ipts.append(np.argmax(iprob))
    g = iux[ipts]
    if ax is not None:
        # useful for debugging or by-eye checking of interpolation
        ax.plot(iux, iprob, color='r')
        ax.plot(iux[ipts], iprob[ipts], 'o', color='r')
    print('quantiles {} are actually at {}'.format(qs, fac[ipts]))
    return g


def fitgauss1D(ux, prob):
    """
    Fit a 1D Gaussian to a marginalized probability.

    Parameters

    xattr : str
        column name (and will be attribute name)
    ux :
        unique xattr values
    prob :
        marginalized probability at ux.

    Returns

    g : astropy.models.Gaussian1D object
    sets g as attribute 'xattr'g
    """
    assert ux is not None, \
        'need to supply values and probability to fitgauss1D'
    from astropy.modeling import models, fitting
    weights = np.ones(len(ux))
    fit_g = fitting.LevMarLSQFitter()
    # g_init = models.Gaussian1D(amplitude=1.,
    #                           mean=np.mean(ux),
    #                           stddev=np.diff(ux)[0])
    deg = 4
    if deg == len(ux):
        deg = 2
    g_init = models.Polynomial1D(degree=deg)
    noweight, = np.nonzero(prob == 2 * np.log(1e-323))
    weights[noweight] = 0.
    return fit_g(g_init, ux, prob)


def lnprob(prob):
    lenp = len(prob)
    prob[prob == 0] = 1e-323
    prob = 2 * np.log(prob)
    pmin = np.min(np.ravel(prob)[(np.isfinite(np.ravel(prob)))])
    # prob[np.isneginf(prob)] = pmin
    prob -= pmin
    assert len(prob) == lenp, 'lnprob has changed array shape of prob.'
    return prob


def marg(x, z, unx=None, log=True):
    """
    marginalize in 1d.
    Does not normalize probability.
    z should be ssp.absprob or some linear probability
    NOTE: this simple code only works on a grid with equally spaced data
    """
    ux = unx
    if ux is None:
        ux = np.unique(x)
    prob = np.zeros(len(ux))
    for i in range(len(ux)):
        iz, = np.nonzero(x == ux[i])
        prob[i] = np.sum(z.iloc[iz])
    if log:
        prob = lnprob(prob)
    return prob, ux


def marg2d(x, y, z, unx=None, uny=None, log=True):
    """
    marginalize in 2d.
    Does not normalize probability.
    z should be ssp.absprob or some linear probability
    NOTE: this simple code only works on a grid with equally spaced data
    """
    ux = unx
    if ux is None:
        ux = np.unique(x)
    uy = uny
    if uy is None:
        uy = np.unique(y)
    prob = np.zeros(shape=(len(ux), len(uy)))
    models = np.array([])
    for i in range(len(ux)):
        for j in range(len(uy)):
            iz, = np.nonzero((x == ux[i]) & (y == uy[j]))
            models = np.append(models, len(iz))
            if len(iz) > 0:
                prob[i, j] = np.sum(z.iloc[iz])
    unm = np.unique(models)
    if len(unm) != 1:
        print('Warning: Uneven grid')
        print(unm)
    if log:
        prob = lnprob(prob)
    return prob, ux, uy


def centered_meshgrid(x, y, unx=None, uny=None):
    """call meshgrid with bins shifted so x, y will be at bin center"""
    X, Y = np.meshgrid(center_grid(x, unx=unx),
                       center_grid(y, unx=uny),
                       indexing="ij")
    return X, Y


def center_grid(a, unx=None):
    """uniquify and shift a uniform array half a bin maintaining its size"""
    x = unx
    if x is None:
        x = np.unique(a)
    dx = np.diff(x)[0]
    x = np.append(x, x[-1] + dx)
    x -= dx / 2
    return x


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
    # tf[tf == 10.15] = 10.13

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
    for pre in pref.split('_'):
        if pre == 'IR':
            continue
        try:
            # this could be the proposal ID
            int(pre)
        except:
            # a mix of str and int should be the target
            target = pre
    return target, filters


def parse_argrange(strarr):
    """
    Read a comma separated string or float list into np.arange or np.array.

    If the final value of the list is less than the others, it will assume
    min, max, delta.

    strarr : string or list
    e.g.,
    "0.,1,0.1"
    array([ 0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9])

    arg:
    if no comma in strarr, return an array of this value.
    """
    def arangeit(arr):
        """
        return np.arange(*arr) if there are 3 values in arr
        and last value is less than others.
        """
        if len(arr) == 3 and np.sign(np.diff(arr))[-1] < 1:
            arng = np.arange(*arr)
        else:
            arng = arr
        return arng

    if strarr is None:
        return np.array([strarr])

    if isinstance(strarr, str):
        if ',' in strarr:
            try:
                farr = np.array(strarr.split(','), dtype=float)
                arr = arangeit(farr)
            except ValueError:
                # not floats.
                arr = np.array(strarr.split(','))
        else:
            # not a comma separated list
            arr = np.array([strarr])
    else:
        arr = arangeit(np.array(strarr, dtype=float))
    return arr


def convertz(z=None, oh=None, mh=None, feh=None, oh_sun=8.76, z_sun=0.01524,
             y0=.2485, dy_dz=1.80, y=None):
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
        z = z_sun * 10 ** (feh)

    if mh is not None:
        z = (1 - y0) / ((10**(-1. * mh) / 0.0207) + (1. + dy_dz))

    if z is not None:
        feh = np.log10(z / z_sun)

    if feh is not None:
        z = z_sun * 10**feh

    oh = feh + oh_sun
    if y is None:
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


def float2sci(num):
    """mpl has a better way of doing this?"""
    _, exnt = '{0:.0e}'.format(num).split('e')
    exnt = int(exnt)
    if exnt == 0:
        # 10 ** 0 = 1
        retv = ''
    else:
        retv = r'$10^{{{0:d}}}$'.format(exnt)
    return retv


def strip_header(ssp_file, skip_header=10, return_lines=False):
    outfile = ssp_file + '.dat'
    with open(ssp_file, 'r') as infile:
        lines = [l.strip() for l in infile.readlines()]

    try:
        footer, = [i for i, l in enumerate(lines) if 'Best' in l]
    except:
        footer = None
    lines = lines[skip_header:footer]
    if return_lines:
        retv = lines
    else:
        np.savetxt(outfile, np.array(lines, dtype=str), fmt='%s')
        retv = outfile
    return retv


def writeorappend(outfile, line):
    """If outfile exists, append line to it, else write line to outfile"""
    wstr = 'w'
    wrote = 'wrote'
    if os.path.isfile(outfile):
        wstr = 'a'
        wrote = 'appended'
    with open(outfile, wstr) as outp:
        outp.write(line)
    print('{0:s} {1:s}'.format(wrote, outfile))
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
    return splitext(filename)[0] + newext


def splitext(filename):
    '''split the filename from its extension'''
    return '.'.join(filename.split('.')[:-1]), filename.split('.')[-1]


def ensure_file(f, mad=True):
    '''
    input
    f (string): if f is not a file will print "no file"
    optional
    mad (bool)[True]: if mad is True, will exit program.
    '''
    test = os.path.isfile(f)
    if test is False:
        logger.warning('{0:s} not found'.format(f))
        if mad:
            sys.exit()
    return test


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
