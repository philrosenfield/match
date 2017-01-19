import numpy as np
from ..fileio import read_calcsfh_param
from ..config import EXT
import matplotlib.pyplot as plt
try:
    import seaborn
    seaborn.set()
except:
    pass

def match_diagnostic(param, phot, fake=None, save=True, xlim=None, ylim=None,
                     good=True):
    """
    make two panel cmd figure (ymag = mag1, mag2)
    from match photometery and cmd parameter space from the match param drawn
    """
    moff = 0.1
    off = 0.1
    mag1, mag2 = np.loadtxt(phot, unpack=True)
    if good:
        im1, = np.nonzero(mag1 < 40)
        im2, = np.nonzero(mag2 < 40)
        inds = np.intersect1d(im1, im2)
        mag1 = mag1[inds]
        mag2 = mag2[inds]

    color = mag1 - mag2


    params = read_calcsfh_param(param)

    cmin, cmax = params['vimin'], params['vimax']
    m1min, m1max = params['vmin'], params['vmax']
    m2min, m2max = params['imin'], params['imax']

    filters = [params['v'], params['i']]

    verts = [np.array([[cmin, m1min], [cmin, m1max], [cmax, m1max],
                       [cmax, m1min], [cmin, m1min]]),
             np.array([[cmin, m2min], [cmin, m2max], [cmax, m2max],
                       [cmax, m2min], [cmin, m2min]])]

    magcuts, = np.nonzero((mag1 < m1max) & (mag1 > m1min) &
                          (mag2 < m2max) & (mag2 > m2min))

    _, axs = plt.subplots(ncols=2, sharex=True, figsize=(12, 6))

    for i, ymag in enumerate([mag1, mag2]):
        axs[i].plot(color, ymag, '.', label='all phot')
        axs[i].plot(color[magcuts], ymag[magcuts], '.', label='mag limits')
        axs[i].plot(verts[i][:, 0], verts[i][:, 1], label='param limits')
        axs[i].set_ylabel(r'${}$'.format(filters[i]))
        axs[i].set_xlabel(r'${}-{}$'.format(*filters))
        axs[i].invert_yaxis()

    if fake is not None:
        mag1in, mag2in, _, _ = np.loadtxt(fake, unpack=True)
        colin = mag1in - mag2in
        for i, ymag in enumerate([mag1in, mag2in]):
            fverts = np.array([[colin.min(), ymag.min()],
                               [colin.min(), ymag.max()],
                               [colin.max(), ymag.max()],
                               [colin.max(), ymag.min()],
                               [colin.min(), ymag.min()]])
            axs[i].plot(fverts[:, 0], fverts[:, 1], label='fake limits')
            plt.legend(loc='best')

    if xlim is not None:
        try:
            [ax.set_xlim(xlim) for ax in axs]
        except:
            [ax.set_xlim(parse_limits(xlim)) for ax in axs]

    if ylim is not None:
        try:
            [ax.set_ylim(ylim) for ax in axs]
        except:
            [ax.set_ylim(parse_limits(ylim)) for ax in axs]

    if save:
        plt.savefig(param + EXT)
        print('wrote', param + EXT)
        plt.close()
    return axs

def parse_limits(lim):
    return np.array(lim.split(), dtype=float)
