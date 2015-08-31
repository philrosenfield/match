from __future__ import print_function
import os
#import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.offsetbox import AnchoredText
from matplotlib.patheffects import withStroke

from .utils import MatchSFH, parse_pipeline
from .fileio import get_files

__all__ = ['add_inner_title', 'forceAspect', 'match_plot', 'pgcmd', 'sfh_plot',
           'read_match_cmd', 'MatchCMD', 'match_diagnostic', 'call_pgcmd']


def match_diagnostic(param, phot):
    """
    make two panel cmd figure (ymag = mag1, mag2)
    from match photometery and cmd parameter space from the match param drawn
    """
    #mc = 0
    moff = 0.1
    off = 0.1
    mag1, mag2 = np.loadtxt(phot, unpack=True)
    color = mag1 - mag2

    with open(param, 'r') as f:
        paramlines = f.readlines()

    colmin, colmax = map(float, paramlines[4].split()[3: 5])
    mag1max, mag1min = map(float, paramlines[5].split()[0: 2])
    mag2max, mag2min = map(float, paramlines[6].split()[0: 2])

    ylims = [(mag1min + moff, mag1max - moff), (mag2min + moff, mag2max - moff)]

    filters = paramlines[4].split()[-1].split(',')
    excludes = np.array([paramlines[7].split()], dtype=float)

    nregions = excludes[0][0]
    if nregions == 1:
        vs = np.array(paramlines[7].split()[1:-1], dtype=float)
        exgverts = np.append(vs, vs[:2]).reshape(5, 2)

    verts = [np.array([[colmin, mag1min], [colmin, mag1max], [colmax, mag1max],
                       [colmax, mag1min], [colmin, mag1min]]),
             np.array([[colmin, mag2min], [colmin, mag2max], [colmax, mag2max],
                       [colmax, mag2min], [colmin, mag2min]])]

    magcuts = [np.nonzero((mag1 > mag1max) & (mag1 < mag1min)),
               np.nonzero((mag2 > mag2max) & (mag2 < mag2min))]

    fig, axs = plt.subplots(ncols=2, sharex=True, figsize=(12, 6))

    for i, ymag in enumerate([mag1, mag2]):
        axs[i].plot(color, ymag, '.')
        axs[i].plot(color[magcuts[i]], ymag[magcuts[i]], '.')
        axs[i].plot(verts[i][:, 0], verts[i][:, 1])
        axs[i].set_ylabel(r'${}$'.format(filters[i]))
        axs[i].set_xlabel(r'${}-{}$'.format(*filters))
        axs[i].set_ylim(ylims[i])
        axs[i].set_xlim(colmin - off, colmax + off)
    if nregions > 0:
        axs[0].plot(exgverts[:, 0] , exgverts[:, 1])
        # mag2 = mag1 - color
        axs[1].plot(exgverts[:, 0] , exgverts[:, 1] - exgverts[:, 0])

    plt.savefig(param + '.png')
    print('wrote', param + '.png')
    plt.close()
    return axs


def add_inner_title(ax, title, loc, size=None, **kwargs):
    '''
    adds a title to an ax inside to a location loc, which follows plt.legends loc ints.
    '''
    if size is None:
        size = dict(size=plt.rcParams['legend.fontsize'])
    anct = AnchoredText(title, loc=loc, prop=size, pad=0., borderpad=0.5,
                        frameon=False, **kwargs)
    ax.add_artist(anct)
    anct.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])
    return anct


def match_plot(ZS, extent, labels=["Data", "Model", "Diff", "Sig"],
               imagegrid_kw={}, **kwargs):
    '''
    ex ZS = [h[2],sh[2],diff_cmd,resid]
    '''
    max_diff = kwargs.get('max_diff')
    max_counts = kwargs.get('max_counts')
    max_sig = kwargs.get('max_sig')

    fig = plt.figure(figsize=(9, 9))
    defaults = {'nrows_ncols': (2, 2),
                'direction': "row",
                'axes_pad': .7,
                'add_all': True,
                'label_mode': "all",
                'share_all': True,
                'cbar_location': "top",
                'cbar_mode': "each",
                'cbar_size': "7%",
                'cbar_pad': "2%",
                'aspect': 'square'}

    imagegrid_kw = dict(defaults.items() + imagegrid_kw.items())

    grid = ImageGrid(fig, 111, **imagegrid_kw)

    # scale color bar data and model the same

    for i, (ax, z) in enumerate(zip(grid, ZS)):
        if i > 1:
            zz = z[np.isfinite(z)]
            vmin = -1. * np.max(np.abs(zz))
            vmax = np.max(np.abs(zz))
            # second row: make 0 on the color bar white
            if i == 3:
                if max_sig is not None:
                    vmin = -1 * max_sig
                    vmax = max_sig
            if i == 2:
                if max_diff is not None:
                    vmin = -1 * max_diff
                    vmax = max_diff
            colors = cm.RdBu
        else:
            # first row: make white 0, but will be on the left of color bar
            # scale color bar same for data and model.
            vmin = 0
            if max_counts is None:
                vmax = np.nanmax(ZS[0:2])
            else:
                vmax = max_counts
            if i == 0:
                colors = cm.Blues
            if i == 1:
                colors = cm.Reds
        im = ax.imshow(z, origin='upper', extent=extent, interpolation="nearest",
                       cmap=colors, vmin=vmin, vmax=vmax)
        ax.cax.colorbar(im)
        forceAspect(ax, aspect=1)

    # crop limits to possible data boundary
    ylim = (extent[2], extent[3])
    xlim = (extent[0], extent[1])
    for ax in grid:
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.xaxis.label.set_visible(True)
        ax.yaxis.label.set_visible(True)

    for ax, im_title in zip(grid, labels):
        t = add_inner_title(ax, im_title, loc=1)
        t.patch.set_alpha(0.5)

    return grid


def forceAspect(ax, aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1] - extent[0]) /
                      (extent[3] - extent[2])) / aspect)


def pgcmd(filename=None, cmd=None, labels=None, figname=None, out_dir=None,
          axis_labels='default', filter1=None, filter2=None, max_diff=None,
          max_counts=None, max_sig=None):
    '''
    produces the image that pgcmd.pro makes
    '''
    if cmd is None:
        cmd = read_match_cmd(filename)
        nmagbin = len(np.unique(cmd['mag']))
        ncolbin = len(np.unique(cmd['color']))
        data = cmd['Nobs'].reshape(nmagbin, ncolbin)
        model = cmd['Nsim'].reshape(nmagbin, ncolbin)
        diff = cmd['diff'].reshape(nmagbin, ncolbin)
        sig = cmd['sig'].reshape(nmagbin, ncolbin)
        hesses = [data, model, diff, sig]
        extent = [cmd['color'][0], cmd['color'][-1],
                  cmd['mag'][-1], cmd['mag'][0]]
    else:
        hesses = cmd.hesses
        extent = cmd.extent

    if axis_labels.lower() == 'default':
        if filter1 is None or filter2 is None:
            filter1 = ''
            filter2 = ''
        kwargs = {'xlabel': r'$%s-%s$' % (filter1, filter2),
                  'ylabel': r'$%s$' % filter2,
                  'max_diff': max_diff,
                  'max_sig': max_sig,
                  'max_counts': max_counts}

    if labels is None:
        grid = match_plot(hesses, extent, **kwargs)
    else:
        grid = match_plot(hesses, extent, labels=labels, **kwargs)

    [ax.set_xlabel('$%s-%s$' % (filter1, filter2), fontsize=20)
     for ax in grid.axes_row[1]]
    [ax.set_ylabel('$%s$' % filter2, fontsize=20) for ax in grid.axes_column[0]]
    grid.axes_all[0].xaxis.label.set_visible(True)

    if figname is not None:
        if out_dir is not None:
            figname = os.path.join(out_dir, os.path.split(figname)[1])
        plt.savefig(figname, dpi=300)
        plt.close()
        print('{} wrote {}'.format(pgcmd.__name__, figname))
    plt.close()
    return grid

def sfh_plot(MatchSFH):
    from matplotlib.ticker import NullFormatter
    fig, (ax1, ax2) = plt.subplots(nrows=2)
    MatchSFH.age_plot(ax=ax1)
    MatchSFH.age_plot(val='mh', convertz=False, ax=ax2)
    ax1.xaxis.set_major_formatter(NullFormatter())
    plt.subplots_adjust(hspace=0.1)
    figname = os.path.join(MatchSFH.base, MatchSFH.name + '.png')
    plt.savefig(figname)
    plt.close()
    print('{} wrote {}'.format(sfh_plot.__name__, figname))


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

def read_match_cmd(filename):
    '''
    reads MATCH .cmd file
    '''
    #mc = open(filename, 'r').readlines()
    # I don't know what the 7th column is, so I call it lixo.
    names = ['mag', 'color', 'Nobs', 'Nsim', 'diff', 'sig', 'lixo']
    cmd = np.genfromtxt(filename, skip_header=4, names=names, invalid_raise=False)
    return cmd


def call_pgcmd(filenames, filter1=None, filter2=None, labels=[]):

    if type(filenames) is not list:
        filenames = [filenames]

    for filename in filenames:
        mcmd = MatchCMD(filename)
        figname = filename + '.png'
        if filter1 is None:
            try:
                target, [filter1, filter2] = parse_pipeline(filename)
                labels[0] = '${}$'.format(target)
                if not 'W' in filter1:
                    filter1 = 'V'
                    filter2 = 'I'
            except:
                pass

        pgcmd(cmd=mcmd, filter1=filter1, filter2=filter2, labels=labels,
              figname=figname)
