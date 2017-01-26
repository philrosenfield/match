"""Plot match's mod* files
usage:
cd matchX.X/Model/data
e.g.,
python -m match.makemod.plot_mods -p 'mod*'
"""
from __future__ import print_function
import glob
import os
import sys
import argparse
import numpy as np
import matplotlib.pylab as plt

from ..scripts.config import EXT


def plot_mods(sub=None, pref='mod1_*', overwrite=False, lims=True):
    """
    make a plot of Mbol vs Log Te of tracks to go to MATCH or the
    MATCH tracks themselves

    color bar axis is either Nstars (if match tracks) or logAge (if unprocessed
    tracks)

    Parameters
    ----------
    sub : string
        the subdirectory to operate in [optional]

    pref : string
        the prefix search string of the track names [mod1*]
        if plotting unprocessed tracks pref should end with .dat.

    overwrite : bool [False]
        overwrite existing plots

    NB: Unprocessed match track plotting is not tested!
    TODO: either:
    hard code color bar limits or
    set all axes limits as options or
    use axes limits from makemod.cpp
    """
    here = os.getcwd()
    if sub is not None:
        assert os.path.isdir(sub), 'sub directory not found'
        os.chdir(sub)

    # i = LogT, j = Mbol k = Nstars or logAge
    # mod format:Mbol Log_Te Nstars ...
    i = 1
    j = 0
    k = 2
    zstr = 'Nstars'
    # unprocessed format logAge Mass logTe Mbol ...
    if pref.endswith('.dat'):
        i = 2
        j = 3
        k = 1
        zstr = 'Age'

    modfiles = glob.glob(pref)
    if len(modfiles) == 0:
        print('{} not found'.format(pref))
    for fname in modfiles:
        figname = '{}{}'.format(fname, EXT)
        if fname.endswith('.png'):
            continue
        if os.path.isfile(figname) and not overwrite:
            print('found {} and not overwriting'.format(figname))
            continue
        data = np.loadtxt(fname).T
        _, ax = plt.subplots()
        try:
            scr = ax.scatter(data[i], data[j], c=np.log10(data[k]),
                             cmap=plt.cm.Blues, edgecolor='none')
            cbar = plt.colorbar(scr)
            cbar.set_label('log {}'.format(zstr))
        except:
            ax.plot(data[i], data[j], color='k')
        ax.set_xlim(ax.get_xlim()[::-1])
        ax.set_ylim(ax.get_ylim()[::-1])
        if lims:
            ax.set_ylim(13, -14.0)
            ax.set_xlim(5.5, 3.0)

        ax.set_title(fname.replace('_', r'\_'))
        ax.set_xlabel('Log T')
        ax.set_ylabel('Mbol')
        plt.savefig(figname)
        print('wrote {}'.format(figname))
        plt.close()
    os.chdir(here)


def main(argv):
    """Main caller for plot_mods."""
    parser = argparse.ArgumentParser(description="Plot mod* files in data/")

    parser.add_argument('-s', '--sub', type=str,
                        help='subdirectory name')

    parser.add_argument('-r', '--recursive', action='store_true',
                        help='run on all subdirectories')

    parser.add_argument('-f', '--overwrite', action='store_true',
                        help='overwrite plots')

    parser.add_argument('-l', '--lims', action='store_false',
                    help='do not set plot limits automatically')

    parser.add_argument('-p', '--pref', type=str, default='mod1*',
                        help='search string for mod files')

    args = parser.parse_args(argv)

    if not os.getcwd().endswith('data'):
        print('warning, this should run in matchX.X/Model/data')

    subs = [args.sub]
    if args.recursive:
        subs = [s for s in os.listdir('.') if os.path.isdir(s)]

    for sub in subs:
        plot_mods(sub=sub, pref=args.pref, overwrite=args.overwrite,
                  lims=args.lims)

if __name__ == '__main__':
    main(sys.argv[1:])
