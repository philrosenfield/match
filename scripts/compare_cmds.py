import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

from .cmd import CMD
from .config import EXT
from .fileio import filename_data, get_files
from .graphics.graphics import zeroed_cmap, square_aspect
from .graphics.graphics import add_inner_title
from .ssp import SSP

import matplotlib.colors as colors
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)


def comp_cmd(cmd0, cmd, label=None, figname=None, plot_data=True):
    """Make a hess diagram of the difference between two cmd.models"""
    extent = cmd0.extent
    aspect = abs((extent[1] - extent[0]) / (extent[3] - extent[2]))
    cmd0.model[cmd0.model < 1e-6] = 0.0
    cmd.model[cmd.model < 1e-6] = 0.0
    hess = cmd0.model - cmd.model
    extrema = np.max(abs(hess))

    # log spaced cmap norm, if used with diverging cmap, white will be at 0
    norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                             vmin=-1.0, vmax=1.0)

    # log spaced colorbar but with linear numbers.
    sm = plt.cm.ScalarMappable(cmap='RdBu_r',
                               norm=plt.Normalize(vmin=-extrema, vmax=extrema))
    sm._A = []

    fig, ax = plt.subplots(figsize=(8, 8))

    kw = {'extent': extent, 'origin': 'upper', 'interpolation': 'nearest',
          'aspect': aspect}

    if cmd0.data is not None and plot_data:
        ax.imshow(cmd0.data, cmap='Greys', **kw)

    im = ax.imshow(hess, cmap='RdBu_r', norm=norm, vmin=-extrema,
                   vmax=extrema, alpha=0.5, **kw)

    if label is not None:
        add_inner_title(ax, label, loc=4)
    plt.colorbar(sm, extend='both')
    # cb = plt.colorbar(im, extend='both')
    cmd0.set_axis_labels(ax=ax)
    plt.tight_layout()
    if figname is not None:
        plt.savefig(figname, dpi=300)

    return fig, ax


def diff_(cmd0, cmd1, ssp=None):
    """
    create a dictionary of filename values that are different between two files
    """
    d0 = filename_data(cmd0.name, exclude='')
    d1 = filename_data(cmd1.name, exclude='')
    if ssp is not None:
        d0.update(ssp.data.iloc[np.argmin(np.abs(ssp.data.fit
                                                 - cmd0.fit))].to_dict())
        d1.update(ssp.data.iloc[np.argmin(np.abs(ssp.data.fit
                                                 - cmd1.fit))].to_dict())
    a = {o: '{0:g}, {1:g}'.format(d1[o], d0[o])
         for o in set(d0.keys()).intersection(d1.keys())
         if d0[o] != d1[o]}

    return a


def main(argv=None):
    """main function for compare_cmds"""
    # shouldn't this be a 3 panel plot?
    args = parse_args(argv)
    cmd1 = CMD(args.cmd1)
    cmd2 = CMD(args.cmd2)
    fig, ax = comp_cmd(cmd1, cmd2, label=args.label)
    fig.savefig(args.figname, dpi=300)
    ssp = args.ssp
    if args.ssp is not None:
        ssp = SSP(args.ssp)
    print(diff_(cmd1, cmd2, ssp=ssp))
    return


def parse_args(argv=None):
    """parse_args for main"""
    parser = argparse.ArgumentParser(description="plot cmd file")

    parser.add_argument('-o', '--figname', type=str,
                        default='comp_cmd{0:s}'.format(EXT),
                        help='name of figure')

    parser.add_argument('cmd1', type=str, help='cmd file')

    parser.add_argument('cmd2', type=str, help='cmd file')

    parser.add_argument('--ssp', type=str, help='combined ssp_file file')

    parser.add_argument('--label', type=str, nargs=2, help='label')

    return parser.parse_args(argv)


if __name__ == "__main__":
    sys.exit(main())
