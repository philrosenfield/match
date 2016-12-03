import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from .cmd import CMD
from .config import EXT
from .fileio import filename_data, get_files
from .graphics.graphics import zeroed_cmap, square_aspect
from .graphics.graphics import add_inner_title
from .ssp import SSP

plt.style.use('presentation')
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.5)


def comp_cmd(cmd0, cmd, label=None):
    """Make a hess diagram of the difference between two cmd.models"""
    extent = cmd0.extent
    aspect = abs((extent[1] - extent[0]) / (extent[3] - extent[2]))

    hess = cmd0.model - cmd.model
    colors = zeroed_cmap(hess, transparent=True)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    kw = {'extent': extent, 'origin': 'upper', 'interpolation': 'nearest',
          'aspect': aspect}
    if cmd0.data is not None:
        ax.imshow(cmd0.data, **kw)
    im = ax.imshow(hess, cmap=colors, alpha=0.5, **kw)
    if label is not None:
        add_inner_title(ax, label, loc=4)
    cb = plt.colorbar(im)
    cmd0.set_axis_labels(ax=ax)
    return fig, ax


def diff_(cmd0, cmd1, ssp=None):
    """create a dictionary of filename values that are different between two files"""
    d0 = filename_data(cmd0.name, exclude='')
    d1 = filename_data(cmd1.name, exclude='')
    if ssp is not None:
        d0.update(ssp.data.iloc[np.argmin(np.abs(ssp.data.fit - cmd0.fit))].to_dict())
        d1.update(ssp.data.iloc[np.argmin(np.abs(ssp.data.fit - cmd1.fit))].to_dict())
    a = {o : (d1[o], d0[o]) for o in set(d0.keys()).intersection(d1.keys())
         if d0[o] != d1[o]}
    return a


def main(argv=None):
    """main function for compare_cmds"""
    args = parse_args(argv)
    cmd1 = CMD(args.cmd1)
    cmd2 = CMD(args.cmd2)
    fig, ax = comp_cmd(cmd1, cmd2)
    fig.savefig(args.figname)
    print(diff_(cmd1, cmd2))
    return

def parse_args(argv=None):
    """parse_args for main"""
    parser = argparse.ArgumentParser(description="plot cmd file")

    parser.add_argument('-o', '--figname', type=str, default='comp_cmd{}'.format(EXT),
                        help='name of figure')

    parser.add_argument('cmd1', type=str, help='cmd file')

    parser.add_argument('cmd2', type=str, help='cmd file')

    return parser.parse_args(argv)


if __name__ == "__main__":
    sys.exit(main())
