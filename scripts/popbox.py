"""Everything Population Box"""
import matplotlib.pyplot as plt
import numpy as np
import os

from ..graphics.graphics import add_inner_title

class PopBox(object):
    def __init__(self, filename=None):
        self.base, self.name = os.path.split(filename)
        if filename is not None:
            self.load_popbox(filename)

    def load_popbox(self, filename):
        """
        Parse MATCH's popbox file

        Parameters
        ----------
        filename : str
            population box file

        Returns
        -------
        sets as attributes:
            mh : array
                metalicity (line 2 of popbox file)
            lagei : array
                log age bin 0 (first column after line 2)
            lagef : array
                log age bin 1 (second column after line 2)
            sfr : array
                star formation rate grid shape Nages x Nmhs
        """
        with open(filename) as inp:
            lines = [np.float_(next(inp).strip().split()) for l in range(2)]

        self.mh = lines[1]
        data_ = np.genfromtxt(filename, skip_header=2)
        self.lagei = data_.T[0]
        self.lagef = data_.T[1]
        self.sfr =  data_.T[2:].T

    def plot_popbox(self, z=None, ax=None, cmap=plt.cm.Reds, figname=None,
                    clab=None, title=None, cbar=None, clim=None, vmin=None,
                    vmax=None):
        """
        Plot population box.
        Paramters
        ---------
        z : ndarray shape self.lagei, self.mh
            array to plot if not self.sfr.T

        ax : Axes instance
            axes to plot to (will create if not passed)

        cmap : mpl colormap
            color map for z

        figname : string
            filename to save figure to

        clab :
            color bar label (seems to be broken with mpl2.0)

        title : string
            add a title to the plot

        cbar : axes instance
            axes to plot colorbar (ax by default)
            useful when using AxisGrid

        clim : tuple or list
            color bar limits

        vmin, vmax : float, float
            passed to plt.pcolor

        Returns
        -------
        ax : axes instance
        (if figname was passed, wrote a figure.)
        """
        cbar = cbar or ax
        if z is None:
            z = self.sfr.T / 1e-3
            clab = r'$SFR (M_\odot / 10^{-3} {\rm yr})$'

        if ax is None:
            fig, ax = plt.subplots()

        c = ax.pcolor(self.lagei, self.mh, z, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xlabel(r'$\log\ {\rm Age}$')
        ax.set_ylabel(r'${\rm[M/H]}$')

        if clab is not None:
            c.set_label(clab)

        if cbar is not None:
            cbar.colorbar(c)

        if clim is not None:
            cbar.set_clim(clim)

        if title is not None:
            add_inner_title(ax, title, 4)

        if figname is not None:
            plt.savefig(figname)
        return ax


def compare_popboxes(pb1, pb2, titles=None, outfig=None):
    gridkw = {'nrows_ncols': (1, 3),
              'axes_pad': 0.7,
              'label_mode': "all",
              'share_all': True,
              'cbar_location': "top",
              'cbar_mode': "each",
              'cbar_size': "7%",
              'cbar_pad': "2%"}
    if titles is None:
        titles = [pb1.name, pb2.name, 'difference']

    difgrid = (pb1.sfr - pb2.sfr) / (pb1.sfr + pb2.sfr)
    fig = plt.figure(figsize=(12, 3))
    grid = AxesGrid(fig, 111, **gridkw)
    pb1.plot_popbox(ax=grid[0], title=titles[0], cbar=grid.cbar_axes[0])
    pb2.plot_popbox(ax=grid[1], title=titles[1], cbar=grid.cbar_axes[1])
    pb1.plot_popbox(ax=grid[2], title=titles[2], cbar=grid.cbar_axes[2],
                    z=difgrid.T)
    if outfig is not None:
        plt.savefig(outfig)
    return grid

def main(argv=None):
    """Main function for popbox.py for population boxes from MATCH"""
    parser = argparse.ArgumentParser(description="Plot population boxes")

    parser.add_argument('popboxes', nargs='*', type=str,
                        help='population box files')

    args = parser.parse_args(argv)
    for popbox in args.popboxes:
        pb = PopBox(popbox)
        pb.plot_popbox()


if __name__ == "__main__":
    sys.exit(main())
