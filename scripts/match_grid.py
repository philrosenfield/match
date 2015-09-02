import os

import matplotlib.pylab as plt
import numpy as np

from .utils import grab_val
from .fileio import savetxt

__all__ = ['MatchGrid']

class MatchGrid(object):
    def __init__(self, filenames, ssp=False, base=os.getcwd()):
        self.ssp = ssp
        if type(filenames) is str:
            if os.path.isfile(filenames) is True:
                self.data = fileio.readfile(filenames)
            elif '*' in filenames:
                filenames = fileio.get_files(base, filenames)

        if type(filenames) is list:
            self.fname_info(filenames)
            self.load_grid(filenames)

    def fname_info(self, filenames):
        self.base = os.path.split(filenames[0])[0]
        self.names = [os.path.split(f)[1] for f in filenames]
        self.zs = [grab_val(f, 'lz') - 4. for f in self.names]
        self.covs = [grab_val(f, 'ov') for f in self.names]

    def check_complete(self, liness):
        inds = []
        for i, lines in enumerate(liness):
            if len(lines) <= 1:
                print 'skipping %s' % self.names[i]
                del liness[i]
                inds.append(i)
        if len(inds) > 0:
            for i in inds:
                self.names.pop(i)
                self.zs.pop(i)
                self.covs.pop(i)

    def load_grid(self, filenames):
        """
        file names is a set of ssp screen outputs for various COV values.
        Will make an array with the data attaching a column of COV.
        """
        liness = []
        for filename in filenames:
            with open(filename, 'r') as infile:
                lines = infile.readlines()
                foot = -2
                if not self.ssp:
                    lines = lines[10:]
                    foot = -1
            liness.append(lines[:foot])

        self.check_complete(liness)

        if self.ssp is True:
            col_keys = ['Av', 'IMF', 'dmod', 'logAge', 'mh', 'fit', 'bg1', 'bg2',
                        'COV']
        else:
            col_keys = ['Av', 'IMF', 'dmod', 'fit', 'Z', 'COV']

        dtype = [(c, float) for c in col_keys]
        nrows = len(np.concatenate(liness))
        self.data = np.ndarray(shape=(nrows,), dtype=dtype)
        row = 0
        for i, lines in enumerate(liness):
            for j in range(len(lines)):
                if self.ssp is True:
                    line = lines[j]
                else:
                    line = re.sub('[Avimfdmodfit=,:]', '', lines[j])
                datum = np.array(line.strip().split(), dtype=float)
                if not self.ssp:
                    datum = np.append(datum, self.zs[i])
                datum = np.append(datum, self.covs[i])
                self.data[row] = datum
                row += 1

    def write_grid(self, outfile, overwrite=True):
        savetxt(outfile, self.data, fmt='%.3g',
                header='# %s\n' % ' '.join(self.data.dtype.names),
                overwrite=overwrite)

    def pdf_plots(self, xcol, ycol, zcol, stat='median', bins='uniq'):
        import matplotlib.gridspec as gridspec

        if self.ssp is True:
            met_col = 'mh'
            unzs = np.unique(self.data[met_col])
        else:
            met_col = 'Z'
            unzs = np.unique(self.zs)
        uncovs = np.unique(self.covs)
        gs = gridspec.GridSpec(len(unzs), len(uncovs))
        k = 0
        vmin = np.min(self.data[zcol])
        vmax = np.max(self.data[zcol])
        fig = plt.figure()
        for i in range(len(unzs)):
            izs, = np.nonzero(self.data[met_col] == unzs[i])
            for j in range(len(uncovs)):
                icovs, = np.nonzero(self.data['COV'] == uncovs[j])

                inds = list(set(izs) & set(icovs))
                ax = plt.subplot(gs[k])
                #ax.set_adjustable('box-forced')
                #ax.autoscale(False)
                ax.set_title('%g %g' % (uncovs[j], unzs[i]))
                ax = self.pdf_plot(xcol, ycol, zcol, stat=stat,
                                   bins=bins, ax=ax, cbar=False,
                                   inds=inds, vmin=vmin, vmax=vmax)
                k += 1
        return fig, gs

    def pdf_plot(self, xcol, ycol, zcol, stat='median', bins='uniq',
                 log=False, cbar=True, inds=None, ax=None, vmin=None,
                 vmax=None):
        try:
            from astroML.stats import binned_statistic_2d
        except ImportError:
            print('need astroML.stats.binned_statistic_2d for this function')
            return -1

        if inds is None:
            inds = np.arange(len(self.data[xcol]))

        if bins == 'uniq':
            bins=[np.unique(self.data[xcol][inds]),
                  np.unique(self.data[ycol][inds])]

        N, xe, ye = binned_statistic_2d(self.data[xcol][inds],
                                        self.data[ycol][inds],
                                        self.data[zcol][inds],
                                        stat, bins=bins)
        if log is True:
            n = np.log10(N.T)
        else:
            n = N.T

        aspect = (xe[-1] - xe[0]) / (ye[-1] - ye[0])

        if ax is None:
            fig, ax = plt.subplots()

        im = ax.imshow(n, extent=[xe[0], xe[-1], ye[0], ye[-1]], vmin=vmin,
                       vmax=vmax, cmap=plt.cm.Blues_r, aspect=aspect,
                       interpolation='nearest')

        ax.set_xlabel(key2label(xcol), fontsize=16)
        ax.set_ylabel(key2label(ycol), fontsize=16)

        if cbar:
            cb = plt.colorbar(im)
            if callable(stat):
                stat = stat.__name__
            cb.set_label(r'$\rm{%s}$\ %s' % (stat, key2label(zcol)))
            ax = (ax, cb)

        return ax

def key2label(string):
    def_fmt = r'$%s$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu$',
               'logAge': r'$\log\ \rm{Age\ (yr)}$',
               'mh': r'$\rm{[M/H]}$',
               'fit': r'$\rm{Fit\ Parameter}$',
               'COV': r'$\Lambda_c$',
               'chi2': r'$\chi^2$'}
    if not string in convert.keys():
        return def_fmt % string
    return convert[string]
