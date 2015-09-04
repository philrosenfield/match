import numpy as np
import os
import matplotlib.pyplot as plt
import sys

try:
    plt.style.use('presentation')
except:
    pass

class SSP(object):
    def __init__(self, filename):
        self.base, self.name = os.path.split(filename)
        self.data = read_ssp_output(filename)
        self.load_best()

    def load_best(self):
        self.best = self.data[np.argmin(self.data['fit'])]

    def marginalize(self, marg_attr, value='best'):
        vals = ['Av', 'dmod', 'lage', 'logZ']
        assert marg_attr in vals, 'must choose between {}'.format(vals)

        if value == 'best':
            value = self.best[marg_attr]

        inds, = np.nonzero(self.data[marg_attr]==value)

        new_attr = '{}_inds'.format(marg_attr)
        self.__setattr__(new_attr, inds)
        return inds

    def age_plot(self):
        [self.marginalize(v) for v in ['Av', 'dmod', 'logZ']]
        inds = list(set(self.Av_inds) & set(self.dmod_inds) & set(self.logZ_inds))

        mfit = np.max(self.data['fit'])
        fit = 1 - self.data['fit'] / mfit
        fig, ax = plt.subplots()

        ax.set_yticks([])
        ax.plot(self.data['lage'], fit, 'o', alpha=0.3, color='gray', mec='none')
        bestage = self.data['lage'][inds]
        bestfit = fit[inds]
        isort = np.argsort(bestage)
        ax.plot(bestage[isort], bestfit[isort], 'o', ms=8, color='k')
        ax.set_title(translate('Av={}\ dmod={}\ logZ={}'.format(self.best['Av'],
                                                             self.best['dmod'],
                                                             self.best['logZ'])))
        ax.set_xlabel(translate('lage'))
        return fig, ax

    def ssp_plot0(self, marg_attr):
        # not sure this was worth doing...
        vals = ['Av', 'dmod', 'lage', 'logZ', 'sfr']

        try:
            vals.pop(vals.index(marg_attr))
        except ValueError:
            print('must choose between {}'.format(vals))
            return

        # marginalize
        inds = list(set(self.marginalize(marg_attr)) & set(self.marginalize('Av')))
        data = self.data[inds]

        fig, axs = plot_config()
        for i, j in itertools.product(range(4), range(4)):
            if j > i:
                continue
            y = vals[i]
            x = vals[j]
            ax = axs[i, j]
            ax.set_ylabel(translate(vals[i]))
            ax.set_xlabel(translate(vals[j]))
            if i == j:
                ax.plot(data[x], data['fit'])
            # plot x vs y
            ax.scatter(data[x], data[y], c=data['fit'], marker='s', s=40,
                       cmap=plt.cm.Spectral, alpha=.5, edgecolors='none')
            ax.axhline(self.best[y])
            ax.axvline(self.best[x])

        axs[0, -1].annotate(translate('{}={}'.format(marg_attr, self.best[marg_attr])),
                            (1, 1), ha='right', fontsize=24, color='k')

        return fig, axs

def translate(string):
    if 'log' in string:
        string = string.replace('log', '\log\ ')
    if 'lage' in string:
        string = string.replace('la', '\log\ A')
    if 'dmod' in string:
        string = string.replace('dmod', '\mu')
    if 'Av' in string:
        string = string.replace('Av', 'A_V')
    if 'sfr' in string:
        string = string.replace('sfr', 'SFR')
    return r'${}$'.format(string)

def read_ssp_output(filename):
    colnames = ['Av', 'IMF', 'dmod', 'lage', 'logZ', 'fit', 'sfr', 'sfrperr',
                'sfrmerr']
    return np.genfromtxt(filename, skip_header=10, skip_footer=1, names=colnames)


def plot_config():
    """
    Make an 4x4 figure with top corner frames off and ticks on the left
    and bottom edges
    """
    fig, axs = plt.subplots(4, 4, figsize=(9.7, 9.7))
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.96, top=0.96,
                        wspace=0.05, hspace=0.05)
    for i in range(4):
        for j in range(4):
            ax = axs[i, j]
            # only have y ticks on left axes (axs[:, 0])
            ax.set_yticks([])
            # turn off top right frames
            if j > i:
                ax.set_frame_on(False)
            # only have x ticks on bottom axes (axs[-1, 0])
            if i != 4:
                ax.set_xticks([])
    return fig, axs


if __name__ == '__main__':
    ssp = SSP(sys.argv[1])
    fig, ax = ssp.age_plot()
    plt.savefig(sys.argv[1] + '.png')
