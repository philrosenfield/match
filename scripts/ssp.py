import os

import matplotlib.pylab as plt
import numpy as np
import itertools

from .fileio import read_ssp_output

__all__ = ['SSP']

def add_filename_info_to_file(fname):
    import numpy.lib.recfunctions as nlr
    oname = fname.replace('.dat', '.fdat')
    with open(fname, 'r') as inp:
        lines = inp.readlines()

    header = lines[:10]
    h = 'Av IMF dmod lage logZ fit sfr sfrperr sfrmerr '
    footer = lines[-1].strip()

    old_data, av, dmod, fit = read_ssp_output(fname)
    npts = len(old_data)
    new_stuff = filename_data(fname)
    names = new_stuff.keys()
    h += ' '.join(names)
    header.append(h)
    new_data = np.array([np.repeat(v, npts) for v in new_stuff.values()])

    data = nlr.append_fields(np.asarray(old_data), names, new_data).data
    np.savetxt(oname, data, fmt='%g', header=''.join(header), footer=footer.strip())
    print 'wrote', oname
    return data

def filename_data(fname, ext='.dat', skip=1, delimiter='_'):
    import re
    keyvals = fname.replace(ext, '').split(delimiter)[skip:]
    d = {}
    for keyval in keyvals:
        kv = re.findall(r'\d+|[a-z]+', keyval)
        neg = ''
        if '-' in keyval:
            neg = '-'
        if kv[0] == 'imf':
            continue
        d[kv[0]] = float(neg + '.'.join(kv[1:]))
    return d


class SSP(object):
    def __init__(self, filenames, base=os.getcwd()):
        """
        filenames are the calcsfh -ssp terminal or console output.
        They do not need to be stripped of their header or footer or
        be concatenated as is typical in MATCH useage.
        See .fileio.read_ssp_output.
        """
        self.data = []
        self.fit = []
        self.name = []
        self.base = []

        [self.loaddata(f) for f in filenames]

    def loaddata(self, filename, i=0):
        """Load the data and other book keeping information"""
        data, _, _, fit = read_ssp_output(filename)
        base, name = os.path.split(filename)

        # so far marginalize needs this concatenated. Not sure if
        # a list of arrays is useful.
        self.data.append(data)

        # best fit -- which is just the min fit of each data, not
        # really needed.
        self.fit.append(fit)

        # book keeping
        self.name.append(name)
        self.base.append(base)

    def _haskey(self, key):
        """make sure the key requested is available."""
        possible = self.data[0].dtype.names
        assert key in possible, 'choose from {}'.format(possible)

    def _hasalldata(self):
        if not hasattr(self, 'all_data'):
            self.all_data = np.concatenate(self.data)

    def marginalize(self, attr):
        """Find the best fit for each unique value of attr"""
        self._hasalldata()
        self._haskey(attr)
        vals = np.unique(self.all_data[attr])
        dbins = np.digitize(self.all_data[attr], vals, right=True)
        probs = [np.min(self.all_data['fit'][dbins==i]) for i in np.unique(dbins)]
        return vals, probs

    def marginalize_2d(self, attr, attr2):
        """Find the best fit for each set of unique values of attr and attr2"""
        self._hasalldata()
        self._haskey(attr)
        self._haskey(attr2)
        vals = np.unique(self.all_data[attr])
        dbins = np.digitize(self.all_data[attr], vals, right=True)

        vals2 = np.unique(self.all_data[attr2])
        dbins2 = np.digitize(self.all_data[attr2], vals2, right=True)
        probs = []
        vs = []
        vs2 = []
        for i, j in itertools.product(np.unique(dbins), np.unique(dbins2)):
            inds = list(set(np.nonzero(dbins == i)[0]) & set(np.nonzero(dbins2 == j)[0]))
            probs.append(np.min(self.all_data['fit'][inds]))
            vs.append(vals[i])
            vs2.append(vals2[j])

        return vs, vs2, probs

    def pdf_plot(self, attr, ax=None):
        """Plot prob vs marginalized attr"""
        if ax is None:
            fig, ax = plt.subplots()

        vals, probs = self.marginalize(attr)
        best = np.argmin(probs)

        c = ax.plot(vals, 1 - np.array(probs), linestyle='steps-mid')
        ax.axhline(vals[best])
        ax.axvline(1 - np.array(probs)[best])

        ax.set_xlabel(key2label(attr))
        ax.set_ylabel('Probability')

        return ax


    def pdf_plot2d(self, attr, attr2, ax=None):
        """Plot prob vs marginalized attr and attr2"""
        if ax is None:
            fig, ax = plt.subplots()

        vs, vs2, probs = self.marginalize_2d(attr, attr2)
        best = np.argmin(probs)
        vals = np.unique(vs)
        vals2 = np.unique(vs2)


        #P = np.array(1. - probs/max(probs)).reshape(len(vals), len(vals2))
        # I don't know why this fails sometimes e.g, IMF vs Av or logZ
        #ax.contour(vals, vals2, P, 20, cmap=plt.cm.Blues_r)

        c = ax.scatter(vs, vs2, c=1.-np.array(probs)/np.max(probs),
                       cmap=plt.cm.Blues, s=1000, marker='s')
        ax.axhline(vs2[best], color='w')
        ax.axvline(vs[best], color='w')
        ax.plot(vs[best], vs2[best], 'o', color='w')

        l = plt.colorbar(c)
        ax.set_xlabel(key2label(attr))
        ax.set_ylabel(key2label(attr2))
        l.set_label('Probability')
        return ax

    def pdf_plots2d(self, marginals='default'):
        """Call pdf_plot2d for a list of attr and attr2"""
        if marginals=='default':
            marg = ['Av', 'IMF', 'dmod', 'lage', 'logZ']
        else:
            marg = marginals

        for i, j in itertools.product(marg, marg):
            if i >= j:
                continue
            self.pdf_plot2d(i, j)
        return

    def pdf_plots(self, marginals='default'):
        """Call pdf_plot for a list of attr"""
        if marginals=='default':
            marg = ['Av', 'IMF', 'dmod', 'lage', 'logZ']
        else:
            marg = marginals

        for i, j in itertools.product(marg, marg):
            if i >= j:
                continue
            self.pdf_plot(i, j)
        return


def key2label(string):
    """latex labels for different strings"""
    def_fmt = r'${}$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu$',
               'lage': r'$\log\ \rm{Age\ (yr)}$',
               'logZ': r'$\log\ \rm{Z}$',
               'fit': r'$\rm{Fit\ Parameter}$',
               'COV': r'$\Lambda_c$',
               'chi2': r'$\chi^2$',
               'bf': r'$\rm{Binary\ Fraction}$',
               'dav': r'$dA_V$'}
    if not string in convert.keys():
        return def_fmt.format(string)
    return convert[string]
