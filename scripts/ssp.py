"""Stats and visualization of calcsfh -ssp runs"""
import argparse
import itertools
import os
import sys

import matplotlib.pylab as plt
import numpy as np

from .fileio import read_ssp_output

__all__ = ['SSP']
try:
    plt.style.use('presentation')
except:
    pass

def add_filename_info_to_file(fname, ofile=None):
    """
    add filename info to the data.
    E.g, ssp_imf4.85_bf0.3_dav0.0.dat
    will add two columns, bf, and dav. See filename_data.
    Parameters
    ----------
    fname : str
        name of the file

    ofile : str
        output file or  will write to file substiting fname's .dat with .fdat

    Returns
    -------
    data : np.array
        data with new columns attached

    """
    import numpy.lib.recfunctions as nlr
    if ofile is None:
        ofile = fname.replace('.dat', '.fdat')
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
    np.savetxt(ofile, data, fmt='%g', header=''.join(header), footer=footer.strip())
    return data

def filename_data(fname, ext='.dat', skip=1, delimiter='_'):
    """
    return a dictionary of key and values from a filename.
    E.g, ssp_imf4.85_bf0.3_dav0.0.fdat
    returns bf: 0.3, dav: 0.0
    imf is excluded because it's already included in the file.

    Parameters
    ----------
    fname : str
        filename

    ext : str
        extension (sub string to remove from the tail)

    delimiter : str
        how the keyvals are separated '_' in example above

    skip : int
        skip n items (skip=1 skips ssp in the above example)

    Returns
    -------
    dict of key and values from filename
    """
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
        self.name = []
        self.base = []

        [self.loaddata(f) for f in filenames]
        self.all_data = np.concatenate(self.data)
        # via Dan Weisz:
        self.absprob = np.exp(0.5 * (self.all_data['fit'].min() \
                                            - self.all_data['fit']))

    def loaddata(self, filename):
        """Load the data and other book keeping information"""
        data, _, _, fit = read_ssp_output(filename)
        base, name = os.path.split(filename)

        # so far marginalize needs this concatenated. Not sure if
        # a list of arrays is useful.
        self.data.append(data)

        # book keeping
        self.name.append(name)
        self.base.append(base)

    def _getmarginals(self):
        """get the values to marginalize over that exist in the data"""
        marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav', 'cov', 'bf'])
        inds = [i for i, m in enumerate(marg) if self._haskey(m)]
        return marg[inds]


    def _haskey(self, key):
        """test if the key requested is available"""
        possible = self.data[0].dtype.names
        ecode = True
        if not key in possible:
            ecode = False
        return ecode

    def marginalize(self, attr):
        """Find the best fit for each unique value of attr"""
        assert self._haskey(attr), '{} not found'.format(attr)
        vals = np.unique(self.all_data[attr])
        dbins = np.digitize(self.all_data[attr], vals, right=True)
        probs = np.array([np.sum(self.absprob[dbins==i]) for i in np.unique(dbins)])
        probs /= probs.sum()
        return vals, probs

    def marginalize_2d(self, attr, attr2):
        """Find the best fit for each set of unique values of attr and attr2"""
        assert self._haskey(attr), '{} not found'.format(attr)
        assert self._haskey(attr2), '{} not found'.format(attr)
        vals = np.unique(self.all_data[attr])
        dbins = np.digitize(self.all_data[attr], vals, right=True)

        vals2 = np.unique(self.all_data[attr2])
        dbins2 = np.digitize(self.all_data[attr2], vals2, right=True)
        probs = []
        vs = []
        vs2 = []
        for i, j in itertools.product(np.unique(dbins), np.unique(dbins2)):
            inds = list(set(np.nonzero(dbins == i)[0]) & set(np.nonzero(dbins2 == j)[0]))
            probs.append(np.sum(self.absprob[inds]))
            vs.append(vals[i])
            vs2.append(vals2[j])
        probs = np.array(probs) / np.sum(probs)
        return vs, vs2, probs

    def pdf_plot(self, attr, ax=None):
        """Plot prob vs marginalized attr"""
        if ax is None:
            fig, ax = plt.subplots()

        vals, probs = self.marginalize(attr)
        best = np.argmax(probs)

        c = ax.plot(vals, np.array(probs), linestyle='steps-mid')

        ax.set_xlabel(key2label(attr))
        ax.set_ylabel('Probability')
        plt.savefig('{}_{}.png'.format(self.name[0].split('_')[0], attr))
        return ax


    def pdf_plot2d(self, attr, attr2, ax=None):
        """Plot prob vs marginalized attr and attr2"""
        if ax is None:
            fig, ax = plt.subplots()

        vs, vs2, probs = self.marginalize_2d(attr, attr2)
        best = np.argmax(probs)
        vals = np.unique(vs)
        vals2 = np.unique(vs2)

        #P = np.array(probs).reshape(len(vals), len(vals2))
        # I don't know why this fails sometimes e.g, IMF vs Av or logZ
        #ax.contour(vals, vals2, P, 20, cmap=plt.cm.Blues_r)
        # really ugly, slow, and silly way to do this:
        c = ax.scatter(vs, vs2, c=probs, cmap=plt.cm.Blues, s=1000, marker='s')
        ax.axhline(vs2[best], color='w')
        ax.axvline(vs[best], color='w')
        ax.plot(vs[best], vs2[best], 'o', color='w')

        l = plt.colorbar(c)
        ax.set_xlabel(key2label(attr))
        ax.set_ylabel(key2label(attr2))
        l.set_label('Probability')
        plt.savefig('{}_{}_{}.png'.format(self.name[0].split('_')[0], attr, attr2))
        return ax

    def pdf_plots2d(self, marginals='default'):
        """Call pdf_plot2d for a list of attr and attr2"""
        if marginals=='default':
            marg = self._getmarginals()
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
            marg = self._getmarginals()
        else:
            marg = marginals

        [self.pdf_plot(i) for i in marg]
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


def main(argv):
    parser = argparse.ArgumentParser(description="Plot match ssp output")

    parser.add_argument('-r', '--format', action='store_true',
                        help='reformat the files to have the data in the filename part of the file')

    parser.add_argument('-o', '--oned', action='store_true',
                        help='make val vs prob plots')

    parser.add_argument('-t', '--twod', action='store_true',
                        help='make val vs val vs prob plots')

    parser.add_argument('fnames', nargs='*', type=str,
                        help='ssp output(s) or formated output(s)')

    args = parser.parse_args(argv)

    if args.format:
        [add_filename_info_to_file(fname) for fname in args.fnames]
    else:
        ssp = SSP(args.fnames)

        if args.oned:
            ssp.pdf_plots()

        if args.twod:
            ssp.pdf_plots2d()


if __name__ == "__main__":
    main(sys.argv[1:])
