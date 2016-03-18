"""Stats and visualization of calcsfh -ssp runs"""
import argparse
import itertools
import os
import sys

import matplotlib.pylab as plt
import numpy as np

from .config import EXT, match_base
from .fileio import read_ssp_output
from .utils import strip_header

__all__ = ['SSP']

try:
    plt.style.use('presentation')
except:
    pass

def add_filename_info_to_file(fname, ofile=None, ext='.dat'):
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
    if ofile is None:
        ofile = fname.replace(ext, '.f{}'.format(ext.replace('.','')))
    with open(fname, 'r') as inp:
        lines = inp.readlines()

    header = lines[:10]
    #h = 'Av IMF dmod lage logZ fit sfr sfrperr sfrmerr '
    h = 'Av IMF dmod lage logZ fit sfr '
    footer = lines[-1].strip()

    old_data, av, dmod, fit = read_ssp_output(fname)
    npts = len(old_data)
    new_stuff = filename_data(fname)
    names = new_stuff.keys()
    h += ' '.join(names)
    header.append(h)
    new_data = np.array([np.repeat(v, npts) for v in new_stuff.values()])

    # these three lines could be their own function in fileio
    import numpy.lib.recfunctions as nlr
    data = nlr.append_fields(np.asarray(old_data), names, new_data).data
    np.savetxt(ofile, data, fmt='%g', header=''.join(header),
               footer=footer.strip())
    return ofile, data

# this function may need to go in fileio
def filename_data(fname, ext='.dat', skip=2, delimiter='_', exclude='imf'):
    """
    return a dictionary of key and values from a filename.
    E.g, ssp_imf4.85_bf0.3_dav0.0.fdat
    returns bf: 0.3, dav: 0.0
    NB: imf is excluded because it's already included in the file.

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

    exclude : str
        do not include this key/value in the file (default: 'imf')

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
        if kv[0].lower() == exclude.lower():
            continue
        try:
            d[kv[0]] = float(neg + '.'.join(kv[1:]))
        except ValueError:
            #print e
            print(sys.exc_info()[0])
            pass
    return d

def sspcombine(fname, dry_run=True, outfile=None):
    sspname = strip_header(fname)
    if outfile is None:
        outfile = '> {}.stats'.format(fname)
    else:
        outfile = '>> {}'.format(outfile)
    cmd = '{} {} {}'.format(os.path.join(match_base, 'bin/sspcombine'), sspname, outfile)
    if not dry_run:
        print('excecuting: {}'.format(cmd))
        os.system(cmd)
    return cmd


class SSP(object):
    def __init__(self, filenames=None, base=os.getcwd(), data=None):
        """
        filenames are the calcsfh -ssp terminal or console output.
        They do not need to be stripped of their header or footer or
        be concatenated as is typical in MATCH useage.
        See .fileio.read_ssp_output.
        """
        self.data = []
        self.name = []
        self.base = []
        if data is None:
            [self.loaddata(f) for f in filenames]
        else:
            [self.loaddata(filenames[i], data=data[i]) for i in range(len(filenames))]

        self.all_data = np.concatenate(self.data)
        self.ibest = np.argmin(self.all_data['fit'])
        self.bestfit = np.array([self.all_data[np.argmin(self.all_data['fit'])]],
                                dtype=self.all_data.dtype)
        # via Dan Weisz:
        self.absprob = np.exp(0.5 * (self.all_data['fit'].min() \
                                            - self.all_data['fit']))

    def bestline(self):
        line = ''
        for k, v in dict(zip(self.all_data.dtype.names, *self.bestfit)).items():
            line += '{}: {} '.format(k, v)
        return line


    def boundries(self, stitched_col=None):
        """
        Print the min, max and best fit values
        Parameters
        stitched_col : str
            split the output for unique values of stitched_col (e.g., 'ov')
        Returns
        msg : str
            multiline message with format value, min, best, max.
        """
        if not hasattr(self, 'bestfit'):
            self.load_bestfits(stitched_col=stitched_col)
        bestfits = [self.bestfits]

        if stitched_col is not None:
            bestfits = self.bestfits

        msg = '# value min best max\n'
        fmt = '{} {} {} {}\n'

        for bestfit in bestfits:
            # not necessarily self.ibest because of stitched_col
            bestofbest = np.array(bestfit[np.argmin(bestfit['fit'])],
                                  dtype=self.all_data.dtype)
            if stitched_col is not None:
                msg += '# {}\n'.format(np.unique(bestfit[stitched_col]))

            for val in self.all_data.dtype.names:
                msg += fmt.format(val, np.min(bestfit[val]), bestofbest[val], np.max(bestfit[val]))
        return msg


    def bestfile(self, ind=None):
        """
        Not sure why I did this much work...
        if fmt == 'default':
            fmt = 'ssp_imf{IMF:.2f}_bf{bf:.1f}_dav{dav:.1f}_ov{ov:.2f}.fdat'
        if line == 'best':
            line = self.bestfit
        elif type(line) == int:
            line = self.all_data[line]
        else:
            assert len(line) == len(self.all_data.dtype.names), \
                'Either line=best, some int, or an actual line of self.data or self.all_data'
        return fmt.format(**dict(zip(self.all_data.dtype.names, *line)))
        """
        ind = ind or self.ibest
        return os.path.join(np.array(self.base)[ind], np.array(self.name)[ind])

    def load_bestfits(self, stitched_col=None):
        """
        cull the best fit in each ssp
        Parameters
        stitched_col : str
            an attribute like OV that was not varied during a MATCH run
            will separate bestfits on unique values of stitched_col

        Returns
        self.bestfits : array
            array of the min fit from each self.data. If stitch_col is not None,
            array is split for each unique value of stitched_col
        """
        self.bestfits = np.array([self.data[i][np.argmin(self.data[i]['fit'])]
                        for i in range(len(self.data))], dtype=self.data[0].dtype)

        if stitched_col is not None:
            usc, inds = np.unique(self.bestfit, return_inverse=True)
            self.bestfits = np.array([self.bestfits[inds==i]
                                      for i in np.unique(inds)])

        return

    def loaddata(self, filename, data=None):
        """Load the data and other book keeping information"""
        if data is None:
            data = read_ssp_output(filename)[0]

        base, name = os.path.split(filename)

        self.data.append(data)

        self.name.append(name)
        self.base.append(base)

    def _getmarginals(self):
        """get the values to marginalize over that exist in the data"""
        marg = np.array(['Av', 'IMF', 'dmod', 'lage', 'logZ', 'dav', 'ov',
                         'bf'])
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
        probs = np.array([np.sum(self.absprob[dbins==i])
                          for i in np.unique(dbins)])
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
            inds = list(set(np.nonzero(dbins == i)[0]) &
                        set(np.nonzero(dbins2 == j)[0]))
            probs.append(np.sum(self.absprob[inds]))
            vs.append(vals[i])
            vs2.append(vals2[j])
        probs = np.array(probs) / np.sum(probs)
        return vs, vs2, probs

    def pdf_plot(self, attr, ax=None, sub=''):
        """Plot prob vs marginalized attr"""
        if len(sub) > 0:
            sub = '_' + sub

        if ax is None:
            fig, ax = plt.subplots()

        vals, probs = self.marginalize(attr)
        best = np.argmax(probs)

        c = ax.plot(vals, np.array(probs), linestyle='steps-mid')

        ax.set_xlabel(key2label(attr))
        ax.set_ylabel('Probability')
        # probably need to make this a conditional if ax is passed...
        plt.savefig('{}_{}{}{}'.format(self.name[0].split('_')[0], attr, sub, EXT))
        plt.close()
        return ax


    def pdf_plot2d(self, attr, attr2, ax=None, sub=''):
        """Plot prob vs marginalized attr and attr2"""
        if len(sub) > 0:
            sub = '_' + sub

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
        # probably need to make this a conditional if ax is passed...
        plt.savefig('{}_{}_{}{}{}'.format(self.name[0].split('_')[0],
                                          attr, attr2, sub, EXT))
        plt.close()
        return ax

    def pdf_plots2d(self, marginals='default', sub=''):
        """Call pdf_plot2d for a list of attr and attr2"""
        if marginals=='default':
            marg = self._getmarginals()
        else:
            marg = marginals

        for i, j in itertools.product(marg, marg):
            # e.g., skip Av vs Av and Av vs IMF if already plotted IMF vs Av
            if i >= j:
                continue
            self.pdf_plot2d(i, j, sub=sub)
        return

    def pdf_plots(self, marginals='default', sub=''):
        """Call pdf_plot for a list of attr"""
        if marginals=='default':
            marg = self._getmarginals()
        else:
            marg = marginals

        [self.pdf_plot(i, sub=sub) for i in marg]
        return


def key2label(string):
    """latex labels for different strings"""
    def_fmt = r'${}$'
    convert = {'Av': r'$A_V$',
               'dmod': r'$\mu$',
               'lage': r'$\log\ \rm{Age\ (yr)}$',
               'logZ': r'$\log\ \rm{Z}$',
               'fit': r'$\rm{Fit\ Parameter}$',
               'ov': r'$\Lambda_c$',
               'chi2': r'$\chi^2$',
               'bf': r'$\rm{Binary\ Fraction}$',
               'dav': r'$dA_V$'}
    if not string in convert.keys():
        return def_fmt.format(string)
    return convert[string]


def main(argv):
    """
    Main function for ssp.py plot or reformat ssp output.

    e.g., Reformat and then plot a OV=0.30 run:
    python -m dweisz.match.scripts.ssp -fot --sub=ov0.30 *ov0.30.dat
    """
    parser = argparse.ArgumentParser(description="Plot or reformat MATCH.calcsfh -ssp output")

    parser.add_argument('-f', '--format', action='store_true',
                        help='reformat the files so data in the filename are columns in the file')

    parser.add_argument('-o', '--oned', action='store_true',
                        help='make val vs prob plots')

    parser.add_argument('-t', '--twod', action='store_true',
                        help='make val vs val vs prob plots')

    parser.add_argument('-s', '--sub', type=str, default='',
                        help='add substring to figure names')

    parser.add_argument('-l', '--list', action='store_true',
                        help='fnames is a file with a list of files to read')

    parser.add_argument('-c', '--sspcombine', action='store_true',
                        help='run sspcombine on the file(s) (and exit)')

    parser.add_argument('fnames', nargs='*', type=str,
                        help='ssp output(s) or formated output(s)')

    args = parser.parse_args(argv)

    if args.list:
        args.fnames = map(str.strip, open(args.fnames[0], 'r').readlines())

    if args.format:
        fnames, data = zip(*[add_filename_info_to_file(fname, ext='.scrn') for fname in args.fnames])
        if args.oned or args.twod:
            data = list(data)
            fnames = list(fnames)
            ssp = SSP(filenames=fnames, data=data)
    elif args.sspcombine:
        [sspcombine(f, dry_run=False) for f in args.fnames]
        sys.exit(0)
    else:
        import pdb; pdb.set_trace()
        ssp = SSP(filenames=args.fnames)
        print(ssp.bestline())

    if args.oned:
        ssp.pdf_plots(sub=args.sub)

    if args.twod:
        ssp.pdf_plots2d(sub=args.sub)


if __name__ == "__main__":
    main(sys.argv[1:])
