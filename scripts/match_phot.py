"""Make 2-column photometry file for CALCSFH"""
from __future__ import print_function
import argparse
import sys
import numpy as np
from astropy.io import fits
from .utils import splitext
from .config import PHOTEXT


def match_fmt(data, filter1, filter2, crowd=None):
    """CALCSFH input format"""
    if crowd is not None:
        filterext = filter1.split('_')[1]
        filter1c = '{}CROWD'.format(filter1.replace(filterext, ''))
        filter2c = '{}CROWD'.format(filter2.replace(filterext, ''))
        inds, = np.nonzero((data[filter1c] < crowd) & (data[filter2c] < crowd))
        data = data[inds]

    return np.column_stack([data[filter1], data[filter2]])


def asteca2matchphot(filename):
    """Convert asteca output to match photometic input"""
    # Columns: ID x y mag e_mag col1 e_col1 memb_prob sel
    # ASSUMING: mag == mag1 since asteca_fmt and match assumes that.
    mag1, color = np.genfromtxt(filename, usecols=(3,5), unpack=True)
    outfile = filename.replace('.dat', PHOTEXT)
    np.savetxt(outfile, np.column_stack([mag1, mag1 - color]), fmt='%.3f')
    print('wrote {0:s}'.format(outfile))
    return


def asteca_fmt(data, filter1, filter2, filterext='VEGA', crowd=None):
    """ASteCA input format"""
    filter1e = '{}ERR'.format(filter1.replace(filterext, ''))
    filter2e = '{}ERR'.format(filter2.replace(filterext, ''))
    if crowd is not None:
        filter1c = '{}CROWD'.format(filter1.replace(filterext, ''))
        filter2c = '{}CROWD'.format(filter2.replace(filterext, ''))
        inds, = np.nonzero((data[filter1c] < crowd) & (data[filter2c] < crowd))
        data = data[inds]
    color = data[filter1] - data[filter2]
    colore = np.sqrt(data[filter1e]**2 + data[filter2e]**2)
    idx = np.arange(len(data)) + 1.
    return np.column_stack([idx, data['RA'], data['DEC'],
                            data[filter1], data[filter1e],
                            color, colore])

def outputfmt(data, filter1, filter2, asteca=False, filterext='VEGA',
              crowd=None):
    """Choose either ASteCA or MATCH format photometry"""
    if asteca:
        retdat = asteca_fmt(data, filter1, filter2, filterext=filterext,
                            crowd=None)
    else:
        retdat = match_fmt(data, filter1, filter2, crowd=crowd)
    return retdat


def make_phot(fitsfile, filterext='VEGA', nexts=2, precision='%.6f',
              filter2=None, dryrun=False, asteca=False, crowd=None):
    """Make mag1, mag2 files from fits file"""
    if nexts == 2:
        # e.g., .gst.fits is the extension keep .gst.
        pref0, ext0 = splitext(fitsfile)
        pref, ext = splitext(pref0)
        if ext == pref0:
            ext = ''
        else:
            ext = '.{}'.format(ext)
        if len(pref) == 0:
            pref = pref0
    else:
        pref, _ = splitext(fitsfile)
        ext = ''
    fts = [p for p in pref.split('_') if p.startswith('F') and p.endswith('W')]
    pref = pref.replace(''.join(fts), '').replace('__', '_')
    fnames = []
    #if asteca:
    #    data = load_asteca(fitsfile)
    #else:
    data = fits.getdata(fitsfile)
    # All filters, but no accidental narrow bands
    filters = [i for i in data.dtype.names
               if filterext.upper() in i.upper() and 'N_' not in i.upper()]
    if filter2 is None:
        filter2 = filters[-1]
        if '814' not in filter2:
            print("Warning: Assuming {} is mag2.".format(filter2))

    filters.pop(filters.index(filter2))

    for filter1 in filters:
        filts = '-'.join([f.replace('_{}'.format(filterext), '')
                          for f in [filter1, filter2]])
        final_ext = PHOTEXT
        if asteca:
            final_ext = '.asteca'
        fname = '{}_{}{}{}'.format(pref, filts, ext, final_ext)
        fname = fname.replace('__', '_')
        if not dryrun:
            wrote = 'wrote'
            np.savetxt(fname, outputfmt(data, filter1, filter2, asteca=asteca,
                                        filterext=filterext, crowd=crowd),
                       fmt=precision)
        else:
            wrote = 'would write'
        print('{} {}'.format(wrote, fname))
        fnames.append(fname)
    return fnames


def main(argv):
    """Main function for make_phot."""
    parser = argparse.ArgumentParser(
        description="Make match photometry formatted files from fits")

    parser.add_argument('fitsfiles', nargs='*', type=str,
                        help='fitsfiles')

    parser.add_argument('--asteca', action='store_true',
                        help='format photometry for ASteCA')

    parser.add_argument('--a2m', action='store_true',
                        help='Convert ASteCA to matchphot')

    parser.add_argument('--dryrun', action='store_true',
                        help='only print filename that would be written')

    parser.add_argument('--crowd', type=float,
                        help='cull fitsfile to this level of crowding')

    args = parser.parse_args(argv)

    if args.a2m:
        [asteca2matchphot(f) for f in args.fitsfiles]
    else:
        _ = [make_phot(fitsfile, dryrun=args.dryrun, asteca=args.asteca,
                       crowd=args.crowd)
             for fitsfile in args.fitsfiles]

if __name__ == "__main__":
    main(sys.argv[1:])
