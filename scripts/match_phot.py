"""Make 2-column photometry file for CALCSFH"""
from __future__ import print_function
import argparse
import sys
import numpy as np
from astropy.io import fits
from .utils import splitext
from .config import PHOTEXT

def make_phot(fitsfile, filterext='VEGA', nexts=2, precision='%.6f',
              filter2=None, dryrun=False):
    """Make mag1, mag2 files from fits file"""
    if nexts == 2:
        # e.g., .gst.fits is the extension keep .gst.
        pref0, _ = splitext(fitsfile)
        pref, ext = splitext(pref0)
        ext = '.{}'.format(ext)
    else:
        pref, _ = splitext(fitsfile)
        ext = ''

    fnames = []
    data = fits.getdata(fitsfile)
    # All filters, but no accidental narrow bands
    filters = [i for i in data.dtype.names if filterext.upper() in i.upper()
               and not 'N_' in i.upper()]
    if filter2 is None:
        filter2 = filters[-1]
        if not '814' in filter2:
            print("Warning: Assuming {} is mag2.".format(filter2))

    filters.pop(filters.index(filter2))

    for filter1 in filters:
        filts = '-'.join([f.replace('_{}'.format(filterext), '')
                          for f in [filter1, filter2]])
        fname = '{}_{}{}{}'.format(pref, filts, ext, PHOTEXT)
        if not dryrun:
            wrote = 'wrote'
            np.savetxt(fname, np.column_stack([data[filter1], data[filter2]]),
                       fmt=precision)
        else:
            wrote = 'would write'
        print('{} {}'.format(wrote, fname))
        fnames.append(fname)
    return fnames


def main(argv):
    """Main function for make_phot."""
    parser = argparse.ArgumentParser(description="Make match photometry formatted files from fits")

    parser.add_argument('fitsfiles', nargs='*', type=str,
                        help='fitsfiles')

    parser.add_argument('--dryrun', action='store_true',
                        help='only print filename that would be written')

    args = parser.parse_args(argv)

    _ = [make_phot(fitsfile, dryrun=args.dryrun) for fitsfile in args.fitsfiles]

if __name__ == "__main__":
    main(sys.argv[1:])
