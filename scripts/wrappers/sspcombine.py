import os
from ..config import match_base
from ..utils import strip_header

def sspcombine(fname, dry_run=True, outfile=None):
    """
    call bin/sspcombine
    fname : string
        input filename
    outfile : None default: [fname].stats
        output filename
    dry_run : bool
        do not actually run sspcombine
    return string command to run sspcombine
    """
    sspname = strip_header(fname)
    if outfile is None:
        outfile = '> {}.stats'.format(fname)
    else:
        outfile = '>> {}'.format(outfile)
    cmd = '{} {} {}'.format(os.path.join(match_base, 'bin/sspcombine'),
                            sspname, outfile)
    if not dry_run:
        print('excecuting: {}'.format(cmd))
        os.system(cmd)
    return cmd
