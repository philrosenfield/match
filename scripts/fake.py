import argparse
import os
import sys
from .fileio import read_calcsfh_param, fake_param_fmt
from .sfh import SFH


def make_fakeparam(param, sfhfile, outfile=None, overwrite=False):
    """
    Convert a calcsfh solution and parameter file into a fake parameter file
    """
    if outfile is None:
        ext = os.path.splitext(sfhfile)[1]
        outfile = sfhfile.replace(ext, '.fparam')

    pdict = read_calcsfh_param(param)
    sfh = SFH(sfhfile)

    pdict['dmod'] = sfh.dmod
    pdict['av'] = sfh.Av
    pdict['Zspread'] = sfh.dlogZ
    pdict['dmag_min'] = -1.50
    ntbins = len(sfh.data.sfr)

    line = '\n'.join([f % pdict for f in fake_param_fmt()])
    line += '\n{}\n'.format(ntbins)
    line += '\n'.join(['{:.2f} {:.2f} {:e} {:.2f}'.format(sfh.data.lagei[i],
                                                          sfh.data.lagef[i],
                                                          sfh.data.sfr[i],
                                                          sfh.data.mh[i])
                       for i in range(ntbins)])
    if not os.path.isfile(outfile) or overwrite:
        with open(outfile, 'w') as f:
            f.write(line)
        print('wrote {:s}'.format(outfile))
    return


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=make_fakeparam.__doc__)

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite existing fake param file')

    parser.add_argument('--outfile', default=None,
                        help='filename to write to')

    parser.add_argument('param', default=None,
                        help='calcsfh input parameter file')

    parser.add_argument('sfh', default=None,
                        help='SFH file (output of zcombine)')

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    make_fakeparam(args.param, args.sfh, outfile=args.outfile,
                   overwrite=args.overwrite)
    return


if __name__ == "__main__":
    sys.exit(main())
