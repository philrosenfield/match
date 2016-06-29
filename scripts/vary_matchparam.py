"""Vary the IMF, BF, -dAv, -sub calls to calcsfh"""
from __future__ import print_function
import argparse
import itertools
import os
import sys

from config import calcsfh, calcsfh_flag, OUTEXT, SCRNEXT
from utils import splitext, writeorappend, parse_argrange

def getflags(dav=0.0, sub=None, imf=None):
    """Add -dAv, -sub, and/or -kroupa or -chabrier to config.calcsfh_flag"""
    flag = calcsfh_flag
    if dav is not None:
        flag += " -dAv={:.2f}".format(dav)
    if sub is not None:
        flag += "  -sub={}".format(sub)
    if imf is not None:
        if imf == 'kroupa' or imf == 'chabrier':
            flag += " -{}".format(imf)
    return flag

#generalize this to use param values... should actually read in values of param...
def vary_matchparam(param_file, imfarr, bfarr):
    """
    Vary parameters from a match param template file.
    param_file : string
        calcsfh input (aka parameter) file
    imfarr : array
        imf values to vary
    bfarr : array
        binary fraction values to vary

    Returns
    -------
    new_names : list
        list of string new parameter file names (with new paramters in the
        filename)
    """
    new_names = []
    lines = open(param_file).readlines()

    for imf, bfrac in itertools.product(imfarr, bfarr):
        imfline = lines[0].strip().split()
        try:
            # imf is a power law
            float(imf)
            if len(imfline) == 6:
                # imf was not in the template param file
                imfline.insert(0, '{}'.format(imf))
            elif len(imfline) > 6:
                # new power law
                imfline[0] = '{:.2f}'.format(imf)
        except ValueError:
            # imf is -kroupa or -chabrier called from command line.
            if len(imfline) > 6:
                print('First line of param file formatted for powerlaw IMF')

        newimfline = ' '.join(imfline) + '\n'
        lines[0] = newimfline

        bfline = lines[2].split()
        # new binary fraction
        bfline[0] = '{:.2f}'.format(bfrac)

        newbfline = ' '.join(bfline) + '\n'
        lines[2] = newbfline

        # place new params
        pname, ext = splitext(param_file)
        new_name = '{}_imf{}_bf{}.{}'.format(pname, imf, bfrac, ext)

        with open(new_name, 'w') as outp:
            outp.write(''.join(lines))

        new_names.append(new_name)
    return new_names


def main(argv):
    """
    With an input match parameter file, replace values and create new
    parameter files over a range of IMF, BF, dAv.
    """
    parser = argparse.ArgumentParser(description="vary BF, IMF, and dAv")

    parser.add_argument('-n', '--nproc', type=int, default=12,
                        help='number of simultaneous calls to calcsfh')

    parser.add_argument('-o', '--outfile', type=str, default='calcsfh_ssp.sh',
                        help='file to save the script')

    parser.add_argument('-i', '--imf', type=str, default='-1,3,1',
                        help='IMF min, max, dIMF')

    parser.add_argument('-b', '--bf', type=str, default='0,1.00,0.3',
                        help='BF min, max, dBF')

    parser.add_argument('-a', '--dav', type=str, default='0,1.1,0.5',
                        help='dAv min, max, ddAv')

    parser.add_argument('-s', '--sub', type=str,
                        help='track sub directory')

    parser.add_argument('-e', '--extra', type=str, default='',
                        help='add an extra string to output filenames.')

    parser.add_argument('-d', '--destination', type=str, default=None,
                        help='destination directory for calcsfh output')

    parser.add_argument('-c', '--check', action='store_true',
                        help=('check if calcsfh output file exists, '
                              'useful for completeing interrupted runs.'))

    parser.add_argument('param_file', type=str,
                        help='template parameter file')

    parser.add_argument('phot', type=str,
                        help='match photometry file')

    parser.add_argument('fake', type=str,
                        help='match ast file')

    args = parser.parse_args(argv)

    if args.destination is not None:
        assert (os.path.isdir(args.destination)), \
            'destination must be an existing directory'

    extra = ''
    if len(args.extra) > 0:
        extra = '_{}'.format(args.extra)

    imfarr = parse_argrange(args.imf, args.imf)
    imf = args.imf
    bfarr = parse_argrange(args.bf, 0.0)
    davarr = parse_argrange(args.bf, 0.0)
    subs = parse_argrange(args.sub, None)

    # write the parameter files
    params = vary_matchparam(args.param_file, imfarr, bfarr)

    # loop over all to create output filenames and calcsfh calls
    line = ''
    nproc = 0
    for sub in subs:
        subfmt = ''
        if sub is not None:
            subfmt = '_{}'.format(sub)
        for dav in davarr:
            for param in params:
                parfile = param
                if args.destination is not None:
                    parfile = os.path.join(args.destination,
                                           os.path.split(param)[1])
                prefx, _ = splitext(parfile)
                suffx = 'dav{}{}{}_ssp'.format(dav, subfmt, extra)
                name = '_'.join([prefx, suffx])

                out = '{}{}'.format(name, OUTEXT)
                scrn = '{}{}'.format(name, SCRNEXT)

                if args.check and os.path.isfile(out):
                    print('{} exists, not overwriting'.format(out))
                else:
                    nproc += 1
                    flags = getflags(dav, sub=sub, imf=imf)
                    line += '{}\n'.format(' '.join([calcsfh, param, args.phot,
                                                    args.fake, out, flags, '>',
                                                    scrn, '&']))
                    if nproc == args.nproc:
                        line += 'wait \n'
                        nproc = 0

    writeorappend(args.outfile, line)
    return


if __name__ == "__main__":
    main(sys.argv[1:])
