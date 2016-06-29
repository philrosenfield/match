"""Vary the IMF, BF, -dAv, -sub calls to calcsfh"""
from __future__ import print_function
import argparse
import itertools
import os
import sys

from config import calcsfh, calcsfh_flag, OUTEXT, SCRNEXT
from utils import splitext, writeorappend, parse_argrange
from fileio import read_calcsfh_param, calcsfh_input_parameter

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

def vary_matchparam(param_file, varyarrs=None, power_law_imf=True,
                    params=None):
    """
    Vary parameters from a match param template file.
    param_file : string
        calcsfh input (aka parameter) file
    varyarrs : dict
        a dictionary of array values where each key is XXXarr where
        XXX is a key in calcsfh_input_parameter

    Returns
    -------
    new_names : list
        list of string new parameter file names (with new parameters in the
        filename)
    """
    new_names = []
    varyarrs = {} or varyarrs
    params = {} or params

    pname, ext = splitext(param_file)
    template = dict(read_calcsfh_param(param_file).items() + params.items())
    # using tbin, tmin, tmax:
    del template['ntbins']

    for vals in itertools.product(*varyarrs.values()):
        name = []
        for i, val in enumerate(vals):
            key = varyarrs.keys()[i].replace('arr', '')
            template[key] = val
            name.append('{}{}'.format(key, val))
        new_param = calcsfh_input_parameter(power_law_imf=power_law_imf,
                                            **template)
        new_name = '{}_{}.{}'.format(pname, '_'.join(name), ext)

        with open(new_name, 'w') as outp:
            outp.write(new_param)
        print('wrote {}'.format(new_name))
        new_names.append(new_name)
    return new_names

def vary_matchparam1(param_file, imfarr, bfarr):
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

    parser.add_argument('--outfile', type=str, default='calcsfh_ssp.sh',
                        help='file to save the script')

    parser.add_argument('--imf', nargs='?', default=[-1, 3, 1],
                        help='IMF min, max, dIMF or one value/string')

    parser.add_argument('--bf', type=float, nargs='*', default=[0, 1., 0.3],
                        help='BF min, max, dBF or list')

    parser.add_argument('--dav', type=float, nargs='*', default=[0, 1.1, 0.5],
                        help='dAv min, max, ddAv or list')

    parser.add_argument('--tbin', type=float, nargs='*', default=[0.5],
                        help='dAv min, max, ddAv')

    parser.add_argument('--vstep', type=float, nargs='*', default=[0.1],
                        help='dAv min, max, ddAv')

    parser.add_argument('--vistep', type=float, nargs='*', default=[0.05],
                        help='dAv min, max, ddAv')

    parser.add_argument('--tmin', type=float, default=6.6,
                        help='min log age')

    parser.add_argument('--tmax', type=float, default=10.25,
                        help='max log age')

    parser.add_argument('--sub', type=str,
                        help='track sub directory')

    parser.add_argument('-e', '--extra', type=str, default='',
                        help='add an extra string to output filenames.')

    parser.add_argument('--destination', type=str, default=None,
                        help='destination directory for calcsfh output')

    parser.add_argument('--check', action='store_true',
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

    import pdb; pdb.set_trace()
    subs = parse_argrange(args.sub)
    davs = parse_argrange(args.dav)

    cparams = {'tmin': args.tmin, 'tmax': args.tmax}

    # write the parameter files
    varyarrs = {'bfarr': parse_argrange(args.bf),
                'tbinarr': parse_argrange(args.tbin), # np.array([0.01, 0.05, 0.1, 0.5]),
                'v-isteparr': parse_argrange(args.vistep), # np.array([0.01, 0.05, 0.1, 0.15]),
                'vsteparr': parse_argrange(args.vstep)} #np.array([0.01, 0.05, 0.1, 0.15])}

    power_law_imf = False
    imf = args.imf
    if not isinstance(imf, str):
        varyarrs['imfarr'] = parse_argrange(args.imf)
        power_law_imf = True

    params = vary_matchparam(args.param_file, varyarrs=varyarrs, params=cparams,
                             power_law_imf=power_law_imf)

    # loop over all to create output filenames and calcsfh calls
    line = ''
    nproc = 0
    for sub in subs:
        subfmt = ''
        if sub is not None:
            subfmt = '_{}'.format(sub)
        for dav in davs:
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
