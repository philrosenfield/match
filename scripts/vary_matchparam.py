import argparse
import itertools
import os
import sys

import numpy as np

from config import calcsfh, calcsfh_flag

def getflags(dav=0.0, sub=None, imf=None):
    flag = calcsfh_flag
    if dav is not None:
        flag += " -dAv={:.2f}".format(dav)
    if sub is not None:
        flag += "  -sub={}".format(sub)
    if imf is not None:
        if imf == 'kroupa' or imf == 'chabrier':
            flag += " -{}".format(imf)
    return flag

def vary_matchparam(param_file, imfarr, bfarr):
    new_names = []
    lines = open(param_file).readlines()

    for imf, bf in itertools.product(imfarr, bfarr):
        imfline = lines[0].split()
        try:
            float(imf)
            if len(imfline) == 6:
                imfline.insert(0, '{}'.format(imf))
            elif len(imfline) > 6 :
                imfline[0] = '{:.2f}'.format(imf)
        except:
            pass

        newimfline = ' '.join(imfline) + '\n'
        lines[0] = newimfline

        bfline = lines[2].split()
        bfline[0] = '{:.2f}'.format(bf)

        newbfline = ' '.join(bfline) + '\n'
        lines[2] = newbfline

        pname, ext = '.'.join(param_file.split('.')[:-1]), param_file.split('.')[-1]
        new_name = '{}_imf{}_bf{}.{}'.format(pname, imf, bf, ext)

        with open(new_name, 'w') as outp:
            outp.write(''.join(lines))
            #print 'wrote {}'.format(new_name)
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

    parser.add_argument('param_file', type=str,
                        help='template parameter file')

    parser.add_argument('phot', type=str,
                        help='match photometry file')

    parser.add_argument('fake', type=str,
                        help='match ast file')

    args = parser.parse_args(argv)

    extra = ''
    if len(args.extra) > 0:
        extra = '_{}'.format(args.extra)


    if ',' in args.imf:
        # vary the IMF slope
        imfarr = np.arange(*map(float, args.imf.split(',')))
    else:
        imf = args.imf
        imfarr = [imf]

    bfarr = [0.0]
    if ',' in args.bf:
        bfarr = np.arange(*map(float, args.bf.split(',')))

    davarr = [0.0]
    if ',' in args.dav:
        davarr = np.arange(*map(float, args.dav.split(',')))

    subs = [None]
    if args.sub is not None:
        subs = args.sub.replace(' ','').split(',')

    # write the parameter files
    params = vary_matchparam(args.param_file, imfarr, bfarr)

    # loop over all to create output filenames and calcsfh calls
    line = ''
    n = 0
    for sub in subs:
        subfmt = ''
        if sub is not None:
            subfmt = '_{}'.format(sub)
        for dav in davarr:
            for param in params:
                n += 1
                outdir = os.path.split(param)[0]
                if args.destination is not None:
                    outdir = args.destination
                pname = os.path.split(param)[1]
                oname = os.path.join(outdir, pname)
                name = '_'.join(np.concatenate([['.'.join(oname.split('.')[:-1])],
                                               ['dav{}{}{}_ssp'.format(dav, subfmt, extra)]]))
                out = '{}.out'.format(name)
                scrn = '{}.scrn'.format(name)
                flags = getflags(dav, sub=sub, imf=imf)
                line += ' '.join([calcsfh, param, args.phot, args.fake, out, flags, '>', scrn, '&']) + '\n'
                if n == args.nproc:
                    line += 'wait \n'
                    n = 0

    wstr = 'w'
    wrote = 'wrote'
    if os.path.isfile(args.outfile):
        wstr = 'a'
        wrote = 'appended'
    with open(args.outfile, wstr) as outp:
        outp.write(line)
    print('{} {}'.format(wrote, args.outfile))
    return

if __name__ == "__main__":
    main(sys.argv[1:])
