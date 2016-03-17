import argparse
import itertools
import sys

import numpy as np

calcsfh="$HOME/match2.6/bin/calcsfh"

def getflags(dav=0.0, sub=None, imf=None):
    # flag = "-PARSEC -ssp -full -dAvy=0.0 -dAv={:.2f}".format(dav)
    flag = "-PARSEC -ssp -dAvy=0.0 -dAv={:.2f}".format(dav)
    if sub is not None:
        flag += "  -sub={}".format(sub)
    if imf is not None:
        print('IMF flag not implemented since IMF array is hard coded and match param need proper formatting')
        #if imf == 'kroupa' or imf == 'chabrier':
        #    flag += " -{}".format(imf)
    return flag

def vary_matchparam(param_file, imfarr, bfarr):
    new_names = []
    lines = open(param_file).readlines()

    for imf, bf in itertools.product(imfarr, bfarr):
       imfline = lines[0].split()
       if len(imfline) == 6:
          imfline.insert(0, '{}'.format(imf))
       elif len(imfline) > 6 :
          imfline[0] = '{:.1f}'.format(imf)

       newimfline = ' '.join(imfline) + '\n'
       lines[0] = newimfline

       bfline = lines[2].split()
       bfline[0] = '{:.1f}'.format(bf)

       newbfline = ' '.join(bfline) + '\n'
       lines[2] = newbfline

       pname, ext = '.'.join(param_file.split('.')[:-1]), param_file.split('.')[-1]
       new_name = '{}_imf{}_bf{}.{}'.format(pname, imf, bf, ext)

       with open(new_name, 'w') as outp:
          outp.write(''.join(lines))
          print 'wrote {}'.format(new_name)
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

    parser.add_argument('-b', '--bf', type=str, default='0,1.05,0.3',
                        help='BF min, max, dBF')

    parser.add_argument('-a', '--dav', type=str, default='0,1.1,0.5',
                        help='dAv min, max, ddAv')

    parser.add_argument('-s', '--sub', type=str,
                        help='track sub directory')

    parser.add_argument('param_file', type=str,
                        help='template parameter file')

    parser.add_argument('phot', type=str,
                        help='match photometry file')

    parser.add_argument('fake', type=str,
                        help='match ast file')

    args = parser.parse_args(argv)

    imfarr = np.arange(*map(float, args.imf.split(',')))
    bfarr = np.arange(*map(float, args.bf.split(',')))
    davarr = np.arange(*map(float, args.dav.split(',')))
    line = ''

    subs = args.sub.replace(' ','').split(',')
    n = 0
    for sub in subs:
        for dav in davarr:
            params = vary_matchparam(args.param_file, imfarr, bfarr)
            for param in params:
                n += 1
	        out = '.'.join(np.concatenate([param.split('.')[:-1], ['{}_ssp.out'.format(sub)]]))
                scrn = '.'.join(np.concatenate([param.split('.')[:-1], ['{}_ssp.scrn'.format(sub)]]))
                flags = getflags(dav, sub)
                line += ' '.join([calcsfh, param, args.phot, args.fake, out, flags, '>', scrn, '&']) + '\n'
                if n == args.nproc:
                    line += 'wait \n'
                    n = 0

    with open(args.outfile, 'w') as outp:
        outp.write(line)
    print('wrote {}'.format(args.outfile))
    return

if __name__ == "__main__":
    main(sys.argv[1:])
