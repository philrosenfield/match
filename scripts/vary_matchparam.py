import argparse
import itertools
import sys

import numpy as np

def vary_matchparam(param_file, imfarr, bfarr, davarr):
    lines = open(param_file).readlines()

    for imf, bf, dav in itertools.product(imfarr, bfarr, davarr):
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
       new_name = '{}_imf{}_bf{}_dav{}.{}'.format(pname, imf, bf, dav, ext)

       with open(new_name, 'w') as outp:
          outp.write(''.join(lines))
          print 'wrote {}'.format(new_name)
    return


def main(argv):
    """
    With an input match parameter file, replace values and create new
    parameter files over a range of IMF, BF, dAv.

    With dAv: will only add dav to the new file name. The calcsfh caller
    can be used for that....

    Maybe a better solution is to combine the calsfh call here too...
    """
    parser = argparse.ArgumentParser(description="vary BF, IMF, or dAv")

    parser.add_argument('param_file', type=str,
                        help='template parameter file')

    args = parser.parse_args(argv)

    # search ranges
    imfarr = np.arange(-1., 3., 1)
    bfarr = np.arange(0.0, 1.05, 0.3)
    davarr = np.arange(0, 1.1, 0.5)

    vary_matchparam(args.param_file, imfarr, bfarr, davarr)

if __name__ == "__main__":
    main(sys.argv[1:])
