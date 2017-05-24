import argparse
import os
import sys
import numpy as np


def check_boundaries(param, scrn):
    """
    check if the best fit file hits the Av or dmod edge of parameter search
    space.

    print information to terminal, nothing is printed if dmod and Av are
    within bounds.

    Parameters
    ----------
    param : match parameter file
    scrn : match console output (saved as a file)
    """
    def betweenie(val, upper, lower, retval=0, msg=''):
        if upper == lower:
            msg += 'value unchanging'
        if val >= upper:
            msg += 'extend boundary higher, %f >= %f\n' % (val, upper)
            retval += 1
        if val <= lower:
            msg += 'extend boundary lower, %f <= %f\n' % (val, lower)
            retval += 1
        return retval, msg

    msg = '{0:s} / {1:s}\n'.format(os.path.split(param)[1],
                                   os.path.split(scrn)[1])
    # parse scrn
    bfit = open(scrn).readlines()[-1]
    if 'Best' not in bfit:
        msg += 'error calcsfh not finished'
        retval = 1
    else:
        av, dmod, _ = bfit.split(',')
        dmod = float(dmod.replace(' dmod=', ''))
        av = float(av.replace('Best fit: Av=', ''))

        # parse param
        pars = open(param).readline()
        try:
            dmod0, dmod1, ddmod, av0, av1, dav = \
                np.array(pars.split(), dtype=float)
        except ValueError:
            imf, dmod0, dmod1, ddmod, av0, av1, dav = \
                np.array(pars.split(), dtype=float)
            # print(sys.exc_info()[1], param)
            # raise
        retval, msg = betweenie(dmod, dmod1, dmod0, msg=msg)
        retval, msg = betweenie(av, av1, av0, retval=retval, msg=msg)

    if retval > 0:
        print(msg)
    else:
        print('{}: Looks good.'.format(msg))
    return


def main(argv=None):
    parser = argparse.ArgumentParser(description="check simulation bounds")

    parser.add_argument('param', type=str, help='calcsfh input parameter file')

    parser.add_argument('scrn', type=str, help='calcsfh terminal output file')

    args = parser.parse_args(argv)

    return check_boundaries(args.param, args.scrn)


if __name__ == "__main__":
    sys.exit(main())
