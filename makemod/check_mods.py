""" something """
from __future__ import print_function
import sys
import os
import argparse
import glob
import numpy as np
try:
    from ..scripts.config import match_base
except:
    from scripts.config import match_base


def zlimits_from_makemod(model='PARSEC', dz=0.01):
    """Get zmin zmax from makemod"""
    def parse_line(string):
        """Parse the static const int modelIZm... line to get value"""
        return int(string.split('=')[-1].replace(';', ''))

    here = os.getcwd()
    modeldir = os.path.join(match_base, model)
    os.chdir(modeldir)

    with open('makemod.cpp') as inp:
        lines = inp.readlines()

    zminstr, = [l for l in lines if 'static const int modelIZmin' in l]
    zmaxstr, = [l for l in lines if 'static const int modelIZmax' in l]
    zmin = parse_line(zminstr) / 100. + 4. + dz / 2
    zmax = parse_line(zmaxstr) / 100. + 4. - dz / 2

    os.chdir(here)
    return zmin, zmax


def filenames_loop(logt0, logt1, dlogt, z0, z1, dz, mod='mod1'):
    """Loop through the data directory counting all mod? files"""
    ishere = 0
    missing = 0
    for logt in np.arange(logt0, logt1, dlogt):
        for z in np.arange(z0, z1, dz):
            fname = '{}_{:6.4f}_{:6.4f}_{:5.3f}_0.005'.format(mod, logt,
                                                              (logt + dlogt), z)
            res = glob.glob(fname)
            if len(res) == 0:
                print(fname)
                missing += 1
            else:
                ishere += 1
    print('{} {} files, {} missing'.format(ishere, mod, missing))


def check_mods(model='PARSEC', dz=0.01, dlogt=0.001, sub=None,
               logt0=6.6, logt1=10.25, mod='mod1'):
    """check the mod files completeness"""
    here = os.getcwd()
    datadir = os.path.join(match_base, model, 'data')
    if sub is not None:
        datadir = os.path.join(datadir, sub)
    assert os.path.isdir(datadir), 'Data directory not found {}'.format(datadir)

    os.chdir(datadir)
    z0, z1 = zlimits_from_makemod(model=model, dz=dz)
    filenames_loop(logt0, logt1, dlogt, z0, z1, dz, mod=mod)
    os.chdir(here)


desc = """
Check for complete makemod run
Usage: e.g., python -m match.makemod.check_mods -p -s ov0.50 -m mod2
will check all mod2 files from agelimits (-a) [0,1] in
[match_base]/[model]/data/[sub]

(given dz and dlogt)

match_base set in config.py
model is set to MIST unless -p will use PARSEC
sub is omitted unless -s flag.
"""


def main(argv):
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-z', '--dz', type=float, default=0.01,
                        help='metallicity resolution')

    parser.add_argument('-t', '--dlogt', type=float, default=0.001,
                        help='time resolution')

    parser.add_argument('-a', '--agelimits', type=float, nargs=2,
                        default=[7.001, 10.25],
                        help='makemod log age limits')

    parser.add_argument('-s', '--sub', type=str,
                        help='subdirectory of data/')

    parser.add_argument('-m', '--mod', type=str, default='mod1',
                        help='mod file prefix (mod1, mod2 etc)')

    parser.add_argument('-p', '--parsec', action='store_true',
                        help='PARSEC models (otherwise MIST)')

    args = parser.parse_args(argv)
    model = 'MIST'
    if args.parsec:
        model = 'PARSEC'

    check_mods(model=model, dz=args.dz, dlogt=args.dlogt, sub=args.sub,
               logt0=args.agelimits[0], logt1=args.agelimits[1],
               mod=args.mod)


if __name__ == "__main__":
    main(sys.argv[1:])
