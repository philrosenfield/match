"""Vary the IMF, BF, -dAv, -sub calls to calcsfh"""
from __future__ import print_function
import argparse
import itertools
import os
import sys

import numpy as np

try:
    from .config import calcsfh, calcsfh_flag, OUTEXT, SCRNEXT
    from .utils import splitext, writeorappend, parse_argrange
    from .fileio import read_calcsfh_param, calcsfh_input_parameter
except SystemError:
    from config import calcsfh, calcsfh_flag, OUTEXT, SCRNEXT
    from utils import splitext, writeorappend, parse_argrange
    from fileio import read_calcsfh_param, calcsfh_input_parameter


def write_slurm(cmdscript, outdir='slurm', istart=1, use_bg=False):
    """
    Read a calcsfh_ssp script (output of varyparams) and split calls into
    separate calcsfh_*.sh scripts to be called in a job array.

    All outputs of this function will be put in [outdir].
    The match parameter, phot, and fake files referenced in the cmdscript
    will also be copied into [outdir] with the idea it will be tar'd up and
    shipped to odyssey.

    calcsfh.slurm is also written (all hardcoded but job array size)
    for sbatch submission.

    All the ssp.scrn and ssp.out file names are changed to have the same job
    array index as the calcsfh_*.sh file for fast debugging if a job gets
    killed or fails.
    """
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    with open(cmdscript, 'r') as inp:
        cmds = inp.readlines()
    cmds = [c.replace('&', '') for c in cmds if 'wait' not in c and len(c) > 0]

    iend = istart + len(cmds) - 1

    copys = []
    num = istart
    for i, cmd in enumerate(cmds):
        # num = ('{0:d}'.format(i + 1)).zfill(4)
        _, param, phot, fake = cmd.split()[:4]
        cmdfn = os.path.join(outdir, 'calcsfh_{0:d}.sh'.format(num))
        cmd = cmd.replace('ssp.', 'ssp{0:d}.'.format(num))
        with open(cmdfn, 'w') as outp:
            outp.write(cmd)
        num += 1
        copys.append(param)
        copys.append(phot)
        copys.append(fake)

    if use_bg:
        copys.append('bg.dat')

    for c in np.unique(copys):
        os.system('cp {0:s} {1:s}/{0:s}'.format(c, outdir))

    line = "#!/bin/bash\n\n"
    line += "#SBATCH -n 2\n"
    line += "#SBATCH -N 1\n"
    line += "#SBATCH -t 36:00:00\n"
    line += "#SBATCH --mem 60000\n"
    line += "#SBATCH -p conroy\n"
    line += "#SBATCH --array={0:d}-{1:d}\n".format(istart, iend)
    line += "#SBATCH -o calcsfh_%a.o\n"
    line += "#SBATCH -e calcsfh_%a.e\n"
    line += "#SBATCH --mail-type=END,FAIL\n"
    line += "#SBATCH --mail-user=philip.rosenfield@cfa.harvard.edu\n\n"
    line += "bash calcsfh_${SLURM_ARRAY_TASK_ID}.sh\n"
    line += "exit\n"

    with open(os.path.join(outdir, 'calcsfh.slurm'), 'w') as outp:
        outp.write(line)


def getflags(dav=None, sub=None, imf=None):
    """Add -dAv, -sub, and/or -kroupa or -chabrier to config.calcsfh_flag"""
    flag = calcsfh_flag

    if dav is not None:
        flag += " -dAv={0:.2f}".format(dav)

    if sub is not None:
        flag += "  -sub={0:s}".format(sub)

    if isinstance(imf, str) and imf == 'kroupa' or imf == 'chabrier':
        flag += " -{0:s}".format(imf)

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

    power_law_imf : bool
        passed to calcsfh_input_parameter

    params : dict
        parameters to overwite param_file with but not vary. (probably
        tmin, tmax)

    Returns
    -------
    new_names : list
        list of string new parameter file names (with new parameters in the
        filename) that were written
    """
    new_names = []
    varyarrs = {} or varyarrs
    params = {} or params

    pname, ext = splitext(param_file)
    template = read_calcsfh_param(param_file)
    template.update(params)

    # force using tbin, tmin, tmax:
    del template['ntbins']

    for vals in itertools.product(*varyarrs.values()):
        # Add the varied parameters to the filename
        name = []
        for i, val in enumerate(vals):
            key = list(varyarrs.keys())[i].replace('arr', '')
            template[key] = val
            name.append('{0:s}{1:g}'.format(key, val))

        new_param = calcsfh_input_parameter(power_law_imf=power_law_imf,
                                            **template)
        new_name = '{0:s}_{1:s}.{2:s}'.format(pname, '_'.join(np.sort(name)), ext)

        with open(new_name, 'w') as outp:
            outp.write(new_param)
        # print('wrote {0:s}'.format(new_name))
        new_names.append(new_name)
    return new_names


def vary_calcsfh_calls(phot, fake, params, outfile, subs, davs,
                       calcsfh=calcsfh, destination=None, check=False,
                       nproc=12, extra='', imf=None):
    """
    loop over param files to create output filenames and calcsfh calls for
    parameters that vary that are not in the calcsfh parameter file.
    """
    runtot = len(params) * len(davs) * len(subs)
    print('Requested {0:d} calcsfh calls'.format(runtot))

    line = ''
    inproc = 0
    for sub in subs:
        subfmt = ''
        if sub is not None:
            subfmt = '_{0:s}'.format(sub)
        for dav in davs:
            for param in params:
                parfile = param
                if destination is not None:
                    parfile = os.path.join(destination,
                                           os.path.split(param)[1])
                prefx, _ = splitext(parfile)
                suffx = 'dav{0:g}{1:s}{2:s}_ssp'.format(dav, subfmt, extra)

                name = '_'.join([prefx, suffx])

                out = '{0:s}{1:s}'.format(name, OUTEXT)
                scrn = '{0:s}{1:s}'.format(name, SCRNEXT)

                if check and os.path.isfile(out):
                    print('{0:s} exists, not overwriting'.format(out))
                else:
                    inproc += 1
                    flags = getflags(dav, sub=sub, imf=imf)
                    line += ' '.join([calcsfh, param, phot, fake,
                                      out, flags, '>', scrn, '&\n'])
                    if nproc == inproc:
                        line += 'wait \n'
                        inproc = 0

    writeorappend(outfile, line)


def main(argv):
    """
    With an input match parameter file, replace values and create new
    parameter files over a range of IMF, BF, dAv.
    """
    parser = argparse.ArgumentParser(description="vary BF, IMF, and dAv",
                                     fromfile_prefix_chars='@')

    parser.add_argument('-n', '--nproc', type=int, default=12,
                        help='number of simultaneous calls to calcsfh')

    parser.add_argument('--outfile', type=str, default='calcsfh_ssp.sh',
                        help='file to save the script')

    parser.add_argument('--imf', nargs='*', default=[-1, 3, 1],
                        help='IMF min, max, dIMF or one value/string')

    parser.add_argument('--bf', type=float, nargs='*', default=[0, 1., 0.3],
                        help='BF min, max, dBF or list')

    parser.add_argument('--dav', type=float, nargs='*', default=[0, 1.1, 0.5],
                        help='dAv min, max, ddAv or list')

    parser.add_argument('--tbin', type=float, nargs='*', default=[0.01],
                        help='time min, max, ddAv')

    parser.add_argument('--vstep', type=float, nargs='*', default=[0.15],
                        help='vstep min, max, ddAv')

    parser.add_argument('--vistep', type=float, nargs='*', default=[0.05],
                        help='v-i step min, max, ddAv')

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

    parser.add_argument('--param_file', type=str,
                        help='template parameter file (see .match_param)')

    parser.add_argument('--phot', type=str,
                        help='match photometry file (see .match_phot)')

    parser.add_argument('--fake', type=str,
                        help='match ast file (see .asts)')

    parser.add_argument('--calcsfh', type=str, default=calcsfh,
                        help='over ride default calcsfh base [config.calcsfh]')

    parser.add_argument('--slurm', action='store_true',
                        help='convert output script to slurm job array.')

    parser.add_argument('--slurmstart', type=int, default=1,
                        help='start the slurm numbering [1]')

    parser.add_argument('--use_bg', action='store_true',
                        help='use backgound file named bg.dat with smoothing=5')

    args = parser.parse_args(argv)

    if args.destination is not None:
        assert (os.path.isdir(args.destination)), \
            'destination must be an existing directory'

    extra = ''
    if len(args.extra) > 0:
        extra = '_{0:s}'.format(args.extra)

    # get the array or calculate arange
    subs = parse_argrange(args.sub)
    davs = parse_argrange(args.dav)
    varyarrs = {'bfarr': parse_argrange(args.bf),
                'tbinarr': parse_argrange(args.tbin),
                'visteparr': parse_argrange(args.vistep),
                'vsteparr': parse_argrange(args.vstep)}
    power_law_imf = False
    imf = args.imf
    if not isinstance(imf, str):
        # i.e, not kroupa or chabrier
        varyarrs['imfarr'] = parse_argrange(imf)
        power_law_imf = True

    cparams = {'tmin': args.tmin, 'tmax': args.tmax}

    if args.use_bg:
        cparams.update({'use_bg': True, 'bg_file': 'bg.dat', 'bg_smooth': 5})

    # do the "internal" matchparam varying
    params = vary_matchparam(args.param_file, varyarrs=varyarrs,
                             params=cparams, power_law_imf=power_law_imf)

    vary_calcsfh_calls(args.phot, args.fake, params, args.outfile, subs, davs,
                       calcsfh=args.calcsfh, destination=args.destination,
                       nproc=args.nproc, extra=extra, imf=imf)

    if args.slurm:
        if args.calcsfh == calcsfh:
            print('Warning: calcsfh path is default -- {0:s}'.format(calcsfh))
        write_slurm(args.outfile, istart=args.slurmstart, use_bg=args.use_bg)
    return


if __name__ == "__main__":
    main(sys.argv[1:])
