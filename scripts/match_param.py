"""
Make a match param file
Note -- will have to edit the param file by hand to insert dmod and Av.
"""
#!/usr/bin/env python
from __future__ import print_function
import argparse
import os
import sys
import time

import numpy as np

import matplotlib.pylab as plt

from .config import PARAMEXT
from .fileio import read_calcsfh_param, calcsfh_input_parameter, match_filters
from .utils import replaceext, parse_pipeline
from .match_phot import make_phot
from .graphics import match_diagnostic

def move_on(okay, msg='0 to move on: '):
    """read and return raw input"""
    okay = int(raw_input(msg))
    time.sleep(1)
    return okay


def within_limits(params, fakefile, offset=1.):
    """
    Cull param cmd limits to that of the fake file
    params : dict
        calcsfh_input_parameter dictionary (only need CMD limits)
    fakefile : string
        match AST file
    offset : float
        mag below
    """
    vimin = params['vimin']
    vimax = params['vimax']
    vmin = params['vmin']
    imin = params['imin']
    vmax = params['vmax']
    imax = params['imax']

    mag1in, mag2in, _, _ = np.loadtxt(fakefile, unpack=True)
    colin = mag1in - mag2in
    msg = 'Overwrote'
    if vimin < colin.min():
        vimin = colin.min()
        msg += ' vimin'
    if vimax > colin.max():
        vimax = colin.max()
        msg += ' vimax'
    if vmin < mag1in.min():
        vmin = mag1in.min()
        msg += ' vmin'
    if vmax > mag1in.max() - 1.:
        vmax = mag1in.max() - 1.
        msg += ' vmax'
    if imin < mag2in.min():
        imin = mag2in.min()
        msg += ' imin'
    if imax > mag2in.max() - 1:
        imax = mag2in.max() - 1.
        msg += ' imax'
    msg += ' with values from matchfake'
    print(msg)
    params['vimin'] = vimin
    params['vimax'] = vimax
    params['vmin'] = vmin
    params['imin'] = imin
    params['vmax'] = vmax
    params['imax'] = imax
    return params

def find_match_limits(mag1, mag2, comp1=90., comp2=90., color_only=False,
                      xlim=None, ylim=None):
    """
    click color limits on a cmd and mag1 mag2 limits on a plot of mag1 vs mag2
    """
    col = mag1 - mag2

    _, ax = plt.subplots()
    ax.plot(col, mag2, 'o', color='k', ms=3, alpha=0.3, mec='none')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ax.get_ylim()[::-1])

    if comp1 < 90.:
        ax.hlines(comp2, *ax.get_xlim())

    okay = 1
    while okay == 1:
        print('click color extrema')
        pts = plt.ginput(2, timeout=-1)
        colmin, colmax = [pts[i][0] for i in range(2)]
        if colmin > colmax:
            colmin, colmax = colmax, colmin
        ax.vlines(colmin, *ax.get_ylim())
        ax.vlines(colmax, *ax.get_ylim())
        plt.draw()
        okay = move_on(0)

    plt.close()

    inds, = np.nonzero((col < colmax) & (col > colmin))
    data = (colmin, colmax)
    if not color_only:
        _, ax = plt.subplots()
        ax.plot(mag1, mag2, '.', color='k')
        okay = 1
        while okay == 1:
            print('click mag extrema')
            pts = plt.ginput(2, timeout=-1)
            mag1max, mag2max = pts[0]
            mag1min, mag2min = pts[1]
            if mag1min > mag1max:
                mag1min, mag1max = mag1max, mag1min
            if mag2min > mag2max:
                mag2min, mag2max = mag2max, mag2min

            ax.plot(mag1max, mag2max, 'o', color='r')
            ax.plot(mag1min, mag2min, 'o', color='r')
            plt.draw()
            okay = move_on(okay)

        plt.close()

        inds, = np.nonzero((mag1 < mag1min) & (mag1 > mag1max) &
                           (mag2 < mag2min) & (mag2 > mag2max) &
                           (col < colmax) & (col > colmin))

    _, ax = plt.subplots()
    ax.plot(col, mag2, '.', color='k')
    ax.plot(col[inds], mag2[inds], '.', color='r')
    ax.set_ylim(ax.get_ylim()[::-1])
    if comp2 < 90.:
        ax.hlines(comp2, *ax.get_xlim(), lw=2)
    ax.vlines(colmin, *ax.get_ylim(), lw=2)
    ax.vlines(colmax, *ax.get_ylim(), lw=2)
    if not color_only:
        ax.hlines(mag2max, *ax.get_xlim(), lw=2)
        data = (colmin, colmax, mag1min, mag1max, mag2min, mag2max)

    plt.draw()

    print(data)
    return data


def find_gates(mag1, mag2, param):
    """Click 4 points to make an exclude gate -- does not work in calcsfh!"""
    print('not supported')
    sys.exit()
    col = mag1 - mag2

    lines = open(param, 'r').readlines()
    colmin, colmax = map(float, lines[4].split()[3:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])
    #mag2min, mag2max = map(float, lines[5].split()[:-1])
    # click around
    _, ax = plt.subplots()
    ax.plot(col, mag2, ',', color='k', alpha=0.2)
    ax.set_ylim(mag1max, mag1min)
    ax.set_xlim(colmin, colmax)

    okay = 1
    while okay != 0:
        print('click ')
        pts = np.asarray(plt.ginput(n=4, timeout=-1))
        exclude_gate = '1 {} 0 \n'.format(' '.join(['%.4f' % p for p in pts.flatten()]))
        pts = np.append(pts, pts[0]).reshape(5, 2)
        ax.plot(pts[:, 0], pts[:, 1], color='r', lw=3, alpha=0.3)
        plt.draw()
        okay = move_on(0)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!

    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        _ = [outp.write(l) for l in lines]
    print('wrote %s' % param)


def match_param(mag1, mag2, filters, param, interactive=False, fake=None,
                comp_frac=0.5, param_kw=None, power_law_imf=False, zinc=False,
                bright_lim=20., clobber=False):
    """
    Make match param file

    Note:
    Will check filter list against templates/match_filters.json
    Will check the CMD limits against the AST limits (if fake is passed)

    mag1, mag2 : array, array
        v mag and i mag (extrema used for CMD limits)
    filters : list of strings
        v and i filter names
    param : string
        template parameter file or if clobber parameter file name
    interactive : bool
        choose cmd limits interactively
    fake : string
        matchfake filename if using comp_frac or want to check CMD limits
        are within fake limits (see FK overflow in MATCH README)
    comp_frac : float
        completeness fraction to set faint mag limit
    param_kw : dict
        parameters of template/calcsfh_input_parameter.json to write
    power_law_imf : bool
        passed to fileio.calcsfh_input_parameter
    zinc : bool
        passed to fileio.calcsfh_input_parameter
    bright_lim : float
        passed to asts.ast.get_completeness_fraction
    clobber : bool
        overwrite param file if exisiting
    """
    param_kw = param_kw or {}

    print('Using filters {}, {}'.format(*filters))
    if os.path.isfile(param) and not clobber:
        template = read_calcsfh_param(param)
        template.update(param_kw)
    else:
        template = param_kw

    if interactive:
        vimin, vimax, vmin, vmax, imin, imax = find_match_limits(mag1, mag2)
    else:
        color = mag1 - mag2
        vimin = np.min(color)
        vimax = np.max(color)
        vmin = np.min(mag1)
        imin = np.min(mag2)
        if fake is None:
            vmax = np.max(mag1)
            imax = np.max(mag2)
            dove = 'data'
        else:
            from .asts import ASTs
            ast = ASTs(fake)
            ast.completeness(combined_filters=True, interpolate=True)
            print('Using {} completeness fraction from {}'.format(comp_frac,
                                                                  fake))
            vmax, imax = ast.get_completeness_fraction(comp_frac,
                                                       bright_lim=bright_lim)
            dove = 'completeness'
        print('From {}: vmax={} imax={}'.format(dove, vmax, imax))

    template['vimin'] = vimin
    template['vimax'] = vimax
    template['vmin'] = vmin
    template['imin'] = imin
    template['vmax'] = vmax
    template['imax'] = imax
    if fake is not None:
        template = within_limits(template, fake)
    template['v'] = filters[0]
    template['i'] = filters[1]

    # HACK!! Is it WFC3 or UVIS?
    # (calcsfh_input_parameter will throw assertion for unrecognized filters)
    # First: UW Pipeline says FXXXW, MATCH would think that's WFPC2
    # So this won't work with WFPC2 data.
    #  If it's not WFC3, both are UVIS.
    # However, this is a bad way to do it, should have a -wfc or -uvis flag.
    itsuvis = False
    possible_filters = match_filters()['filters']
    if filters[0].startswith('F') and filters[0].endswith('W'):
        template['v'] = filters[0].replace('F', 'WFC')
    if not template['v'] in possible_filters:
        template['v'] = filters[0].replace('F', 'UVIS')
        itsuvis = True

    if filters[1].startswith('F') and filters[1].endswith('W'):
        template['i'] = filters[1].replace('F', 'WFC')
    if not template['i'] in possible_filters or itsuvis:
        template['i'] = filters[1].replace('F', 'UVIS')

    with open(param, 'w') as out:
        # see fileio.match_param_fmt for dictionary defaults
        out.write(calcsfh_input_parameter(power_law_imf=power_law_imf,
                                          zinc=zinc,
                                          **template))
    print('wrote {}'.format(param))
    return param


def main(argv):
    """main function for match_param"""
    parser = argparse.ArgumentParser(description="make calcsfh param file")

    parser.add_argument('--imf', default=None,
                        help='IMF power law value (None if using calcsfh flag)')

    parser.add_argument('--bf', type=float, default=0.0,
                        help='Binary fraction')

    parser.add_argument('--tbin', type=float, default=0.05,
                        help='age bin width(s)')

    parser.add_argument('--vstep', type=float, default=0.15,
                        help='mag step size')

    parser.add_argument('--vistep', type=float, default=0.05,
                        help='color step size')

    parser.add_argument('--tmin', type=float, default=6.6,
                        help='min log age')

    parser.add_argument('--tmax', type=float, default=10.24,
                        help='max log age')

    parser.add_argument('--zinc', action='store_true',
                        help='use zinc')

    parser.add_argument('--dmod', type=float, nargs=2,
                        help='dmod0, dmod1')

    parser.add_argument('--av', type=float, nargs=2,
                        help='av0, av1')

    parser.add_argument('--dav', type=float, default=0.05,
                        help='Av step -- NOT -dAv flag')

    parser.add_argument('--ddmod', type=float, default=0.10,
                        help='dmod step')

    parser.add_argument('-s', '--slice', type=float, default=99.,
                        help='cut out mags outside of this value')

    parser.add_argument('-i', '--interactive', action='store_true',
                        help='find cmd limits interactively')

    parser.add_argument('-f', '--filters', type=str, default=None,
                        help=('comma separated filter names (if filename does '
                              'not follow: PID_TARGET_FILTER1_FILTER2.ext)'))

    parser.add_argument('--fake', type=str, default=None,
                        help=('match fake file, used to calculate completeness '
                              'to set faint mag limits of param file'))

    parser.add_argument('-c', '--comp_frac', type=float, default=0.50,
                        help=('completeness fraction as faint mag limit '
                              '(use with --fake=)'))

    parser.add_argument('-b', '--bright_lim', type=float, default=20,
                        help=('bright limit to consider for calculating '
                              'completeness (use with --fake=)'))

    parser.add_argument('-p', '--param', type=str,
                        help='template match param file')

    parser.add_argument('--clobber', action='store_true',
                        help='overwrite')


    parser.add_argument('phot', type=str, help='photometry file match or fits')

    args = parser.parse_args(argv)
    assert(not isinstance(args.imf, str)), 'Only set IMF if it is a power law.'
    if args.phot.endswith('fits'):
        args.phot = make_phot(args.phot)[0]

    mag1, mag2 = np.loadtxt(args.phot, unpack=True)
    inds, = np.nonzero((np.abs(mag1) < args.slice) &
                       (np.abs(mag2) < args.slice))
    mag1 = mag1[inds]
    mag2 = mag2[inds]

    args.param = args.param or replaceext(args.phot, PARAMEXT)

    filters = args.filters
    if args.filters is None:
        _, filters = parse_pipeline(args.phot)

    if not os.path.isfile(args.param) or args.clobber:
        print('Making param file')
        param_kw = {'imf': args.imf,
                    'bf': args.bf,
                    'tbin': args.tbin,
                    'vstep': args.vstep,
                    'vistep': args.vistep,
                    'tmin': args.tmin,
                    'tmax': args.tmax,
                    'dmod0': args.dmod[0],
                    'dmod1': args.dmod[1],
                    'ddmod': args.ddmod,
                    'av0': args.av[0],
                    'av1': args.av[1],
                    'dav': args.dav}

        power_law_imf = True
        if args.imf is None:
            power_law_imf = False

        match_param(mag1, mag2, filters, args.param, fake=args.fake,
                    interactive=args.interactive, comp_frac=args.comp_frac,
                    param_kw=param_kw, power_law_imf=power_law_imf,
                    zinc=args.zinc, bright_lim=args.bright_lim,
                    clobber=args.clobber)

        match_diagnostic(args.param, args.phot, fake=args.fake)
    else:
        print('{} file found, not overwriting'.format(args.param))

if __name__ == "__main__":
    main(sys.argv[1:])
