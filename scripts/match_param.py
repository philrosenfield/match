"""
Make a match param file
Note -- will have to edit the param file by hand to insert dmod and Av.
"""
#!/usr/bin/env python
import argparse
import matplotlib.pylab as plt
import numpy as np
import os
import sys
import time

from astropy.io import fits
from astropy.table import Table

from .diagnostics import match_diagnostic
from .fileio import match_param_default_dict, match_param_fmt
from .fileio import replace_ext, parse_pipeline

def move_on(ok, msg='0 to move on: '):
    ok = int(raw_input(msg))
    time.sleep(1)
    return ok

def exclude_tpagb(phot, param, mtrgb):
    from scipy.interpolate import interp1d

    lines = open(param, 'r').readlines()
    dmag, dcol, _, colmin, colmax = map(float, lines[4].split()[:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])

    maxcol = colmax - dcol
    brightlim = mag1min + dmag
    mag1, mag2 = np.loadtxt(phot, unpack=True)
    # choose color cut
    mincol = find_match_limits(mag1, mag2, color_only=True)[0]

    xarr, yarr = put_a_line_on_it(mtrgb, consty=False, filter1=True)
    f = interp1d(xarr, yarr)
    faintlim1 = f(mincol)
    faintlim2 = f(maxcol)
    exclude_gate = '1 {0} {2} {0} {3} {1} {3} {1} {4} 0 \n'.format(mincol,
                                                                   maxcol,
                                                                   faintlim1,
                                                                   brightlim,
                                                                   faintlim2)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!
    print(exclude_gate)
    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        [outp.write(l) for l in lines]
    print('wrote %s' % param)

def put_a_line_on_it(val, npts=100, consty=True, ax=None,
                     ls='--', annotate=True, filter1=None,
                     annotate_fmt='$TRGB=%.2f$', xmin=-3, xmax=6,
                     ymin=14, ymax=32, ann_kwargs={}, pltkw={}):
    """
    if consty is True: plots a constant y value across ax.xlims().
    if consty is False: plots a constant x on a plot of y vs x-y
    """
    xarr = np.linspace(xmin, xmax, npts)
    yarr = np.linspace(ymin, ymax, npts)
    if consty:
        # just a contsant y value over the plot range of x.
        new_xarr = xarr
    else:
        # a plot of y vs x-y and we want to mark
        # where a constant value of x is
        # e.g, f814w vs f555-f814; val is f555
        new_xarr = val - yarr
        # e.g, f555w vs f555-f814; val is f814
        if filter1 is True:
            yarr = xarr + val
            new_xarr = xarr

    if ax is not None:
        ax.plot(new_xarr, yarr, ls, **pltkw)
        if annotate:
            xy = (new_xarr[-1] - 0.1, yarr[-1] - 0.2)
            ax.annotate(annotate_fmt % val, xy=xy, ha='right', fontsize=16,
                        **ann_kwargs)
    return new_xarr, yarr


def find_match_limits(mag1, mag2, comp1=90., comp2=90., color_only=False,
                      xlim=None, ylim=None):
    """
    click color limits on a cmd and mag1 mag2 limits on a plot of mag1 vs mag2
    """
    col = mag1 - mag2

    fig, ax = plt.subplots()
    ax.plot(col, mag2, 'o', color='k', ms=3, alpha=0.3, mec='none')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ax.get_ylim()[::-1])

    if comp1 < 90.:
        ax.hlines(comp2, *ax.get_xlim())

    ok = 1
    while ok == 1:
        print 'click color extrema'
        pts = plt.ginput(2, timeout=-1)
        colmin, colmax = [pts[i][0] for i in range(2)]
        if colmin > colmax:
            colmin, colmax = colmax, colmin
        ax.vlines(colmin, *ax.get_ylim())
        ax.vlines(colmax, *ax.get_ylim())
        plt.draw()
        ok = move_on(0)

    plt.close()

    inds, = np.nonzero((col < colmax) & (col > colmin))
    data = (colmin, colmax)
    if not color_only:
        fig, ax = plt.subplots()
        ax.plot(mag1, mag2, '.', color='k')
        ok = 1
        while ok == 1:
            print 'click mag extrema'
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
            ok = move_on(ok)

        plt.close()

        inds, = np.nonzero((mag1 < mag1min) & (mag1 > mag1max) &
                           (mag2 < mag2min) & (mag2 > mag2max) &
                           (col < colmax) & (col > colmin))

    fig, ax = plt.subplots()
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

    print data
    return data


def find_gates(mag1, mag2, param):
    col = mag1 - mag2

    lines = open(param, 'r').readlines()
    colmin, colmax = map(float, lines[4].split()[3:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])
    #mag2min, mag2max = map(float, lines[5].split()[:-1])
    # click around
    fig, ax = plt.subplots()
    ax.plot(col, mag2, ',', color='k', alpha=0.2)
    ax.set_ylim(mag1max, mag1min)
    ax.set_xlim(colmin, colmax)

    ok = 1
    while ok == 1:
        print 'click '
        pts = np.asarray(plt.ginput(n=4, timeout=-1))
        exclude_gate = '1 {} 0 \n'.format(' '.join(['%.4f' % p for p in pts.flatten()]))
        pts = np.append(pts, pts[0]).reshape(5,2)
        ax.plot(pts[:,0], pts[:,1], color='r', lw=3, alpha=0.3)
        plt.draw()
        ok = move_on(0)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!

    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        [outp.write(l) for l in lines]
    print('wrote %s' % param)


def match_param(mag1, mag2, filters, phot, param, interactive=False, fake=None,
                comp_frac=0.5, param_kw={}):
    print filters
    p = dict(match_param_default_dict().items() + param_kw.items())

    if interactive:
        p['V-Imin'], p['V-Imax'], p['Vmin'], p['Vmax'], p['Imin'], p['Imax'] = \
            find_match_limits(mag1, mag2)
    else:
        color = mag1 - mag2
        p['V-Imin'] = np.min(color)
        p['V-Imax'] = np.max(color)
        p['Vmin'] = np.min(mag1)
        p['Imin'] = np.min(mag2)
        p['Vmax'] = np.max(mag1)
        p['Imax'] = np.max(mag2)
        print p['Vmax'], p['Imax']
        if fake is not None:
            from .asts import ASTs
            ast = ASTs(fake)
            ast.completeness(combined_filters=True, interpolate=True)
            print('Using {} completeness fraction from {}'.format(comp_frac, fake))
            p['Vmax'], p['Imax'] = ast.get_completeness_fraction(comp_frac)
            print p['Vmax'], p['Imax']

    p['V'] = filters[0]
    p['I'] = filters[1]
    with open(param, 'w') as out:
        out.write(match_param_fmt() % p)
    print('wrote {}'.format(param))
    return param

def match_limits(mag1, mag2, color_only=False, comp1=99., comp2=99.):
    plt.ion()
    ok = 1
    while ok == 1:
        data = find_match_limits(mag1, mag2, comp2=comp2, comp1=comp1,
                                 color_only=color_only)
        ok = move_on(0)

    plt.close()

    if color_only:
        colmin, colmax = data
        data_str = '%.2f %.2f' % (colmin, colmax)
    else:
        colmin, colmax, mag1max, mag2max = data
        data_str = '%.2f %.2f %.2f %.2f' % (colmin, colmax, mag1max, mag2max)

    print data_str


def make_phot(fitsfile, filters):
    assert len(filters) > 0, 'need filters'

    tab = Table.read(fitsfile)

    test = [f in tab.colnames for f in filters]
    if False in test:
        print(tab.colnames)
        sys.exit()
    else:
        pref, ext = fitsfile.split('.')[0], '.'.join(fitsfile.split('.')[1:])
        ext = '.{}'.format(ext.replace('.fits', ''))
        photname = '{0}_{1}_{2}{3}.match'.format(pref, filters[0], filters[1], ext)
        photname = photname.replace('_VEGA', '')
        np.savetxt(photname, np.column_stack((tab[filters[0]], tab[filters[1]])),
                   fmt='%.3f')
        print('wrote {}'.format(photname))
    return photname


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find color mag limits of CMDs interactively")

    parser.add_argument('-x', '--color_only', action='store_true',
                        help='with -i skip the magnitude finding')

    parser.add_argument('-s', '--slice', type=float, default=99.,
                        help='cut out mags outside of this value')

    parser.add_argument('-m', '--mtrgb', type=float, default=None,
                        help='trgb magnitude (run with -t)')

    parser.add_argument('-e', '--exgates', action='store_true',
                        help='find exclude gates (will force -i')

    parser.add_argument('-t', '--extpgates', action='store_true',
                        help='with -m exclude tp-agb (will force -i)')

    parser.add_argument('-i', '--interactive', action='store_true',
                        help='find limits interactively')

    parser.add_argument('-f', '--filters', type=str, default=None,
                        help='comma separated fits mag column names to make match phot')

    parser.add_argument('-a', '--fake', type=str, default=None,
                        help='match fake file')

    parser.add_argument('-c', '--comp_frac', type=float, default=None,
                        help='completeness fraction to use for faint limit (run with --fake=)')

    parser.add_argument('-p', '--param', type=str, help='match param file if existing')

    parser.add_argument('phot', type=str, help='match phot file or fits file')

    args = parser.parse_args(sys.argv[1:])

    if args.phot.endswith('fits'):
        filters = args.filters.split(',')
        args.phot = make_phot(args.phot, filters)

    mag1, mag2 = np.loadtxt(args.phot, unpack=True)
    inds, = np.nonzero((np.abs(mag1) < args.slice) & (np.abs(mag2) < args.slice))
    mag1 = mag1[inds]
    mag2 = mag2[inds]

    args.param = args.param or replace_ext(args.phot, '.param')
    if args.filters is None:
        target, filters = parse_pipeline(args.phot)

    if not os.path.isfile(args.param):
        print('Making param file')
        match_param(mag1, mag2, filters, args.phot, args.param, fake=args.fake,
                    interactive=args.interactive, comp_frac=args.comp_frac)

    if args.exgates:
        find_gates(mag1, mag2, args.param)

    if args.extpgates:
        exclude_tpagb(args.phot, args.param, args.mtrgb)

    match_diagnostic(args.param, args.phot)
