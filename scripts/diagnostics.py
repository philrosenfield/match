from __future__ import print_function
import argparse
import os
import sys

from .fileio import get_files
from .graphics import call_pgcmd, match_diagnostic, sfh_plot
from .utils import MatchSFH, check_boundaries


def main(argv):
    parser = argparse.ArgumentParser(description="Plot match diagnostics")

    parser.add_argument('-f', '--filters', type=str, default=None,
                        help='comma separated filter names')

    parser.add_argument('-d', '--directory', type=str, default=os.getcwd(),
                        help='specify directory')

    parser.add_argument('-n', '--name', nargs='*', type=str, help='match cmd, sfh, zc\
                        file(s)')

    args = parser.parse_args(argv)

    if args.name is None:
        cmd_names = get_files(args.directory, '*cmd')
        sfh_files = get_files(args.directory, '*sfh')
        sfh_files.extend(get_files(args.directory, '*zc'))
        params = get_files(args.directory, '*.param')
        phots = get_files(args.directory, '*match')
        scrns = get_files(args.directory, '*scrn')
        scrns = [s for s in scrns if not 'mcmc' in s]
    else:
        cmd_names = [n for n in args.name if n.endswith('cmd')]
        sfh_files = [n for n in args.name if n.endswith('sfh')]
        sfh_files.extend([n for n in args.name if n.endswith('zc')])
        params = [n for n in args.name if n.endswith('param')]
        phots = [n for n in args.name if n.endswith('match')]
        scrns = [n for n in args.name if n.endswith('scrn')]
        scrns = [s for s in scrns if not 'mcmc' in s]

    [check_boundaries(p, s) for p, s in zip(params, scrns)]

    try:
        filter1, filter2 = args.filters.split(',')
    except AttributeError:
        filter1 = 'V'
        filter2 = 'I'

    labels = ['${\\rm %s}$' % i for i in ('data', 'model', 'diff', 'sig')]

    call_pgcmd(cmd_names, filter1, filter2, labels=labels)

    if len(sfh_files) > 0:
        for sfh_file in sfh_files:
            msfh = MatchSFH(sfh_file)
            if len(msfh.data) != 0:
                sfh_plot(msfh)
                msfh.plot_csfr()

    [match_diagnostic(params[i], phots[i]) for i in range(len(phots))]


if __name__ == "__main__":
    main(sys.argv[1:])
