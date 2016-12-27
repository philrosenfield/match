import argparse
import sys
from match.scripts.graphics.match_diagnostic import match_diagnostic

def parse_args(argv):
    parser = argparse.ArgumentParser(description="make a diagnostic plot of the calcsfh search space")

    parser.add_argument('phot', type=str,
                        help='photometry')

    parser.add_argument('param', type=str,
                        help='calcsfh input parameter file')

    parser.add_argument('--fake', type=str,
                        help='fake file')

    parser.add_argument('--xlim', type=str,
                        help=', separated xlimits')

    parser.add_argument('--ylim', type=str,
                        help=', separated ylimits')

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    match_diagnostic(args.param, args.phot, fake=args.fake, xlim=args.xlim,
                     ylim=args.ylim)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
