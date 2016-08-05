""" Likelihood used in MATCH """
import numpy as np
from .utils import read_binned_sfh
import sys
import argparse


def stellar_prob(obs, model, normalize=False):
    '''
    FROM MATCH README
    The quality of the fit is calculated using a Poisson maximum likelihood
    statistic, based on the Poisson equivalent of chi^2.
      2 m                                if (n=0)
      2 [ 0.001 + n * ln(n/0.001) - n ]  if (m<0.001)
      2 [ m + n * ln(n/m) - n ]          otherwise
    m=number of model points; n=number of observed points

    This statistic is based on the Poisson probability function:
       P =  (e ** -m) (m ** n) / (n!),
    Recalling that chi^2 is defined as -2lnP for a Gaussian distribution and
    equals zero where m=n, we treat the Poisson probability in the same
    manner to get the above formula.

    '''
    n = obs
    m = model

    if normalize is True:
        n /= np.sum(n)
        m /= np.sum(m)

    d = 2. * (m + n * np.log(n / m) - n)

    smalln = np.abs(n) < 1e-10
    d[smalln] = 2. * m[smalln]

    smallm = (m < 0.001) & (n != 0)
    d[smallm] = 2. * (0.001 + n[smallm] *
                      np.log(n[smallm] / 0.001) - n[smallm])

    sig = np.sqrt(d) * np.sign(n - m)
    # fit = np.sum(sig)
    # pct_dif = (m - n) / n
    # sum(d) = -2 ln P
    prob = np.exp(-1 * np.sum(d) / 2)
    # prob = np.sum(d) / float(len(np.concatenate(n)) - 1)
    return d, prob  # , pct_dif, sig


def match_stats(sfh_file, match_cmd_file, nfp_nonsfr=5, nmc_runs=10000,
                outfile='cmd_stats.dat', dry_run=False, extra=''):
    '''
    NFP = # of non-zero time bins
          + dmod + av + 1 for metallicity (zinc) + 2 for background.

    run match/bin/stats on a match_cmd_file. Calculates the non-zero sfr bins
    in sfh_file.
    '''
    stats_exe = '$HOME/match2.5/bin/stats'
    sfr_data = read_binned_sfh(sfh_file)
    inds, = np.nonzero(sfr_data.sfr)

    nonzero_bins = len(inds)
    nfp = nonzero_bins + nfp_nonsfr
    cmd = '%s %s %s %i %i >> %s \n' % (extra, stats_exe, match_cmd_file,
                                       nmc_runs, nfp, outfile)

    if nmc_runs > 0:
        perr_frac = sfr_data.sfr_errp[inds] / sfr_data.sfr[inds]
        merr_frac = sfr_data.sfr_errm[inds] / sfr_data.sfr[inds]
        line = ('# min_sfr_merr max_sfr_perr med_sfr_merr med_sfr_perr ',
                'max_sfr_merr max_sfr_perr \n')
        line += '%.3f %.3f %.3f %.3f %.3f %.3f\n' % \
            (np.min(perr_frac), np.min(merr_frac), np.median(perr_frac),
             np.median(merr_frac), np.max(perr_frac), np.max(merr_frac))
        line += '# %s' % cmd
        with open(outfile, 'a') as out:
            out.write(line)

        print('wrote %s' % outfile)

    if not dry_run:
        import os
        os.system(cmd)
    return cmd


def read_match_stats(statsfile):
    with open(statsfile, 'r') as inp:
        lines = inp.readlines()
    stats = {}
    for line in lines:
        if ':' not in line:
            continue
        key, val = line.split(':')
        stats[''.join(key.replace('^', '').split())] = float(val)
    return stats


def main(argv):
    parser = argparse.ArgumentParser(description="run match stats")

    parser.add_argument('-n', '--nmc_runs', type=int, default=10000,
                        help='number of MC runs')

    parser.add_argument('-f', '--nfp_nonsfr', type=int, default=5,
                        help='non-sfr free parameters')

    parser.add_argument('-o', '--outfile', type=str, default='cmd_stats.dat',
                        help='output file')

    parser.add_argument('sfh_file', type=str, help='match sfh file')

    parser.add_argument('match_cmd_file', type=str,
                        help='match cmd grid')

    args = parser.parse_args(argv)

    match_stats(args.sfh_file, args.match_cmd_file, nfp_nonsfr=args.nfp_nonsfr,
                nmc_runs=args.nmc_runs, outfile=args.outfile)

if __name__ == '__main__':
    main(sys.argv[1:])
