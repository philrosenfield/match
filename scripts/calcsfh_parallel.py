"""Run calcsfh or hybridMC in Parallel (using subprocess)"""
import argparse
import logging
import os
import subprocess
import sys

from glob import glob1
import numpy as np

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Could be in a config or environ
calcsfh = '$HOME/research/match2.5/bin/calcsfh'
zcombine = '$HOME/research/match2.5/bin/zcombine'
hybridmc = '$HOME/research/match2.5/bin/hybridMC'

def test_files(prefs, run_calcsfh=True):
    """make sure match input files exist"""
    return_code = 0
    for pref in prefs:
        if run_calcsfh:
            pfiles = calcsfh_existing_files(pref)
        else:
            pfiles = [hybridmc_existing_files(pref)]
        test = [os.path.isfile(f) for f in pfiles]
        if False in test:
            logger.error('missing a file in {}'.format(pref))
            logger.error(pfiles)
            return_code += 1
    if return_code > 0:
        sys.exit(2)
    return

def uniform_filenames(prefs, dry_run=False):
    """
    make all fake match and par files in a directory follow the format
    target_filter1_filter2.gst.suffix all lower case
    use dry_run to print the mv command, or will call os.system.
    """
    from glob import glob1
    for pref in prefs:
        dirname, p = os.path.split(pref)
        filters = '_'.join(p.split('_')[1:])
        print dirname, p, filters
        fake, = glob1(dirname, '*{}*fake'.format(filters))
        match, = glob1(dirname, '*{}*match'.format(filters))
        param, = glob1(dirname, '*{}*param'.format(filters))
        ufake = '_'.join(fake.split('_')[1:]).replace('_gst.fake1',
                                                      '.gst').lower()
        umatch = '_'.join(match.split('_')[1:]).lower()
        uparam = param.replace('.param', '.gst.param').lower()
        for old, new in zip([fake, match, param],[ufake, umatch, uparam]):
            cmd = 'mv {dir}/{old} {dir}/{new}'.format(dir=dirname, old=old,
                                                      new=new)
            logger.info(cmd)
            if not dry_run:
                os.system(cmd)

def calcsfh_existing_files(pref, optfilter1=''):
    """file formats for param match and matchfake"""
    param = pref + '.param'
    match = pref + '.match'
    fake = pref + '.matchfake'
    return (param, match, fake)


def calcsfh_new_files(pref):
    """file formats for match grid, sdout, and sfh file"""
    out =  pref + '.out'
    scrn = pref + '.scrn'
    sfh = pref + '.sfh'
    return (out, scrn, sfh)


def hybridmc_existing_files(pref):
    """file formats for the HMC, based off of calcsfh_new_files"""
    mcin = pref + '.out.dat'
    return mcin


def hybridmc_new_files(pref):
    """file formats for HybridMC output and the following zcombine output"""
    pref = pref.strip()
    mcmc = pref + '.mcmc'
    mcscrn = mcmc + '.scrn'
    mczc = mcmc + '.zc'
    return (mcmc, mcscrn, mczc)


def run_parallel(prefs, dry_run=False, nproc=8, run_calcsfh=True):
    """run calcsfh and zcombine in parallel, flags are hardcoded."""
    test_files(prefs, run_calcsfh)

    rdict = {'calcsfh': calcsfh, 'zcombine': zcombine,'hybridmc': hybridmc}
    # calcsfh
    # calcsfh, param, match, fake, out, scrn
    cmd1 = '{calcsfh} {param} {match} {fake} {out} -PARSEC -mcdata -kroupa -zinc -sub=v2 > {scrn}'
    # zcombine
    #zcombine, out, sfh
    cmd2 = '{zcombine} {out} -bestonly > {sfh}'
    # hybridmc
    #hybridmc, mcin, mcmc, mcscrn
    cmd3 = '{hybridmc} {mcin} {mcmc} -tint=2.0 -nmc=10000 -dt=0.015 > {mcscrn}'
    # zcombine w/ hybrid mc
    #zcombine, mcmc, mczc
    cmd4 = '{zcombine} {mcmc} -unweighted -medbest -jeffreys -best={mczc}'

    niters = np.ceil(len(prefs) / float(nproc))
    sets = np.arange(niters * nproc, dtype=int).reshape(niters, nproc)
    logging.debug('{} prefs, {} niters'.format(len(prefs), niters))

    for j, iset in enumerate(sets):
        # don't use not needed procs
        iset = iset[iset < len(prefs)]
        
        # run calcsfh
        procs = []
        for i in iset:
            if run_calcsfh:
                rdict['param'], rdict['match'], rdict['fake'] = calcsfh_existing_files(prefs[i])
                rdict['out'], rdict['scrn'], rdict['sfh'] = calcsfh_new_files(prefs[i])
                cmd = cmd1.format(**rdict)
            else:
                rdict['mcin'] = hybridmc_existing_files(prefs[i])
                rdict['mcmc'], rdict['mcscrn'], rdict['mczc'] = hybridmc_new_files(prefs[i])
                cmd = cmd3.format(**rdict)
            if not dry_run:
                procs.append(subprocess.Popen(cmd, shell=True))
            logger.info(cmd)
        
        # wait for calcsfh
        if not dry_run:
            [p.wait() for p in procs]
            logger.debug('calcsfh or hybridMC set {} complete'.format(j))
        
        # run zcombine
        procs = []
        for i in iset:
            if run_calcsfh:
                rdict['out'], rdict['scrn'], rdict['sfh'] = calcsfh_new_files(prefs[i])
                zcom = cmd2.format(**rdict)
            else:
                zcom = cmd4.format(**rdict)
            if not dry_run:
                procs.append(subprocess.Popen(zcom, shell=True))
            logger.info(zcom)
        
        # wait for zcombine
        if not dry_run:
            [p.wait() for p in procs]
            logger.debug('zcombine set {} complete'.format(j))


def main(argv):
    """parse in put args, setup logger, and call run_parallel"""
    parser = argparse.ArgumentParser(description="Run calcsfh in parallel. Note: bg cmd, if in use, need to be in the current folder")

    parser.add_argument('-d', '--dry_run', action='store_true',
                        help='only print commands')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='set logging to debug')

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('-m', '--hmc',  action='store_false',
                        help='run hybridMC (must be after a calcsfh run)')

    parser.add_argument('-f', '--logfile', type=str, default='calcsfh_parallel.log',
                        help='log file name')

    parser.add_argument('-s', '--simplify', action='store_true',
                        help='make filename uniform and exit (before calcsfh run)')

    parser.add_argument('pref_list', type=argparse.FileType('r'),
                        help="list of prefixs to run on. E.g., ls */*.match | sed 's/.match//' > pref_list")

    args = parser.parse_args(argv)
    prefs = [l.strip() for l in args.pref_list.readlines()]
    
    handler = logging.FileHandler(args.logfile)
    if args.verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if args.simplify:
        uniform_filenames(prefs, dry_run=args.dry_run)
    else:
        logger.info('running on {}'.format(', '.join([p.strip() for p in prefs])))
        run_parallel(prefs, dry_run=args.dry_run, nproc=args.nproc, run_calcsfh=args.hmc)


if __name__ == '__main__':
    main(sys.argv[1:])
