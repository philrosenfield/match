"""Utilities for running match"""
from .fileio import readfile
import numpy as np

MAX_PROC = 8
MIN_PROC = 0

__all__ = ['write_script', 'insert_ext', 'check_boundaries', 'check_proc',
           'reset_proc', 'call_stats', 'call_sspcombine', 'call_calcsfh',
           'call_zcombine', 'diag_calls', 'call_hmc', 'mcmc_run']

def write_script(filename, cmd):
    """add diag_calls to cmd and write to filename"""
    cmd = diag_calls(cmd)
    with open(filename, 'w') as out:
        out.write(cmd)
    print('wrote %s' % filename)


def insert_ext(ext, string):
    """instert text before the last '.xxx' """
    strings = string.split('.')
    return '.'.join(np.insert(strings, -1, ext))


def check_boundaries(mparam, res_file='setz_results.dat'):
    def grab_line(mparam):
        return map(float, open(mparam).readline().strip().split())
    try:
        # powerlaw imf
        _, dmod0, dmod1, _, av0, av1, _ = grab_line(mparam)
    except ValueError:
        try:
            # kroupa imf
            dmod0, dmod1, _, av0, av1, _ = grab_line(mparam)
        except ValueError:
            # chabrier imf
            _, _, _, dmod0, dmod1, _, av0, av1, _ = grab_line(mparam)

    data = readfile(res_file)

    if len(np.nonzero(data['dmod'] == dmod1)[0]) > 0:
        print 'error need to increase dmod1 past %.2f' % dmod1
    if len(np.nonzero(data['dmod'] == dmod0)[0]) > 0:
        print 'error need to decrease dmod0 past %.2f' % dmod0
    if len(np.nonzero(data['Av'] == av1)[0]) > 0:
        print 'error need to increase av1 past %.2f' % av1
    if av0 > 0:
        if len(np.nonzero(data['Av'] == av0)[0]) > 0:
            print 'error need to decrease av0 past %.2f' % av0


def check_proc(nproc, cmd, max_proc=MAX_PROC):
    """if nproc >= max_proc will call reset_proc"""
    if nproc >= max_proc:
        nproc, cmd = reset_proc(nproc, cmd)
    return nproc, cmd


def reset_proc(nproc, cmd, min_proc=MIN_PROC):
    """add a wait signal to cmd and set nproc to 0"""
    #cmd += 'wait\n'
    cmd += "\nfor job in `jobs -p`\ndo\n    echo $job\n    wait $job\ndone\n\n"
    nproc = min_proc
    return nproc, cmd


def call_stats(outfile, cmd='', nproc=MIN_PROC, nfp_nonsfr=0, max_proc=MAX_PROC):
    """ add a line to call my call to match stats program """
    sfh_file = outfile.replace('out', 'sfh')
    cmd_file = outfile + '.cmd'
    stats_file = cmd_file + '.stats'

    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)
    extra = 'taskset -c %i' % nproc
    nproc += 1
    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)

    cmd += 'taskset -c %i python -c \"from ResolvedStellarPops.match.likelihood import match_stats; ' % nproc
    cmd += 'match_stats(\'%s\', \'%s\', nfp_nonsfr=%i, nmc_runs=0, outfile=\'%s\', dry_run=False, extra=\'%s\')\" & \n' % \
            (sfh_file, cmd_file, nfp_nonsfr, stats_file, extra)
    nproc += 1
    return cmd, nproc, stats_file


def call_sspcombine(outfile, nproc=0, cmd='', max_proc=MAX_PROC):
    """add a line calling sspcombine"""
    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)
    sspcombine = 'taskset -c %i $HOME/research/match2.5/bin/sspcombine' % nproc
    stats_file = outfile.replace('scrn', 'stats')

    cmd += '%s %s.dat > %s\n' % (sspcombine, outfile, stats_file)
    return cmd, nproc, stats_file


def call_calcsfh(phot, fake, mparam, flag0='-setz', flag='', nproc=0, cmd='',
                 max_proc=MAX_PROC):
    """add a line to a script to run calcsfh"""
    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)

    if flag == '':
        if 'ov' in mparam:
            flag = 'ov%s' % mparam.split('ov')[1][:4]
            flag = '-sub=%s' % flag
        else:
            flag = ''
    if 's12' in flag and not '_' in flag:
        flag = ''
    calcsfh = 'taskset -c %i $HOME/research/match2.5/bin/calcsfh' % nproc
    scrn = mparam.replace('matchparam', 'scrn')
    outfile = mparam.replace('matchparam', 'out')
    cmd += '%s %s %s %s %s -PARSEC %s %s > %s &\n' % (calcsfh, mparam, phot,
                                                      fake, outfile, flag0,
                                                      flag, scrn)
    return cmd, nproc, outfile


def call_zcombine(outfile, nproc=MIN_PROC, cmd='', flag='-bestonly',
                  max_proc=MAX_PROC):
    """add a line to a script to run zcombine"""
    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)

    zcombine = 'taskset -c %i $HOME/research/match2.5/bin/zcombine' % nproc
    if 'bestonly' in flag:
        sfh_file =  outfile.replace('out', 'sfh')
    else:
        sfh_file = outfile + '.zc'
    cmd += '%s %s %s > %s & \n' % (zcombine, outfile, flag, sfh_file)
    return cmd, nproc, sfh_file


def diag_calls(cmd):
    """call match plotting routine and grep the best fits to file"""
    cmd += 'python ~/research/python/ResolvedStellarPops/match/graphics.py F555W F814W\n'
    cmd += 'grep Best *.scrn | sed \'s/=/ /g\' | sort -g -r -k 8 > sorted_best_fits.dat\n'
    return cmd


def call_hmc(hmcinp, nproc=MIN_PROC, cmd='', flag='-tint=2.0 -nmc=10000 -dt=0.015',
             max_proc=MAX_PROC):
    """add a line to a script to run hybricMC"""
    nproc, cmd = check_proc(nproc, cmd, max_proc=max_proc)
    hmcout = hmcinp + '.mcmc'
    hmcscrn = hmcinp + '.scrn'
    hmc = 'taskset -c %i $HOME/research/match2.5/bin/hybridMC' % nproc
    cmd += '%s %s.dat %s %s > %s & \n' % (hmc, hmcinp, hmcout, flag, hmcscrn)
    return cmd, nproc, hmcout


def mcmc_run(phot, fake, match_params, cmd='', nproc=MIN_PROC, flags=None,
             flag0='-setz -mcdata', max_proc=MAX_PROC):
    outfiles = []
    hmc_files = []
    sfh_files = []

    # calcsfh
    for mparam in match_params:
        cmd, nproc, outfile = call_calcsfh(phot, fake, mparam, flag0=flag0,
                                           nproc=nproc, cmd=cmd, max_proc=max_proc)
        outfiles.append(outfile)
        nproc += 1
    nproc, cmd = reset_proc(nproc, cmd)

    # zcombine from calcsfh
    for outfile in outfiles:
        cmd, nproc, sfh_file = call_zcombine(outfile, nproc=nproc, cmd=cmd,
                                             max_proc=max_proc)
        sfh_files.append(sfh_file)
        nproc += 1
    nproc, cmd = reset_proc(nproc, cmd)

    # Hybrid MonteCarlo
    for outfile in outfiles:
        cmd, nproc, hmc_file = call_hmc(outfile, nproc=nproc, cmd=cmd,
                                        max_proc=max_proc)
        hmc_files.append(hmc_file)
        nproc += 1
    nproc, cmd = reset_proc(nproc, cmd)

    # zcombine on HMC
    flag = ' -unweighted -medbest -jeffreys -best=%s '
    for i in range(len(hmc_files)):
        cmd, nproc, zc_file = call_zcombine(hmc_files[i], flag=flag % sfh_files[i],
                                           cmd=cmd, nproc=nproc, max_proc=max_proc)
        nproc += 1
    nproc, cmd = reset_proc(nproc, cmd)

    write_script('mcdata_script.sh', cmd)
    return
