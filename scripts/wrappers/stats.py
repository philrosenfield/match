from ..config import match_base
import os

def call_stats(cmdfiles, outdir=None, nfp=3, dryrun=False):
    """ call match/bin/stats on a list of .cmd files"""
    if type(cmdfiles) is not list:
        cmdfiles = [cmdfiles]

    stats = os.path.join(match_base, 'bin', 'stats')
    assert os.path.isfile(stats), 'stats program not found. {}'.format(stats)

    for cmdfile in cmdfiles:
        outfile = cmdfile + '.stats'
        if outdir is not None:
            assert os.path.isdir(outdir), \
                '{} directory not found'.format(outdir)
            outfile = os.path.join(outdir, outfile)
        cmd = '{:s} {:s} 0 {:d} > {:s}'.format(stats, cmdfile, nfp, outfile)
        print(cmd)
        if dryrun:
            print(cmd)
        else:
            os.system(cmd)
    return
