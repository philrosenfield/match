"""
A hacky way to run match in parallel. Search for the best Z, IMF, COV etc for
a population or ssp.

The first reason to do this was to get SFH and HMC uncertainties for the best
value of one fixed Z. Then it was expanded to do the same but for each COV.
Then, IMF...

These functions don't run anything. They make bash script files to be run on
their own. I think it's better to separate creating files to do the runs and
doing the runs themselves.
"""
from __future__ import absolute_import
import argparse
import glob
import sys

import numpy as np

from .fileio import readfile, savetxt
from .likelihood import read_match_stats
from .utils import grab_val
from . import match_scripts


def vary_mpars(template_match_param, gmin=-1., gmax=-.4, dg=0.05, zdisp=0.05,
               flag='', vary='setz', gs=None, imf_flag=''):
    """
    take a match parameter file set up with the -setz option and make many
    at a given z and zdisp.
    """
    with open(template_match_param, 'r') as infile:
        lines = [l.strip() for l in infile.readlines()]

    if vary == 'setz':
        if gs is None:
            gs = np.arange(gmin, gmax + dg / 2, dg)
        mdisp = zdisp

    if vary == 'imf':
        if gs is None:
            gs = np.arange(gmin, gmax, dg)
        if imf_flag != '':
            if imf_flag == '-kroupa':
                # it is a string because of sext definition below
                gs = ['']
            elif imf_flag == '-chabrier':
               gs = ['1.3 0.08 0.69']
            else:
                print 'imf_flag not recognized'
                return []
    fnames = []
    for i in range(len(gs)):
        if vary == 'setz':
            lines[1] = '%.3f' % mdisp
            for j in np.arange(9, len(lines)-1):
                to, tf, mh = lines[j].strip().split()
                lines[j] = '     %s %s %.3f' % (to, tf, gs[i])
                sext = 'lz%g' % (gs[i] + 4.)

        if vary == 'imf':
            data = lines[0].split()
            data[0] = str(gs[i])
            lines[0] = ' '.join(data)
            try:
                sext = 'imf%g' % np.abs(gs[i])
            except:
                sext = 'imf%s' % imf_flag.replace('-', '_')

        outfile = match_scripts.insert_ext(sext, template_match_param)
        if flag != '':
            sflag = flag.split('=')[1]
            if not sflag in outfile:
                outfile = match_scripts.insert_ext(sflag, outfile)
        with open(outfile, 'w') as out:
            out.write('\n'.join(lines))

        fnames.append(outfile)
    return fnames


def run_grid(phot, fake, mparams, gmin=-1., gmax=-.4, dg=0.05, zdisp=0.05,
             flags=[''], flag0='-setz', cmd='', nproc=0, vary='setz',
             gs=None, max_proc=8):
    """
    create a bash script to run calsfh, zcombine, and make plots in parallel
    """
    imf_flag = ''
    if 'kroupa' in flag0:
        imf_flag = '-kroupa'
    elif 'chabrier' in flag0:
        imf_flag = '-chabrier'

    match_params = []
    outfiles = []

    for flag in flags:
        for template_match_param in mparams:
            fnames = vary_mpars(template_match_param, gmin, gmax, dg,
                                zdisp=zdisp, flag=flag, vary=vary, gs=gs,
                                imf_flag=imf_flag)
            for fname in fnames:
                cmd, nproc, outfile = \
                    match_scripts.call_calcsfh(phot, fake, fname, flag0=flag0,
                                               flag=flag, nproc=nproc, cmd=cmd,
                                               max_proc=max_proc)
                outfiles.append(outfile)
                match_params.append(fname)
                nproc += 1

    nproc, cmd = match_scripts.reset_proc(nproc, cmd)

    if '-ssp' in flag0:
        scrns = []
        # creates a outfile.dat with no header
        for mparam in match_params:
            scrn = mparam.replace('matchparam', 'scrn')
            nproc, cmd = match_scripts.check_proc(nproc, cmd, max_proc=max_proc)
            cmd += 'python -c \"from match.utils import strip_header; strip_header(\'%s\')\" & \n' % scrn
            nproc += 1
            scrns.append(scrn)
        nproc, cmd = match_scripts.reset_proc(nproc, cmd)

    for i, outfile in enumerate(outfiles):
        if '-ssp' in flag0:
            cmd, nproc, stats_file = \
                match_scripts.call_sspcombine(scrns[i], nproc=nproc, cmd=cmd,
                                              max_proc=max_proc)
        else:
            cmd, nproc, sfh_file = \
                match_scripts.call_zcombine(outfile, nproc=nproc, cmd=cmd,
                                            max_proc=max_proc)
        nproc += 1
        nfp_nonsfr = 2  # dmod, Av
        cmd, nproc, stats_file = \
            match_scripts.call_stats(outfile, cmd=cmd, nproc=nproc,
                                     nfp_nonsfr=nfp_nonsfr, max_proc=max_proc)

    proc, cmd = match_scripts.reset_proc(nproc, cmd)
    if vary not in flag0:
        flag0 = '%s %s' % (vary, flag0)

    match_scripts.write_script('%s_script.sh' % '_'.join(flag0.replace('-', '').split()), cmd)
    return match_params


def setz_results(fname='setz_results.dat'):
    import matplotlib.pylab as plt
    from palettable import sequential
    from ResolvedStellarPops.match.match_grid import MatchGrid

    def cov_plot():
        unc = np.unique(data['COV'])
        colors = sequential.Blues[len(unc)].mpl_colors
        #colors = ['black', 'orange', 'green', 'purple', 'darkred']
        fig, ax = plt.subplots()
        for i, c in enumerate(unc):
            inds, = np.nonzero(data['COV'] == c)
            isort = np.argsort(data['Z'][inds])
            xarr = data['Z'][inds][isort]
            yarr = data['fit'][inds][isort]
            ax.plot(xarr, yarr, color='k', lw=3)
            ax.plot(xarr, yarr, color=colors[i], lw=2,
                    label='$\Lambda_c=%.2f$' % c)
            ax.plot(xarr, yarr, 'o', color=colors[i], mec='white')

        for i in range(len(data)):
            ax.annotate('$%.2f, %.1f$' % (data['dmod'][i], data['Av'][i]),
                        (data['Z'][i], data['fit'][i]), fontsize=7)
        ax.legend(loc=0, frameon=False)
        ax.set_ylabel('Fit Parameter', fontsize=20)
        ax.set_xlabel('$Z$', fontsize=20)

        plt.savefig(fname.replace('.dat', '.png'))
        print 'wrote %s' % fname.replace('.dat', '.png')
        plt.close()

    def pdf_plots():
        mg = MatchGrid(fname, ssp=False)
        ycols = ['COV', 'COV', 'dmod']
        xcols = ['Z', 'IMF', 'Av']
        pops = []
        #import pdb; pdb.set_trace()
        # don't plot if nothing varies
        for i in range(len(ycols)):
            y = mg.data[ycols[i]]
            x = mg.data[xcols[i]]
            nuiqs = np.min([len(np.unique(x)), len(np.unique(y))])
            if nuiqs <= 1:
                print '%s or %s unchanging' % (xcols[i], ycols[i])
                pops.append(i)
            elif np.sum(np.isnan(x)) == len(x) or np.sum(np.isnan(y)) == len(y):
                print '%s or %s nans' % (xcols[i], ycols[i])
                pops.append(i)

        if len(pops) > 0:
            for p in pops:
                xcols.pop(p)
                ycols.pop(p)

        zcols = ['fit', 'fit', 'fit']
        zcols = ['chi2', 'chi2', 'chi2']

        stat = np.min
        for i in range(len(xcols)):
            ax, cb = mg.pdf_plot(xcols[i], ycols[i], zcols[i], stat=stat)
            ffmt = '_'.join([xcols[i], ycols[i], zcols[i]])
            figname = fname.replace('.dat', '_%s.png' % ffmt)
            plt.savefig(figname)
            print 'wrote %s' % figname
            plt.close()

    data = readfile(fname)
    if len(np.unique(data['IMF'])) <= 1:
        # not so interesting if varying IMF....
        cov_plot()
    pdf_plots()
    return


def find_best(infile='sorted_best_fits.dat', outfile='setz_results.dat', ssp=False):
    """

    """
    def get_imf(filename):
        """
        imf is the second line, second value of the screen output from match
        """
        with open(filename, 'r') as f:
            f.readline()
            line = f.readline().replace(',', '')
        return float(line.split()[1])

    outfile = infile.replace('sorted', '').replace('.dat', '_results.dat')
    with open(infile, 'r') as inp:
        lines = inp.readlines()
    fnames, data = zip(*[l.strip().split(':Best fit:') for l in lines])
    # values from filename
    z = np.array([grab_val(f, 'lz') for f in fnames]) - 4.

    try:
        ov = np.array([grab_val(f, 'ov', v2='s', v3='v') for f in fnames])
    except:
        ov = np.zeros(len(z)) + 5

    # values after the filename
    _, av, _, dmod, _, fit = zip(*[d.replace(',', '').split() for d in data])
    fit = np.array(fit, dtype=float)
    av = np.array(av, dtype=float)
    dmod = np.array(dmod, dtype=float)

    # values from within the file
    imf = np.array([get_imf(f) for f in fnames])

    if ssp:
        savetxt(outfile,
                np.column_stack([z, ov, av, dmod, imf, fit]),
                fmt='%.4f %.3f %.2f %.2f %.2f %6f',
                header='# Z COV Av dmod IMF fit \n',
                overwrite=True)
    else:
        # best_fit, expect, sigma, variance, fitz, effchi2
        best_fit, expect, sigma, variance, fitz, effchi2 = load_match_stats(best_fit_file=infile)

        savetxt(outfile,
                np.column_stack([z, ov, av, dmod, imf, fit, best_fit, expect,
                                 sigma, variance, fitz, effchi2]),
                fmt='%.4f %.3f %.2f %.2f %.2f %6f %6f %.3f %.2f %.2f %.3f %.3f',
                header='# Z COV Av dmod IMF fit best_fit expect sig variance fit_z chi2\n',
                overwrite=True)

    inds = np.digitize(ov, bins=np.unique(ov))
    mparams = []
    for iov in range(len(np.unique(ov))):
        best = np.min(fit[np.nonzero(inds==iov+1)])
        ind, = np.nonzero(fit == best)
        #print lines[ind]
        mparam = fnames[ind].replace('scrn', 'matchparam')
        mparams.append(mparam)
    return mparams, outfile


def load_match_stats(best_fit_file='sorted_best_fits.dat'):
    with open(best_fit_file, 'r') as inp:
        lines = inp.readlines()
    fnames, data = zip(*[l.strip().split('.scrn:Best fit:') for l in lines])
    best_fit = np.array([])
    expect = np.array([])
    sigma = np.array([])
    variance = np.array([])
    fitz = np.array([])
    effchi2 = np.array([])
    for fname in fnames:
        fname += '.out.cmd.stats'
        sdict = read_match_stats(fname)
        best_fit = np.append(best_fit, sdict['Bestfit'])
        expect = np.append(expect, sdict['expectation'])
        sigma = np.append(sigma, sdict['sigma'])
        variance = np.append(variance, sdict['variance'])
        fitz = np.append(fitz, sdict['fitz'])
        effchi2 = np.append(effchi2, sdict['effchi2'])
    return best_fit, expect, sigma, variance, fitz, effchi2


def main(argv):
    """
    usage e.g,:
    setzsearch.py template_match_param setz imf kroupa
    setzsearch.py template_match_param setz imf chabrier
    setzsearch.py template_match_param ssp imf
    then:
    setzsearch.py zres
    then:
    setzsearch.py mcmc
    """
    #flags = ['-sub=s12_hb2']
    parser = argparse.ArgumentParser(description="Create bash script to run MATCH many times")

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('-z', '--setz', action='store_true',
                        help='use setz flag in MATCH calls')

    parser.add_argument('-s', '--ssp', action='store_true',
                        help='use ssp flag in MATCH calls')

    parser.add_argument('-i', '--imf', action='store_true',
                        help='set or vary the imf')

    parser.add_argument('-k', '--kroupa', action='store_true',
                        help='set imf to kroupa')

    parser.add_argument('-c', '--chabrier', action='store_true',
                        help='set imf to chabrier')

    parser.add_argument('-m', '--mcmc', action='store_true',
                        help='do MCMC run')

    parser.add_argument('-b', '--bestfits', type=str, default=None,
                        help='Best fit file to use with -r flag')

    parser.add_argument('-r', '--res', type=str, default=None,
                        help='Cull results')

    parser.add_argument('-x', '--imfrange', type=str, default='0.0,0.14,0.02',
                        help='comma separated imfmin, imfmax, dimf')

    parser.add_argument('-g', '--setz_at', type=str, default='-0.70',
                        help='set metallicity to use with -z flag')

    parser.add_argument('-o', '--ovrange', type=str, default='0.3,0.75,0.05',
                        help='core overshoot range. Must to correspond with MATCH sub directories')

    parser.add_argument('matchparam', type=str,
                        help='match input file template')

    args = parser.parse_args(argv)

    ovmin, ovmax, dov = map(float, args.ovrange.split(','))
    flags = ['-sub=ov%.2f' % o for o in np.arange(ovmin, ovmax, dov)]
    flag0 = '-full'

    phot, = glob.glob1('.', '*match')
    fake, = glob.glob1('.', '*matchfake')
    template_match_param = args.matchparam

    if args.setz:
        flag0 += ' -setz'
        gs = [map(float(args.setzat))]
        match_params = run_grid(phot, fake, [template_match_param], gs=gs,
                                vary='setz', flags=flags, flag0=flag0,
                                max_proc=args.nproc)
    elif args.res:
        # write the best values
        best_mpars, res_file = find_best(infile=args.bestfits)
        # check the solution for boundary issues
        match_scripts.check_boundaries(template_match_param, res_file=res_file)
        #if 'z' in func:
        # make a plot
        setz_results(fname=res_file)
    elif args.mcmc:
        # run mcmc on the best mparams
        match_params = [template_match_param]
        match_scripts.mcmc_run(phot, fake, match_params)

    if args.imf:
        # make a grid of varying IMF values based on either the best
        # mparams or a given mparam
        mparams = [template_match_param]
        if args.kroupa:
            flag0 += ' -kroupa'
        elif args.chabrier:
            flag0 += ' -chabrier'
        if args.setz:
            # use the grid already made, don't need to loop over the flags
            # again.
            mparams = match_params
            flags = ['']
        elif args.ssp:
            flag0 += ' -ssp'
        imf_min, imf_max, dimf = map(float, args.imfrange.split(','))

        mparams = run_grid(phot, fake, mparams, gmin=imf_min, gmax=imf_max,
                           dg=dimf, vary='imf', flag0=flag0, flags=flags,
                           max_proc=args.nproc)

if __name__ == '__main__':
    main(sys.argv[1:])
