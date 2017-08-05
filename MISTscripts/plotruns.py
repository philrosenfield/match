#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from match.scripts import ssp
from match.scripts import cmd
import numpy as np
import glob
import sys
import os
import shutil
import matplotlib.image as mpimg
import seaborn as sns
# no seaborn gridlines on plot, and white bg:
sns.set_context('paper')
sns.set(font='serif')
sns.set_style("white", {'font.family' : 'serif', 'font.serif': ['Times', 'Palatino', 'serif']})
plt.axes(frameon=False)
from MIST_codes.scripts import read_mist_models as rmm


# Add something to perform diagnostic checks...
# + show that the results are accurate via e.g. plot errors and values as function of star number.
# + show sensitivity to completeness?
# + show effects of alteration of binsize
# + show effects of altering binary fraction? (maybe useful)

# General idea: show that match is accurate (or not) in smallish survey size regime that these clusters possess.


def marginalize(clustername, photobase, jobidstr, outname, params = ['lage', 'vvcrit', 'logZ', 'dmod', 'Av'], local=False):

    """
        This function calls certain methods from Philip Rosenfield/Daniel Weisz's MATCH scripts. Particularly,

        + ssp.combine_files(); to compile various v/vcrit MATCH runs into a single run where all best fits will be judged.
        + ssp.SSP(); to create a useable SSP object from the MATCH .scrn files.
            > methods therein of the SSP object are also called.

        + cmd.CMD(); to create a CMD object from the MATCH .cmd files.
            > mehthods therein of the CMD object are called as well.

        Arguments:
        ----------
        + clustername (str): Name of the cluster whose data will be utilized.
        
    """

    global home_dir
    global match_dir
    global pathtofiles
    global scrn_dir
    global outpath
    global plotpath

    # All paths defined below are system dependent:
    home_dir = os.environ['HOME_DIR']

    # dir where e.g. Hyades/output is located:
    # There are two directories here because on my local machine, I keep supercomp
    # output in a separate dir that I DL to, and local runs w/ MATCH are stored in
    # their own location.
    if not local:
        # downloaded runs (wouldn't be used on remote machine):
        match_dir = os.path.join(os.environ['MATCH_ODYOUT'], clustername)
    else:
        # local runs (should be used on remote machine, or your local machine):
        match_dir = os.path.join(os.environ['MATCH_DIR'], 'MIST', clustername)

    # general path to output files:
    pathtofiles = os.path.join(photobase, 'SLURM{:s}'.format(jobidstr))
    # path to MATCH .scrn files:
    scrn_dir = os.path.join(match_dir, 'scrns/', pathtofiles)
    # path to MATCH .cmd and .out files:
    outpath = os.path.join(match_dir, 'output/', pathtofiles)
    # path to place output plots in:
    plotpath = os.path.join(outpath, 'plots/')

    # list of .scrn files:
    flist = glob.glob(os.path.join(scrn_dir, '*'))
    print(scrn_dir) # debug
    # print .scrn files found:
    for afile in flist:
        print("Reading {:s}...".format(afile))

    # calling MATCH scripts' combine_files() to combine various v/vcrit runs:
    ssp.combine_files(flist, outfile=outname)

    # Create a match SSP object from the combined scrns of the given runs:
    combineddata = ssp.SSP(outname)

    # Corner plots:
    # 'lage', 'vvcrit', 'logZ', 'dmod', 'Av'
    fig, pdfax, quant = combineddata.pdf_plots(['lage', 'vvcrit', 'logZ'], twod=False, quantile=True, cmap=plt.cm.Reds)
    plt.close()
    fig, pdfax, quant = combineddata.pdf_plots(['lage', 'vvcrit', 'logZ'], twod=True, quantile=True, cmap=plt.cm.Reds)
    plt.savefig(os.path.join(plotpath, 'pdfcornerplot.png'))
    plt.close()

    # quantile dictionary (the quantiles will be used to draw bounding isochrones for the best-fits at ~+/- 1 sigma parameter values):
    qdict = {'lage': quant[0], 'vvcrit': quant[1], 'logZ':quant[2]}
    print(qdict) # debug

    # remove the combined .csv file since it's no longer needed:
    os.remove(os.path.join(os.getcwd(), outname))

    # A plot of the best fit ages found vs. the corresp. v/vcrit values...
    # This is also getting the name of the saved figure for the plot, the matplotlib figure for the plot,
    # the axis object too, the best parameters for each v/vcrit run and the corresp. v/vcrits:
    vvcritagename, vvcritagefig, vvcritageax, bestdict, vvcrits = agevsvvcrit(flist, jobidstr, outpath = plotpath) # vbests, vvcrits
    # Move the best fit ages vs v/vcrit plot:
    shutil.move(os.path.join(os.getcwd(), vvcritagename), os.path.join(outpath, vvcritagename))

    # Write a .txt file containing the best fit parameters found, along with their assoc. logP:
    combo_bestdict = {param: getbest(param, combineddata) for param in params}
    with open(os.path.join(plotpath, 'marginalized_params.txt'), 'w+') as f:
        for i, param in enumerate(params):
            #if param == 'lage':
            #    f.write('Best age: {:.2e}, logP = {:.4f}\n'.format(10**bestlist[i][0], bestlist[i][1]))
            #else:
            #    f.write('Best {:s}: {:.2f}, logP = {:.4f}\n'.format(param, bestlist[i][0], bestlist[i][1]))
            if param == 'Av':
                # Converting Av to E(B-V):
                f.write('Best E(B-V): {:.2f}, logP = {:.4f}\n'.format(combo_bestdict[param][0]/3.1, combo_bestdict[param][1]))
            else:
                f.write('Best {:s}: {:.2f}, logP = {:.4f}\n'.format(param, combo_bestdict[param][0], combo_bestdict[param][1]))


    #return vbests, vvcrits, qdict, bestdict#, vvcritagefig
    return vvcrits, bestdict, qdict, combo_bestdict

def agevsvvcrit(flist, jobarr_id, outpath, ext='.png'):

    """
        Creates best-fit age vs. v/vcrit plot.

        Arguments:
        ----------
             + flist (list, str): A list of file paths to various MATCH .scrn files to be used.
             + jobarr_id (str): String for the SLURM job array ID (maybe you don't use SLURM); used to distringuish
                                runs that otherwise have the same name & is the numerical sequence assigned to that run's
                                job array ID. Specific to how I do things & perhaps should be generalized...
             + outpath (str): Filepath to which the output plots will be saved,
             + ext (str): String for the saved plot's file extension, e.g. '.png' or '.eps', etc.

        Returns:
        --------
            + savename (str): The save name of the saved best age vs. v/vcrit plot.
            + fig (matplotlib Figure): Figure object of the saved best age vs. v/vcrit plot.
            + ax (matplotlib Axes): Axes object of      '                                 '.
            + bestlst (list, list, float): A list of lists of the best-fit values for log10 age, [Fe/H], distance modulus, and Av.
                                           Elements of each sublist correspond to a particular v/vcrit value.
            + vvcrits (list, float): A list of the v/vcrit values that each element in the aforementioned sublists corresponds to.
    """

    outnames = flist#[fpath.split('/')[-1].split('.scrn')[0] + '.csv' for fpath in flist]
    vvcrits=[]
    data_arr=[]
    bestages=[]
    bestfehs=[]
    bestdmods=[]
    bestavs=[]
    best = {}
    
    # For each v/vcrit value, get the best fit age found:
    for i, name in enumerate(outnames):

        #ssp.combine_files([flist[i]], outfile=name)
        data_arr.append(ssp.SSP(name))
        # Get the v/vcrit value:
        vvcrits.append(float(name.split('vvcrit')[-1].split('_')[0]))
        # For this v/vcrit value, get the best-fit age, logZ, and dist. modulus:
        val, lagelnP = getbest('lage', data_arr[i])
        bestages.append(val)
        val, fehlnP = getbest('logZ', data_arr[i])
        bestfehs.append(val)
        val, dmodlnP = getbest('dmod', data_arr[i])
        bestdmods.append(val)
        val, AvlnP = getbest('Av', data_arr[i])
        bestavs.append(val)

        # writes a .txt file storing the best parameters for a particular v/vcrit run:
        with open(os.path.join(match_dir, outpath, 'bestparams_vvcrit{:.1f}.txt'.format(vvcrits[i])), 'w+') as f:
            #for j, param in enumerate(params):
                #if param == 'lage':
            f.write('Best age: {:.2e}, logP = {:.4f}\n'.format(10**bestages[i], lagelnP))
                #else:
            f.write('Best [Fe/H]: {:.2f}, logP = {:.4f}\n'.format(bestfehs[i], fehlnP))
            f.write('Best dmod (fixed): {:.2f}, logP = {:.4f}\n'.format(bestdmods[i], dmodlnP))
            f.write('Best Av (fixed): {:.2f}, logP = {:.4f}\n'.format(bestavs[i], AvlnP))

        # Move the combined file:
        shutil.move(os.path.join(os.getcwd(), name), os.path.join(plotpath, name))

    # sort v/vcrits in ascending order and their best params too:
    # log10 ages found:
    bestages = [lage_val for (vvcrit_val, lage_val) in sorted(zip(vvcrits, bestages))]
    # logZ values found:
    bestfehs = [logz_val for (vvcrit_val, logz_val) in sorted(zip(vvcrits, bestfehs))]
    # dmod (fixed) values found:
    bestdmods = [dmod_val for (vvcrit_val, dmod_val) in sorted(zip(vvcrits, bestdmods))]
    # E(B-V)s (fixed) found; note: converting from Av to E(B-V):
    bestebvs = [av_val/3.1 for (vvcrit_val, av_val) in sorted(zip(vvcrits, bestavs))]
    # Sort the v/vcrits in ascending order:
    vvcrits = sorted(vvcrits)

    best_dict = {}
    for i, vvcrit in enumerate(vvcrits):
        best_dict[vvcrit] = {'lage': bestages[i], 'feh': bestfehs[i], 'dmod': bestdmods[i], 'ebv': bestebvs[i]}

    # Make a scatter plot of the best fit ages vs. v/vcrit:

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(vvcrits, bestages)
    ax.set_ylabel('Best-fit Age')
    ax.set_xlabel('v/v_crit')

    savename = jobarr_id + '_BestAgevsVVCrit' + ext
    plt.savefig(savename)

    #bestlst = [bestages, bestfehs, bestdmods, bestavs]

    return savename, fig, ax, best_dict, vvcrits#bestlst, vvcrits

def getbest(param, data):

    """
        This method gets the best-fit (max likelihood) parameter value from an ssp object (SSP class defined in Phil R. and Dan W.'s -- RW's --  match scripts).

        Argument:
        ---------
            + param (str): String for the parameter of interest (would be e.g.: lage, logZ, dmod, see marginalize method of SSP class in RW's scripts)
            + data (SSP): RW script's SSP class instance.

       Returns:
       --------
            + best (float): best-fit value of the parameter corresponding to the string designated by the 'param' argument.
            + lnP (float): Associated log probability of the best-fit parameter.

    """
    # Index where the marginalized dist. for the given param is at a maximum probability.
    bestindex = np.where(data.marginalize(param)[1] == np.amax(data.marginalize(param)[1]))[0][0]
    # 0 indexes the parameters values (1 would index the corresponding lnP)
    # This is the max lnP parameter value.
    best, lnP = data.marginalize(param)
    best = best[bestindex]
    lnP = lnP[bestindex]

    return float(best), float(lnP)

if __name__ == "__main__":

    """
        This is a script runnable as e.g.

        >> ./plotruns.py Hyades hyades_debruTIGS_TychoBV.phot 12345678 somename.csv Tycho_B Tycho_V 
     
        The intent is to plot the best-fit models found via a MATCH run; this makes use of Phil R./Dan W.'s (RW) match scripts and also
    Jieun Choi's MIST codes (available on her Github site). The former is mainly used to conglomerate MATCH output across separate
    runs made at various v/vcrit values (as this parameter is not explored in the usual way by MATCH, RW's scripts concatenate the output
    of the runs and marginalize to acheive a best-fit over any new parameters not used within MATCH itself, e.g. v/vcrit). The latter is
    used to plot the best-fit MIST models found by the run(s).

    """

# NOTE:
# sys.argv[x] = cluster name, photometry file, job array slurm id, outputfilename, bluer band filter, redder band filter.
    photbase = os.path.splitext(sys.argv[2])[0]

    try:
        local = sys.argv[7]
    except IndexError:
        local = False

    #vbests, vvcrits, qdict, bfdict = marginalize(sys.argv[1], photbase, sys.argv[3], sys.argv[4], local=local)
    vvcrits, best_dict, qdict, bfdict = marginalize(sys.argv[1], photbase, sys.argv[3], sys.argv[4], local=local)
    #pdf.savefig(vvcage_fig)
    plt.clf()
    # clearing out pre-existing output files that this script may have produced in prev. runs:
    for f in glob.glob(os.path.join(outpath,'*.png')):
        os.remove(f)
    for f in glob.glob(os.path.join(plotpath,'*.png')):
        os.remove(f)
    for f in glob.glob(os.path.join(plotpath,'*.txt')):
        os.remove(f)
    #Av = float(sys.argv[7])

    v_filter = sys.argv[5]
    i_filter = sys.argv[6]
    filters = [v_filter, i_filter]
    for n, band in enumerate(filters):
        if band == 'B':
            filters[n] = 'Bessell_B'#1
        elif band == 'V':
            filters[n] = 'Bessell_V'#2
        elif band == 'J':
            filters[n] = '2MASS_J'#5
        elif band == 'H':
            filters[n] = '2MASS_H'#6
        elif band == 'Ks':
            filters[n] = '2MASS_Ks'#7
        elif band == 'Tycho_B':
            filters[n] = 'Tycho_B'#11
        elif band == 'Tycho_V':
            filters[n] = 'Tycho_V'#12


    # Retrieve the magnitude limits used in the fit (from the calcsfh .param file used):
    param_fname = glob.glob(os.path.join(plotpath, 'calcsfh_{:s}.param'.format(photbase)))[0]
    with open(param_fname, 'r') as f:
        lines = f.readlines()
        # Color limits:
        color_lims = lines[4].split(' ')
        minvmi = float(color_lims[3])
        maxvmi = float(color_lims[4])
        ibin = float(color_lims[0])
        vibin = float(color_lims[1])
        ilims = lines[6].split(' ')
        imin = float(ilims[1])
        imax = float(ilims[0])
        ntbins= int(lines[8])
        tbins=lines[9:8+ntbins]
        for i, aline in enumerate(tbins):
            tbins[i] = map(float, aline.split('\n')[0].split(' ')[-2:])
        print(tbins)

    # sort v/vcrits in ascending order and their best params too:
    # log10 ages found:
    #vbests[0] = [lage_val for (vvcrit_val, lage_val) in sorted(zip(vvcrits, vbests[0]))]
    # logZ values found:
    #vbests[1] = [logz_val for (vvcrit_val, logz_val) in sorted(zip(vvcrits, vbests[1]))]
    # dmod (fixed) values found:
    #vbests[2] = [dmod_val for (vvcrit_val, dmod_val) in sorted(zip(vvcrits, vbests[2]))]
    # E(B-V)s (fixed) found:
    #vbests[3] = [av_val/3.1 for (vvcrit_val, av_val) in sorted(zip(vvcrits, vbests[3]))]
    # Sort the v/vcrits in ascending order:
    #vvcrits = sorted(vvcrits)

    #best_dict = {}
    #for i, vvcrit in enumerate(vvcrits):
    #    best_dict[vvcrit] = {'lage': vbests[0][i], 'feh': vbests[1][i], 'dmod': vbests[2][i], 'ebv': vbests[3][i]}

    # Also get the residuals:
    # For overplotting data points, get the data points CMD x, y values from the photometry file:
    photf_v = np.genfromtxt(os.path.join(match_dir, 'photometry', '{:s}.phot'.format(photbase)), usecols=(0,))
    photf_i = np.genfromtxt(os.path.join(match_dir, 'photometry', '{:s}.phot'.format(photbase)), usecols=(1,))
    photf_vmi = photf_v - photf_i
    # Tuple storing the x, y values for plotting the data on a CMD:
    photf_pts = (photf_vmi, photf_v)

    plt.clf()
    bestfit = 99999.9
    # Plot best isochrone for each v/vcrit individually:
    #=================================================================
    allfig = plt.figure()
    allax = allfig.add_subplot(111)
    for j, vvcrit in enumerate(vvcrits):

        fig = plt.figure()
        ax = fig.add_subplot(111)

        fname = glob.glob(os.path.join(outpath,'*vvcrit{:.1f}*.cmd'.format(vvcrit)))[0]
        a_cmd = cmd.CMD(fname)
        extent = a_cmd.extent
        
        # best isochrone:
        iso = rmm.ISOCMD(best_dict[vvcrit]['feh'], vvcrit, ebv=best_dict[vvcrit]['ebv'])
        iso.set_isodata(best_dict[vvcrit]['lage'],filters[1], filters[0], dmod=best_dict[vvcrit]['dmod'])

        # older age isochrone:
        isou = rmm.ISOCMD(best_dict[vvcrit]['feh'], vvcrit, ebv=best_dict[vvcrit]['ebv']) #qdict['logZ'][1], float('{:.1f}'.format(qdict['vvcrit'][1]))
        isou.set_isodata(qdict['lage'][1],filters[1], filters[0], dmod=best_dict[vvcrit]['dmod'])
        # younger age:
        isol = rmm.ISOCMD(best_dict[vvcrit]['feh'], vvcrit, ebv=best_dict[vvcrit]['ebv']) #qdict['logZ'][0], float('{:.1f}'.format(qdict['vvcrit'][0]))
        isol.set_isodata(qdict['lage'][0],filters[1], filters[0], dmod=best_dict[vvcrit]['dmod'])

        iso.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, label=True)
        isou.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')
        isol.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')
        #ax.fill_between(iso.x, isol.y, isou.y)

        # all best isochrones
        iso.isoplot(allax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, legloc='upper right')

        ax.scatter(*photf_pts, lw=0, s=8, c='r')
        #if j == 0:
        #    allax.scatter(*photf_pts, lw=0, c='r')
            #isou.isoplot(allax, xlim=extent[0:2], ylim=extent[2:], shade = vvcrit, alpha=0.6, ls='--', c='k')
            #isol.isoplot(allax, xlim=extent[0:2], ylim=extent[2:], shade = vvcrit, alpha=0.6, ls='--', c='k')

        savename = os.path.join(plotpath, 'cmd_vvcrit{:.1f}_m2lnP{:.2f}.png'.format(vvcrit, a_cmd.fit))  
        fig.savefig(savename, dpi=300)
        mist_pts = [(isol.x, isol.y), (iso.x, iso.y), (isou.x, isou.y)]

        a_cmd.pgcmd(photf_pts=photf_pts, mist_pts=mist_pts, best_list=best_dict[vvcrit].values())
        fig.clf()

        if a_cmd.fit < bestfit:
            bestfit = a_cmd.fit
            bestvvc = vvcrit
    #=======================================================

    # For plotting all best fits together:
    #=======================================================
    allax.scatter(*photf_pts, lw=0, s=8, c='r')

    # older age isochrone:
    isou = rmm.ISOCMD(qdict['logZ'][1], bestvvc, ebv=best_dict[bestvvc]['ebv']) #qdict['logZ'][1], float('{:.1f}'.format(qdict['vvcrit'][1]))
    isou.set_isodata(qdict['lage'][1],filters[1], filters[0], dmod=best_dict[bestvvc]['dmod'])
    # younger age:
    isol = rmm.ISOCMD(qdict['logZ'][0], bestvvc, ebv=best_dict[bestvvc]['ebv']) #qdict['logZ'][0], float('{:.1f}'.format(qdict['vvcrit'][0]))
    isol.set_isodata(qdict['lage'][0],filters[1], filters[0], dmod=best_dict[bestvvc]['dmod'])

    isou.isoplot(allax, xlims=extent[0:2], ylims=extent[2:], shade = bestvvc, alpha=0.6, ls='--', c='k')
    isol.isoplot(allax, xlims=extent[0:2], ylims=extent[2:], shade = bestvvc, alpha=0.6, ls='--', c='k')

    savename = os.path.join(plotpath, 'cmdall_bestvvcrit{:.1f}.png'.format(bestvvc))
    allfig.savefig(savename, dpi=300)
    allfig.clf()
    #=========================================================
