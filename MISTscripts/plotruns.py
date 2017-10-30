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
from matplotlib.patches import Rectangle
import seaborn as sns
# no seaborn gridlines on plot, and white bg:
sns.set_context('paper')
sns.set(font='serif')
sns.set_style("white", {'font.family' : 'serif', 'font.serif': ['Times', 'Palatino', 'serif']})
plt.axes(frameon=False)
from MIST_codes.scripts import read_mist_models as rmm
from read_fake import plt_truth
from fileio import *


# Add something to perform diagnostic checks...
# + show that the results are accurate via e.g. plot errors and values as function of star number.
# + show sensitivity to completeness?
# + show effects of alteration of binsize
# + show effects of altering binary fraction? (maybe useful)

# General idea: show that match is accurate (or not) in smallish survey size regime that these clusters possess.


def marginalize(clustername, photobase, jobidstr, outname, params = ['lage', 'vvcrit', 'logZ'], 
                local=False, saveagevsvvc=True,debug=False, interp=True, incvvc='all'):

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

    #global match_dir
    #global scrn_dir
    #global outpath
    #global plotpath

    # There are two directories here because on my local machine, I keep supercomp
    # output in a separate dir that I DL to, and local runs w/ MATCH are stored in
    # their own location.
    #if not local:
        # downloaded runs (wouldn't be used on remote machine):
    #    match_dir = os.path.join(os.environ['MATCH_ODYOUT'], clustername)
    #else:
        # local runs (should be used on remote machine, or your local machine):
    #    match_dir = get_workdir(data=clustername)#os.path.join(os.environ['MATCH_DIR'], 'MIST', clustername)

    # path to MATCH .scrn files:
    #scrn_dir = os.path.join(get_scrndir(data=clustername), photobase, 'SLURM{:s}'.format(jobidstr))
    # path to MATCH .cmd and .out files:
    #outpath = os.path.join(get_outdir(data=clustername), photobase, 'SLURM{:s}'.format(jobidstr))
    # path to place output plots in:
    #plotpath = os.path.join(outpath, 'plots')

    # list of .scrn files:
    if incvvc != 'all':
        flist = glob.glob(os.path.join(scrn_dir, '*vvcrit{:.1f}*.scrn'.format(incvvc)))
        params = ['lage', 'logZ']
    else:
        flist = glob.glob(os.path.join(scrn_dir, '*.scrn'))

    #print(scrn_dir)
    #print(flist)

    # print .scrn files found:
    #for afile in flist:
        #print("Reading {:s}...".format(afile))

    # calling MATCH scripts' combine_files() to combine various v/vcrit runs:
    ssp.combine_files(flist, outfile=outname)

    # Create a match SSP object from the combined scrns of the given runs:
    combineddata = ssp.SSP(outname)

    # Corner plots:
    # 'lage', 'vvcrit', 'logZ', 'dmod', 'Av'
    pdffig, pdfax, quant = combineddata.pdf_plots(params, twod=False, quantile=True, cmap=plt.cm.Reds, interpolateq=interp, debug=debug)
    plt.close()
    pdffig, pdfax, quant = combineddata.pdf_plots(params, twod=True, quantile=True, cmap=plt.cm.Reds, interpolateq=interp, debug=debug)
    plt.savefig(os.path.join(plotpath, 'pdfcornerplot.png'))
    plt.close()

    # quantile dictionary (the quantiles will be used to draw bounding isochrones for the best-fits at ~+/- 1 sigma parameter values):
    qdict = {key: value for (key, value) in zip(params, quant)}#{'lage': quant[0], 'vvcrit': quant[1], 'logZ':quant[2]}
    #print(qdict)

    # remove the combined .csv file since it's no longer needed:
    os.remove(os.path.join(os.getcwd(), outname))

    # A plot of the best fit ages found vs. the corresp. v/vcrit values...
    # This is also getting the name of the saved figure for the plot, the matplotlib figure for the plot,
    # the axis object too, the best parameters for each v/vcrit run and the corresp. v/vcrits:
    vvcritagename, bestdict, vvcrits = agevsvvcrit(flist, jobidstr, outpath = plotpath, save=saveagevsvvc) # vbests, vvcrits
    # Move the best fit ages vs v/vcrit plot:
    if saveagevsvvc:
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

def agevsvvcrit(flist, jobarr_id, outpath, ext='.png', save=True):

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

    savename = jobarr_id + '_BestAgevsVVCrit' + ext
    if save:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(vvcrits, bestages)
        ax.set_ylabel('Best-fit Age')
        ax.set_xlabel('v/v_crit')
        plt.savefig(savename)

    #bestlst = [bestages, bestfehs, bestdmods, bestavs]

    return savename, best_dict, vvcrits#bestlst, vvcrits

def getbest(param, data):

    """
       This method gets the best-fit (max likelihood) parameter value from an ssp object (SSP class defined 
       in Phil R. and Dan W.'s -- RW's --  match scripts).

       Argument:
       ---------
           + param (str): String for the parameter of interest (would be e.g.: lage, logZ, dmod, see 
              marginalize method of SSP class in RW's scripts)
           + data (SSP): RW script's SSP class instance.

       Returns:
       --------
           + best (float): best-fit value of the parameter corresponding to the string designated by the 
              'param' argument.
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

# auxillary function for main(). Not majorly useful on its own.
def plotisocmd(ax, vvcrit, best_dict, filters, truths, extent):

    color_name = '{:s}-{:s}'.format(filters[0], filters[1])
    redmag_name = filters[1]
    bluemag_name = filters[0]

    # best isochrone (for current v/vcrit):
    iso = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'], vvcrit, ebv=best_dict[vvcrit]['ebv'], exttag='TP')
    iso.set_isodata(best_dict[vvcrit]['lage']['val'], color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])

    # older, redder isochrone (1 sigma away in metallicity and age):
    isou = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'] + best_dict[vvcrit]['feh']['err'], vvcrit, 
                      ebv=best_dict[vvcrit]['ebv'], exttag='TP')

    isou.set_isodata(best_dict[vvcrit]['lage']['val'] + best_dict[vvcrit]['lage']['err'], 
                     color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])

    # younger, bluer isochrone (1 sigma away in [Fe/H], age but in other direction):
    isol = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'] - best_dict[vvcrit]['feh']['err'], vvcrit, 
                                ebv=best_dict[vvcrit]['ebv'], exttag='TP')
    isol.set_isodata(best_dict[vvcrit]['lage']['val'] - best_dict[vvcrit]['lage']['err'], 
                     color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])

    # plot best, 1 sig older age/[Fe/H], 1 sig lower age/[Fe/H] isochrones.
    iso.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, label=True)
    isou.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')
    isol.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')

    if truths is not None:
        # also plot "true" isochrones if told to (assumes 'truth' equals a +/- bound for age, met. 
        # and a value for vvc):
        true_isou = rmm.ISOCMD(truths['feh'][1], truths['vvc'], ebv=best_dict[vvcrit]['ebv'], exttag='TP')
        true_isou.set_isodata(truths['age'][1], color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])
        true_isou.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], alpha=0.6, c='k', ls=':')

        true_isol = rmm.ISOCMD(truths['feh'][0], truths['vvc'], ebv=best_dict[vvcrit]['ebv'], exttag='TP')
        true_isol.set_isodata(truths['age'][0], color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])
        true_isol.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], alpha=0.6, c='k', ls=':')                   
  
    #fig.savefig(savename, dpi=300)
    # maintain collection of isochrone pts. used in plot. Passed to Phil's pgcmd code to overlay
    # isochrones on Hess diagrams.
    #mist_pts = [(isol.x, isol.y), (iso.x, iso.y), (isou.x, isou.y)]
    mist_pts = [(isol.get_data([color_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'] - best_dict[vvcrit]['lage']['err'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0], 
                 isol.get_data([bluemag_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'] - best_dict[vvcrit]['lage']['err'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0]),
                (iso.get_data([color_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0], 
                 iso.get_data([bluemag_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0]),
                 (isou.get_data([color_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'] + best_dict[vvcrit]['lage']['err'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0], 
                 isou.get_data([bluemag_name], phasemask=[], 
                               lage=best_dict[vvcrit]['lage']['val'] + best_dict[vvcrit]['lage']['err'], 
                               dmod=best_dict[vvcrit]['dmod']).values()[0])
                 ]

    return ax, iso, mist_pts

def main(cluster_name, photbase, jobid, filters, dmod, ebv, local, star_num, truth, incvvc, savebest=True, bestax=None, justbest=False):

    global match_dir
    global scrn_dir
    global outpath
    global plotpath

    match_dir = get_workdir(data=cluster_name)
    # path to MATCH .scrn files:
    scrn_dir = os.path.join(get_scrndir(data=cluster_name), photbase, 'SLURM{:s}'.format(jobid))
    # path to MATCH .cmd and .out files:
    outpath = os.path.join(get_outdir(data=cluster_name), photbase, 'SLURM{:s}'.format(jobid))
    # path to place output plots in:
    plotpath = os.path.join(outpath, 'plots')
    #print(outpath)

    #vbests, vvcrits, qdict, bfdict = marginalize(cluster_name, photbase, sys.argv[3], sys.argv[4], local=local)
    #vvcrits, best_dict, qdict, bfdict = marginalize(cluster_name, photbase, sys.argv[3], sys.argv[4], local=local, incvvc=incvvc)
    # by default, readssp() will return an AVERAGE (over all runs) of the bestfit parameter value and uncertainty.
    # the original, unaveraged values are derived from the sspcombine solutions.
    if incvvc == 'all':
        vvcrits = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    else:
        vvcrits = [incvvc]
    best_dict = {}
    for vvc in vvcrits:

        # read sspcombine output for age results (subsampled solution being used):
        starnums, sspbest_age, ssp_ageerr = readssp(data=cluster_name, solution='subsampled', 
                                                    fromfit=False, param='age', photname_base=photbase, incvvc=vvc)
        # read sspcombines for metallicity results:
        starnums, sspbest_logz, ssp_logzerr = readssp(data=cluster_name, solution='subsampled', 
                                                      fromfit=False, param='logZ', photname_base=photbase, incvvc=vvc)

        # dictionary whose keys are the corresponding vvc values for the best fit values found above (dmod, Av fixed):
        best_dict[vvc] = {'lage': {'val': float(sspbest_age[0]), 'err': float(ssp_ageerr[0])}, 
                          'feh': {'val': float(sspbest_logz[0]), 'err': float(ssp_logzerr[0])}, 
                          'dmod': dmod, 
                          'ebv': ebv}

    # if told to plot "truth" values:
    if truth:
        if 'Fakedata' in cluster_name:
            # true age bin (list), logZ +/- spread (list), v/vcrit (float) of artificial data fit to:
            true_age = plt_truth(None, cluster_name, photbase, 'lage', justval=True)
            true_feh = plt_truth(None, cluster_name, photbase, 'logZ', justval=True)
            true_vvc = plt_truth(None, cluster_name, photbase, 'vvcrit', justval=True)

        elif 'Hyades' in cluster_name:
            # Perryman et al.1998.
            true_age = [np.log10(575e6), np.log10(675e6)]
            # [Fe/H] -- e.g. deB01 says 0.14 +/- 0.05
            true_feh = [0.15, 0.15]
            # assume canonically non-rotating:
            true_vvc = 0.0
      
        elif 'Praesepe' in cluster_name:
            # source???.
            true_age = [np.log10(575e6), np.log10(675e6)]
            # [Fe/H] source???
            true_feh = [0.15, 0.15]
            # assume canonically non-rotating:
            true_vvc = 0.0

        elif 'Pleiades' in cluster_name:
            # source???
            true_age = [np.log10(100e6), np.log10(125e6)]
            # [Fe/H] e.g. Soderblom+ 2009 says 0.03 +/- 0.02
            true_feh = [0.00, 0.00]
            # assume canonically non-rotating:
            true_vvc = 0.0

        truths = {'age': true_age, 'feh': true_feh, 'vvc': true_vvc}
    else:
        truths = None

    #plt.clf()
    # clearing out pre-existing output files that this script may have produced in prev. runs:
    if not justbest:
        for f in glob.glob(os.path.join(outpath,'cmd*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(plotpath,'cmd*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(plotpath,'*.txt')):
            os.remove(f)

    # if req., correct 2MASS filter names to those recognizable by MIST.
    for i, afilter in enumerate(filters):
        if afilter == 'J':
            filters[i] = '2MASS_J'
        elif afilter == 'Ks':
            filters[i] = '2MASS_Ks'
        elif afilter == 'B':
            filters[i] = 'Bessell_B'
        elif afilter == 'V':
            filters[i] = 'Bessell_V'


    # Retrieve the magnitude limits used in the fit (from the calcsfh .param file used):
    param_fname = glob.glob(os.path.join(match_dir, 'csfh_param', 'calcsfh_{:s}.param'.format(photbase)))[0]
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
        exclusions = lines[7].split('\n')[0]
        ntbins= int(lines[8])
        tbins=lines[9:8+ntbins]
        for i, aline in enumerate(tbins):
            tbins[i] = map(float, aline.split('\n')[0].split(' ')[-2:])

    # Manage exclude gates (right now only works for one region):
    # If there are exclude gates...
    if int(exclusions[0]) > 0:
        # Get the exclude gate points:
        expts = map(float, exclusions.split(' ')[1:-1])
        ex_x = expts[::2]
        ex_y = expts[1::2]
        ex_xy = zip(ex_x, ex_y)

    else:
        expts = None
        

    # Also get the residuals (???):
    # For overplotting data points, get the data points CMD x, y values from the photometry file:
    photf_v = np.genfromtxt(os.path.join(match_dir, 'photometry', '{:s}.phot'.format(photbase)), usecols=(0,))
    photf_i = np.genfromtxt(os.path.join(match_dir, 'photometry', '{:s}.phot'.format(photbase)), usecols=(1,))
    photf_vmi = photf_v - photf_i
    # Tuple storing the x, y values for plotting the data on a CMD (using color on x, red filter on y):
    photf_pts = (photf_vmi, photf_i)
    # match uses blue filter on y axis.
    photf_pts_alt = (photf_vmi, photf_v)

    # store number of data pts if told to:
    if star_num:
        nstars = len(photf_vmi)

    #plt.clf()

    bestfit = 99999.9
    # Plot best isochrone for each v/vcrit individually:
    #=================================================================
    allfig = plt.figure()
    allax = allfig.add_subplot(111)

    # iterate through all v/vcrits used for fitting:
    for j, vvcrit in enumerate(vvcrits):

        #print(outpath)
        # get the .cmd file for making a MATCH pg style plot:
        fname = glob.glob(os.path.join(outpath,'*vvcrit{:.1f}*.cmd'.format(vvcrit)))[0]
        # use Phil's code to make a CMD object from the .cmd file:
        a_cmd = cmd.CMD(fname)
        # store the axis extents of the CMD data
        extent = a_cmd.extent

        if not justbest:   

            fig = plt.figure()
            ax = fig.add_subplot(111)

            # scatter plot of stars from data file (observations):
            # w/ num of data points in the plot legend...
            if star_num:
                ax.scatter(*photf_pts, lw=0, s=8, c='r', label=r"$N_*$ = "+"{:d}".format(nstars))
                ax.legend(loc = 'best')
            # just the scatter plot...
            else: 
                ax.scatter(*photf_pts, lw=0, s=8, c='r')
            # X out masked data points.
            if expts != None:
                for k in range(len(photf_pts[0])):
                    if (min(ex_x) < photf_vmi[k] < max(ex_x)) & (min(ex_y) < photf_v[k] < max(ex_y)):
                        ax.scatter(photf_vmi[k], photf_i[k], s=12, c='k', marker='x')  

            # Plot MIST isochrones corresponding to MATCH best-fit parameters.
            ax, iso, mist_pts = plotisocmd(ax=ax, vvcrit=vvcrit, best_dict=best_dict, filters=filters, truths=truths, extent=extent)

            # saving the plot of data points w/ isochrones overlaid:
            savename = os.path.join(plotpath, 'cmd_vvcrit{:.1f}_m2lnP{:.2f}.png'.format(vvcrit, a_cmd.fit))
            fig.savefig(savename, dpi=300)

            # create a MATCH pg style plot using the .cmd file:
            # using photf_pts_alt because MATCH does its CMDs with V-I vs. V & not sure if this is changeable.
            a_cmd.pgcmd(photf_pts=photf_pts_alt, mist_pts=mist_pts, 
                        best_list=[best_dict[vvcrit]['lage']['val'], best_dict[vvcrit]['feh']['val']])
        
            # closing figure before moving on to next plots.
            fig.clf()

            #==================================================================
            # For plotting all best fits together (on another matplotlib axis):
            #==================================================================
            # Only scatter plot data points once (on axis 2/ all isochrones).
            if j == 0:
                allax.scatter(*photf_pts, lw=0, s=8, c='r')
                # Mark excluded region(s)
                if expts != None:
                    for k in range(len(photf_pts[0])):
                        if (min(ex_x) < photf_vmi[k] < max(ex_x)) & (min(ex_y) < photf_v[k] < max(ex_y)):
                            allax.scatter(photf_vmi[k], photf_i[k], s=12, c='k', marker='x') 

            # all best isochrones are plotted together (no +/- 1 sigma isochrones are plotted):
            iso.isoplot(allax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, legloc='upper right')

        # track which cmd has the overall highest lnP & corresponding v/vcrit.
        if a_cmd.fit < bestfit:

            bestfit = a_cmd.fit
            bestvvc = vvcrit




        
    # write best fit solutions +/- uncertainties to a summary file:
    # -------------------------------------------------------------
    outlines = ['Hyades\n','_______\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n',
                '2MASS J, Ks\n','-------\n', '\n', '\n', 
                '\n','Praesepe\n','________\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n', 
                '2MASS J, Ks\n','-------\n', '\n', '\n',
                '\n','Pleiades\n_','_______\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n', 
                '2MASS J, Ks\n','-------\n', '\n', '\n', '\n']

    with open('/n/home12/sgossage/match_results.txt', 'rw+') as sumfile:

        # get current lines of the summary file:
        inlines  = sumfile.readlines()
        #print(inlines)
        #print(len(inlines))
        #print(len(outlines))

        for i, line in enumerate(inlines):
            # maintain lines that already contain best fits
            if 'Best' in line:
                print(line)
                outlines[i] = line

        # update best fit line for the current run:
        if 'Hyades' in cluster_name:
            n = 4
        elif 'Praesepe' in cluster_name:
            n = 16
        elif 'Pleiades' in cluster_name:
            n = 28

        if '2MASSJK' in photbase:
            n += 5

        # print these lines when updating...
        # ...if using a distribution of roation rates...
        if '_rot' in cluster_name:
            n += 1
            outlines[n] = 'Best (rotation distribution) v/vcrit = {:.1f}, ' \
                          'age = {:.1f} +/- {:.1f} Myrs,' \
                          '[Fe/H] = {:.2f} +/- {:.2f}\n'.format(bestvvc, 
                                                                10**best_dict[bestvvc]['lage']['val']/1e6, 
                                                                (10**(best_dict[bestvvc]['lage']['val'] + \
                                                                best_dict[bestvvc]['lage']['err']) - \
                                                                10**best_dict[bestvvc]['lage']['val'])/1e6, 
                                                                best_dict[bestvvc]['feh']['val'],  
                                                                best_dict[bestvvc]['feh']['err'])

        # ...or else if not using a distribution of rotation rates.
        else:
            outlines[n] = 'Best v/vcrit = {:.1f}, ' \
                          'age = {:.1f} +/- {:.1f} Myrs,' \
                          '[Fe/H] = {:.2f} +/- {:.2f}\n'.format(bestvvc, 
                                                                (10**best_dict[bestvvc]['lage']['val'])/1e6, 
                                                                (10**(best_dict[bestvvc]['lage']['val'] + \
                                                                best_dict[bestvvc]['lage']['err']) - \
                                                                10**best_dict[bestvvc]['lage']['val'])/1e6, 
                                                                best_dict[bestvvc]['feh']['val'], 
                                                                best_dict[bestvvc]['feh']['err'])

        # return to top of file & then write lines out.
        sumfile.seek(0)
        sumfile.writelines(outlines)

    # After all, save image of the plot containing all isochrones plotted side by side.
    savename = os.path.join(plotpath, 'cmdall_bestvvcrit{:.1f}.png'.format(bestvvc))
    allfig.savefig(savename, dpi=300)
    allfig.clf()
    #=========================================================


    if bestax == None:
        bestfig = plt.figure()
        bestax = bestfig.add_subplot(111)

    bestax.scatter(*photf_pts, lw=0, s=8, c='r')
    if expts != None:
        for k in range(len(photf_pts[0])):
            if (min(ex_x) < photf_vmi[k] < max(ex_x)) & (min(ex_y) < photf_v[k] < max(ex_y)):
                bestax.scatter(photf_vmi[k], photf_i[k], s=12, c='k', marker='x') 

    bestax, bestiso, mist_pts = plotisocmd(ax=bestax, vvcrit=bestvvc, best_dict=best_dict, filters=filters, truths=truths, extent=extent)

    if savebest:
        # save the figure.
        bestfig.savefig(os.path.join(plotpath, 'bestcmd_vvcrit{:.1f}_m2lnP{:.2f}.png'.format(vvcrit, a_cmd.fit)), dpi=300)
        bestfig.clf()

    return bestax, allax

if __name__ == "__main__":

    """
        This is a script runnable as e.g.

        >> ./plotruns.py Hyades hyades_debruTIGS_TychoBV.phot 12345678 somename.csv Tycho_B Tycho_V 
     
        The intent is to plot the best-fit models found via a MATCH run; this makes use of Phil R./Dan 
        W.'s (RW) match scripts and also Jieun Choi's MIST codes (available on her Github site). The 
        former is mainly used to conglomerate MATCH output across separate runs made at various v/vcrit 
        values (as this parameter is not explored in the usual way by MATCH, RW's scripts concatenate 
        the output of the runs and marginalize to acheive a best-fit over any new parameters not used 
        within MATCH itself, e.g. v/vcrit). The latter is used to plot the best-fit MIST models found 
        by the run(s).

    """

# NOTE:
# sys.argv[x] = cluster name, photometry file, job array slurm id, outputfilename, bluer band filter, redder band filter, dmod, E(B-V), and then optional flags.
    cluster_name = sys.argv[1]
    photbase = os.path.splitext(sys.argv[2])[0]
    jobid = sys.argv[3]

    # blue and red filters:
    v_filter = sys.argv[5]
    i_filter = sys.argv[6]
    # organize filter inputs from user (blue, red is the expected order):
    filters = [v_filter, i_filter]

    # dmod and E(B-V) used:
    dmod = float(sys.argv[7])
    ebv = float(sys.argv[8])

    # check for script options from user:
    if '-local' in sys.argv[9:]:
        local = True
    else:
        local = False
    if '-nstar' in sys.argv[9:]:
        star_num = True
    else:
        star_num = False
    if '-truth' in sys.argv[9:]:
        truth = True
    else:
        truth = False

    for arg in sys.argv[9:]:
        print('.')
        if '-vvcin=' in arg:
            incvvc = float(arg.split('-vvcin=')[-1])
            break
        else:
            incvvc = 'all'

    ax, allax = main(cluster_name=cluster_name, photbase = photbase, jobid = jobid, filters = filters, dmod = dmod, ebv = ebv,
                     local=local, star_num=star_num, truth = truth, incvvc = incvvc)
