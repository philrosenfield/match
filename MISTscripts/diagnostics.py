import os
import subprocess
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pyl
from plotruns import marginalize
from fileio import *
from read_fake import *
from read_phot import *

params = {'axes.labelsize': 20,'axes.titlesize':20, 'text.fontsize': 14, 'legend.fontsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
mpl.rcParams.update(params)
#from MIST_scripts.scripts.read_mist_models import expand_lims
#from MIST_scripts.scripts.read_mist_models import mkpad

# Count lines in a file:
#def file_len(fname):
#    with open(fname) as f:
#        for i, l in enumerate(f):
#            pass
#    return i + 1

#SFRS=['0.01','0.001','0.0001', '0.00001', '0.000007', '0.000005', '0.000001']
def plotouts(params, data='Fakedata', SFRS=['0.01','0.001'], logZ='0.05', vvcs=[None],
             local=True, solution='single', fromfit=False, sub=None, showplt=True, 
             save_path=None, debug=False, marg=True, ssp=False, incvvc='all', filtertag='TychoBV'):

    """
       X binlst may be unecessary; its useage was meant to showcase how different bin sizes may affect things.

       params given as a dic

       Arguments
       ______________________

           + params (str, list): A list of the parameters (log10 Age, log Z, v/vcrit) desired for plotting.
 
           ========
           Optional
           ========

           + data (str): A string to designate which MATCH output data directory to access.

           + SFRS (str, list): List of strings specifying which star formation rates to access.
        
           + logZ (str): String specifying which log Z is approriate for the desired run.
         
           + vvcs: List of v/vcrit values used to create the artificial data desired for plotting.

           + local (bool): Boolean value designating whether or not the files should access the local directory.
  
           + solution (str): MATCH sspcombine solution type to use for extracting uncertainties and best fit values.

           + fromfit (bool): Boolean saying whether or not to use the values 'from fit' when extracting best fit/uncert. values.

           + sub (bool):

           + showplt (bool): Boolean for whether or not to display plots created.

           + save_path (str): String that may be specified to say where plots should be saved.
 
           + debug (bool): Boolean that will turn on several debugging outputs.

           + marg (bool):

           + ssp (bool):

           + incvvc (str):

           + filtertag (str):


    """
    vvcs = vvcs*len(SFRS)
    # Close pre-existing figures:
    plt.close("all")
    fig = plt.figure(figsize=(16,9))
    subplt_n = int('{:d}11'.format(len(params)))
    axarr = []

    data_dir = data#'Fakedata'

    # mean is initialized to 0; this will be the mean of the best fits across all sfrs eventually.
    mean = 0
    maxbestfit = -np.inf
    minbestfit = np.inf

    # Loop through the specified sfrs:
    for i, sfr in enumerate(SFRS):

        best = dict((ele, 0) for ele in params)
        lerr = dict((ele, 0) for ele in params)
        uerr = dict((ele, 0) for ele in params)

        if vvcs[i] == None:
            photname_base = 'SFR{:s}*_logZ{:s}_dmag0.10_dcol0.05_{:s}'.format(sfr, logZ, filtertag)
        else:
            photname_base = 'SFR{:s}*_logZ{:s}_dmag0.10_dcol0.05_vvc{:s}_{:s}'.format(sfr, logZ, vvcs[i], filtertag)

        jobids = get_jobids(data_dir, photname_base)

        for jobid in jobids:
            # get the number of stars (N), age (or other param) best fits, and associated uncertainties:
            starnum = file_len(get_photof(photname_base, data=data_dir))
            if marg:
                # there are issues with interpolation in match scripts when star numbers get very low (N~200).
                if starnum <= 200:
                    interp = False
                else:
                    interp = True
                
                # this is the v/vcrits, dictionary of best fits, quantile dictionary, and best-fit dict containing lnP values corresponding to ONE jobid (SLURM run) for the given set of population parameters.
                vvcrits, best_dict, qdict, bfdict = marginalize(data_dir, photname_base, jobid, '12345678.csv', local=local, saveagevsvvc=False, debug=debug, interp=interp, incvvc = incvvc)
                #print(qdict)

            starnum = file_len(get_photof(photname_base, data=data_dir))

            # starnum in certain cut...
            phot_vmi, phot_v = get_mags(photname_base, data=data_dir, vcuts=[8, 3])
            #print('NUM CUT: {:d}'.format(len(phot_v)))

            # Convert age (or other parameter) from log to linear:
            #linages = [10**age for age in ages]

            # valid params are: lage, logZ, vvcrit
            for j, param in enumerate(params):

                if ssp:
                    # by default, readssp() will return an AVERAGE (over all runs) of the bestfit parameter value and uncertainty.
                    # the original, unaveraged values are derived from the sspcombine solutions.
                    starnums, sspbest, ssperr = readssp(data_dir, sfr, solution, fromfit, param, photname_base=photname_base, sub=sub, incvvc=incvvc)
                    #print(sspbest)

                # acquire best fits and create an average from all runs at a particular SFR/logZ/input vvcrit.
                if marg:
                    # for the best-fit & uncertainties derived via match scripts.
                    best[param]+=qdict[param][2]
                    lerr[param]+=(qdict[param][2] - qdict[param][0])**2
                    uerr[param]+=(qdict[param][1] - qdict[param][2])**2
                    if jobid == jobids[-1]:
                        best[param] /= len(jobids)
                        lerr[param] = np.sqrt(lerr[param])/len(jobids)
                        uerr[param] = np.sqrt(uerr[param])/len(jobids)

                # create a new axis/subplot for each parameter:
                try:
                    ax = axarr[j]
                except IndexError:
                    ax = fig.add_subplot(subplt_n)
                    axarr.append(ax)
                    subplt_n += 1

                
                # plot average created via Phil's marginalization:
                if marg and jobid == jobids[-1]:
                    ax.errorbar(starnum, best[param], yerr=[[lerr[param]],[uerr[param]]], fmt='o', c = 'k')
                    _ = plt_truth(ax, data_dir, photname_base, param, vvc=vvcs[i])
                    print('match scripts: best = {:f}'.format(best[param]))
                    print('match scripts: lerr = {:f}, uerr = {:f}'.format(lerr[param], uerr[param]))
                # plot average according to match sspcombine solutions:
                if ssp and jobid == jobids[-1]:
                    #print(sspbest)
                    #print(len(np.array([sspbest])))
                    if sspbest[0] is not None and np.isfinite(sspbest[0]):
                        ssp_eb = ax.errorbar(starnums, sspbest, yerr=float(ssperr), fmt='o', c = 'r', ecolor='r')
                        ssp_eb[-1][0].set_linestyle(':')
                        _ = plt_truth(ax, data_dir, photname_base, param, vvc=vvcs[i])
                        print('sspcombine: best {:s} = {:f}'.format(param, sspbest[0]))
                        print('sspcombine: uncerts = +/- {:f}'.format(ssperr[0]))

                # axis labeling
                if param=='lage':
                    ax.set_ylabel(r"$log Age$")
                elif param=='logZ':
                    ax.set_ylabel(r"$[Fe/H]$")
                elif param=='vvcrit':
                    ax.set_ylabel(r"$\frac{\Omega}{\Omega_c}$")

                ax.set_xlabel(u"$N$") 
                ax.set_xscale('log')

            #print "Plotted {:d} points.".format(len(ages))
            print "SFR = {:s} ({:d} stars)".format(sfr, starnum)
            #print "scatter (std.dev. of best fits): {:f}".format(np.std(ages))
            #print "mean uncert: {:f}".format(np.mean(uncerts))
            #print "scatter / mean uncert: {:f}".format(np.std(ages)/np.mean(uncerts))

            # Keep tabs on the mean best fit for each SFR
            #sfr_bestfitmean = np.log10(sum(linages)/len(ages))
            #mean += sfr_bestfitmean

            print "=====================================================" 

    plt.tight_layout()

    output_dir = os.path.join(get_workdir(data), "output")

    # Directory to save plots in:
    if sub == None:
        plot_dir = os.path.join(output_dir, "plots")
    else: 
        plot_dir = os.path.join(os.path.join(output_dir, sub), "plots") 

    # Check if the save path exists and save an image of the plot there:
    if not os.path.isdir(plot_dir):
        # If the plot save directory does not exist, create it and save there:
        print "Creating plot save directory {:s}.".format(plot_dir)
        os.mkdir(plot_dir)

    if save_path != None:
        print "Saving plot to {:s}.".format(save_path)
        plt.savefig(save_path, dpi=300)

    # Display the plot:
    if showplt:
        plt.show()

    return
