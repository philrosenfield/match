import os
import subprocess
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pyl
from plotruns import marginalize
from fileio import *
from read_fake import *
from read_phot import *

def readssp(data='Fakedata', SFR='0.001', solution='single', fromfit=False, param='age', photfn=None, sub=None):

    workdir = get_workdir(data)

    if param=='lage':
        param = 'age'

    # The photometry directory is used to count the number of stars used; the ssp output directory
    # is used to extract the solutions found by MATCH.
    photdir = os.path.join(workdir, "photometry")
    sspoutdir = os.path.join(workdir, "sspcombines")

    # If desired, look in a specified subdirectory for the star num and best fits:
    if sub != None:
        photdir = os.path.join(photdir, sub)
        sspoutdir = os.path.join(sspoutdir, sub)

    # Get the number of stars in each file
    #sfrdir = os.path.join(photdir, "SFR{:s}")
    #os.chdir(sfrdir)
    os.chdir(photdir)
    #if not bins:
    if photfn is None:
        photfn = "SFR{:s}_logZ0.05.phot".format(SFR)
    #else:
    #    photfn = "SFR{:s}_logZ0.05_dmag{:.2f}_dcol{:.2f}.phot".format(SFR, bins[0], bins[1])

    photfn_base = photfn.split('.phot')[0]
    photfile = glob.glob(photfn)[0]

    starnum = file_len(photfile)

    # Now go to the sspcombine output directories and get ages, etc.
    #sfrdir = os.path.join(sspoutdir, "SFR"+SFR)
    sspoutdir = os.path.join(sspoutdir, photfn_base, 'SLURM12345678')
    #os.chdir(sfrdir)
    os.chdir(sspoutdir)
    sspfiles = glob.glob("*.ssp")

    # age solutions (could be adapted to any parameter):
    ages = [] 
    # associated uncertainties:
    uncerts = []

    badi = []
    # Look through all files at the current SFR:
    print(sspfiles) #debug
    for j, f in enumerate(sspfiles):
        # only plotting the v/vcrit = 0.0 solutions with this if statement.
        if 'vvcrit0.0' not in f:
            continue
        print(f)
        # Figure out which file matches the current jobid:
        # Open the file:
        currf = open(f)
        print("Reading {:s}".format(f))
        lines = currf.readlines()
          
        # Search for the desired solution:
        
        if solution == 'single':
            if 'Single solution failed; points not consistent with normal distribution\n' in lines:
                print 'Points inconsistent with normal distribution. No single solution in {:s}.'.format(f)
                # next file:
                continue
            else:
                inblock = False
                for line in lines:
                    if 'Estimate from single solution' in line:
                        inblock = True
                    if (param + ' =' in line) & inblock:
            
                        ages.append(float((line.split("+")[0]).split("=")[-1]))
                        uncerts.append(float(line.split("+/-")[-1].split('(')[0]))
                        break
                    elif ('\n' == line) & inblock:
                        # If there is no solution from fit:
                        print "**No solution from fit in {:s}; it will be excluded.**".format(f)
                        badi.append(j)
                        inblock = False
                        break
            continue

        if solution == 'marginalized':
            if 'Single solution failed; points not consistent with normal distribution\n' in lines:
                print 'Points inconsistent with normal distribution. No single solution in {:s}.'.format(f)
                continue
            else:
                inblock = False
                for line in lines:
                    if 'Estimates from marginalized solutions' in line:
                        inblock = True
                    if (param + ' =' in line) & inblock:
            
                        ages.append(float((line.split("+")[0]).split("=")[-1]))
                        uncerts.append(float(line.split("+/-")[-1].split('(')[0]))
                        break
                    elif ('\n' == line) & inblock:
                        # If there is no solution from fit:
                        print "**No solution from fit in {:s}; it will be excluded.**".format(f)
                        badi.append(j)
                        inblock = False
                        break
            continue

        if solution == 'subsampled':
            #if 'Single solution failed; points not consistent with normal distribution\n' in lines:
                #print 'Points inconsistent with normal distribution. No single solution in {:s}.'.format(f)
                #continue
            #else:
            # Go line by line and look for the block containing the subasmpled marginalized distribution
            #print(lines)
            inblock = False
            for k, line in enumerate(lines):
                if 'Subsampled marginalized distributions' in line:
                    print('Getting subsampled marginalized solution...')
                    inblock = True
                #print(line)
                if param in line and inblock:
                    if fromfit:
                        if line[k+1] != '\n':
                            ages.append(float((line[k+1].split("+")[0]).split("=")[-1]))
                            uncerts.append(float(line[k+1].split("+/-")[-1]).split('from fit')[0])
                            break

                    else:
                        ages.append(float((line.split("+")[0]).split("=")[-1]))
                        uncerts.append(float(line.split("+/-")[-1].split('(')[0]))
                        break

                elif k == len(lines)-1:
                    print "**No solution from fit in {:s}; it will be excluded.**".format(f)
                    badi.append(j)
                    inblock = False
                    break
            continue


    # switch back to the work directory and return the found solutions and uncertainties, along with star number used
    # in finding the respective solutions:
    os.chdir(workdir)
    print('*')
    print(badi)
    #for i, badindex in enumerate(badi):
        #del uncerts[badindex - i]#uncerts[badindex - i]
        #del ages[badindex - i]#ages[badindex - i]
        #del starnums[badindex]#starnums[badindex - i]
    print('**')
    starnums = np.array([starnum]*len(ages))
    return np.array(starnums), np.array(ages), np.array(uncerts)

    # Count lines in a file:
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#SFRS=['0.01','0.001','0.0001', '0.00001', '0.000007', '0.000005', '0.000001']
def plotouts(params, data='Fakedata', SFRS=['0.01','0.001'], vvcs=[None],
             jobid='12345678', local=True, solution='single', fromfit=False, sub=None, showplt=True, save_path=None, debug=False, marg=True):

    """
       X binlst may be unecessary; its useage was meant to showcase how different bin sizes may affect things.

       params given as a dic
    """

    print ""
    vvcs = vvcs*len(SFRS)
    # Close pre-existing figures:
    plt.close("all")
    fig = plt.figure(figsize=(16,9))
    subplt_n = int('1{:d}1'.format(len(params)))
    axes = []
    #ax = fig.add_subplot(111)

    # sys.argv[x] = cluster name, photometry file, job array slurm id, outputfilename, bluer band filter, redder band filter.
    #if photbase != None:
    #    vvcrits, best_dict, qdict, bfdict = marginalize('Fakedata', photbase, jobid, '12345678.csv', local=local)
    #    print(qdict)
    #    print(bfdict)

    # mean is initialized to 0; this will be the mean of the best fits across all sfrs eventually.
    mean = 0
    maxbestfit = -np.inf
    minbestfit = np.inf
    # Loop through the specified sfrs:
    for i, sfr in enumerate(SFRS):

        if vvcs[i] == None:
            photname_base = 'SFR{:s}_logZ0.05_dmag0.10_dcol0.05'.format(sfr)
        else:
            photname_base = 'SFR{:s}_logZ0.05_dmag0.10_dcol0.05_vvc{:s}'.format(sfr, vvcs[i])
        # get the number of stars (N), age (or other param) best fits, and associated uncertainties:
        print(sfr)
        if marg:
            vvcrits, best_dict, qdict, bfdict = marginalize('Fakedata', photname_base, jobid, '12345678.csv', local=local, saveagevsvvc=False, debug=debug)

        starnum = file_len(get_photof(photname_base))

        # starnum in certain cut...
        phot_vmi, phot_v = get_mags(photname_base, vcuts=[7, 3])
        print('NUM CUT: {:d}'.format(len(phot_v)))

        # Convert age (or other parameter) from log to linear:
        #linages = [10**age for age in ages]

        # valid params are: lage, logZ, vvcrit
        for j, param in enumerate(params):

            if not marg:
                starnums, ages, uncerts = readssp(data, sfr, solution, fromfit, param, photfn='{:s}.phot'.format(photname_base), sub=sub)
                print(ages)

            if marg:
                best=qdict[param][2]
                lerr=[best - qdict[param][0]]
                uerr=[qdict[param][1] - best]
            try:
                ax = axes[j]
            except IndexError:
                ax = fig.add_subplot(subplt_n)
                axes.append(ax)
                subplt_n += 1

            if marg:
                ax.errorbar(starnum, best, yerr=[lerr,uerr], fmt='o', c = 'k')
                _ = plt_truth(ax, photname_base, param, vvc=vvcs[i])
            else:
                ax.errorbar(starnums, ages, yerr=uncerts, fmt='o', c = 'k')
                _ = plt_truth(ax, photname_base, param, vvc=vvcs[i])

            if param=='lage':
                ax.set_ylabel(r"$log Age$")
            elif param=='logZ':
                ax.set_ylabel(r"$[Fe/H]$")
                #ax.set_ylim([0.00, 0.09])
            elif param=='vvcrit':
                ax.set_ylabel(r"$\frac{\Omega}{\Omega_c}$")

            ax.set_xlabel(u"$log N$") 
            #ax.set_ylim([minbestfit-1.5*dminbestfit, maxbestfit+1.5*dmaxbestfit])
            ax.set_xscale('log')

        #loy, upy = ax.get_ylim()
        #if any(1.01*(ages+uncerts) > upy):
        #    upy = max(1.01*(ages+uncerts))
        #elif any(0.99*(ages-uncerts) < loy):
        #    loy = min(0.99*(ages-uncerts))
        #ax.set_ylim([loy, upy])

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
    # "True" age value calculated as the average age of the age bin (in which star formation takes place).
    #truth = np.log10((10**agebin[0] + 10**agebin[1])/2.0)

    # Take mean of the best fit ages across all SFRS involved:
    #mean = mean / len(SFRS)

    #x = np.linspace(0, 1e6, 1e2)
    # Plot the mean of the best fits as a horizontal line
    #ax.plot(x, [mean]*len(x), label='Mean of best fits')
    # Plot the "truth" as another horizontal line:
    #ax.plot(x, [truth]*len(x), c='k', label='Truth')
    #ax.plot(x, [agebin[0]]*len(x), c='k', ls='--', label='Age bin limits')
    #ax.plot(x, [agebin[1]]*len(x), c='k', ls='--')
    #ax.fill_between([0, 1e6], agebin[0], agebin[1], alpha=0.4)
    #for binbound in np.arange(8.0, 9.2, 0.02):
    #    ax.axhline(y=binbound, xmax=0.1, ls = '--')

    #print "Mean Fit = {:f}".format(mean)

    # alternate scaling
    #if param == 'age':
    #    ax2 = ax.twinx()
    #    ax2.set_ylabel('age [Myr]')


    # print subprocess.call("pwd", shell=True)
    output_dir = os.path.join(get_workdir(data), "output")

    # Directory to save plots in:
    if sub == None:
        plot_dir = os.path.join(output_dir, "plots")
    else: 
        plot_dir = os.path.join(os.path.join(output_dir, sub), "plots") 

    # Plot save file name:
   # if fromfit:
   #     save_path = os.path.join(plot_dir, "out_{:s}_{:s}_fromfit.eps".format(param, solution))
   # else:
   #     save_path = os.path.join(plot_dir, "out_{:s}_{:s}.eps".format(param, solution))

    # Check if the save path exists and save an image of the plot there:
    if not os.path.isdir(plot_dir):
        # If the plot save directory does not exist, create it and save there:
        print "Creating plot save directory {:s}.".format(plot_dir)
        os.mkdir(plot_dir)

    #ax.set_title(u"{:s} (Best Fit) vs. $log_{{10}}(N)$".format(param))
    #ax.set_title(u"{:s} (Best Fit) vs. $log_{{10}}(N)$".format(param))
    #ax.legend(loc='lower right', prop={'size':10})
    print(save_path)
    if save_path != None:
        print "Saving plot to {:s}.".format(save_path)
        plt.savefig(save_path, dpi=300)

    # Display the plot:
    if showplt:
        plt.show()

    return
