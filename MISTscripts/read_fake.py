import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from fileio import *

def plt_truth(ax, photbase, param, vvc=0.0, justval=False):

    paramf = glob.glob(os.path.join(get_workdir('Fakedata'), 'fake.param_{:s}'.format(photbase)))[0]
    with open(paramf, 'r') as f:
        flines = f.readlines()

    firstline = flines[0]
    lastline = flines[-1]

    if ax != None:
        xlims = [1e-1,1e8]
        x = np.linspace(*xlims, num=1e2)
        ax.set_xlim(xlims)

    if param == 'lage':
        # get the age bin:
        truth = map(float, lastline.split(' ')[-4:-2]) 
        if justval:
            return truth

        # Plot the mean of the best fits as a horizontal line
        #ax.plot(x, [mean]*len(x), label='Mean of best fits')
        # Plot the "truth" as another horizontal line:
        ax.plot(x, [truth[0]]*len(x), c='k', ls='--')
        ax.plot(x, [truth[1]]*len(x), c='k', ls='--')
        ax.fill_between(xlims, truth[0], truth[1], alpha=0.4)

    elif param == 'logZ':
        
        # get the logZ value, and logZ spread:
        logz = float(lastline.split(' ')[-1])
        zspread = float(firstline.split(' ')[3])
        truth = [logz-zspread/2., logz+zspread/2.]
        if justval:
            return truth

        # Plot the mean of the best fits as a horizontal line
        #ax.plot(x, [mean]*len(x), label='Mean of best fits')
        # Plot the "truth" as another horizontal line:
        #ax.plot(x, [truth]*len(x), c='k')
        ax.plot(x, [truth[1]]*len(x), c='k', ls='--')
        ax.plot(x, [truth[0]]*len(x), c='k', ls='--')
        ax.fill_between(xlims, truth[0], truth[1], alpha=0.4)

    elif param == 'vvcrit':
        
        # as of now the truth is 0.0...
        truth = float(vvc)
        if justval:
            return truth

        # Plot the mean of the best fits as a horizontal line
        #ax.plot(x, [mean]*len(x), label='Mean of best fits')
        # Plot the "truth" as another horizontal line:
        ax.plot(x, [truth]*len(x), c='k')

    return truth
