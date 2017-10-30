import glob
import numpy as np
import os
import matplotlib.pyplot as plt
from fileio import *

def plt_truth(ax, data, photbase, param, vvc=0.3, justval=False):

    assert data == 'Fakedata' or data == 'Fakedata_rot', "Attempting to plot truth with data from {:s}/ (need artificial data to define \"truth\").".format(data)

    paramf = glob.glob(os.path.join(get_workdir(data), 'fake_params', 'fake.param_{:s}'.format(photbase)))[0]
    with open(paramf, 'r') as f:
        flines = f.readlines()

    firstline = flines[0]
    lastline = flines[-1]

    # fixing plot limits if an axis is given:
    if ax != None and not justval:
        xlims = [1e1,3e4]
        x = np.linspace(*xlims, num=1e2)
        ax.set_xlim(xlims)

    if param == 'lage':
        # get the age bin (** truth in this case is a list of TWO numbers, bounds of 'true' age bin):
        truth = map(float, lastline.split(' ')[-4:-2]) 
        if justval:
            return truth

        # Plot the mean of the best fits as a horizontal line
        # Plot the MATCH age bin using the upper and lower bounds from the param file:
        ax.plot(x, [truth[0]]*len(x), c='k', ls='--')
        ax.plot(x, [truth[1]]*len(x), c='k', ls='--')

        # fill the region between the bounding lines:
        ax.fill_between(xlims, truth[0], truth[1], alpha=0.4)

    elif param == 'logZ':
        
        # get the logZ value, and logZ spread:
        logz = float(lastline.split(' ')[-1])
        zspread = float(firstline.split(' ')[2])
        # ** truth in this case is a list of numbers also, the upper and lower bounds of the log Z bin
        #    that considers the log Z spread required by MATCH.
        truth = [logz-zspread/2., logz+zspread/2.]
        if justval:
            return truth

        # Plot lines for marking the 'true' log Z bounded region
        # and fill the area in between the lines.
        ax.plot(x, [truth[1]]*len(x), c='k', ls='--')
        ax.plot(x, [truth[0]]*len(x), c='k', ls='--')
        ax.fill_between(xlims, truth[0], truth[1], alpha=0.4)

    elif param == 'vvcrit':
        
        # as of now the truth is 0.3...means I create all my artificial data with v/vcrit = 0.3 as the input.
        # should make this more flexible...
        if not isinstance(vvc, float) and not isinstance(vvc, str):
            vvc = 0.3

        # ** truth in this case is just a number, representing the input v/vcrit value for the artificial data
        #    the models were fit to.
        truth = float(vvc)
        if justval:
            return truth

        # Plot the truth as a horiz. line marking where the input v/vc lies.
        ax.plot(x, [truth]*len(x), c='k')

    return truth
