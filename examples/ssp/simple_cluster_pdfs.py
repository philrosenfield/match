from __future__ import print_function
import numpy as np
import pylab as plt
import pdb
import sys
import glob

''' Takes the output from a MATCH SSP run and produces various margainlized distributions '''
''' This is a simplier version of the make_cluster_pdfs.py code, which requires more columns and doesn't yet check for uniqueness '''


''' 
> python simple_cluster_pdfs.py <sspfilename>

columns in sspfile:
0 - Av
1 - IMF Slope
2 - dmod
3 - log Age
4 - logZ
5 - fit value
6 - SFR
7 - background scale factor
'''


# cat syn02.all | awk '$2==av {x+=exp(0.5*(187-$6))} END {print x}' av=5.05 -- Andy's 1d marginalization using AWK


# make 1D pdfs given a set of values and fit values
def makepdfs(x, weights):
    # get unique grid values
    vals = np.unique(x)
    prob = np.zeros(len(vals))
    # compute linear probabilites
    absprob =  np.exp(0.5*(weights.min() - weights))
    # sum over probabilites for each unique grid value
    for i, j in enumerate(vals):
        prob[i] = np.sum(absprob[x == j])
    prob /= prob.sum()
    return vals, prob, np.log(prob)

# make 2D joint distributions from x,y grids and fit values
def make2dpdfs(x, y, weights):
    vals1, vals2 = np.unique(x), np.unique(y)
    prob = np.zeros(len(vals1) * len(vals2))
    vals_x, vals_y = np.zeros(len(vals1) * len(vals2)), np.zeros(len(vals1) * len(vals2))
    absprob = np.exp(0.5*(weights.min() - weights))
    k=0
    for i, j in enumerate(vals1):
        for m, n in enumerate(vals2):
            prob[k] = np.sum(absprob[np.where((x == j) & (y == n))])
            vals_x[k] = j
            vals_y[k] = n
            k+=1
            #pdb.set_trace()
    prob /= prob.sum()
    return vals_x, vals_y, prob, np.log(prob)

if __name__ == '__main__':

    cm = plt.get_cmap('Blues')

    # name of match console input file
    name = sys.argv[1]
    
    # read in match ssp file
    data = np.loadtxt(name)

    
    plt.close('all')
    
    # 1-d distributions for av, age, and 2d distribution for age-av
    plt.figure()
    res = makepdfs(data[:,0], data[:,5])
    plt.hist(res[0], weights=res[1], bins=len(res[0]), histtype='step', lw=4, color='k', align = 'mid')
    plt.xlabel('Av')
    plt.savefig(name+'av_1d_pdf.png', dpi=300, bbox_inches='tight')

    plt.figure()
    res = makepdfs(data[:,3], data[:,5])
    plt.hist(res[0], weights=res[1], bins=len(res[0]), histtype='step', lw=4, color='k', align = 'mid')
    plt.xlabel('Log(Age)')
    plt.savefig(name+'age_1d_pdf.png', dpi=300, bbox_inches='tight')

    res2 = make2dpdfs(data[:,0], data[:,3], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('Log(Age)')
    plt.xlabel('Av')
    cb.set_label('Probability')
    plt.savefig(name+'_joint_av_age.png', dpi=300, bbox_inches='tight')
    

    # 1-d distribution for dmod and 2d distribution for age-dmod
    plt.figure()
    res = makepdfs(data[:,2], data[:,5])
    plt.hist(res[0], weights=res[1], bins=len(res[0]), histtype='step', lw=4, color='k', align = 'mid')
    plt.xlabel('dmod')
    plt.savefig(name+'dmod_1d_pdf.png', dpi=300, bbox_inches='tight')

    res2 = make2dpdfs(data[:,3], data[:,2], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.xlabel('Log(Age)')
    plt.ylabel('dmod')
    cb.set_label('Probability')
    plt.savefig(name+'_joint_dmod_age.png', dpi=300, bbox_inches='tight')





    plt.show()



    
    #print(res[1])
    plt.show()
    
