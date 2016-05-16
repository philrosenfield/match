''' Takes the output from a MATCH SSP run and produces various margainlized distributions '''

''' 
> python make_cluster_pdfs.py <sspfilename>

columns in sspfile:
0 - Av
1 - IMF Slope
2 - dmod
3 - log Age
4 - logZ
5 - fit value
6 - SFR
7 - background scale factor
8 - binary fraction
9 - dAv
'''

from __future__ import print_function
import numpy as np
import pylab as plt
import pdb
import seaborn as sns
import sys
import glob

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


def main():

    cm = plt.get_cmap('Blues')
    sns.set_style("white")
    #cm.set_gamma(0.2)

    # 
    name = sys.argv[1]
    outname = name.split('.')[0]
    
    #

    data = np.loadfromtxt(name)

    res = makepdfs(data[:,1], data[:,5])
    plt.close('all')
    plt.hist(res[0], weights=res[1], bins=31, histtype='step', lw=4, color='k')
    plt.xlim(0.,3.)
    #plt.plot(res[0], res[1], c='k', lw=2)
    plt.xlabel('IMF Slope', fontsize=20)
    plt.ylabel('Probability', fontsize=20)
    plt.savefig(outname+'_marginal_gamma.png', dpi=300, bbox_inches='tight')

    res2 = make2dpdfs(data[:,1], data[:,3], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('Log(Age)')
    plt.xlabel('IMF Slope')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_gamma_age.png', dpi=300, bbox_inches='tight')

    
    plt.figure()
    plt.hist(res2[0], weights = res2[2])

    res2 = make2dpdfs(data[:,1], data[:,0], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('Av')
    plt.xlabel('IMF Slope')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_gamma_av.png', dpi=300, bbox_inches='tight')


    plt.figure()
    plt.hist(res2[0], weights = res2[2])

    res2 = make2dpdfs(data[:,0], data[:,3], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('Log(Age)')
    plt.xlabel('Av')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_av_age.png', dpi=300, bbox_inches='tight')
    

    res2 = make2dpdfs(data[:,1], data[:,8], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('BF')
    plt.xlabel('IMF Slope')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_bf_gamma.png', dpi=300, bbox_inches='tight')


    res2 = make2dpdfs(data[:,0], data[:,9], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('dAv')
    plt.xlabel('Av')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_av_dav.png', dpi=300, bbox_inches='tight')


    res2 = make2dpdfs(data[:,8], data[:,9], data[:,5])
    plt.figure()
    ax = plt.subplot(111)
    h, xedges, yedges = np.histogram2d(res2[0], res2[1], weights=res2[2])
    plt.imshow(h.T, origin='low', interpolation='None', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=cm, aspect='auto')
    cb = plt.colorbar()
    plt.ylabel('dAv')
    plt.xlabel('BF')
    cb.set_label('Probability')
    plt.savefig(outname+'_joint_bf_dav.png', dpi=300, bbox_inches='tight')
    
    #print(res[1])
    plt.show()
    

main()
