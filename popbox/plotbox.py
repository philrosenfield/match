import numpy as np
import pylab as plt
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1 import host_subplot
import brewer2mpl
from matplotlib.font_manager import FontProperties
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import sys
import pdb

plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica')


galaxy_name = sys.argv[1]


''' plotbox <zcombine file> <popbox file> '''


# plots the cumulative SFH
def plotcum(time, csfr, color, linestyle):
    plt.plot(10**np.append(time, 10.15)/10**9, np.append(csfr, 0), lw=3, color=color, ls=linestyle)
    plt.xlim(14.9,0.9)
    plt.ylim(-0.09, 1.09)


# plot uncertainties on the cumulative SFH
def plotfill(time, csfr, csfrhi, csfrlo, color, alpha=1.0):
    plt.fill_between(10**np.append(time, 10.15)/10**9, np.append(csfr, 0), np.append(csfr+csfrhi, 0), alpha=alpha, color=color)
    plt.fill_between(10**np.append(time, 10.15)/10**9, np.append(csfr, 0), np.append(csfr-csfrlo, 0), alpha=alpha, color=color)

# plots the abolute SFH
def plot_sfh(time0, time1, sfr, color, lw=2):
    for i in range(len(time0)-1):
        plt.plot((10**time1[i]/1e9, 10**time1[i]/1e9), (sfr[i], sfr[i+1]), color=color,lw=2)
        plt.plot((10**time1[i]/1e9, 10**time1[i+1]/1e9), (sfr[i+1], sfr[i+1]), color=color,lw=2)


# assumes all the header information is present for zcombine file
# assumes 2 background files
sfh = np.genfromtxt(galaxy_name+'.sfh', skip_header=6, skip_footer=2)

# reads in the population box
lines = [np.float_(line.strip().split()) for line in open(galaxy_name+'.popbox')]
a= np.zeros((len(lines[2:]), len(lines[2:][0][2:])))    
b=()
t0_1, t1_1 = (), ()        
for i in range(len(a)):
    a[i] = lines[2:][i][2:]
    b = np.append((b), (lines[2:][i][2:]))
    t0_1 = np.append(t0_1, (lines[2:][i][0]))
    t1_1 = np.append(t1_1, (lines[2:][i][1]))


cmap = plt.get_cmap('Reds')
cmap.set_gamma(0.5)
plt.close('all')

metal = lines[1] + 0.05 # assumes a metallcity bin width of 0.1 dex

fig=plt.figure()

plt.subplot(331)
ax = fig.add_subplot(331,frameon=False)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='white', top='off', bottom='off', left='off', right='off')
ax1 = host_subplot(3,3,1)
ax01 = ax1.twin()

# hard coded age-redshift relationship -- should update using astropy
ax01.set_xticks([13.1389, 12.060, 10.263, 7.6327, 4.9217, 1.2513, 0])
ax01.set_xticklabels(["10", "5", "2", "1", "0.5", "0.1", "0"], fontsize=10)


ax01.set_xlabel('Redshift (z)', fontsize=12)
ax01.tick_params(which='both', length=6)
plt.setp(ax01.get_yticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
plotcum(sfh[:,0], sfh[:,12], 'DarkRed', '-')
plt.plot((14.1, 1.), (0,1), '0.7', ls='--')
plt.ylabel('Cumulative SFH')
plt.xticks(np.arange(0,15,3))

ax3 = plt.subplot(334)
scalefactor = 10**np.around(np.log10(sfh[:,3].mean()),0)
plot_sfh(sfh[:,0], sfh[:,1], sfh[:,3]/scalefactor, 'DarkRed')
plt.xlim(14.9, 0.9)
plt.ylim(-1.5,1.05 * sfh[:,3].max()/scalefactor)
plt.xticks(np.arange(0,15,3))
plt.setp(ax3.get_xticklabels(), visible=False)
plt.ylabel('SFR (10$^{{{0}}}$ M$_{{\odot}}$ /yr)'.format(np.int(np.log10(scalefactor))))

plt.subplot(337)
#plt.ylim(0.1, -2.3)
plt.ylim(metal.min()-0.1, metal.max())
plt.xlim(14.9, 0.9)


# plot the population boxes by shading indivudal rectangles
from matplotlib.patches import Rectangle

k=0
print('Plotting the population box, may take a couple of minutes...')
for i in range(len(t1_1)):
    for j in range(len(metal)):
        if a[::-1].flatten()[k] > 0:
            rect = Rectangle((10**t1_1[::-1][i]/1e9, metal[j]), (10**t0_1[::-1][i]-10**t1_1[::-1][i])/1e9, 0.1, facecolor=cmap((a[::-1].flatten()/a.max())[k]*5), angle=0, edgecolor='0.7',lw=0.1)
            plt.gca().add_patch(rect)
        else:
            rect1 = Rectangle((10**t1_1[::-1][i]/1e9, metal[j]), (10**t0_1[::-1][i]-10**t1_1[::-1][i])/1e9, 0.1, facecolor='white', angle=0, edgecolor='0.7', lw=0.1)
            plt.gca().add_patch(rect1)
        k+=1

plt.xticks(np.arange(0,15,3))
        #plt.colorbar(mappable=cm)
plt.subplots_adjust(hspace=0)
plt.ylabel('[M/H]')
plt.xlabel('Lookback Time (Gyr)')

plt.savefig(galaxy_name+'_popbox.png', dpi=300, bbox_inches='tight')
plt.show()




