import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import sys

depwarn = """
Hess_plot is deprecated. To make hess diagram plots (ala pgpro) of the .cmd
output files from MATCH please use scripts/cmd.py from the command line
(with -h for options).

The main plotting function is in scripts/graphics/match_plot.py
"""

import warnings
warnings.warn(depwarn, DeprecationWarning)

cmdfile_name = sys.argv[1]
plot_name = sys.argv[2]
figsize=(8, 8)
subsize = 0.38
subx1, suby1 = 0.08, 0.08
subx2, suby2 = 0.58, 0.58
interpolation = 'hanning'
xlabel, ylabel = 'F606W - F814W', 'F606W'

def makecmap(arr):
    x = np.linspace(arr.min(), arr.max(), 50)
    steps = (x - x.min()) / (x.max() - x.min())
    numerator = np.arcsinh(x) - np.arcsinh(x.min())
    denominator = np.arcsinh(x.max()) - np.arcsinh(x.min())
    mapping = 1 - numerator / denominator
    cdict = {}
    for key in ('red', 'green', 'blue'):
        cdict[key] = np.vstack([steps, mapping, mapping]).transpose()
    return cl.LinearSegmentedColormap('new_colormap', cdict, N=1024)

def plot_HessD(fig, arr, subx1, suby2, subsize, extent, interpolation,
               title, xlabel, ylabel):
    cm = makecmap(arr)
    #cm.set_gamma(0.2)
    ax = fig.add_axes([subx1, suby2, subsize, subsize])
    ax.imshow(arr, cmap=cm, aspect='auto', extent=extent,
              interpolation=interpolation)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xticklabels(ax.get_xticks(), size=11)
    ax.set_yticklabels(ax.get_yticks(), size=11)
    #cax = fig.add_axes([subx1+0.13, suby2+0.35, subsize-0.15, subsize-0.365])
    cax = fig.add_axes([subx1+0.03, suby2+0.35,
                        subsize-0.20, subsize-0.365])
    cb_arr = np.linspace(arr.min(), arr.max(), 1e2).reshape(1, 1e2)
    cax_limits = [np.floor(arr.min()), np.ceil(arr.max()), 0, 1]
    cax.imshow(cb_arr, cmap=cm, aspect='auto', extent=cax_limits, interpolation='nearest')
    cax.plot([arr.min(), arr.min()], [0, 1], 'k-')
    cax.plot([arr.max(), arr.max()], [0, 1], 'w-')
    cax.axis(cax_limits)
    cax.yaxis.set_visible(False)
    cax.xaxis.tick_bottom()
    cax.xaxis.set_major_locator(plt.LinearLocator(5))
    cax.set_xticklabels(cax.get_xticks(), size=8)
    #cax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    return ax, cax

cmdfile = open(cmdfile_name, 'r')
line = cmdfile.readline()
shape = [int(i) for i in cmdfile.readline().split()]
line_col, line_mag = cmdfile.readline(), cmdfile.readline()
magbins, colbins = np.zeros(shape[0]), np.zeros(shape[1])
obs_arr, mod_arr = np.zeros(shape), np.zeros(shape)
res_arr, sig_arr = np.zeros(shape), np.zeros(shape)
for i in range(shape[0]):
    for j in range(shape[1]):
        line = cmdfile.readline().split()
        mag, col, obs, mod, res, sig = [float(n) for n in line[:6]]
        colbins[j] = col
        obs_arr[i, j], mod_arr[i, j] = obs, mod
        res_arr[i, j], sig_arr[i, j] = res, sig
    magbins[i] = mag
cmdfile.close()

xlabel = ' - '.join(line_col.replace('WFC','F').split('-'))
ylabel = line_mag.replace('WFC','F')

fig = plt.figure(figsize=figsize)

extent = [colbins[0], colbins[-1], magbins[-1], magbins[0]]
ax_obs, cax_obs = plot_HessD(fig, obs_arr, subx1, suby2, subsize, extent, interpolation, '(a) Observed CMD', xlabel, ylabel)

ax_mod, cax_mod = plot_HessD(fig, mod_arr, subx2, suby2, subsize, extent, interpolation, '(b) Modeled CMD', xlabel, ylabel)

ax_res, cax_res = plot_HessD(fig, res_arr, subx1, suby1, subsize, extent, interpolation, '(c) Residual (Obs. - Mod.)', xlabel, ylabel)

ax_sig, cax_sig = plot_HessD(fig, sig_arr, subx2, suby1, subsize, extent, interpolation, '(d) Residual Significance', xlabel, ylabel)

fig.savefig(plot_name, dpi=300, bbox_inches='tight')
