import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
import glob
import os

def main(sub='', pref='mod1_*', overwrite=False):
    """
    make a plot of Mbol vs Log Te of tracks to go to MATCH or the
    MATCH tracks themselves

    Parameters
    ----------
    sub : string
        the subdirectory to operate in
    pref : string
        the prefix search string of the track names
    """
    if sub != '':
        os.chdir(sub)

    i = 1  # mod1_ fmt: Mbol Log_Te Nstars ...
    j = 0
    k = 2
    if 'match' in pref:
        i = 2 # match_ fmt: logAge Mass logTe Mbol ...
        j = 3
        k = 1

    #pref += '.dat'
    modfiles = glob.glob(pref)
    for f in modfiles:
        if f.endswith('.png'):
            continue
        if os.path.isfile('%s.png' % f) and not overwrite:
            continue
        data = np.loadtxt(f)
        fig, ax = plt.subplots()
        try:
            l, = ax.scatter(data.T[i], data.T[j], c=np.log10(data.T[k]),
                            cmap=plt.cm.Blues, edgecolor='none')
            cb = plt.colorbar(l)
            cb.set_label('log Nstars')
        except:
            ax.plot(data.T[i], data.T[j], color='k')
        #ax.set_xlim(ax.get_xlim()[::-1])
        #ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylim(13, -14.0)
        ax.set_xlim(5, 3.30)

        ax.set_title(f)
        ax.set_xlabel('Log Te')
        ax.set_ylabel('Mbol')
        plt.savefig('%s.png' % f)
        plt.close()

def find_flag(argv, flag, default=''):
    string = '-%s=' % flag
    try:
        val, = [s.replace(string, '') for s in argv if string in s]
    except:
        val = default
    return val

def find_arg(argv, arg):
    try:
        arg, = [s.replace('-', '') for s in argv if arg in s]
    except:
        arg = None
    return arg

if __name__ == "__main__":
    import sys
    sub = find_flag(sys.argv, 'sub', '')
    pref = find_flag(sys.argv, 'pref', 'mod1_*')
    recursive = find_arg(sys.argv, 'R')

    here = os.getcwd()
    if recursive is not None:
        if sub is not None:
            os.chdir(sub)
        subs = [s for s in os.listdir('.') if os.path.isdir(s)]
        for s in subs:
            print s
            main(sub=s, pref=pref)
            os.chdir('..')
    else:
        main(sub=sub, pref=pref)
