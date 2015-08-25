import errno
import os
import subprocess
import sys
import numpy as np
import pdb

def safe_mkdir(path):
    """Create the directory only if it does not already exist."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def zcombine_hmc_individual(hmcfile, sfhfile=None, zcbfile=None):
    """Get individual .sfh and/or .zcb files for each run in the hybridMC
    output file.

    Note that both `sfhfile` and `zcbfile` are optional and control where
    .sfh and .zcb files are saved. If None (default), then temp files are
    written in the parent directory of `hmcfile` and automatically deleted.

    Parameters
    ----------
    hmcfile : str
        Absolute path to the hybridMC output file.
    sfhfile : str, optional
        Template absolute file path for the .sfh (calcsfh SFH) files.
        Individual files are distinguished by number, so the template path
        must contain a format string "{0:0{1}d}" where the number is to
        appear in each filename. If None (default), then temp files are
        written (of the form "temp_{0:0{1}d}.sfh") in the parent directory
        of `hmcfile` and are automatically deleted.
    zcbfile : str, optional
        Template absolute file path for the .zcb (zcombine SFH) files.
        Individual files are distinguished by number, so the template path
        must contain a format string "{0:0{1}d}" where the number is to
        appear in each filename. No files are written if None (default).

    Returns
    -------
    None

    """
    if sfhfile is None:
        dirname = os.path.dirname(hmcfile)
        sfhfile = os.path.join(dirname, 'temp_{0:0{1}d}.sfh')
        writesfh = False
    else:
        dirname = os.path.dirname(sfhfile)
        safe_mkdir(dirname)
        writesfh = True

    if zcbfile is None:
        writezcb = False
    else:
        dirname = os.path.dirname(zcbfile)
        safe_mkdir(dirname)
        writezcb = True

    if writesfh or writezcb:
        with open(hmcfile, 'r') as f:
            sfhtext_list = f.read().split('\n\n')
        sfhtext_list.pop(-1)  # last element is empty

        nfiles = len(sfhtext_list)
        nzeros = len(str(nfiles))

        for n, sfhtext in enumerate(sfhtext_list):
            sfhtext += '\n\n'

            sfh = sfhfile.format(n+1, nzeros) if sfhfile else sfhfile
            with open(sfh, 'w') as g:
                g.writelines(sfhtext)

            if writezcb:
                zcb = zcbfile.format(n+1, nzeros)
		#pdb.set_trace()
                #cmd = 'zcombine {0:s} -bestonly -jeffreys > {1:s}'.format(sfh, zcb)
		cmd = 'zcombine '+np.str(sfh)+' -bestonly -jeffreys -out=popbox_'+np.str(n)+' > '+np.str(zcb)
                subprocess.call(cmd, shell=True)

            if not writesfh:
                cmd = 'rm {:s}'.format(sfh)
                subprocess.call(cmd, shell=True)

    return None


def main():
    """Set `hmcfile` so that it points to the hybridMC output file, and set
    `sfhfile` and `zcbfile` however you want. The files are numbered, so if
    you wanted myfile_00001.sfh, myfile_00002.sfh, etc., you would use
    'myfile_{0:0{1}d}.sfh' (the actual formatting is done within
    `zcombine_hmc_individual` and it expects this exact format string).
    Note that all paths are absolute.

    """
    galaxy = sys.argv[1]
    name = galaxy.split('_')[0]
    hmcname = sys.argv[2]
    hmcfile = os.getcwd()+'/'+sys.argv[2] # name of HMC output file
    sfhfile = os.getcwd()+'/'+name+'_{0:0{1}d}.sfh'
    zcbfile = os.getcwd()+'/'+name+'_{0:0{1}d}.zcb'

    zcombine_hmc_individual(hmcfile, sfhfile=sfhfile, zcbfile=zcbfile)


if __name__ == '__main__':
    main()
