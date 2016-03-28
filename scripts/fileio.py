from __future__ import print_function
import numpy as np
import os
import re
import glob
import sys
import logging
logger = logging.getLogger()

__all__ = ['ensure_file', 'get_files',  'readfile', 'replace_ext', 'savetxt',
           'match_param_default_dict', 'match_param_fmt', 'process_match_sfh',
           'read_match_cmd', 'calcsfh_dict', 'make_match_param',
           'read_ssp_output', 'read_binned_sfh', 'parse_pipeline']


def read_fake(filename):
    """
    Read in the file produced with fake (or fake -full)
    """
    colnames = ['mag1', 'mag2', 'mass', 'Mbol', 'logTe', 'logg', 'logZ', 'CO', 'Av', 'age']
    try:
        data = np.genfromtxt(filename, names=colnames)
    except ValueError:
        # -full not run, just a 2 column file
        data = np.genfromtxt(filename, names=colnames[:2])
    return data


def read_calcsfh_param(filename):
    """
    Read the calcsfh parameter file into a dictionary.
    NB:
    Exclude gates line is just a string 'gates'
    Anything after the time bins are read in as a string. 'footer'
    """
    lines = open(filename).readlines()
    d = {}
    try:
        d['dmod0'], d['dmod1'], d['ddmod'], d['av0'], d['av1'], d['dav'] = np.array(lines[0].split(), dtype=float)
    except:
        d['imf'], d['dmod0'], d['dmod1'], d['ddmod'], d['av0'], d['av1'], d['dav'] = np.array(lines[0].split(), dtype=float)
    try:
        d['logzmin'], d['logzmax'], d['dlogz'], d['logzmin0'], d['logzmax0'], d['logzmin1'], d['logzmax1'] = np.array(lines[1].split(), dtype=float)
    except:
        d['logzmin'], d['logzmax'], d['dlogz'] = np.array(lines[1].split(), dtype=float)

    d['BF'], d['bad0'], d['bad1'] = np.array(lines[2].split(), dtype=float)
    d['ncmds'] = int(lines[3].strip())
    Vstep, VIstep, fake_sm, VImin, VImax, filters = lines[4].strip().split()
    d['V'], d['I'] = filters.split(',')
    d['Vstep'], d['VIstep'], d['fake_sm'], d['VImin'], d['VImax'] = map(float, [Vstep, VIstep, fake_sm, VImin, VImax])
    Vmin, Vmax, d['V'] = lines[5].strip().split()
    Imin, Imax, d['I'] = lines[6].strip().split()
    d['Vmin'], d['Vmax'], d['Imin'], d['Imax'] = map(float, [Vmin, Vmax, Imin, Imax])
    # fuck it.
    d['gates'] = lines[7].strip()
    d['ntbins'] = int(lines[8].strip())
    d['to'], d['tf'] = np.array([l.strip().split() for l in lines[9:] if not l.startswith('-')], dtype=float).T
    d['footer'] = ''.join([l for l in lines[8:] if l.startswith('-')])
    return d



def read_match_cmd(filename):
    '''
    reads MATCH .cmd file
    '''
    #mc = open(filename, 'r').readlines()
    # I don't know what the 7th column is, so I call it lixo.
    names = ['mag', 'color', 'Nobs', 'Nsim', 'diff', 'sig', 'lixo']
    cmd = np.genfromtxt(filename, skip_header=4, names=names, invalid_raise=False)
    return cmd


def parse_pipeline(filename):
    '''find target and filters from the filename'''
    name = os.path.split(filename)[1].upper()

    # filters are assumed to be F???W
    starts = np.array([m.start() for m in re.finditer('_F', name)])
    starts += 1
    if len(starts) == 1:
        starts = np.append(starts, starts+6)
    filters = [name[s: s+5] for s in starts]

    # the target name is assumed to be before the filters in the filename
    pref = name[:starts[0]-1]
    for t in pref.split('_'):
        if t == 'IR':
            continue
        try:
            # this could be the proposal ID
            int(t)
        except:
            # a mix of str and int should be the target
            target = t
    return target, filters


def match_param_default_dict():
    ''' default params for match param file'''

    dd = {'ddmod': 0.05,
          'dav': 0.05,
          'logzmin': -2.3,
          'logzmax': 0.1,
          'dlogz': 0.1,
          'logzmin0': -2.3,
          'logzmax0': -1.0,
          'logzmin1': -1.3,
          'logzmax1': -0.1,
          'BF': 0.35,
          'bad0': 1e-6,
          'bad1': 1e-6,
          'ncmds': 1,
          'Vstep': 0.1,
          'V-Istep': 0.05,
          'fake_sm': 5,
          'nexclude_gates': 0,
          'exclude_gates': '',
          'ninclude_gates': 0,
          'include_gates': ''}

    therest = ['imf', 'dmod1', 'dmod2', 'av1', 'av2', 'V-Imin', 'V-Imax', 'V',
               'I', 'Vmin', 'Vmax', 'Imin', 'Imax']
    for key in therest:
        dd[key] = np.nan
    return dd


def fake_param_fmt(power_law_imf=False, fake=True):
    """
    I donno... without Ntbins and age binning I think this is stupid.
    Does not allow for -diskav or -mag or anything besides trying to reproduce
    the cmds used in calcsfh.

    IMF (m-M)o Av Zspread BF dmag_min
    Vstep V-Istep fake_sm V-Imin V-Imax V,I  (per CMD)
    Vmin Vmax V                              (per filter)
    Imin Imax I                              (per filter)

    NOT INCLUDED:
    Ntbins
      To Tf SFR Z ^see below^ (for each time bin)

    FROM README :
        dmag_min is the minimum good output-input magnitude in the
        fake star results; usually -0.75 (the recovered star is twice as bright
        as the input star) is a good value (-1.50 would match identically the
        value used by calcsfh).
    """
    return ['%(dmod).3f %(av).3f %(Zspread).3f %(BF).3f %(dmag_min).3f',
            '%(Vstep).2f %(VIstep).2f %(fake_sm)i %(VImin).2f %(VImax).2f %(V)s,%(I)s',
            '%(Vmin).2f %(Vmax).2f %(V)s',
            '%(Imin).2f %(Imax).2f %(I)s']

def match_param_fmt(set_z=False, zinc=True):
    '''
    calcsfh parameter format, set up for dan's runs and parsec M<12.
    NOTE exclude and include gates are strings and must have a space at
    their beginning.
    '''

    return '''%(imf)s %(dmod1).3f %(dmod2).3f %(ddmod).3f %(av1).3f %(av2).3f %(dav).3f
%(logzmin).2f %(logzmax).2f %(dlogz).2f %(logzmin0).2f %(logzmax0).2f %(logzmin1).2f %(logzmax1).2f
%(BF).2f %(bad0).6f %(bad1).6f
%(ncmds)i
%(Vstep).2f %(V-Istep).2f %(fake_sm)i %(V-Imin).2f %(V-Imax).2f %(V)s,%(I)s
%(Vmin).2f %(Vmax).2f %(V)s
%(Imin).2f %(Imax).2f %(I)s
%(nexclude_gates)i%(exclude_gates)s %(ninclude_gates)i%(include_gates)s
50
6.60 6.70
6.70 6.80
6.80 6.90
6.90 7.00
7.00 7.10
7.10 7.20
7.20 7.30
7.30 7.40
7.40 7.50
7.50 7.60
7.60 7.70
7.70 7.80
7.80 7.90
7.90 8.00
8.00 8.10
8.10 8.20
8.20 8.30
8.30 8.40
8.40 8.50
8.50 8.60
8.60 8.70
8.70 8.75
8.75 8.80
8.80 8.85
8.85 8.90
8.90 8.95
8.95 9.00
9.00 9.05
9.05 9.10
9.10 9.15
9.15 9.20
9.20 9.25
9.25 9.30
9.30 9.35
9.35 9.40
9.40 9.45
9.45 9.50
9.50 9.55
9.55 9.60
9.60 9.65
9.65 9.70
9.70 9.75
9.75 9.80
9.80 9.85
9.85 9.90
9.90 9.95
9.95 10.00
10.00 10.05
10.05 10.10
10.10 10.15
-1 5 -1bg.dat
-1  1 -1
'''


def process_match_sfh(sfhfile, outfile='processed_sfh.out', sarah_sim=False,
                      zdisp=0.):
    '''
    turn a match sfh output file into a sfr-z table for trilegal.

    todo: add possibility for z-dispersion.
    '''

    fmt = '%.6g %.6g %.4g %s \n'

    data = read_binned_sfh(sfhfile)
    sfr = data['sfr']
    # Trilegal only needs populated time bins, not fixed age array
    inds, = np.nonzero(sfr > 0)
    sfr = sfr[inds]
    to = data['lagei'][inds]
    tf = data['lagef'][inds]
    dlogz = data['mh'][inds]
    half_bin = np.diff(dlogz[0: 2])[0] / 2.
    if zdisp > 0:
        zdisp = '%.4g' % (0.02 * 10 ** zdisp)
    else:
        zdisp = ''

    # correct age for trilegal isochrones.
    # with PARSEC V1.1 and V1.2 no need!
    #tf[tf == 10.15] = 10.13

    with open(outfile, 'w') as out:
        for i in range(len(to)):
            if sarah_sim is True:
                z1 = dlogz[i] - half_bin
                z2 = dlogz[i] + half_bin
                sfr[i] /= 2.
            else:
                sfr[i] *= 1e3  # sfr is normalized in trilegal
                # MATCH conversion:
                z1 = 0.02 * 10 ** (dlogz[i] - half_bin)
                z2 = 0.02 * 10 ** (dlogz[i] + half_bin)
            age1a = 1.0 * 10 ** to[i]
            age1p = 1.0 * 10 ** (to[i] + 0.0001)
            age2a = 1.0 * 10 ** tf[i]
            age2p = 1.0 * 10 ** (tf[i] + 0.0001)

            out.write(fmt % (age1a, 0.0, z1, zdisp))
            out.write(fmt % (age1p, sfr[i], z1, zdisp))
            out.write(fmt % (age2a, sfr[i], z2, zdisp))
            out.write(fmt % (age2p, 0.0, z2, zdisp))
            out.write(fmt % (age1a, 0.0, z2, zdisp))
            out.write(fmt % (age1p, sfr[i], z2, zdisp))
            out.write(fmt % (age2a, sfr[i], z1, zdisp))
            out.write(fmt % (age2p, 0.0, z1, zdisp))

    print('wrote', outfile)
    return outfile


def make_matchfake(fname):
    """
    make four-column Vin, Iin, Vdiff, Idiff artificial stars

    made to work with pipeline fake.fits files with 3 filters, should work
    for two but not tested
    assumes _F*W_F*W_ etc in the file name label the filters.
    """
    try:
        tbl = Table.read(fname, format='fits')
    except:
        logger.error('problem with {}'.format(fname))
        return
    filters = [f for f in fname.split('_') if f.startswith('F')]
    pref = fname.split('F')[0]
    sufx = fname.split('W')[-1].replace('fits', 'matchfake')
    for i in range(len(filters)-1):
        mag1in_col = 'MAG{}IN'.format(i+1)
        mag2in_col = 'MAG{}IN'.format(len(filters))

        if mag1in_col == mag2in_col:
            continue

        mag1out_col = 'MAG{}OUT'.format(i+1)
        mag2out_col = 'MAG{}OUT'.format(len(filters))

        try:
            mag1in = tbl[mag1in_col]
            mag2in = tbl[mag2in_col]
            mag1diff = tbl[mag1in_col] - tbl[mag1out_col]
            mag2diff = tbl[mag2in_col] - tbl[mag2out_col]
        except:
            logger.error('problem with column formats in {}'.format(fname))
            return

        fout = pref + '{}_{}'.format(filters[i],filters[-1]) + sufx
        np.savetxt(fout, np.column_stack((mag1in, mag2in, mag1diff, mag2diff)),
                   fmt='%.4f')
        logger.info('wrote {}'.format(fout))
    return


def calcsfh_dict():
    '''
    default dictionary for calcsfh.
    '''
    return {'dmod': 10.,
            'Av': 0.,
            'filter1': None,
            'filter2': None,
            'bright1': None,
            'faint1': None,
            'bright2': None,
            'faint2': None,
            'color': None,
            'mag': None,
            'dmod2': None,
            'colmin': None,
            'colmax': None,
            'Av2': None,
            'imf': 1.30,
            'ddmod': 0.050,
            'dAv': 0.050,
            'logzmin': -2.3,
            'logzmax': 0.1,
            'dlogz': 0.1,
            'zinc': True,
            'bf': 0.35,
            'bad0': 1e-6,
            'bad1': 1e-6,
            'Ncmds': 1,
            'dmag': 0.1,
            'dcol': 0.05,
            'fake_sm': 5,
            'nexclude_gates': 0,
            'exclude_poly': None,
            'ncombine_gates': 0,
            'combine_poly': None,
            'ntbins': 0,
            'dobg': -1,
            'bg_hess': .0,   # neg if it's a .CMD, else it's same fmt as match_phot
            'smooth': 1,
            'ilogzmin': -2.3,
            'ilogzmax': -1.3,
            'flogzmin': -1.9,
            'flogzmax': -1.1,
            'match_bg': ''}


def read_ssp_output(filename):
    """
    Read calcsfh -ssp console output.
    """
    if filename.endswith('fdat') or filename.endswith('fscrn'):
        """file with added columns from dAv, COV, etc."""
        data = readfile(filename, commented_header=True)
    else:
        skip_header = 10
        skip_footer = 1
        colnames = ['Av', 'IMF', 'dmod', 'lage', 'logZ', 'fit', 'sfr', 'bg1',
                    'bg2']
        try:
            data = np.genfromtxt(filename, skip_header=skip_header,
                                 skip_footer=skip_footer, names=colnames)
        except:
            # no bg?
            try:
                data = np.genfromtxt(filename, skip_header=skip_header,
                                     skip_footer=skip_footer, names=colnames[:-2])
            except:
                print('can not load file: {}'.format(filename))
                return np.array([]), np.nan, np.nan, np.nan

    bfline, = os.popen('tail -n 1 {}'.format(filename)).readlines()
    Av, dmod, fit = map(float, bfline.strip().translate(None, '#Bestfit:Av=dmod').split(','))
    return data, Av, dmod, fit

def read_binned_sfh(filename, hmc_file=None):
    '''
    reads the file created using zcombine or HybridMC from match
    into a np.recarray.

    If hmc_file is not None, will overwrite error columns in file name with
    those from the hmc_file.

    NOTE
    calls genfromtext up to 3 times. There may be a better way to figure out
    how many background lines/what if there is a header... (it's a small file)
    '''
    dtype = [('lagei', '<f8'),
             ('lagef', '<f8'),
             ('dmod', '<f8'),
             ('sfr', '<f8'),
             ('sfr_errp', '<f8'),
             ('sfr_errm', '<f8'),
             ('mh', '<f8'),
             ('mh_errp', '<f8'),
             ('mh_errm', '<f8'),
             ('mh_disp', '<f8'),
             ('mh_disp_errp', '<f8'),
             ('mh_disp_errm', '<f8'),
             ('csfr', '<f8'),
             ('csfr_errp', '<f8'),
             ('csfr_errm', '<f8')]
    def _loaddata(filename, dtype):
        try:
            data = np.genfromtxt(filename, dtype=dtype)
        except ValueError:
            try:
                # zcmerge file
                data = np.genfromtxt(filename, dtype=dtype, skip_header=1)
            except ValueError:
                # background footers
                try:
                    data = np.genfromtxt(filename, dtype=dtype, skip_header=6,
                                         skip_footer=1)
                except ValueError:
                    data = np.genfromtxt(filename, dtype=dtype, skip_header=6,
                                         skip_footer=2)
        return data

    data = _loaddata(filename, dtype)
    if hmc_file is not None:
        # overwrite errors
        hmc_data = _loaddata(hmc_file, dtype)
        for attr in hmc_data.dtype.names:
            if 'err' in attr:
                data[attr] = hmc_data[attr]

    return data.view(np.recarray)


def savetxt(filename, data, fmt='%.4f', header=None, overwrite=False,
            loud=False):
    '''
    np.savetxt wrapper that adds header. Some versions of savetxt
    already allow this...
    '''
    if overwrite is True or not os.path.isfile(filename):
        with open(filename, 'w') as f:
            if header is not None:
                if not header.endswith('\n'):
                    header += '\n'
                f.write(header)
            np.savetxt(f, data, fmt=fmt)
        if loud:
            print('wrote', filename)
    else:
        logger.error('%s exists, not overwriting' % filename)
    return


def readfile(filename, col_key_line=0, comment_char='#', string_column=None,
             string_length=16, only_keys=None, delimiter=' ', commented_header=False):
    '''
    reads a file as a np array, uses the comment char and col_key_line
    to get the name of the columns.
    '''
    if commented_header:
        with open(filename, 'r') as f:
            lines = f.readlines()
        header = [l for l in lines[:-10] if l.startswith(comment_char)]
        col_keys = header[-1].replace(comment_char, '').strip().translate(None, '/[]-').split()
    elif col_key_line == 0:
        with open(filename, 'r') as f:
            line = f.readline()
        col_keys = line.replace(comment_char, '').strip().translate(None, '/[]-').split()
    else:
        with open(filename, 'r') as f:
            lines = f.readlines()
        col_keys = lines[col_key_line].replace(comment_char, '').strip().translate(None, '/[]').split()
    usecols = range(len(col_keys))

    if only_keys is not None:
        only_keys = [o for o in only_keys if o in col_keys]
        usecols = list(np.sort([col_keys.index(i) for i in only_keys]))
        col_keys = list(np.array(col_keys)[usecols])

    dtype = [(c, '<f8') for c in col_keys]
    if string_column is not None:
        if type(string_column) is list:
            for s in string_column:
                dtype[s] = (col_keys[s], '|S%i' % string_length)
        else:
            dtype[string_column] = (col_keys[string_column], '|S%i' % string_length)
    data = np.genfromtxt(filename, dtype=dtype, invalid_raise=False,
                         usecols=usecols, skip_header=col_key_line + 1)
    return data


def replace_ext(filename, ext):
    '''
    input
    filename string with .ext
    new_ext replace ext with new ext
    eg:
    $ replace_ext('data.02.SSS.v4.dat', '.log')
    data.02.SSS.v4.log
    '''
    return split_on_extention(filename)[0] + ext


def split_on_extention(filename):
    '''
    split the filename from its extension
    '''
    return '.'.join(filename.split('.')[:-1]), filename.split('.')[-1]


def ensure_file(f, mad=True):
    '''
    input
    f (string): if f is not a file will print "no file"
    optional
    mad (bool)[True]: if mad is True, will exit program.
    '''
    test = os.path.isfile(f)
    if test is False:
        logger.warning('{} not found'.format(f))
        if mad:
            sys.exit()
    return test


def get_files(src, search_string):
    '''
    returns a list of files, similar to ls src/search_string
    '''
    if not src.endswith('/'):
        src += '/'
    try:
        files = glob.glob1(src, search_string)
    except IndexError:
        logger.error('Can''t find %s in %s' % (search_string, src))
        sys.exit(2)
    files = [os.path.join(src, f)
             for f in files if ensure_file(os.path.join(src, f), mad=False)]
    return files
