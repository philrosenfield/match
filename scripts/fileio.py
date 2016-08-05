"""Methods to read and write MATCH files"""
from __future__ import print_function
import json
import os
import glob
import sys
import logging

import pandas as pd
import numpy as np

from .config import SCRNEXT

logger = logging.getLogger()

__all__ = ['add_filename_info_to_file', 'add_gates', 'calcsfh_dict',
           'calcsfh_input_parameter', 'fake_param_fmt', 'filename_data',
           'get_files', 'make_matchfake', 'read_binned_sfh',
           'read_calcsfh_param', 'read_fake', 'read_match_cmd',
           'read_ssp_output']


def add_filename_info_to_file(fname, best=False, stats=False):
    """
    add filename info to the data.
    E.g, ssp_imf4.85_bf0.3_dav0.0.dat
    will add two columns, bf, and dav. See filename_data.
    Parameters
    ----------
    fname : str
        name of the file

    ofile : str
        output file or  will write to file substiting fname's .dat with .fdat

    Returns
    -------
    data : np.array
        data with new columns attached

    """
    def getheader(infile):
        """get the length of the header"""
        idx = -1
        with open(infile) as inp:
            while True:
                line = inp.readline()
                try:
                    idx += 1
                    np.array(line.strip().split(), dtype=float)
                    return idx
                except ValueError:
                    pass

    ihead = getheader(fname)
    names = 'Av IMF dmod lage logZ fit sfr'.split()
    df = pd.read_table(fname, names=names, delim_whitespace=True,
                       skiprows=ihead)
    print(fname)
    try:
        ibest, = np.where(df['Av'] == 'Best')[0]
    except ValueError:
        print('Problem in {}'.format(fname))
        print(sys.exc_info()[1])
        raise
    # av, dmod, fit = map(float, [d.replace(',','':.split('=')[1]
    #                             for d in df.iloc[ibest].values
    # d                            if type(d) == str and '=' in d])
    if best:
        # only best fit line
        df = df.iloc[ibest+1].copy(deep=True)
    else:
        # delete best line and below (it's a duplicate entry)
        df = df.iloc[:ibest].copy(deep=True)

    new_stuff = filename_data(fname)
    print('adding columns: {}'.format(new_stuff.keys()))
    for name, val in new_stuff.items():
        df[name] = val

    if stats:
        try:
            cmdfile, = get_files(fname.replace(SCRNEXT, '.out.cmd.stats'))
        except:
            print('{} not found.'.format(cmdfile))
            return df
        data = read_cmd_stats(cmdfile)
    return df

def read_cmd_stats(cmdfile):
    pass


def filename_data(fname, ext=None, skip=2, delimiter='_', exclude='imf'):
    """
    return a dictionary of key and values from a filename.
    E.g, ssp_imf4.85_bf0.3_dav0.0.fdat
    returns bf: 0.3, dav: 0.0
    NB: imf is excluded because it's already included in the file.

    Parameters
    ----------
    fname : str
        filename

    ext : str
        extension (sub string to remove from the tail)

    delimiter : str
        how the keyvals are separated '_' in example above

    skip : int
        skip n items (skip=1 skips ssp in the above example)

    exclude : str
        do not include this key/value in the file (default: 'imf')

    Returns
    -------
    dict of key and values from filename
    """
    import re
    if ext is None:
        ext = '.{}'.format(fname.split('.')[-1])
    keyvals = fname.replace(ext, '').split(delimiter)[skip:]
    d = {}
    for keyval in keyvals:
        kv = re.findall(r'\d+|[a-z]+', keyval)
        neg = ''
        if '-' in keyval:
            neg = '-'
        if kv[0].lower() == exclude.lower():
            continue
        try:
            d[kv[0]] = float(neg + '.'.join(kv[1:]))
        except ValueError:
            # print e
            # print(sys.exc_info()[1])
            pass
    return d


def read_fake(filename):
    """
    Read in the file produced with fake (or fake -full)
    """
    colnames = ['mag1', 'mag2', 'mass', 'Mbol', 'logTe', 'logg', 'logZ', 'CO',
                'Av', 'age']
    try:
        data = np.genfromtxt(filename, names=colnames)
    except ValueError:
        # -full not run, just a 2 column file
        data = np.genfromtxt(filename, names=colnames[:2])
    return data


def read_match_cmd(filename):
    '''read MATCH .cmd file'''
    if not filename.endswith('.cmd'):
        print('Warning: {} might not be a .cmd file'.format(filename))
    names = ['mag', 'color', 'Nobs', 'Nsim', 'diff', 'sig', 'gate']
    cmd = np.genfromtxt(filename, skip_header=4, names=names,
                        invalid_raise=False)
    with open(filename, 'r') as inp:
        header = [next(inp).strip() for _ in range(4)]
    fit = float(header[0].split()[0])
    colors = header[2]
    yfilter = header[-1]
    return cmd, fit, colors, yfilter


def add_gates(ngates):
    """TODO: figure out how to include gates progammatically."""
    if ngates == 0:
        return ''
    else:
        print('Automated in/exclude gates is not implemented!')
    return ''


def calcsfh_input_parameter(zinc=False, power_law_imf=True, **params):
    '''
    Returns a formatted string of the calcsfh input parameter file.
    params is a dictionary to update calcsfh_dict which has most of the same
    keys as are directly needed in the calcsfh input parameter file.

    imf dmod0 dmod1 ddmod av0 av1 dav
    logzmin logzmax dlogz[1]
    bf bad0 bad1
    ncmds[2]
    vstep v-istep fake_sm v-imin v-imax v,i
    vmin vmax v
    imin imax i
    nexclude_gates exclude_gates ninclude_gates include_gates[3]
    ntbins[4]
        ...
    use_bg bg_smooth bg_sample bg_file[5]

    [1] if zinc the following is added:
        logzmin0 logzmin1 logzmax0 logzmax1
    [2] ncmds > 1 not implemented yet
    [3] Not implemented yet.
    [4] The time array is calculated here.
        Supply minimum log (or linear) age as tmin, maximum log age as tmax and
        either the number of time bins as ntbins or the length of a time bin as
        tbin.
    '''

    param_dict = calcsfh_dict()
    param_dict.update(params)

    possible_filters = match_filters()
    ohno = 0
    for filt in [param_dict['v'], param_dict['i']]:
        if filt not in possible_filters['filters']:
            print('{} not in filter list'.format(filt))
            ohno += 1
    assert(ohno == 0), 'Filters need to be in match filter list.'

    # the logZ line changes if using -zinc flag
    zincfmt = '{logzmin:.2f} {logzmax:.2f} {dlogz:.2f}'
    if zinc:
        zincfmt += \
            ' {logzmin0:.2f} {logzmax0:.2f} {logzmin1:.2f} {logzmax1:.2f}\n'
    else:
        zincfmt += '\n'

    # the first line changes if using power_law_imf.
    line0 = ''
    if power_law_imf:
        line0 = '{imf:.2f}'
    line0 += \
        ' {dmod0:.3f} {dmod1:.3f} {ddmod:.3f} {av0:.3f} {av1:.3f} {dav:.3f}\n'

    # Parse the in/exclude gates
    param_dict['exclude_gates'] = add_gates(param_dict['nexclude_gates'])
    param_dict['include_gates'] = add_gates(param_dict['ninclude_gates'])

    # Prepare time bins
    if param_dict['ntbins'] > 0:
        # set tbin size
        dtarr = np.linspace(param_dict['tmin'], param_dict['tmax'],
                            param_dict['ntbins'])
    else:
        # set ntbins
        if len(param_dict['tbreak']) > 0:
            tbreak = np.asarray(param_dict['tbreak'])
            tbin = np.asarray(param_dict['tbin'])
            dtarr = np.array([])
            assert len(param_dict['tbin']) == len(param_dict['tbreak']) + 1, \
                "to use tbreak, tbin must be an array len(tbreak) + 1"
            tmins = np.concatenate([[param_dict['tmin']], tbreak])
            tmaxs = np.concatenate([tbreak, [param_dict['tmax']]])

            for i in range(len(tmaxs)):
                subarr = np.arange(tmins[i], tmaxs[i] + tbin[i], tbin[i])
                dtarr = np.append(dtarr, subarr)
        else:
            dtarr = np.arange(param_dict['tmin'],
                              param_dict['tmax'] + param_dict['tbin'],
                              param_dict['tbin'])
        # linear ages as input converted to log
        if dtarr[0] > 100.:
            dtarr = np.log10(dtarr)

    if param_dict['ntbins'] > 100:
        print('Warning {} time bins'.format(param_dict['ntbins']))

    param_dict['ntbins'] = len(dtarr) - 1

    # Add background information (if a file is supplied)
    # This might not be the correct formatting...
    if param_dict['bg_file'] != '':
        footer = '{use_bg:d} {bg_smooth:d} {bg_sample:d}{bg_file:s}\n'
    else:
        footer = '\n'

    fmt = line0
    fmt += zincfmt
    fmt += '{bf:.2f} {bad0:.6f} {bad1:.6f}\n'
    fmt += '{ncmds:d}\n'
    fmt += '{vstep:.2f} {vistep:.2f} {fake_sm:d} '
    fmt += '{vimin:.2f} {vimax:.2f} {v:s},{i:s}\n'
    fmt += '{vmin:.2f} {vmax:.2f} {v:s}\n'
    fmt += '{imin:.2f} {vmax:.2f} {i:s}\n'
    fmt += '{nexclude_gates:d} {exclude_gates:s} '
    fmt += '{ninclude_gates:d} {include_gates:s} \n'
    #    Metallicity information not yet supported
    fmt += '{ntbins:d}\n'
    fmt += ''.join(['   {:.6f} {:.6f}\n'.format(i, j) for i, j in
                    zip(dtarr[:], dtarr[1:])
                    if np.round(i, 4) != np.round(j, 4)])
    fmt += footer
    return fmt.format(**param_dict)


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

        fout = pref + '{}_{}'.format(filters[i], filters[-1]) + sufx
        np.savetxt(fout, np.column_stack((mag1in, mag2in, mag1diff, mag2diff)),
                   fmt='%.4f')
        logger.info('wrote {}'.format(fout))
    return


def calcsfh_dict():
    '''return default dictionary for calcsfh.'''
    base = os.path.split(__file__)[0]
    inp_par = os.path.join(base, 'templates/calcsfh_input_parameter.json')
    with open(inp_par, 'r') as inp:
        cdict = json.load(inp)
    return cdict


def match_filters():
    '''return dictionary of possible filters in match.'''
    base = os.path.split(__file__)[0]
    with open(os.path.join(base, 'templates/match_filters.json')) as inp:
        fdict = json.load(inp)
    return fdict


def read_ssp_output(filename):
    """
    Read calcsfh -ssp console output.
    """
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
    Av, dmod, fit = \
        np.array(bfline.strip().translate(None, '#Bestfit:Av=dmod').split(','),
                 dtype=float)
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
             for f in files if os.path.isfile(os.path.join(src, f), mad=False)]
    return files


# Needs revision:
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
        d['dmod0'], d['dmod1'], d['ddmod'], d['av0'], d['av1'], d['dav'] \
            = np.array(lines[0].split(), dtype=float)
    except ValueError:
        d['imf'], d['dmod0'], d['dmod1'], d['ddmod'], d['av0'], d['av1'], \
            d['dav'] = np.array(lines[0].split(), dtype=float)
    try:
        d['logzmin'], d['logzmax'], d['dlogz'], d['logzmin0'], d['logzmax0'], \
            d['logzmin1'], d['logzmax1'] = np.array(lines[1].split(),
                                                    dtype=float)
    except ValueError:
        d['logzmin'], d['logzmax'], d['dlogz'] \
            = np.array(lines[1].split(), dtype=float)

    d['bf'], d['bad0'], d['bad1'] = np.array(lines[2].split(), dtype=float)
    d['ncmds'] = int(lines[3].strip())
    vstep, vistep, fake_sm, vimin, vimax, filters = lines[4].strip().split()
    d['v'], d['i'] = filters.split(',')
    d['vstep'], d['vistep'], d['vimin'], d['vimax'] = \
        np.array([vstep, vistep, vimin, vimax], dtype=float)
    d['fake_sm'] = int(fake_sm)
    vmin, vmax, _ = lines[5].strip().split()
    imin, imax, _ = lines[6].strip().split()
    d['vmin'], d['vmax'], d['imin'], d['imax'] \
        = map(float, [vmin, vmax, imin, imax])
    # d['gates'] = lines[7].strip()
    d['ntbins'] = int(lines[8].strip())
    d['to'], d['tf'] = np.array([l.strip().split() for l in lines[9:]
                                 if not l.startswith('-') and
                                 len(l.strip().split()) > 0], dtype=float).T
    # d['footer'] = ''.join([l for l in lines[8:] if l.startswith('-')])
    return d
