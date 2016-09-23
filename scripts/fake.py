from .fileio import read_calcsfh_param, fake_param_fmt
from .sfh import SFH


def make_fakeparam(param, sfhfile):
    """
    Convert a calcsfh solution and parameter file into a fake parameter file
    """
    outfile = param.replace('param', 'fparam')
    pdict = read_calcsfh_param(param)
    sfh = SFH(sfhfile)

    pdict['dmod'] = sfh.dmod
    pdict['av'] = sfh.Av
    pdict['Zspread'] = sfh.dlogZ
    pdict['dmag_min'] = -1.50
    ntbins = len(sfh.data.sfr)

    line = '\n'.join([f % pdict for f in fake_param_fmt()])
    line += '\n{}\n'.format(ntbins)
    line += '\n'.join(['{:.2f} {:.2f} {:e} {:.2f}'.format(sfh.data.lagei[i],
                                                          sfh.data.lagef[i],
                                                          sfh.data.sfr[i],
                                                          sfh.data.mh[i])
                       for i in range(ntbins)])
    with open(outfile, 'w') as f:
        f.write(line)
