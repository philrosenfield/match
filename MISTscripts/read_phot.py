import numpy as np
from fileio import *

def get_mags(photbase, data='Fakedata', vcuts = [999.0, -999.0]):

    photf_v = np.genfromtxt(get_photof(photbase, data=data), usecols=(0,))
    photf_i = np.genfromtxt(get_photof(photbase, data=data), usecols=(1,))

    badi = np.where( (abs(photf_v) >= 99) | (abs(photf_i) >= 99) )

    photf_v = np.delete(photf_v, badi)
    photf_i = np.delete(photf_i, badi)

    photf_vmi = photf_v - photf_i

    # cut mags:
    photf_v_mask = photf_v[(photf_v < max(vcuts)) & (photf_v > min(vcuts))]
    photf_vmi_mask = photf_vmi[(photf_v < max(vcuts)) & (photf_v > min(vcuts))]

    # Tuple storing the x, y values for plotting the data on a CMD:
    photf_pts = (photf_vmi_mask, photf_v_mask)

    return photf_pts
