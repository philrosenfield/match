import numpy as np
from fileio import *

def get_mags(photbase, vcuts = [999.0, -999.0]):

    photf_v = np.genfromtxt(get_photof(photbase), usecols=(0,))
    photf_i = np.genfromtxt(get_photof(photbase), usecols=(1,))
    photf_vmi = photf_v - photf_i

    # cut mags:
    photf_v_mask = photf_v[(photf_v < max(vcuts)) & (photf_v > min(vcuts))]
    photf_vmi_mask = photf_vmi[(photf_v < max(vcuts)) & (photf_v > min(vcuts))]

    # Tuple storing the x, y values for plotting the data on a CMD:
    photf_pts = (photf_vmi_mask, photf_v_mask)

    return photf_pts
