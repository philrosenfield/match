import os
import glob

def get_workdir(data='Fakedata'):

    # Environment path to /n/home12/sgossage/match2.6
    matchdir_envkey = 'MATCH_DIR'
    matchdir = os.environ.get(matchdir_envkey)

    if matchdir != None:

        mistdir = os.path.join(matchdir, "MIST")
    
        # returns the path to the work directory...e.g. match2.6/MIST/Hyades:
        return  os.path.join(mistdir, data) 
    
    else:
        print "MATCH_DIR environment variable (path to e.g. ../matchx.xx, where x.xx = version #) not set.)"
        return

def get_photodir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'photometry')

def get_photof(photbase):

    return glob.glob(os.path.join(get_photodir('Fakedata'), photbase+'.phot'))[0]
