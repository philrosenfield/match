import os
import glob
import numpy as np

def get_workdir(data='Fakedata'):

    # Environment path to /n/home12/sgossage/match2.6
    matchdir_envkey = 'MATCH_DIR'
    matchdir = os.environ.get(matchdir_envkey)

    if matchdir != None:

        # append '_rot' to a data directory name to tell this function to access MISTrot (rotating models), rather than MIST.
        if '_rot' == data[-4::]:
            mistdir = os.path.join(matchdir, "MISTrot")
            data = data.split('_rot')[0]
        else:
            mistdir = os.path.join(matchdir, "MISTrot_nodist")

        workdir = os.path.join(mistdir, data)
        #print(workdir)
    
        # returns the path to the work directory...e.g. match2.6/MIST/Hyades:
        return  os.path.join(mistdir, data) 
    
    else:
        print "MATCH_DIR environment variable (path to e.g. ../matchx.xx, where x.xx = version #) not set.)"
        return

def get_photodir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'photometry')

def get_scrndir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'scrns')

def get_outdir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'output')

def get_photof(photbase, data='Fakedata'):

    return glob.glob(os.path.join(get_photodir(data), photbase+'.phot'))[0]

def get_jobids(data, photbase):
    
    # searches in the scrn directory for SLURM jobid subdirectories and extracts the jobids present: 
    return [direc.split('SLURM')[-1] for direc in glob.glob(os.path.join(get_scrndir(data), photbase, 'SLURM*'))]

def readssp(data='Fakedata', SFR='0.001', solution='single', fromfit=False, param='age', photname_base=None, sub=None, incvvc='all', jobidstr=None):

    workdir = get_workdir(data)

    if param=='lage':
        param = 'age'

    # The photometry directory is used to count the number of stars used; the ssp output directory
    # is used to extract the solutions found by MATCH.
    photdir = os.path.join(workdir, "photometry")
    sspoutdir = os.path.join(workdir, "sspcombines")

    # If desired, look in a specified subdirectory for the star num and best fits:
    if sub != None:
        photdir = os.path.join(photdir, sub)
        sspoutdir = os.path.join(sspoutdir, sub)

    # Get the number of stars in each file
    photfile = get_photof(photname_base, data=data)
    starnum = file_len(photfile)

    # Now go to the sspcombine output directories and get ages, etc.
    if jobidstr == None:
        # use all SLURM* subdirs if a jobid isn't specified:
        sspoutdir = os.path.join(sspoutdir, photname_base, 'SLURM*')
    else:
        sspoutdir = os.path.join(sspoutdir, photname_base, 'SLURM{:s}'.format(jobidstr))

    if incvvc != 'all':
        # if told to, just get the vvc specified by incvvc:
        sspfiles = glob.glob(os.path.join(sspoutdir, "*vvcrit{:.1f}*.ssp".format(incvvc)))
    else:
        sspfiles = glob.glob(os.path.join(sspoutdir, "*.ssp"))

    # age solutions (could be adapted to any parameter):
    bestval = np.array([]) #[] 
    # associated uncertainties:
    uncerts = np.array([])

    # Look through all files at the current SFR:
    print(workdir)
    print(sspoutdir)
    print(sspfiles)
    for j, f in enumerate(sspfiles):
        # Open the file, read lines, close file:
        currf = open(f)
        lines = currf.readlines()
        currf.close()
        
        solns = {}
        inblock = False
        for l, line in enumerate(lines):
            if 'Estimate from single solution' in line:
                inblock = True
                m = 1
                # search the block for solutions:
                while inblock:
                    if "{:s} = ".format(param) in lines[l+m]:
                        solns['single'] = [float((lines[l+m].split("+")[0]).split("=")[-1]), float(lines[l+m].split("+/-")[-1].split('(')[0])]
                        break
                    # here's where the end of the block is, so break:
                    elif lines[l+m] == '\n':
                        break
                    # keep searching until solution is found, or break out.
                    else:
                        m += 1

                inblock = False

            if 'Estimates from marginalized solutions' in line:
                inblock = True
                m = 1
                # search the block for solutions:
                while inblock:
                    if ("{:s} = ".format(param) in lines[l+m]) | ("{:s} < ".format(param) in lines[l+m]) | ("{:s} > ".format(param) in lines[l+m]):
                        if fromfit:
                            # the solution 'from fit' is on the line following param = value +/- uncert
                            if lines[l+m+1] != '\n':
                                solns['marginalized'] = [float((lines[l+m+1].split("+")[0]).split("=")[-1]), float(lines[l+m+1].split("+/-")[-1].split('(')[0])]
                            else:
                                break
                        else:
                            #solns['marginalized'] = [float((lines[l+m].split("+")[0]).split("=")[-1]), float(lines[l+m].split("+/-")[-1].split('(')[0])]
                            try:
                                solns['direct'] = [float((lines[l+m].split("+")[0]).split("=")[-1]), float(lines[l+m].split("+/-")[-1].split('(')[0])]
                            except ValueError:
                                try:
                                    solns['direct'] = [float((lines[l+m].split("<")[-1]).split('(')[0]), 0]
                                except ValueError:
                                    try:
                                        solns['direct'] = [float((lines[l+m].split(">")[-1]).split('(')[0]), 0]
                                    except ValueError:
                                        pass
                        break
                    # here's where the end of the block is, so break:
                    elif lines[l+m] == '\n':
                        break
                    # keep searching until solution is found, or break out.
                    else:
                        m += 1

                inblock = False

            if 'Direct marginalized distributions' in line:
                inblock = True
                m = 1
                # search the block for solutions:
                while inblock:
                    if ("{:s} = ".format(param) in lines[l+m]) |  ("{:s} < ".format(param) in lines[l+m]) |  ("{:s} > ".format(param) in lines[l+m]):
                        try:
                            solns['direct'] = [float((lines[l+m].split("+")[0]).split("=")[-1]), float(lines[l+m].split("+/-")[-1].split('(')[0])]
                        except ValueError:
                            try:
                                solns['direct'] = [float((lines[l+m].split("<")[-1]).split('(')[0]), 0]
                            except ValueError:
                                try:
                                    solns['direct'] = [float((lines[l+m].split(">")[-1]).split('(')[0]), 0]
                                except ValueError:
                                    pass
                        break
                    # here's where the end of the block is, so break:
                    elif lines[l+m] == '\n':
                        break
                    # keep searching until solution is found, or break out.
                    else:
                        m += 1

                inblock = False

            if 'Subsampled marginalized distributions' in line:
                inblock = True
                m = 1
                # search the block for solutions:
                while inblock:
                    if "{:s} = ".format(param) in lines[l+m]:
                        if fromfit:
                            # the solution 'from fit' is on the line following param = value +/- uncert
                            if lines[l+m+1] != '\n':
                                solns['subsampled'] = [float((lines[l+m+1].split("+")[0]).split("=")[-1]), float(lines[l+m+1].split("+/-")[-1].split('(')[0])]
                            else:
                                break
                        else:
                            solns['subsampled'] = [float((lines[l+m].split("+")[0]).split("=")[-1]), float(lines[l+m].split("+/-")[-1].split('(')[0])]
                        break
                    # here's where the end of the block is, so break:
                    elif lines[l+m] == lines[-1] or lines[l+m]=='Grid not fully populated\n':
                        break
                    # keep searching until solution is found, or break out.
                    else:
                        m += 1

                inblock = False

        try:
            bestval = np.append(bestval, solns[solution][0])
            uncerts = np.append(uncerts, solns[solution][1])
        except KeyError:
            # solution was not found.
            if solution == 'subsampled':
                try:
                    bestval = np.append(bestval, solns['direct'][0])
                    uncerts = np.append(uncerts, solns['direct'][1])
                except KeyError:
                    print('No solution found!')
                    pass
            else:
                print('No solution found!')
                             

    # switch back to the work directory and return the found solutions and uncertainties, along with star number used
    # in finding the respective solutions:
    starnums = np.array([starnum])
    
    # model resolution
    dT = 0.02
    dZ = 0.02

    if param == 'age':
        res = dT
    elif param == 'logZ':
        res = dZ

    if param == 'age':
        # add half dT to best-fit age, since match reports w/ bottom edge of age bin.
        bestval = bestval + (res/2.0)

    # also add in quadrature the unertainty due to parameter resolution. From e-mail between D. Weisz and A.E. Dolphin:
    #
    # "For reference, a uniform draw from 0 to 1 has a mean of 0.5 and standard deviation of sqrt(1/12) = 0.29.  So, 
    #  one-sigma resolution for any variable sampled over uniform spacing is the step size / sqrt(12)." - AED
    uncerts = np.sqrt(uncerts**2 + (res/np.sqrt(12))**2)

    print(solns)
    print(bestval)
    print(uncerts)

    return starnums, np.array([np.mean(bestval)]), np.array([np.sqrt(np.sum(uncerts**2))/len(uncerts)])

# Count lines in a file:
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
