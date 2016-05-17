# What to do to add MIST models to MATCH:
## "Cannonical set"
1. put the new set of tracks in <code>MIST/raw/</code>
1. run bash masses.sh and bash limits.sh to get the model boundaries
   (if you don't know them already)
1. edit the first 20 lines of <code>makemod.cpp</code> with the new track information
    * if the EEPs have changed, edit lines 26-30 of <code>makemod.cpp</code>.
    * If the log Age, M/H, log L, or log Te space has changed, make sure it's the
      the same in <code>../models_common.cpp, ../models.cpp, ../calcbinmod.cpp</code>
1. the following command will run makemod (If this is a *rerun* make sure the data directory is empty first.): 
    
<code>cd ../; make MIST; cd MIST; ./makemod </code> (this makes <code>data/mod1*</code> files)

1. After makemod finishes, cd .. (go to the main match directory) and make
   binaries (this makes <code>data/mod2*</code> files):
      
      <code>./bin/calcbinmod -MIST</code>

1. for faster computing turn all of match's mod* files from ascii to binary (this makes data/mod?b* files):
      
      <code>./bin/ascii2binmod -MIST</code>

## If you are working in a subdirectory (e.g., an overshooting grid) a "non-connonical set”
(I’m not positive this is implemented in <code>MIST/makemod.cpp</code>, but I can add it easily from <code>PARSEC/makemod.cpp</code>)
for each set in the grid (e.g, ov0.30, ov0.40. ov0.50, etc):
1. put the new set of tracks in <code>MIST/raw/subdirname/</code>
1. get model boundaries (may need to edit <code>masses.sh</code> and <code>limits.sh</code> to read subdirname/)
1. same as 3 above — but be careful the entire grid is using the same values here!
1. add <code>-sub=subdirname</code> to every command.
1. add <code>-sub=subdirname</code> to every command.
1. add <code>-sub=subdirname</code> to every command.
1. add the flag <code>-sub=subdirname</code> in <code>run_test.sh</code> to the calcsfh lines.

## For both:
1. Go into the test directory: <code>MIST/test/</code> and <code>bash run_tests.sh</code> to verify all
   is working. (do a <code>tail *scrn</code> after it finishes should see Best fit lines and not
   “Model needs to be created” or anything else)
1. tar ball up the mod*b files for others to use in their match distributions. Or if they want to use the ascii files, mod?_* files.
1. Run <code>calcsfh ... -MIST </code> or if grid: <code>calcsfh ... -MIST -sub=subdirname </code>

# All the details

## Files you should care about/may need to adjust:

### In <code>match2.6/MIST/</code>
The directory structure that matters is <code>raw/ZXXX_YXXX/filename_MXXX.dat</code> and <code>data/mod*</code>
I set up the <code>raw/</code> structure. If you want something else, it's all hardcoded in
<code>makemod.cpp</code> (see below).

inside <code>raw/</code> are the "unprocessed" match tracks.
inside <code>data/</code> are the "processed" match tracks: mod1_... mod2_... etc

<code>./makemod</code> will not overwrite files in <code>data/</code>. If you want to replace the files in
<code>data/</code> you'll need to <code>rm -rf data</code>. Often I will have a safe/ or bkup/ directory
that I move data/ into instead of rm -rf temporarily.

If you want several options within data, (e.g., an overshooting grid
then you'll need subdirectories in raw/.

e.g.,
1. <code>raw/ov0.30/ZXXX_YXXX</code>
1. run <code>./makemod -sub=ov0.30</code>
1. makemod will populate <code>data/ov0.30/mod1...</code>
1. run calcsfh: <code>calcsfh ... -MIST -sub=ov0.30</code>

#### File: makemod.cpp
**Purpose:** convert "unprocessed" tracks for match to models that match reads.

**Useage e.g.:** cd ../; make MIST; cd MIST; ./makemod

**Pay attention to:** You should only need to worry about the first 30 lines.

    const double Zsun = 0.0142; // MIST's Zsun

The Z array below should match directly with the metallicities in the directory structure. I prefer yo keep the directories raw/ZXXX_YXXX so there is no confusion. Here however, one could also add alpha, etc.

    const double Z[] = {0.0000045,0.0000142,0.0000451,0.0001428,0.0002540,0.0004517,0.0008033,0.0014285,0.0025403,0.0045175,0.0080334,0.0142857,0.0254039,0.0451753};
    static const int NZ = sizeof(Z)/sizeof(double);

The directory names in raw or raw/[subdirectory]/ note the max number of characters in the directory name is hard coded to 18, change it if you need to:

    const char FNZ[NZ][18] = {"Z0.0000045_Y0.249","Z0.0000142_Y0.249","Z0.0000451_Y0.249","Z0.0001428_Y0.249","Z0.0002540_Y0.249","Z0.0004517_Y0.250","Z0.0008033_Y0.250","Z0.0014285_Y0.251","Z0.0025403_Y0.253","Z0.0045175_Y0.256","Z0.0080334_Y0.261","Z0.0142857_Y0.270","Z0.0254039_Y0.287","Z0.0451753_Y0.316"};

The array of masses to use. Use masses.sh to grab a unique set from the files in raw/Z*/ otherwise write code to supply this some other way. If there is a track with a mass not listed below in a directory it will be skipped. If there is not a track with a mass listed below in a directory it will be created by interpolation.

    const double M[] = {0.1,0.15,0.2,0.25,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68,1.7,1.72,1.74,1.76,1.78,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2.0,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.22,2.24,2.26,2.28,2.3,2.32,2.34,2.36,2.38,2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,9.0,10.0,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,175,200,225,250,275,300};
    static const int NM = sizeof(M)/sizeof(double);

Calculate this from min Z and max Z of your track set.

    // limits of age and metallicity coverage =  10 * [M/H]
    static const int modelIZmin = -40;
    static const int modelIZmax = 5;

How many columns are in the file besides logTe and Mbol (will need to +=1 for surface Z later)

    // number of values along isochrone (in addition to logTe and Mbol)
    static const int NHRD=4;

What are the range of Mbol and Teff in the track set:

    // range of Mbol and logTeff
    static const double MOD_L0 = -13.8;
    static const double MOD_LF = 13.1;
    static const double MOD_T0 = 3.3;
    static const double MOD_TF = 5.35;

Don't mess with this yet ... I don't know how Andy is treating RGB mass loss for MIST.

    //static const int ML0 = 9; // number of mass loss steps
    static const int ML0 = 0; // number of mass loss steps
    static const double ACC = 3.0; // CMD subsampling

Easiest with an example?

    001-202: PMS-ZAMS
    202-353: ZAMS-IAMS
    353-454: IAMS-TAMS
    454-605: SGB-TRGB
    605-616: TRGB-ZAHB transition
    616-677: ZACHeB-TACHeB

(Hopefully, that is enough information in case the EEP numbers change. 
If not, I’ll add lots more information here later.)
Note: 

    IAMS = 454 (NPT_LOW) 
    TRGB = 605 (NPT_MS)
    616 - 605 = 11 (TRGB-ZAHB transition = NPT_TR)
    677 - 616 = 61 (HB+AGB points)

Currently - if you add AGB/TPAGB it should add only to NPT_HB. I’m not sure how he deals with M_ini > 3 Msun TP-AGB. (no He-flash)!! This means there are two file lengths in <code>raw/ZXXX_YXXX/*dat</code>: 455 (including the header line) NPT_LOW and 678 (not low)

     static const int NPT_LOW = 454; // low-mass tracks points
     static const int NPT_MS = 605; // MS+RGB tracks points
     static const int NPT_TR = 11;   // transition RGB->HB points
     static const int NPT_HB = 61;  // HB+AGB points


Other things in this file (line 30 and beyond!)

    lines 38-40: max filename length, max file lengths.
    line 49,76,77: directory structure setting
    line 56: only reading .dat files (cutting out the .dat part)
    line 65,70: finding the Mass string from the filename (assumes _Mxxx.dat)
    line 85-92: reading in the MIST grid. Will need to be edited when surface composition is added.


### In match2.6/:

#### File: models.cpp

**Purpose:** Used with makemod (perhaps also with calcsfh?) to load models

**Usage e.g.:** [never directly called by user]

**Pay attention to:** Size of grid should be bigger or as big as in <code>MIST/makemod.cpp</code> Specifically look at <code>MOD_L0</code>, <code>MOD_T0</code>. <code>MOD_NL</code>, <code>MOD_NT</code> only need to change if <code>MOD_LF</code> and <code>MOD_TF</code> in makemod.cpp are much larger than initial values of 13.1, 5.35.


    static const double MOD_L0 = -13.8;  // minimum Mbol to be modeled
    static const int MOD_NL = 1400;
    static const double MOD_T0 = 3.30;   // minimum logTeff to be modeled
    static const int MOD_NT = 1250;


#### File: models_common.h

**Purpose:** sets the model grid size in time and metallicity as well as dL, dT

**Usage e.g.:** [never directly called by user]

**Pay attention to:** Primarily make sure the age in Z boundaries are equal to those available in MIST tracks. Down the road, the <code>MOD_dL</code> and <code>MOD_dT</code> may need to be changed depending on the stellar evolution grid. Because it may change the speed of the code, Andy should to be looped in for any changes outside of <code>ITmin</code>, <code>ITmax</code>, <code>IZmin</code>, and <code>IZmax</code>.


    static const int ITmin = 640; // log Age (Gyr) * 100
    static const int ITmax = 1030; // log Age (Gyr) * 100
    static const int dIT = 5;
    static const double Tstep = 0.01;

    static const int IZmin = -40; // [M/H] * 10
    static const int IZmax = 5; // [M/H] * 10
    static const int dIZ = 1;
    static const double Zstep = 0.1;

    static const double MOD_dL = 0.02;
    static const double MOD_dT = 0.002;


#### File: calcbinmod.cpp
**Purpose:** Makes binary (more than one star) models

**Usage e.g.:** <code>./bin/calcbinmod -MIST</code>

**Pay attention to:** Size of the underlying grid should match models.cpp

    // copied from models.cpp
    static const double MOD_L0 = -13.8;  // minimum Mbol to be modeled
    static const int MOD_NL = 1400;
    static const double MOD_T0 = 3.3;   // minimum logTeff to be modeled
    static const int MOD_NT = 1250;


#### File: makemodFromTracks.h
**Purpose:** Make models from tracks...

**Usage e.g.:** [never directly called by user]

**Pay attention to:** Typically nothing. But this file will throw assertion errors when tracks have MOD_L0/F MOD_T0/F values outside those specified in makemod.cpp. One way around this is to put in a print statement, another is to run MIST/limits.sh or have a previous code output the limits of the computed model grid.


