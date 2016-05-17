const double Zsun = 0.0142;
const double Z[] = {0.0000045,0.0000142,0.0000451,0.0001428,0.0002540,0.0004517,0.0008033,0.0014285,0.0025403,0.0045175,0.0080334,0.0142857,0.0254039,0.0451753};
static const int NZ = sizeof(Z)/sizeof(double);
const char FNZ[NZ][18] = {"Z0.0000045_Y0.249","Z0.0000142_Y0.249","Z0.0000451_Y0.249","Z0.0001428_Y0.249","Z0.0002540_Y0.249","Z0.0004517_Y0.250","Z0.0008033_Y0.250","Z0.0014285_Y0.251","Z0.0025403_Y0.253","Z0.0045175_Y0.256","Z0.0080334_Y0.261","Z0.0142857_Y0.270","Z0.0254039_Y0.287","Z0.0451753_Y0.316"};

const double M[] = {0.1,0.15,0.2,0.25,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.0,1.02,1.04,1.06,1.08,1.1,1.12,1.14,1.16,1.18,1.2,1.22,1.24,1.26,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54,1.56,1.58,1.6,1.62,1.64,1.66,1.68,1.7,1.72,1.74,1.76,1.78,1.8,1.82,1.84,1.86,1.88,1.9,1.92,1.94,1.96,1.98,2.0,2.02,2.04,2.06,2.08,2.1,2.12,2.14,2.16,2.18,2.2,2.22,2.24,2.26,2.28,2.3,2.32,2.34,2.36,2.38,2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.0,9.0,10.0,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,32,34,36,38,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,175,200,225,250,275,300};
static const int NM = sizeof(M)/sizeof(double);

// limits of age and metallicity coverage =  10 * [M/H]
static const int modelIZmin = -40;
static const int modelIZmax = 5;

// number of values along isochrone (in addition to logTe and Mbol)
static const int NHRD=4;

// range of Mbol and logTeff
static const double MOD_L0 = -13.8;
static const double MOD_LF = 13.1;
static const double MOD_T0 = 3.3;
static const double MOD_TF = 5.35;

//static const int ML0 = 9; // number of mass loss steps
static const int ML0 = 0; // number of mass loss steps
static const double ACC = 3.0; // CMD subsampling

static const int NPT_LOW = 454; // low-mass tracks points
static const int NPT_MS = 605; // MS+RGB tracks points
static const int NPT_TR = 11;   // transition RGB->HB points
static const int NPT_HB = 61;  // HB+AGB points

inline bool isDefined(const double*d) {return (d[2]!=0);}
inline void setUndefined(double*d) {d[2]=0;}

#include "../makemodFromTracks.h"

void init_mods(void) {
   FILE *f;
   char fn[161],str[161];
   char fn_ms[1000][161];
   double data[NDATA];

   setup();

   for (int i=0;i<NZ;i++) {
      NM_HB = NM;
      allocmods(i);
      // scan directory
      struct dirent *ep;
      sprintf(fn,"%s/MIST/raw/%s",BASEDIR,FNZ[i]);
      int nfn_ms=0;
      DIR *dp = opendir(fn);
      if (dp==0) ferr(fn);
      while ((ep=readdir(dp))!=0) {
         strcpy(str,ep->d_name);
         if (!strcmp(str,".") || !strcmp(str,".."));
         else if (!strcmp(str+strlen(str)-4,".dat")) strcpy(fn_ms[nfn_ms++],str);
         else {
            fprintf(stderr,"File name not understood\n%s\n",str);
	 }
      }
      printf("Z=%g: %d models\n",Z[i],nfn_ms);
      closedir(dp);
      // read in MS models
      for (int jj=0;jj<nfn_ms;jj++) {
         char *ptr=strstr(fn_ms[jj],"_M");
         if (ptr==0) {
            fprintf(stderr,"Cannot parse filename\n%s\n",fn_ms[jj]);
            endrun(-1);
         }
         double m0 = atof(ptr+2);
         int j=0;
         for (j=0;j<NM && fabs(m0-M[j])>0.0001;j++);
         if (j==NM) fprintf(stderr,"Unmatched file %s; skipping\n",fn_ms[jj]);
         else {
            int nline=0;
            if (modsubdirname[0]) sprintf(fn,"%s/MIST/raw%s/%s/%s",BASEDIR,modsubdirname,FNZ[i],fn_ms[jj]);
            else sprintf(fn,"%s/MIST/raw/%s/%s",BASEDIR,FNZ[i],fn_ms[jj]);
            if ((f=fopen(fn,"r"))==NULL) ferr(fn);
            fgets(str,161,f);
            if (str[0]!='#') {
               fprintf(stderr,"Illegal format: no header\n%s\n",fn);
               endrun(-1);
            }
            //printf("%d %s %d\n",jj,fn_ms[jj],NDATA);
            while (fscanf(f,"%lf %lf %lf %lf %lf %lf",data,data+1,data+2,data+3,data+4,data+5)==NDATA) {
               for (int k=0;k<NML;k++) {
                  grid[i][j][k][nline][0] = data[0]; // log10(age)
                  grid[i][j][k][nline][1] = log(data[1]); // ln(m) -- will be replaced by initial mass later
                  grid[i][j][k][nline][2] = data[2]; // log10(Teff)
                  grid[i][j][k][nline][3] = data[3]; // Mbol
                  grid[i][j][k][nline][4] = data[4]; // log10(g)
                  grid[i][j][k][nline][5] = data[5]; // C/O
                  }
               nline++;
               fgets(str,161,f);
            }
            fclose(f);
            if (nline!=NPT && nline!=NPT_LOW) {
               fprintf(stderr,"Only %d points in file\n%s\n",nline,fn);
            }
         }
      }

      // separate HB models from rest
      NM_HB=0;
      if (NML>1) {
	 for (int j=0;j<NM;j++) if (isDefined(grid[i][j][0][NPT_MS+NPT_TR])) { // If there are enough points such that HeB occurs
	    if (exp(grid[i][j][0][NPT_MS][1])>=M[j]*0.98) j=NM;
	    else {
	       M_HB[NM_HB] = exp(grid[i][j][0][NPT_MS][1]);  // TRGB mass
	       for (int k=0;k<NPT_HB;k++) {
		  dointerp(1.0,grid[i][j][0][k+NPT_MS+NPT_TR],grid[i][j][0][k+NPT_MS+NPT_TR],gridhb[NM_HB][k]);
		  gridhb[NM_HB][k][0] = log10( pow(10,grid[i][j][0][k+NPT_MS+NPT_TR][0]) - pow(10,grid[i][j][0][NPT_MS][0])); // Age relative to TRGB
	       }
	       // delete points from main tracks
	       for (int k=0;k<NML;k++) for (int l=NPT_MS+1;l<NPT;l++) setUndefined(grid[i][j][k][l]);
	       NM_HB++;
	    }
	 }
      }

      if (FIXMODELS) fixmodsinit(i);
      fillMSmods(i);

      // apply mass loss
      for (int j=0;j<NM;j++) if (isDefined(grid[i][j][0][NPT_MS]) && !isDefined(grid[i][j][0][NPT_MS+1])) {
	 const double m0 = M[j];
	 const double mean_loss = m0-exp(grid[i][j][0][NPT_MS][1]);
	 double sigma_loss = sigma_loss=mean_loss*0.3;
	 sigma_loss = std::min(sigma_loss,m0/3.0);
	 sigma_loss = std::min(sigma_loss,0.04);
	 addHBtrack(i,j,mean_loss,3.0*sigma_loss);
      }

      // fix mass defined in tracks
      for (int j=0;j<NM;j++) for (int k=0;k<NML;k++) for (int l=0;l<NPT;l++) if (isDefined(grid[i][j][k][l])) grid[i][j][k][l][1] = log(M[j]);

      freemods();
   }
   return;
}

int main(int argc,char**argv) {
   for (int i=1;i<argc;i++) {
      if (!strcasecmp(argv[i],"-fix")) FIXMODELS = true;
      else if (!strcasecmp(argv[i],"-errok")) IGNOREERRS = true;
      else {
	 fprintf(stderr,"Usage: %s <options>\n",argv[0]);
	 fprintf(stderr,"  -fix   fixes trends of age with mass and point\n");
	 return -1;
      }
   }
   strcpy(modname,"MIST");
   init_mods();
   checkmods();
#ifdef PTHREADS
   init_threads();
   //call_thread(750,-18); call_thread(900,-18); call_thread(1010,-18);
   for (int it=modelITmin;it<modelITmax;it+=dIT) for (int iz=modelIZmin;iz<modelIZmax;iz+=dIZ) call_thread(it,iz);
   for (int it=0;it<PTHREADS;it++) while (thread_done[it]!=1) usleep(10);
#else
   for (int it=modelITmin;it<modelITmax;it+=dIT) for (int iz=modelIZmin;iz<modelIZmax;iz+=dIZ) addmod(0,it,iz);
#endif
   return 0;
}
