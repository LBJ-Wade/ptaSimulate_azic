#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_STRLEN 1024
#define MAX_CUTS 10 // Number of cuts that can be made to a data set
#define MAX_OUTPUT 10 // Maximum number of output definitions
#define MAX_LABELS 50 // Maximum number of labels for output definitions
#define MAX_PARAMS 40 // Maximum number of parameters in an ephemeris
#define MAX_PSRS 50 // Maximum number of pulsars
#define MAX_LINE_PARAMS 20 // Maximum number of parameters on a script line
#define MAX_OBSRUN 30 // Maximum number of observation runs
#define MAX_SCHED 100 // Maximum number of schedules
#define MAX_OBS_SCHED 100 // Maximum number of observations in schedule
#define MAX_TOAS 2000 // Maximum number of ToAs per pulsar
#define MAX_TNOISE 50 // Maximum number of timing noise definitions
#define MAX_PLANETS 50 // Maximum number of planets
#define MAX_RCVR 20 // Maximum number of receivers
#define MAX_BE 20 // Maximum number of backends
#define MAX_FLUX 10 // Maximum flux density values per pulsar
#define MAX_TSKY 10 // Maximum tsky values per pulsar
#define MAX_PROF_FILE 10 // Maximum number of profiles per pulsar
#define MAX_DMVAR 50 // Maximum number of DM variations definitions
#define MAX_DMCOVAR 50 // Maximum number of DM covariance definitions
#define MAX_DMFUNC 50 // Maximum number of DM function definitions
#define MAX_JITTER 50 // Maximum number of jitter definitions
#define MAX_SYS 5 // Maximum number of simultaneous systems in a given obsSys
#define MAX_BEOFFSETS 10 // Maximum number of backend offsets
#define MAX_OBS_SYS 10 
#define MAX_FDIST_VALS 5000
#define MAX_T2TIM 50 // Maximum number of tempo2 tim files
#define MAX_GLITCHES 10 // Maximum number of glitches in total
#define MAX_OBSERVATORIES 100 // Maximum number of observatories in observatories.dat file
#define DM_CONST    2.41e-4
#define SECDAY 86400
#define MAX_GWS 10 // Maximum number of GWs
#define MAX_CWS 100000 // Maximum number of individual SMBH sources to simulate

typedef struct paramStruct {
  char l[MAX_STRLEN]; // Label
  char v[MAX_STRLEN]; // Value
  int type; // 1 = just value, 2 = label and value
} paramStruct;

typedef struct valStruct {
  char inVal[MAX_STRLEN];
  double dval;
  int set; // 1 = set dval
  int constant; // 1 = Remain constant for different realisations, 0 = change
} valStruct;


typedef struct t2TimStruct {
  char psrName[MAX_STRLEN];
  char fileName[MAX_STRLEN];
  int  psrNum;
  valStruct toaErr;
  valStruct efac;
  valStruct equad;
} t2TimStruct;

typedef struct rcvrStruct {
  char name[MAX_STRLEN];
  valStruct flo;
  valStruct fhi;
  valStruct tsys;
  valStruct gain; // In K/Jy
} rcvrStruct;

typedef struct beStruct {
  char name[MAX_STRLEN];
  valStruct bw;
  valStruct nbin;
  int nOffset;
  valStruct offsetMJD[MAX_BEOFFSETS];
  valStruct offsetVal[MAX_BEOFFSETS];
} beStruct;

typedef struct glitchStruct {
  int psrNum;
  valStruct glep;
  valStruct glph;
  valStruct glf0;
  valStruct glf1;
  valStruct glf0d;
  valStruct gltd;
} glitchStruct;

typedef struct obsSysStruct {
  char name[MAX_STRLEN];
  int nSys;
  int rcvrNum[MAX_SYS];
  int beNum[MAX_SYS];
  char rcvrName[MAX_SYS][128];
  char beName[MAX_SYS][128];
  valStruct freq[MAX_SYS];
} obsSysStruct;

typedef struct obsStruct {
  char tel[MAX_STRLEN];
  char sched[MAX_STRLEN];
  char or[MAX_STRLEN];
  int psrNum;
  valStruct start;
  valStruct finish;
  valStruct toaErr;
  valStruct efac;
  valStruct equad;
  valStruct freq;
  valStruct tobs;
  valStruct ha;
  valStruct outlierAmp;
  valStruct outlierProb;
  long double sat;
  int satSet;
  int rcvrNum;
  int beNum;
  int obsSysNum;
} obsStruct;

typedef struct psrStruct {
  char name[MAX_STRLEN]; int setName;
  char label[MAX_STRLEN]; int setLabel;
  char ephem[MAX_STRLEN]; int setEphem;
  char rajStr[MAX_STRLEN];
  char decjStr[MAX_STRLEN];
  int  requireCatRead;
  int  requireEphemRead;
  int  nSetParam;  // Number of PSRCAT parameters set
  char setParamName[MAX_PARAMS][MAX_STRLEN];
  valStruct paramVal[MAX_PARAMS];
  int nToAs;
  obsStruct obs[MAX_TOAS];
  double rajd; // Position in degrees
  double decjd; // Position in degrees
  double dm; // dispersion measure
  double dist; // Distance in metres
  double p0; // Period
  double f0; // frequency
  char profileEqn[MAX_STRLEN]; int setProfileEqn;
  char profileFile[MAX_PROF_FILE][MAX_STRLEN]; valStruct freqProfileFile[MAX_PROF_FILE]; 
  int nProfileFile;
  valStruct flux[MAX_FLUX]; valStruct freqFlux[MAX_FLUX]; int nFlux;
  valStruct tsky[MAX_FLUX]; valStruct freqTsky[MAX_FLUX]; int nTsky;
  valStruct diff_df; valStruct diff_dfFreq; int setDiff_df;
  valStruct diff_ts; valStruct diff_tsFreq; int setDiff_ts;
} psrStruct;

typedef struct obsrunStruct {
  char name[MAX_STRLEN];
  char tel[MAX_STRLEN];
  char sched[MAX_STRLEN];
  int  setSched;
  valStruct start;
  valStruct finish;
  valStruct cadence;
  valStruct probFailure;
  t2TimStruct T2Tim[MAX_T2TIM];
  int nT2Tim;
} obsrunStruct;


typedef struct scheduleStruct {
  char name[MAX_STRLEN];
  obsStruct obs[MAX_OBS_SCHED];
  int nObsSched;
} scheduleStruct;

typedef struct observatoryStruct {
  char name1[MAX_STRLEN];
  char name2[MAX_STRLEN];
  double posX;
  double posY;
  double posZ;
} observatoryStruct;

typedef struct gwStruct {
  valStruct amp;
  valStruct alpha;
  valStruct ra;
  valStruct dec;
  valStruct ap;
  valStruct ac;
  valStruct gwmAmp;
  valStruct gwmEpoch;
  valStruct gwmPhi;

  valStruct gwcsAmp1;
  valStruct gwcsAmp2;
  valStruct gwcsEpoch;
  valStruct gwcsWidth;
  
  valStruct cgw_freq;
  valStruct cgw_h0;
  valStruct cgw_epoch;
  valStruct cgw_cosinc;
  valStruct cgw_angpol;
  valStruct cgw_mc;
  int type; // 1 = GWB, 2 = single defined by A+ and Ax, 3 = GWM, 4 = CW source
  char fname[1024]; // File name for GW source listing
} gwStruct;

typedef struct tnoiseStruct {
  int psrNum;
  valStruct alpha;
  valStruct beta;
  valStruct p0;
  valStruct fc;
  char predictMode[MAX_STRLEN];
  char label[MAX_STRLEN];
} tnoiseStruct;

typedef struct planetStruct {
  int psrNum;
  valStruct pb;
  valStruct ecc;
  valStruct a1;
  valStruct t0;
  valStruct om;
  char label[MAX_STRLEN];
} planetStruct;

typedef struct clkNoiseStruct {
  valStruct alpha;
  valStruct p0;
  valStruct fc;
  valStruct gwAmp;
} clkNoiseStruct;

typedef struct ephemNoiseStruct {
  valStruct alpha;
  valStruct p0;
  valStruct fc;
  valStruct gwAmp;
} ephemNoiseStruct;

typedef struct dmVarStruct {
  int psrNum;
  valStruct d_tscale;
  valStruct dVal;
  valStruct refFreq;
  int type;
} dmVarStruct;

typedef struct dmCovarStruct {
  int psrNum;
  valStruct alpha;
  valStruct a;
  valStruct b;
  int type;
} dmCovarStruct;

typedef struct dmFuncStruct {
  int psrNum;
  valStruct ddm;
} dmFuncStruct;

typedef struct jitterStruct {
  int psrNum;
  valStruct t0;
  valStruct sigma_j;
  valStruct refFreq;
  int type;
} jitterStruct;


typedef struct outputStruct {
  int nAdd;
  char label[MAX_LABELS][MAX_STRLEN];
  char fname[MAX_STRLEN];
} outputStruct;

typedef struct controlStruct {
  char name[MAX_STRLEN];
  int nproc;

  int nCut;
  float mjdCut[MAX_CUTS];
  char  cutName[MAX_CUTS][512];
  
  long seed; // Random number seed
  psrStruct psr[MAX_PSRS];
  int npsr; // Number of pulsars
  char inputScript[MAX_STRLEN];
  int nreal; // Number of realisations

  obsrunStruct obsRun[MAX_OBSRUN];
  int nObsRun; // Number of observation runs

  scheduleStruct sched[MAX_SCHED];
  int nSched; // Number of available schedules

  tnoiseStruct tnoise[MAX_TNOISE];
  int nTnoise;

  planetStruct planets[MAX_PLANETS];
  int nPlanets;

  clkNoiseStruct clkNoise;
  int nClkNoise;

  ephemNoiseStruct ephemNoise;
  int nEphemNoise;

  dmVarStruct dmVar[MAX_DMVAR];
  int nDMvar;

  dmCovarStruct dmCovar[MAX_DMCOVAR];
  int nDMcovar;

  dmFuncStruct dmFunc[MAX_DMFUNC];
  int nDMfunc;

  jitterStruct jitter[MAX_DMVAR];
  int nJitter;
  
  long double minT;
  long double maxT;

  gwStruct gw[MAX_GWS];
  int nGW;

  rcvrStruct rcvr[MAX_RCVR];
  int nRCVR;

  beStruct be[MAX_BE];
  int nBE;

  char simEphem[MAX_STRLEN]; int simTypeEphem;
  char simClock[MAX_STRLEN];
  char useEphem[MAX_STRLEN]; int useTypeEphem;
  char useClock[MAX_STRLEN];
  char useEOP[MAX_STRLEN];
  char simEOP[MAX_STRLEN];
  char useSWM[MAX_STRLEN];
  char simSWM[MAX_STRLEN];
  char useNE_SW[MAX_STRLEN];
  char simNE_SW[MAX_STRLEN];

  obsSysStruct obsSys[MAX_OBS_SYS];
  int nObsSys;

  glitchStruct glitches[MAX_GLITCHES];
  int nGlitches;

  char t2exe[MAX_STRLEN];
  char shell[MAX_STRLEN];
  char shellPth[MAX_STRLEN];

  observatoryStruct observatory[MAX_OBSERVATORIES];
  int nObservatory;

  outputStruct output[MAX_OUTPUT];
  int nOutput;

} controlStruct;

int runEvaluateExpression(char *expression,controlStruct *control);
int changeRandomOnce(char *expression,controlStruct *control);
int checkProbability(valStruct in,controlStruct *control);
void createOutliers(controlStruct *control,int r);
void processEphemNoise(controlStruct *control,int r);
void createEphemNoise(controlStruct *control,int r);
