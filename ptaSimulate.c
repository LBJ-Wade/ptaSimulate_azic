
// gcc -lm -o ptaSimulate *.c -lfftw3 /usr/lib/x86_64-linux-gnu/libfftw3f.so
// gcc -lm -o ptaSimulate *.c -L/usr/local/lib/ -lfftw3f
// gcc -lm -o ptaSimulate *.c -L/usr/local/lib/ -L../../fftw-3.3.4/.libs/  -I../../fftw-3.3.4/api/ -lfftw3f

// psrpop file
// Deal correctly with pulsar positions for different realisations
// Add in searches on ecliptic coordinates or equatorial coordinates
// check for unique names
// probMissed

// Update check that PSRCAT_RUNDIR set

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "ptaSimulate.h"
#include "T2toolkit.h"
#include "toasim.h"
#include "makeRedNoise.h"
#include "TKfit.h"
#include "GWsim.h"
#include "ptimeLib.h"
#include <fftw3.h>
#include <complex.h>
#include "evaldefs.h"

void createDirectoryStructure(controlStruct *control);
void createRunScript(controlStruct *control,char *dir0);
void initialiseControl(controlStruct *control);
void writeTimFiles(controlStruct *control,int r);
void fillDval(valStruct *param,controlStruct *control);
void loadInputs(controlStruct *control,int argc,char *argv[]);
void finishOff(controlStruct *control);
void readScript(controlStruct *control);
void createPulsarName(controlStruct *control,int p);
size_t trimwhitespace(char *out, size_t len, const char *str);
void readRCVRFromScript(controlStruct *control,FILE *fin);
void readBEFromScript(controlStruct *control,FILE *fin);
void readObsSysFromScript(controlStruct *control,FILE *fin);
double BTmodel(double pb,double ecc,double a1,double t0,double om,long double t);
void readPulsarsFromScript(controlStruct *control,FILE *fin);
void readObsRunFromScript(controlStruct *control,FILE *fin);
void readScheduleFromScript(controlStruct *control,FILE *fin);
void readDefineFromScript(controlStruct *control,FILE *fin);
void readT2filesFromScript(controlStruct *control,FILE *fin);
void readAdditionsFromScript(controlStruct *control,FILE *fin);
void readT2TimFile(controlStruct *control,int or,int t2Num,int r);
int getParams(char *line,char *label,paramStruct *p);
void processObsRun(controlStruct *control,int r);
void processSched(controlStruct *control,int r);
void processObsSys(controlStruct *control,int r);
void initialisePulsar(psrStruct *psr);
void processPulsars(controlStruct *control,int r);
void createIdealArrivalTimes(controlStruct *control,int r);
void loadEphemeris(psrStruct *psr,char *fname,controlStruct *control);
void createParSimulate(controlStruct *control,int r);
void convertToUpper(char *str);
int turn_hms(double turn, char *hms);
int turn_dms(double turn, char *dms);
double dms_turn(char *line);
double hms_turn(char *line);
double calculateToaErrRadiometer(controlStruct *control,int a,int b,int sys,double scale,int r);
double getToaErr(double *prof,double *templ,int nbin);
void makeRealScript(controlStruct *control,int r,char *dir0);
void createRadiometerNoise(controlStruct *control, int r);
void processTnoise(controlStruct *control,int r);
void createTnoise(controlStruct *control,int r);
void processPlanets(controlStruct *control,int r);
void createPlanets(controlStruct *control,int r);
void processClkNoise(controlStruct *control,int r);
void createClkNoise(controlStruct *control,int r);
void processGlitches(controlStruct *control,int r);
void processDMvar(controlStruct *control,int r);
void createDMvar(controlStruct *control,int r);
void createBEoffsets(controlStruct *control,int r);
void processDMcovar(controlStruct *control,int r);
void createDMcovar(controlStruct *control,int r);
void createDMfunc(controlStruct *control,int r);
void processJitter(controlStruct *control,int r);
void createJitter(controlStruct *control,int r);
void processGW(controlStruct *control,int r);
void createGW(controlStruct *control,int r);
void processBE(controlStruct *control,int r);
void processRCVR(controlStruct *control,int r);
void loadProfileFromFile(double *prof,double *templ,int nbin,char *fname);
double calcDiffractiveScint(controlStruct *control,int s0,int j,int sys);

void fftfit(double *prof,double *standard,int nmax,double *shift,double *eshift,
	    double *snr,double *esnr,double *b,double *errb,int *ngood);
void fft(double *y,int nmax,double *amp,double *pha);
void four1(double data[], unsigned long nn, int isign);
void fccf(double *amp,double *pha,double *shift,int nprof);
double dchisqr(double tau,double *tmp,double *r,int num);
double zbrent(double x1,double x2,double f1,double f2,double tol,double *tmp,
	      double *pha,int nsum);
double min(double a,double b);
double sign(double a,double b);
void readObservatoryPositions(controlStruct *control);
long double getTimeHA(controlStruct *control,double ha,long double sat,int telID,long double ra);
long double fortran_mod(long double a,long double p);

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

int main(int argc,char *argv[])
{
  controlStruct *control;
  char dir0[MAX_STRLEN];
  int r,p;
  printf("At the start\n");
  getcwd(dir0,MAX_STRLEN);

  control = (controlStruct *)malloc(sizeof(controlStruct));
  printf("init control\n");
  initialiseControl(control);
  printf("Load inputs\n");
  loadInputs(control,argc,argv);
  printf("Reading script\n");
  readScript(control);
  printf("Starting\n");
  createDirectoryStructure(control);
  printf("Complete directory structure\n");
  // Setup pulsars

  readObservatoryPositions(control);

  // Note that this should run for every iteration if the user requests that pulsar positions etc. change
  // for each iteration. Should only run once if kept constant
  for (r=0;r<control->nreal;r++)
    {
      // Reinitialise number of observations
      for (p=0;p<control->npsr;p++)
	control->psr[p].nToAs=0;

      printf("Creating realisation %d\n",r);
      processBE(control,r);
      processRCVR(control,r);

      processPulsars(control,r);

      processGlitches(control,r);
      // Create parameter files used in the simulation
      createParSimulate(control,r);
      printf("process ObsRun\n");
      processObsRun(control,r);
      printf("process Sched\n");
      processSched(control,r);
      printf("process ObsSys\n");
      processObsSys(control,r);
      printf("createIdealArrivaltimes\n");

      createIdealArrivalTimes(control,r);
      printf("writeArrivalTimes %d\n",r);
      writeTimFiles(control,r);
      printf("Create radiometer noise %d\n",r);
      createRadiometerNoise(control,r);
      printf("ProcessTnoise %d\n",r);
      processTnoise(control,r);
      printf("createTnoise %d\n",r);
      createTnoise(control,r);

      printf("ProcessPlanets %d\n",r);
      processPlanets(control,r);
      printf("createPlanets %d\n",r);
      createPlanets(control,r);

      printf("processCLKNoise\n");
      processClkNoise(control,r);
      printf("createCLKnoise\n");
      createClkNoise(control,r);

      printf("processEphemNoise\n");
      processEphemNoise(control,r);
      printf("createEphemNoise\n");
      createEphemNoise(control,r);

      printf("processDMvar\n");
      processDMvar(control,r);
      printf("CreateDMvar\n");
      createDMvar(control,r);

      printf("ProcessDMcovar %d\n",r);
      processDMcovar(control,r);
      createDMcovar(control,r);

      printf("CreateDMfunc %d\n",r);
      createDMfunc(control,r);

      printf("processJitter %d\n",r);
      processJitter(control,r);
      createJitter(control,r);

      printf("ProcessGW %d\n",r);
      processGW(control,r);
      createGW(control,r);
      
      //      printf("CreateBEoffsets %d\n",r);
      //      createBEoffsets(control,r);
      //      printf("CreateOutliers %d\n",r);
      createOutliers(control,r);

      printf("meanRealScript %d\n",r);
      makeRealScript(control,r,dir0);
    }
  createRunScript(control,dir0);

  finishOff(control);
}

void createRunScript(controlStruct *control,char *dir0)
{
  FILE *fout;
  char foutName[128];
  int i,j;
  int num;
  int proc;

  num = control->nreal/control->nproc;
  for (i=0;i<control->nproc;i++)
    {
      sprintf(foutName,"%s/scripts/runScripts_proc_%d",control->name,i);
      fout = fopen(foutName,"w");

      fprintf(fout,"#!%s/%s\n",control->shellPth,control->shell);
      for (j=0;j<num;j++)
	{
	  proc = i*num+j;
	  fprintf(fout,"source %s/%s/scripts/process_real_%d\n",dir0,control->name,proc);

	  if (strcmp(control->shell,"tcsh")==0)
	    {
	      fprintf(fout,"if (-e %s/%s/scripts/status/stopScript) then\n",dir0,control->name);
	      fprintf(fout," echo 'Stopping as stopScript exists'\n");
	      fprintf(fout," exit\n");
	      fprintf(fout,"endif\n");
	    }
	  else if (strcmp(control->shell,"bash")==0)
	    {
	      fprintf(fout,"if [ -a %s/%s/scripts/status/stopScript ]\n",dir0,control->name);
	      fprintf(fout," then\n");
	      fprintf(fout,"  echo 'Stopping as stopScript exists'\n");
	      fprintf(fout,"  exit\n");
	      fprintf(fout," fi\n");
	    }
	}
      // are we missing one?
      if (i==control->nproc-1)
	{
	  for (j=(control->nproc-1)*num+num;j<control->nreal;j++)
	    {
	      fprintf(fout,"source %s/%s/scripts/process_real_%d\n",dir0,control->name,j);
	      if (strcmp(control->shell,"tcsh")==0)
		{
		  fprintf(fout,"if (-e %s/%s/scripts/status/stopScript) then\n",dir0,control->name);
		  fprintf(fout," echo 'Stopping as stopScript exists'\n");
		  fprintf(fout," exit\n");
		  fprintf(fout,"endif\n");
		}
	      else if (strcmp(control->shell,"bash")==0)
		{
		  fprintf(fout,"if [ -a %s/%s/scripts/status/stopScript ]\n",dir0,control->name);
		  fprintf(fout," then\n");
		  fprintf(fout,"  echo 'Stopping as stopScript exists'\n");
		  fprintf(fout,"  exit\n");
		  fprintf(fout," fi\n");
		}
	    }
	}
      fclose(fout);
    }

  // Now create the master script
  printf("Making master script\n");
  sprintf(foutName,"%s/scripts/runScripts_master",control->name);
  fout = fopen(foutName,"w");
  fprintf(fout,"#!%s/%s\n",control->shellPth,control->shell);


  fprintf(fout,"unalias rm\n");
  fprintf(fout,"unalias mv\n");
  fprintf(fout,"unalias cp\n");
  fprintf(fout,"rm %s/%s/scripts/status/runStat\n",dir0,control->name);
  fprintf(fout,"set dte = `date`\n");
  fprintf(fout,"echo \"Processing start: $dte\" > %s/%s/scripts/status/runStat\n",dir0,control->name);
  for (i=0;i<control->nproc;i++)
    {
      fprintf(fout,"source %s/%s/scripts/runScripts_proc_%d &\n",dir0,control->name,i);
    }
  // Now check if the processing has finished
  fprintf(fout,"set endit=0\n"); 
  fprintf(fout,"while ($endit =~ \"0\")\n");
  fprintf(fout," set nf = `grep Complete %s/%s/scripts/status/runStat | wc -l | tail -1`\n",dir0,control->name);
  fprintf(fout,"  if ($nf =~ %d) then\n",control->nproc);
  fprintf(fout,"   set endit=1\n");
  fprintf(fout,"  else\n");
  fprintf(fout,"   echo \"master script for %s sleeping 1 second (processed $nf realisations)\"\n",control->name);
  fprintf(fout,"   sleep 1\n");
  fprintf(fout,"  endif\n");
  fprintf(fout,"end\n");
  fprintf(fout,"set dte = `date`\n");
  fprintf(fout,"echo \"Processing complete: $dte\" >> %s/%s/scripts/status/runStat\n",dir0,control->name);
  //  if (strcmp(control->email,"UNSET")==0)
  //    printf("No email set\n");
  //  else
  //    fprintf(fout,"mailx -s \"ptaSimulate: Complete processing for %s\" %s < %s/%s/scripts/status/runStat\n",control->name,control->email,dir0,control->name);
  fclose(fout);
  printf("Complete making master script\n");
}


void makeRealScript(controlStruct *control,int r,char *dir0)
{
  FILE *fout;
  char foutName[1024];
  char runStr[4096];
  char add[1024];
  int i,j,k,l,m;
  int beOff;
  int beNums[1024];
  int nBE,found;
  int include;

  printf("Making script\n");
  sprintf(foutName,"%s/scripts/process_real_%d",control->name,r);
  fout = fopen(foutName,"w");
  fprintf(fout,"#!%s/%s\n",control->shellPth,control->shell);
  fprintf(fout,"echo 'Processing realisation %d'\n",r);
  fprintf(fout,"set host = `hostname`\n");
  fprintf(fout,"set usr = `whoami`\n");
  fprintf(fout,"set dte = `date`\n");
  fprintf(fout,"set pid = `echo $$`\n");
  fprintf(fout,"echo \"[$dte] [$host] [$usr] [$pid] Processing realisation %d\" >> %s/%s/scripts/status/runStat\n",r,dir0,control->name);
  fprintf(fout,"cd %s/%s/workFiles/real_%d\n",dir0,control->name,r);

  for (i=0;i<control->npsr;i++)
    {
      // Consider stopping the script
      if (strcmp(control->shell,"tcsh")==0)
	{
	  fprintf(fout,"if (-e %s/%s/scripts/status/stopScript) then\n",dir0,control->name);
	  fprintf(fout," echo 'Stopping as stopScript exists'\n");
	  fprintf(fout," exit\n");
	  fprintf(fout,"endif\n");
	}
      else if (strcmp(control->shell,"bash")==0)
	{
	  fprintf(fout,"if [ -a %s/%s/scripts/status/stopScript ]\n",dir0,control->name);
	  fprintf(fout," then\n");
	  fprintf(fout,"  echo 'Stopping as stopScript exists'\n");
	  fprintf(fout,"  exit\n");
	  fprintf(fout," fi\n");
	}
      fprintf(fout,"%s -gr formIdeal -f %s.par.sim %s.itim\n",control->t2exe,control->psr[i].name,control->psr[i].name);
      fprintf(fout,"%s -output add_pulseNumber -f %s.par.sim %s.itim.sim\n",control->t2exe,control->psr[i].name,control->psr[i].name);

      //      fprintf(fout,"mv %s.itim.sim %s.sim\n",control->psr[i].name,control->psr[i].name);
      fprintf(fout,"mv withpn.tim %s.sim\n",control->psr[i].name,control->psr[i].name);

      for (l=0;l<control->nOutput;l++)
	{
	  
	  sprintf(runStr,"%s -gr createRealisation -f %s.sim -corr %s.addGauss",control->t2exe,control->psr[i].name,control->psr[i].name);
	  
	  // Shall we add backend jumps?
	  //	  sprintf(add," -corr %s.addBEoffsets",control->psr[i].name);
	  //	  strcat(runStr,add);
	  
	  // Outliers?
	  sprintf(add," -corr %s.addOutliers",control->psr[i].name);
	  strcat(runStr,add);
	  
	  // Does this pulsar have red noise model?
	  for (j=0;j<control->nTnoise;j++)
	    {
	      include=0;
	      if (l==0) include=1;
	      else {
		for (m=0;m<control->output[l].nAdd;m++)
		  {
		    if (strcmp(control->output[l].label[m],control->tnoise[j].label)==0)
		      include=1;
		  }
	      }

	      if (control->tnoise[j].psrNum == i && include==1)
		{
		  sprintf(add," -corr %s.tnoise.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  // Does this pulsar have planets?
	  for (j=0;j<control->nPlanets;j++)
	    {
	      include=0;
	      if (l==0) include=1;
	      else {
		for (m=0;m<control->output[l].nAdd;m++)
		  {
		    if (strcmp(control->output[l].label[m],control->planets[j].label)==0)
		      include=1;
		  }
	      }
	      
	      if (control->planets[j].psrNum == i && include==1)
		{
		  sprintf(add," -corr %s.planets.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  if (control->nClkNoise == 1)
	    {
	      sprintf(add," -corr %s.clknoise.%d",control->psr[i].name,0);
	      strcat(runStr,add);	      
	    }
	  
	  if (control->nEphemNoise == 1)
	    {
	      sprintf(add," -corr %s.ephemnoise.%d",control->psr[i].name,0);
	      strcat(runStr,add);	      
	    }
	  
	  
	  // Does this pulsar have red noise model?
	  for (j=0;j<control->nDMvar;j++)
	    {
	      if (control->dmVar[j].psrNum == i)
		{
		  sprintf(add," -corr %s.dmvar.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  // Does this pulsar have covariance model file?
	  for (j=0;j<control->nDMcovar;j++)
	    {
	      if (control->dmCovar[j].psrNum == i)
		{
		  sprintf(add," -corr %s.dmcovar.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  // Does this pulsar have DM function?
	  for (j=0;j<control->nDMfunc;j++)
	    {
	      if (control->dmFunc[j].psrNum == i)
		{
		  sprintf(add," -corr %s.dmfunc.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  // Does this pulsar have red noise model?
	  for (j=0;j<control->nJitter;j++)
	    {
	      if (control->jitter[j].psrNum == i)
		{
		  sprintf(add," -corr %s.jitter.%d",control->psr[i].name,j);
		  strcat(runStr,add);
		}
	    }
	  
	  for (j=0;j<control->nGW;j++)
	    {
	      sprintf(add," -corr %s.addGW.%d",control->psr[i].name,j);
	      strcat(runStr,add);
	    }
	  
	  //sprintf(runStr,"tempo2 -gr createRealisation -f %s.sim -corr %s.addGauss -corr %s.addBackendOffsets",control->psr[i].name,control->psr[i].name,control->psr[i].name);
	  //      sprintf(runStr,"tempo2 -gr createRealisation -f %s.sim -corr %s.addGauss",control->psr[i].name,control->psr[i].name);
	  /*      if (control->gwb_amp > 0)
		  {
		  sprintf(add," -corr %s.addGWB",control->psr[i].name);
		  strcat(runStr,add);
		  }
		  if (control->nRedNoise > 0)
		  {
		  sprintf(add," -corr %s.addRedNoise_0",control->psr[i].name);
	      strcat(runStr,add);
	      }*/
	  //      printf("Running: >%s<\n",runStr);
	  fprintf(fout,"%s\n",runStr);
      
	  fprintf(fout,"mv %s.sim.real %s.tim\n",control->psr[i].name,control->psr[i].name);
	  if (l==0)
	    fprintf(fout,"cp %s.tim %s/%s/output/real_%d/.\n",control->psr[i].name,dir0,control->name,r);
	  else
	    fprintf(fout,"cp %s.tim %s/%s/output/real_%d/%s/.\n",control->psr[i].name,dir0,control->name,r,control->output[l].fname);
	  fprintf(fout,"%s -f %s.par %s.tim -newpar\n",control->t2exe,control->psr[i].name,control->psr[i].name);
	  if (l==0)
	    fprintf(fout,"cp new.par %s/%s/output/real_%d/%s.par\n",dir0,control->name,r,control->psr[i].name);
	  else
	    fprintf(fout,"cp new.par %s/%s/output/real_%d/%s/%s.par\n",dir0,control->name,r,control->output[l].fname,control->psr[i].name);

	  // Now cope with the cuts
	  for (j=0;j<control->nCut;j++)
	    {
	      if (l==0)
		{
		  fprintf(fout,"mkdir %s/%s/output/real_%d/%s\n",dir0,control->name,r,control->cutName[j]);
		}
	      else
		{
		  fprintf(fout,"mkdir %s/%s/output/real_%d/%s/%s\n",dir0,control->name,r,control->output[l].fname,control->cutName[j]);
		}
	      fprintf(fout,"\\rm cut.%s.tim; awk '{if ($3 < %f) {print $0}}' %s.tim > cut.%s.tim \n",control->psr[i].name,control->mjdCut[j],control->psr[i].name,control->psr[i].name);
	      fprintf(fout,"%s -f %s.par cut.%s.tim -newpar\n",control->t2exe,control->psr[i].name,control->psr[i].name);
	      if (l==0)
		{
		  fprintf(fout,"cp new.par %s/%s/output/real_%d/%s/%s.par\n",dir0,control->name,r,control->cutName[j],control->psr[i].name);
		  fprintf(fout,"cp cut.%s.tim %s/%s/output/real_%d/%s/%s.tim\n",control->psr[i].name,dir0,control->name,r,control->cutName[j],control->psr[i].name);
		}
	      else
		{
		  fprintf(fout,"cp new.par %s/%s/output/real_%d/%s/%s/%s.par\n",dir0,control->name,r,control->output[l].fname,control->cutName[j],control->psr[i].name);
		  fprintf(fout,"cp cut.%s.tim %s/%s/output/real_%d/%s/%s/%s.tim\n",control->psr[i].name,dir0,control->name,r,control->output[l].fname,control->cutName[j],control->psr[i].name);
		}
	    }

	}
    }
  fprintf(fout,"set dte = `date`\n");
  fprintf(fout,"echo \"[$dte] [$host] [$usr] [$pid] Complete processing realisation %d\" >> %s/%s/scripts/status/runStat\n",r,dir0,control->name);
  fclose(fout);
  printf("End making script\n");
}

void writeTimFiles(controlStruct *control,int r)
{
  int i,j,k,p;
  FILE *fout;
  char fname[1024];

  for (p=0;p<control->npsr;p++)
    {
      sprintf(fname,"%s/workFiles/real_%d/%s.itim",control->name,r,control->psr[p].name);
      fout = fopen(fname,"w");
      fprintf(fout,"FORMAT 1\n");
      for (i=0;i<control->psr[p].nToAs;i++)
	{
	  fprintf(fout,"%d %.5f %15.15Lf %.5f %s -or %s -sched %s -tobs %g\n",i,
		  control->psr[p].obs[i].freq.dval,control->psr[p].obs[i].sat,
		  control->psr[p].obs[i].toaErr.dval*1e6,control->psr[p].obs[i].tel,
		  control->psr[p].obs[i].or,control->psr[p].obs[i].sched,control->psr[p].obs[i].tobs.dval);
	}
      fclose(fout);
    }
}

void createIdealArrivalTimes(controlStruct *control,int r)
{
  int i,j,s0,k,p0,sys;
  char tel[128];
  long double t0,addT,sat;
  double err,efac,equad,freq,scale=1;
  int telID=-1;
  int v=0;
  int ntoa;
  int fail;

  for (i=0;i<control->nObsRun;i++)
    {
      if (control->obsRun[i].setSched==1)
	{
	  printf("Processing obsRun: %d\n",i);
	  strcpy(tel,control->obsRun[i].tel);

	  // Find coordinates for this telescope
	  for (k=0;k<control->nObservatory;k++)
	    {
	      if (strcasecmp(tel,control->observatory[k].name1)==0 ||
		  strcasecmp(tel,control->observatory[k].name2)==0)
		{telID = k; break;}
	    }

	  // Find the correct schedule
	  printf("Schedule: %s\n",control->obsRun[i].sched);
	  s0=-1;
	  for (j=0;j<control->nSched;j++)
	    {
	      if (strcmp(control->obsRun[i].sched,control->sched[j].name)==0)
		{
		  s0=j;
		  break;
		}
	    }
	  if (s0==-1)
	    {
	      printf("ERROR: Cannot find schedule named %s\n",control->obsRun[i].sched);
	      finishOff(control);
	    }
	  printf("Found sched %d\n",s0);
	  // Now start calculating times
	  t0=control->obsRun[i].start.dval;
	  do {
	    // t0 is the start time of the observing session
	    // should add the time to the actual observation
	    for (j=0;j<control->sched[s0].nObsSched;j++)
	      {
		sys=0;
		do {
		  //		  printf("Processing system: %d %d\n",sys,control->sched[s0].obs[j].obsSysNum);
		  // Now add this observation to the correct pulsar
		  p0=control->sched[s0].obs[j].psrNum;
		  ntoa = control->psr[p0].nToAs;
		  
		  //		  printf("Setting tobs for %s\n",control->psr[p0].name);
		  fillDval(&(control->sched[s0].obs[j].tobs),control);
		  control->psr[p0].obs[ntoa].tobs.dval = control->sched[s0].obs[j].tobs.dval;
		  control->psr[p0].obs[ntoa].beNum = control->sched[s0].obs[j].beNum;
		  printf("Setting Result = %g %d %d\n", control->psr[p0].obs[ntoa].tobs.dval,p0,ntoa);
		  
		  // Update the error bar size
		  fillDval(&(control->sched[s0].obs[j].toaErr),control);
		  fillDval(&(control->sched[s0].obs[j].efac),control);
		  fillDval(&(control->sched[s0].obs[j].equad),control);

		  
		  if (control->psr[control->sched[s0].obs[j].psrNum].setDiff_df==1
		      && control->psr[control->sched[s0].obs[j].psrNum].setDiff_ts==1)
		    {
		      printf("Into diff scale %g %d %d\n",control->psr[p0].obs[ntoa].tobs.dval,p0,ntoa);
		      scale = calcDiffractiveScint(control,s0,j,sys);
		      printf("Out diff scale %g\n",control->psr[p0].obs[ntoa].tobs.dval);
		    }
		  else
		    scale=1;
		  
		  if (strcmp(control->sched[s0].obs[j].toaErr.inVal,"radiometer")==0)
		    err = calculateToaErrRadiometer(control,s0,j,sys,scale,r);
		  else
		    err=control->sched[s0].obs[j].toaErr.dval;
		  efac=control->sched[s0].obs[j].efac.dval;
		  equad=control->sched[s0].obs[j].equad.dval;
		  


		  if (control->sched[s0].obs[j].obsSysNum != -1)
		    freq=control->obsSys[control->sched[s0].obs[j].obsSysNum].freq[sys].dval;
		  else
		    {
		      fillDval(&(control->sched[s0].obs[j].freq),control);
		      freq=control->sched[s0].obs[j].freq.dval;
		    }

		  if ((control->sched[s0].obs[j].start.set==0 ||
		       t0 > control->sched[s0].obs[j].start.dval) &&
		      (control->sched[s0].obs[j].finish.set==0 ||
		       t0 < control->sched[s0].obs[j].finish.dval))
		    {
		      
		      sat=t0;
		      // Check if we should check for rise and set times
		      if (control->sched[s0].obs[j].ha.set==1)
			{
			  double ha;

			  if (telID == -1)
			    {
			      printf("Request use of hour angle range, but telescope %s not known in observatories.dat file\n",tel);
			      finishOff(control);
			    }
			  printf("Using ha for telescope %d\n",telID);
			  fillDval(&(control->sched[s0].obs[j].ha),control);
			  ha = control->sched[s0].obs[j].ha.dval;
			  sat = getTimeHA(control,ha,sat,telID,control->psr[p0].rajd);
			  //			  printf("haRange = %g\n",ha);

			  //			  exit(1);

			}

		      //	    printf("Process %s\n",control->psr[p0].name);
		      control->psr[p0].obs[ntoa].sat = sat;
		      if (sat > control->maxT) control->maxT = sat;
		      if (sat < control->minT) control->minT = sat;
		      control->psr[p0].obs[ntoa].freq.dval = freq;
		      control->psr[p0].obs[ntoa].toaErr.dval = err;
		      control->psr[p0].obs[ntoa].efac.dval = efac;
		      control->psr[p0].obs[ntoa].equad.dval = equad;
		      
		      // Look for outliers
		      printf("Outlier here: %s %s\n",control->sched[s0].obs[j].outlierAmp.inVal,control->sched[s0].obs[j].outlierProb.inVal);
		      control->psr[p0].obs[ntoa].outlierAmp.set = control->sched[s0].obs[j].outlierAmp.set;
		      control->psr[p0].obs[ntoa].outlierProb.set = control->sched[s0].obs[j].outlierProb.set;
		      strcpy(control->psr[p0].obs[ntoa].outlierAmp.inVal,control->sched[s0].obs[j].outlierAmp.inVal);
		      strcpy(control->psr[p0].obs[ntoa].outlierProb.inVal,control->sched[s0].obs[j].outlierProb.inVal);
		      strcpy(control->psr[p0].obs[ntoa].tel,tel);
		      strcpy(control->psr[p0].obs[ntoa].sched,control->sched[s0].name);
		      strcpy(control->psr[p0].obs[ntoa].or,control->obsRun[i].name);
		      (control->psr[p0].nToAs)++;
		    }
		  
		  sys++;
		} while (control->sched[s0].obs[j].obsSysNum!=-1 && sys < control->obsSys[control->sched[s0].obs[j].obsSysNum].nSys);
	      }
	    fail=0;
	    do {
	      fillDval(&(control->obsRun[i].cadence),control);
	      addT = control->obsRun[i].cadence.dval; 
	      t0 += addT;
	      if (control->obsRun[i].probFailure.set==1)
		fail = checkProbability(control->obsRun[i].probFailure,control);
	    } while (fail==1);
	  } while (t0 < control->obsRun[i].finish.dval);
	}
    }
}


void createParSimulate(controlStruct *control,int r)
{
  int i,p,j;
  FILE *fout,*fout2;
  char fname[MAX_STRLEN];
  int glnum;

  for (p=0;p<control->npsr;p++)
    {
      sprintf(fname,"%s/workFiles/real_%d/%s.par.sim",control->name,r,control->psr[p].name);
      fout = fopen(fname,"w");
      sprintf(fname,"%s/workFiles/real_%d/%s.par",control->name,r,control->psr[p].name);
      fout2 = fopen(fname,"w");

      fprintf(fout,"%-15.15s %s\n","PSRJ",control->psr[p].name);
      fprintf(fout2,"%-15.15s %s\n","PSRJ",control->psr[p].name);
      fprintf(fout,"CORRECT_TROPOSPHERE N\n");
      fprintf(fout2,"CORRECT_TROPOSPHERE N\n");
      fprintf(fout,"MODE 1\n");
      fprintf(fout2,"MODE 1\n");
      fprintf(fout,"CLK %s\n",control->simClock);
      fprintf(fout2,"CLK %s\n",control->useClock);
      if (control->simTypeEphem==1)
	fprintf(fout,"EPHEM %s\n",control->simEphem);
      else if (control->simTypeEphem==2)
	fprintf(fout,"EPH_FILE %s\n",control->simEphem);

      if (control->useTypeEphem==1)
	fprintf(fout2,"EPHEM %s\n",control->useEphem);
      else
	fprintf(fout2,"EPH_FILE %s\n",control->useEphem);

      fprintf(fout,"EOP_FILE %s\n",control->simEOP);
      fprintf(fout2,"EOP_FILE %s\n",control->useEOP);
      fprintf(fout2,"TRACK -2\n");
      if (strcasecmp(control->simSWM,"wso")==0)
	fprintf(fout,"SWM 1\n");
      if (strcasecmp(control->useSWM,"wso")==0)
	fprintf(fout2,"SWM 1\n");

      fprintf(fout,"NE_SW %s\n",control->simNE_SW);
      fprintf(fout2,"NE_SW %s\n",control->useNE_SW);

      for (i=0;i<control->psr[p].nSetParam;i++)
	{
	  if (strcasecmp(control->psr[p].setParamName[i],"RAJ")==0 || strcasecmp(control->psr[p].setParamName[i],"DECJ")==0 ||
	      strcasecmp(control->psr[p].setParamName[i],"BINARY")==0)
	    {
	      fprintf(fout,"%-15.15s %s\n",control->psr[p].setParamName[i],control->psr[p].paramVal[i].inVal);
	      // Change the 1 here to turn off fitting
	      fprintf(fout2,"%-15.15s %s\n",control->psr[p].setParamName[i],control->psr[p].paramVal[i].inVal);
	    }
	  else
	    {
	      fprintf(fout,"%-15.15s %.15g\n",control->psr[p].setParamName[i],control->psr[p].paramVal[i].dval);
	      // Turn on fitting
	      if (strcmp(control->psr[p].setParamName[i],"F0")==0 ||
		  strcmp(control->psr[p].setParamName[i],"F1")==0 ||
		  strcmp(control->psr[p].setParamName[i],"P0")==0 ||
		  strcmp(control->psr[p].setParamName[i],"P1")==0)
		fprintf(fout2,"%-15.15s %.15g 1\n",control->psr[p].setParamName[i],control->psr[p].paramVal[i].dval);
	      else
		fprintf(fout2,"%-15.15s %.15g\n",control->psr[p].setParamName[i],control->psr[p].paramVal[i].dval);
	    }
	}

      // Glitch events

      glnum=1;
      for (j=0;j<control->nGlitches;j++)
	{
	  printf("glitch: %d %d %d\n",p,j,control->glitches[j].psrNum);
	  if (control->glitches[j].psrNum==p)
	    {
	      fprintf(fout,"GLEP_%d %g\n",glnum,control->glitches[j].glep.dval);
	      if (control->glitches[j].glph.set==1)
		fprintf(fout,"GLPH_%d %g\n",glnum,control->glitches[j].glph.dval);
	      if (control->glitches[j].glf0.set==1)
		fprintf(fout,"GLF0_%d %g\n",glnum,control->glitches[j].glf0.dval);
	      if (control->glitches[j].glf1.set==1)
		fprintf(fout,"GLF1_%d %g\n",glnum,control->glitches[j].glf1.dval);
	      if (control->glitches[j].glf0d.set==1)
		fprintf(fout,"GLF0D_%d %g\n",glnum,control->glitches[j].glf0d.dval);
	      if (control->glitches[j].gltd.set==1)
		fprintf(fout,"GLTD_%d %g\n",glnum,control->glitches[j].gltd.dval);
	      glnum++;
	    }
	}


      fclose(fout);
      fclose(fout2);
    }
  
}

void convertToUpper(char *str)
{
  int i;
  for (i=0;i<strlen(str);i++)
    str[i] = toupper(str[i]);
}

void processPulsars(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
  char expression[1024];
  int errorFlag;
  double dist,dm;
  //
  // 1. Check if these parameters should remain constant, change or be set for the first time
  //
  // 2. Fill parameters into pulsar structure
  //      - run PSRCAT
  //      - load from ephemeris or PSRCAT file
  //      - load parameters from population file
  //      - evaluate the expressions given in script

  // 3. Make name as required (and ensure that the name is unique and warn if the name doesn't match the position)
  // 

  for (p=0;p<control->npsr;p++)
    {
      if (control->psr[p].setDiff_df==1)
	{
	  fillDval(&(control->psr[p].diff_df),control);
	  fillDval(&(control->psr[p].diff_dfFreq),control);
	}

      if (control->psr[p].setDiff_ts==1)
	{
	  fillDval(&(control->psr[p].diff_ts),control);
	  fillDval(&(control->psr[p].diff_tsFreq),control);
	}

      for (i=0;i<control->psr[p].nFlux;i++)
	{
	  fillDval(&(control->psr[p].flux[i]),control);
	  fillDval(&(control->psr[p].freqFlux[i]),control);
	}

      for (i=0;i<control->psr[p].nTsky;i++)
	{
	  fillDval(&(control->psr[p].tsky[i]),control);
	  fillDval(&(control->psr[p].freqTsky[i]),control);
	}
      for (i=0;i<control->psr[p].nProfileFile;i++)
	fillDval(&(control->psr[p].freqProfileFile[i]),control);
      
      if (control->psr[p].requireCatRead==1 && r==0) {
	printf("Reading catalogue for %s\n",control->psr[p].name);
	sprintf(str,"psrcat -all -e %s > psrcat.res\n",control->psr[p].name);
	system(str);

	loadEphemeris(&(control->psr[p]),"psrcat.res",control);
	control->psr[p].requireCatRead=2; // So it doesn't get reloaded for the next iteration
	printf("Done loading ephemeris\n");
      }
      if (control->psr[p].requireEphemRead==1 && r==0) {
	printf("Read ephemeris: %s\n",control->psr[p].ephem);
	loadEphemeris(&(control->psr[p]),control->psr[p].ephem,control);
	control->psr[p].requireEphemRead=2;
      }
      // Evaluate expressions
      //      printf("Pulsar: %d\n",p);
      printf("Looking for parameter set %d\n",control->psr[p].nSetParam);
      for (i=0;i<control->psr[p].nSetParam;i++)
	{
	  printf("str = %s %s\n",control->psr[p].setParamName[i],
		 control->psr[p].paramVal[i].inVal);
	  sprintf(expression,"v = %s",control->psr[p].paramVal[i].inVal);
	  nVariables = 0;
	  //	  printf("evaluating: %s\n",expression);
	  if (r==0 || control->psr[p].paramVal[i].constant == 0)
	    {
	      errorFlag = runEvaluateExpression(expression,control);
	      convertToUpper(control->psr[p].setParamName[i]);
	      //	  printf("[%s] [%d] [%d] result = %g\n",expression,errorFlag,nVariables,variable[0].value);
	      control->psr[p].paramVal[i].dval = variable[0].value;
	      control->psr[p].paramVal[i].set = 1;

	      if (strcmp(control->psr[p].setParamName[i],"RAJ")==0)
		strcpy(control->psr[p].rajStr,control->psr[p].paramVal[i].inVal);
	      else if (strcmp(control->psr[p].setParamName[i],"DECJ")==0)
		strcpy(control->psr[p].decjStr,control->psr[p].paramVal[i].inVal);
	      else if (strcmp(control->psr[p].setParamName[i],"F0")==0)
		{
		  double f0,p0;
		  sscanf(control->psr[p].paramVal[i].inVal,"%lf",&f0);
		  p0 = 1.0/f0;
		  control->psr[p].f0 = f0;
		  control->psr[p].p0 = p0;
		}
	      else if (strcmp(control->psr[p].setParamName[i],"P0")==0 || strcmp(control->psr[p].setParamName[i],"P")==0)
		{
		  double f0,p0;
		  sscanf(control->psr[p].paramVal[i].inVal,"%lf",&p0);
		  f0 = 1.0/p0;
		  control->psr[p].f0 = f0;
		  control->psr[p].p0 = p0;
		}
	    }
	  if (strcmp(control->psr[p].setParamName[i],"DM")==0)
	    control->psr[p].dm = control->psr[p].paramVal[i].dval;
	}

      // Create a name for this pulsar
      if (control->psr[p].setName == 0)
	{
	  printf("Must set a name for pulsar %d %s %s\n",p,control->psr[p].rajStr,control->psr[p].decjStr);
	  createPulsarName(control,p);
	}

      // Create pulsar distance
      // Should check to see if a distance parameter has been set
      dm = control->psr[p].dm;
      dist = (dm/0.03)*3.08567758e19; // Very simple electron density model!
      control->psr[p].dist = dist;
      printf("Dist = %g %g\n",dist,control->psr[p].dist);
      // Obtain rajd and decjd
      control->psr[p].rajd = hms_turn(control->psr[p].rajStr)*360.0;
      control->psr[p].decjd = dms_turn(control->psr[p].decjStr)*360.0;
    }

}
void processObsSys(controlStruct *control,int r)
{
  int i,k,j,l;
  char str[1024];
  int r0,b0;

  for (i=0;i<control->nObsSys;i++)
    {

      for (j=0;j<control->obsSys[i].nSys;j++)
	{
	  fillDval(&(control->obsSys[i].freq[j]),control);      

	  r0=-1;
	  for (k=0;k<control->nRCVR;k++)
	    {
	      if (strcmp(control->rcvr[k].name,control->obsSys[i].rcvrName[j])==0)
		{r0 = k; break;}
	    }

	  b0=-1;
	  for (k=0;k<control->nBE;k++)
	    {
	      if (strcmp(control->be[k].name,control->obsSys[i].beName[j])==0)
		{b0 = k; break;}
	    }
 	  if (r0==-1)
	    {
	      printf("Cannot find recevier in obsSys: %s\n",control->obsSys[i].rcvrName[j]);
	      finishOff(control);
	    }
 	  if (b0==-1)
	    {
	      printf("Cannot find backend in obsSys: %s\n",control->obsSys[i].beName[j]);
	      finishOff(control);
	    }
	  control->obsSys[i].rcvrNum[j] = r0;
	  control->obsSys[i].beNum[j] = b0;

	}
    }
}

void processObsRun(controlStruct *control,int r)
{
  int i,p,k;
  char str[1024];
 
  for (p=0;p<control->nObsRun;p++)
    {
      // Now look for Tim file
      for (k=0;k<control->obsRun[p].nT2Tim;k++)
	{
	  printf("Processing psr %s\n",control->obsRun[p].T2Tim[k].psrName);
	  readT2TimFile(control,p,k,r);
	}

      // Evaluate expressions
      fillDval(&(control->obsRun[p].start),control);
      fillDval(&(control->obsRun[p].finish),control);
      //      fillDval(&(control->obsRun[p].cadence),control);
      //      printf("obsRun %g %g %g\n",control->obsRun[p].start.dval,control->obsRun[p].finish.dval,control->obsRun[p].cadence.dval);
    }

}

void processSched(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (p=0;p<control->nSched;p++)
    {
      for (i=0;i<control->sched[p].nObsSched;i++)
	{
	  // Evaluate expressions
	  if (strcmp(control->sched[p].obs[i].toaErr.inVal,"radiometer")!=0)
	    fillDval(&(control->sched[p].obs[i].toaErr),control);
	  fillDval(&(control->sched[p].obs[i].efac),control);
	  fillDval(&(control->sched[p].obs[i].equad),control);

	  fillDval(&(control->sched[p].obs[i].freq),control);
	  fillDval(&(control->sched[p].obs[i].start),control);
	  fillDval(&(control->sched[p].obs[i].tobs),control);
	  fillDval(&(control->sched[p].obs[i].finish),control);
	}
    }

}

void processTnoise(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nTnoise;i++)
    {
	  // Evaluate expressions
      if (strcmp(control->tnoise[i].alpha.inVal,"gwamp_auto")==0)
	{
	  double fc,alpha,p0;
	  double secperyear = 86400.0*365.25;

	  alpha = -13.0/3.0;
	  fc = 1.0/(((control->maxT-control->minT)/365.25)*2);
	  fillDval(&(control->tnoise[i].p0),control);
	  p0 = control->tnoise[i].p0.dval;

	  control->tnoise[i].alpha.dval = alpha;
	  control->tnoise[i].fc.dval = fc;
	  control->tnoise[i].p0.dval = pow(p0,2)/12.0/M_PI/M_PI*pow(fc,alpha);//*secperyear*secperyear;
	}
      else
	{
	  fillDval(&(control->tnoise[i].alpha),control);
	  fillDval(&(control->tnoise[i].beta),control);
	  fillDval(&(control->tnoise[i].p0),control);
	  fillDval(&(control->tnoise[i].fc),control);
	}
    }

}

void processPlanets(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nPlanets;i++)
    {
      fillDval(&(control->planets[i].pb),control);
      fillDval(&(control->planets[i].ecc),control);
      fillDval(&(control->planets[i].a1),control);
      fillDval(&(control->planets[i].t0),control);
      fillDval(&(control->planets[i].om),control);
    }

}

void processClkNoise(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  // Evaluate expressions
  if (control->clkNoise.gwAmp.set==1)
    {
      double fc,alpha,p0;
      double secperyear = 86400.0*365.25;
      
      alpha = -13.0/3.0;
      fc = 1.0/(((control->maxT-control->minT)/365.25)*2);
      fillDval(&(control->clkNoise.gwAmp),control);
      p0 = control->clkNoise.gwAmp.dval;

      control->clkNoise.alpha.dval = alpha;
      control->clkNoise.fc.dval = fc;
      control->clkNoise.p0.dval = pow(p0,2)/12.0/M_PI/M_PI*pow(fc,alpha);//*secperyear*secperyear;
      
    }
  else
    {
      fillDval(&(control->clkNoise.alpha),control);
      fillDval(&(control->clkNoise.p0),control);
      fillDval(&(control->clkNoise.fc),control);
    }
}


void processDMvar(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nDMvar;i++)
    {
	  // Evaluate expressions
      fillDval(&(control->dmVar[i].d_tscale),control);
      fillDval(&(control->dmVar[i].dVal),control);
      fillDval(&(control->dmVar[i].refFreq),control);
    }
}

void processGlitches(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nGlitches;i++)
    {
      // Evaluate expressions
      fillDval(&(control->glitches[i].glep),control);
      fillDval(&(control->glitches[i].glph),control);
      fillDval(&(control->glitches[i].glf0),control);
      fillDval(&(control->glitches[i].glf1),control);
      fillDval(&(control->glitches[i].glf0d),control);
      fillDval(&(control->glitches[i].gltd),control);
    }
}

void processDMcovar(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nDMcovar;i++)
    {
	  // Evaluate expressions
      fillDval(&(control->dmCovar[i].alpha),control);
      fillDval(&(control->dmCovar[i].a),control);
      fillDval(&(control->dmCovar[i].b),control);
    }
}

void processJitter(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nJitter;i++)
    {
	  // Evaluate expressions
      fillDval(&(control->jitter[i].t0),control);
      fillDval(&(control->jitter[i].refFreq),control);
      fillDval(&(control->jitter[i].sigma_j),control);
    }
}

void processGW(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nGW;i++)
    {
	  // Evaluate expressions
      if (control->gw[i].type==1)
	{
	  fillDval(&(control->gw[i].alpha),control);
	  fillDval(&(control->gw[i].amp),control);
	}
      else if (control->gw[i].type==2)
	{
	  fillDval(&(control->gw[i].ra),control);
	  fillDval(&(control->gw[i].dec),control);
	}
      else if (control->gw[i].type==3)
	{
	  fillDval(&(control->gw[i].ra),control);
	  fillDval(&(control->gw[i].dec),control);
	  fillDval(&(control->gw[i].gwmAmp),control);
	  fillDval(&(control->gw[i].gwmEpoch),control);
	  fillDval(&(control->gw[i].gwmPhi),control);
	}
      else if (control->gw[i].type==4)
	{
	  fillDval(&(control->gw[i].ra),control);
	  fillDval(&(control->gw[i].dec),control);
	  fillDval(&(control->gw[i].cgw_freq),control);
	  fillDval(&(control->gw[i].cgw_h0),control);
	  fillDval(&(control->gw[i].cgw_epoch),control);
	  fillDval(&(control->gw[i].cgw_cosinc),control);
	  fillDval(&(control->gw[i].cgw_angpol),control);
	  fillDval(&(control->gw[i].cgw_mc),control);
	}
      else if (control->gw[i].type==6)
	{
	  fillDval(&(control->gw[i].ra),control);
	  fillDval(&(control->gw[i].dec),control);
	  fillDval(&(control->gw[i].gwcsAmp1),control);
	  fillDval(&(control->gw[i].gwcsAmp2),control);
	  fillDval(&(control->gw[i].gwcsEpoch),control);
	  fillDval(&(control->gw[i].gwcsWidth),control);
	}
    }

}

void processBE(controlStruct *control,int r)
{
  int i,p,j;
  char str[1024];
 
  for (i=0;i<control->nBE;i++)
    {
      // Evaluate expressions
      fillDval(&(control->be[i].bw),control);
      fillDval(&(control->be[i].nbin),control);
      for (j=0;j<control->be[i].nOffset;j++)
	{
	  fillDval(&(control->be[i].offsetMJD[j]),control);
	  fillDval(&(control->be[i].offsetVal[j]),control);
	}
    }

}

void processRCVR(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  for (i=0;i<control->nRCVR;i++)
    {
      // Evaluate expressions
      fillDval(&(control->rcvr[i].flo),control);
      fillDval(&(control->rcvr[i].fhi),control);
      fillDval(&(control->rcvr[i].tsys),control);
      fillDval(&(control->rcvr[i].gain),control);
      printf("gain = %g\n",control->rcvr[i].gain.dval);
    }

}

void fillDval(valStruct *param,controlStruct *control)
{
 char expression[1024];
 int errorFlag;
 
 sprintf(expression,"v = %s",param->inVal);

 nVariables = 0;
 errorFlag = runEvaluateExpression(expression,control);
 param->dval = variable[0].value;
}


void createPulsarName(controlStruct *control,int p)
{
  char *tok;
  int hr,min,deg;
  double sec;
  char psrName1[128],psrName2[128];

  tok = strtok(control->psr[p].rajStr,":");
  sscanf(tok,"%d",&hr);
  tok = strtok(NULL,":");
  sscanf(tok,"%d",&min);
  sprintf(psrName1,"J%02d%02d",hr,min);

  tok = strtok(control->psr[p].decjStr,":");
  sscanf(tok,"%d",&deg);
  tok = strtok(NULL,":");
  sscanf(tok,"%d",&min);
  //  printf("Here with %d %d %02d%02d\n",deg,min,deg,min);
  if (deg < 0)
    strcat(psrName1,"-");
  else
    strcat(psrName1,"+");
  sprintf(psrName2,"%02d%02d",abs(deg),min);
  strcat(psrName1,psrName2);
  printf("Name: %s\n",psrName1);
  strcpy(control->psr[p].name,psrName1);
}

void loadEphemeris(psrStruct *psr,char *fname,controlStruct *control)
{
  FILE *fin;
  char label[1024],val[1024];
  char line[1024];
  int np;
  int exist=0;
  int i;

  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to open ephemeris file: %s\n",fname);
      finishOff(control);
    }
  while (!feof(fin))
    {
      fgets(line,1024,fin);
      sscanf(line,"%s %s",label,val);
      if (strcmp(label,"PSRJ")!=0 && strcmp(label,"RM")!=0 && strcmp(label,"EPHVER")!=0 &&
	  strcmp(label,"UNITS")!=0 && strcmp(label,"CLK")!=0 && strcmp(label,"EPHEM")!=0 &&
	  strcmp(label,"TZRMJD")!=0 && strcmp(label,"TZRSITE")!=0 && strcmp(label,"TZRFRQ")!=0 &&
	  strcmp(label,"NTOA")!=0 && strcmp(label,"TRES")!=0 && strcmp(label,"START")!=0 &&
	  strcmp(label,"FINISH")!=0)
	{
	  np = psr->nSetParam;
	  // Check if already exists
	  exist=0;
	  for (i=0;i<np;i++)
	    {
	      if (strcasecmp(label,psr->setParamName[i])==0)
		{
		  exist=1;
		  break;
		}
	    }
	  if (exist==0)
	    {
	      strcpy(psr->setParamName[np],label);
	      strcpy(psr->paramVal[np].inVal,val);
	      (psr->nSetParam)++;
	    }
	}
      // Should check reading pulsar name here
      else if (strcmp(label,"PSRJ")==0)
	{
	  // Set the pulsar name
	  strcpy(psr->name,val);
	  psr->setName = 1;
	}
    }
  fclose(fin);
  printf("HEre with Name %d\n",psr->nSetParam);
  if (psr->nSetParam < 2)
    {
      printf("Have a problem reading the ephemeris\n");
      printf("Have you set $PSRCAT_RUNDIR?");
      exit(1);
    }
}


void readScript(controlStruct *control)
{
  FILE *fin;
  char line[1024],trimLine[1024];

  if (!(fin = fopen(control->inputScript,"r")))
    {
      printf("Unable to open input script: %s\n",control->inputScript);
      finishOff(control);
    }
  while (!feof(fin))
    {
      if (fgets(line,1024,fin)!=NULL)
	{
	  trimwhitespace(trimLine,1024,line);
	  printf("Reading line %s\n",line);

	  if (strcmp(trimLine,"<pulsars>")==0)
	    {
	      readPulsarsFromScript(control,fin);
	      processPulsars(control,0); // Make sure that pulsar names are set
	    }
	  else if (strcmp(trimLine,"<obsRun>")==0)
	    readObsRunFromScript(control,fin);
	  else if (strcmp(trimLine,"<schedule>")==0)
	    readScheduleFromScript(control,fin);
	  else if (strcmp(trimLine,"<define>")==0)
	    readDefineFromScript(control,fin);
	  else if (strcmp(trimLine,"<t2files>")==0)
	    readT2filesFromScript(control,fin);
	  else if (strcmp(trimLine,"<add>")==0)
	    readAdditionsFromScript(control,fin);
	  else if (strcmp(trimLine,"<rcvr>")==0)
	    readRCVRFromScript(control,fin);
	  else if (strcmp(trimLine,"<obsSys>")==0)
	    readObsSysFromScript(control,fin);
	  else if (strcmp(trimLine,"<be>")==0)
	    readBEFromScript(control,fin);
	  else
	    printf("%s\n",trimLine);

	}
    }
  fclose(fin);
}

void readObsRunFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr;
  int or=control->nObsRun;
  control->obsRun[or].nT2Tim=0;
  control->obsRun[or].setSched=0;
  control->obsRun[or].probFailure.set=0;
  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the obsRun information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</obsRun>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->obsRun[or].name,p[0].v);
	else if (strcmp(label,"tel:")==0)
	  strcpy(control->obsRun[or].tel,p[0].v);
	else if (strcmp(label,"sched:")==0)
	  {strcpy(control->obsRun[or].sched,p[0].v);  control->obsRun[or].setSched=1;}
	else if (strcmp(label,"start:")==0)
	  strcpy(control->obsRun[or].start.inVal,p[0].v);
	else if (strcmp(label,"finish:")==0)
	  strcpy(control->obsRun[or].finish.inVal,p[0].v);
	else if (strcmp(label,"sampling:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"cadence")==0)
		  strcpy(control->obsRun[or].cadence.inVal,p[i].v);
		if (strcmp(p[i].l,"probFailure")==0)
		  {
		    strcpy(control->obsRun[or].probFailure.inVal,p[i].v);
		    control->obsRun[or].probFailure.set=1;
		  }
	      }
	  }
	else if (strcmp(label,"t2tim:")==0)
	  {
	    int nt = control->obsRun[or].nT2Tim;
	    control->obsRun[or].T2Tim[nt].toaErr.set=0;
	    control->obsRun[or].T2Tim[nt].efac.set=1;
	    control->obsRun[or].T2Tim[nt].equad.set=1;
	    strcpy(control->obsRun[or].T2Tim[nt].efac.inVal,"1");
	    strcpy(control->obsRun[or].T2Tim[nt].equad.inVal,"0");

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(control->obsRun[or].T2Tim[nt].psrName,p[i].v);
		else if (strcmp(p[i].l,"file")==0)
		  strcpy(control->obsRun[or].T2Tim[nt].fileName,p[i].v);
		else if (strcmp(p[i].l,"toaerr")==0)
		  {
		    strcpy(control->obsRun[or].T2Tim[nt].toaErr.inVal,p[i].v);
		    control->obsRun[or].T2Tim[nt].toaErr.set=1;
		  }
		else if (strcmp(p[i].l,"efac")==0)
		  strcpy(control->obsRun[or].T2Tim[nt].efac.inVal,p[i].v);
		else if (strcmp(p[i].l,"equad")==0)
		  strcpy(control->obsRun[or].T2Tim[nt].equad.inVal,p[i].v);
	      }
	    (control->obsRun[or].nT2Tim)++;
	  }

	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
  (control->nObsRun)++;
}

void readRCVRFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr;
  int or=control->nRCVR;
  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the rcvr information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</rcvr>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->rcvr[or].name,p[0].v);
	else if (strcmp(label,"flo:")==0)
	  strcpy(control->rcvr[or].flo.inVal,p[0].v);
	else if (strcmp(label,"fhi:")==0)
	  strcpy(control->rcvr[or].fhi.inVal,p[0].v);
	else if (strcmp(label,"tsys:")==0)
	  strcpy(control->rcvr[or].tsys.inVal,p[0].v);
	else if (strcmp(label,"gain:")==0)
	  strcpy(control->rcvr[or].gain.inVal,p[0].v);
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
  (control->nRCVR)++;
}

void readBEFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr;
  int or=control->nBE;

  control->be[or].nOffset=0;

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the backend information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</be>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->be[or].name,p[0].v);
	else if (strcmp(label,"bw:")==0)
	  strcpy(control->be[or].bw.inVal,p[0].v);
	else if (strcmp(label,"nbin:")==0)
	  strcpy(control->be[or].nbin.inVal,p[0].v);
	else if (strcmp(label,"offset:")==0)
	  {
	    int nos = control->be[or].nOffset;
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"mjd")==0)
		  strcpy(control->be[or].offsetMJD[nos].inVal,p[i].v);
		else if (strcmp(p[i].l,"size")==0)
		  strcpy(control->be[or].offsetVal[nos].inVal,p[i].v);
	      }	   
	    (control->be[or].nOffset)++;

	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
  (control->nBE)++;
}

void readObsSysFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr;
  int or=control->nObsSys;
  control->obsSys[or].nSys=0;

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the obssys information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</obsSys>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->obsSys[or].name,p[0].v);
	else if (strcmp(label,"system:")==0)
	  {
	    int nos = control->obsSys[or].nSys;
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"rcvr")==0)
		  strcpy(control->obsSys[or].rcvrName[nos],p[i].v);
		else if (strcmp(p[i].l,"be")==0)
		  strcpy(control->obsSys[or].beName[nos],p[i].v);
		else if (strcmp(p[i].l,"freq")==0)
		  strcpy(control->obsSys[or].freq[nos].inVal,p[i].v);
	      }	   
	    (control->obsSys[or].nSys)++;
	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
  (control->nObsSys)++;
}

void readScheduleFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr,j;
  int ns=control->nSched;
  control->sched[ns].nObsSched=0;
  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the schedule information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</schedule>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->sched[ns].name,p[0].v);
	else if (strcmp(label,"observe:")==0)
	  {
	    int no = control->sched[ns].nObsSched;
	    char pname[1024];
	    char toaErr[1024],freq[1024],start[1024],finish[1024],tobs[1024],ha[1024];
	    char efac[1024],equad[1024];
	    char outlierAmp[1024],outlierProb[1024];
	    int setlabel=0,setStart=0,setFinish=0,setTobs=0;
	    char label[1024];
	    int setRcvr=-1;
	    int setBE=-1;
	    int setObsSys=-1;
	    int setHa=-1;
	    int setOutlier=-1;

	    strcpy(efac,"1");
	    strcpy(equad,"0");

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"toaerr")==0)
		  strcpy(toaErr,p[i].v);
		else if (strcmp(p[i].l,"efac")==0)
		  strcpy(efac,p[i].v);
		else if (strcmp(p[i].l,"equad")==0)
		  strcpy(equad,p[i].v);
		else if (strcmp(p[i].l,"tobs")==0)
		  {strcpy(tobs,p[i].v); setTobs=1;}
		else if (strcmp(p[i].l,"start")==0)
		  {strcpy(start,p[i].v); setStart=1;}
		else if (strcmp(p[i].l,"finish")==0)
		  {strcpy(finish,p[i].v); setFinish=1;}
		else if (strcmp(p[i].l,"freq")==0)
		  strcpy(freq,p[i].v);
		else if (strcmp(p[i].l,"ha")==0)
		  {strcpy(ha,p[i].v); setHa=1;}
		else if (strcmp(p[i].l,"obsSys")==0)
		  {
		    printf("Checking obsSys %d\n",control->nObsSys);
		    for (j=0;j<control->nObsSys;j++)
		      {
			if (strcmp(p[i].v,control->obsSys[j].name)==0)
			  {setObsSys=j; break;}
		      }
		    printf("Got Obssys %d\n",setObsSys);
		  }
		else if (strcmp(p[i].l,"rcvr")==0)
		  {
		    for (j=0;j<control->nRCVR;j++)
		      {
			if (strcmp(p[i].v,control->rcvr[j].name)==0)
			  {setRcvr=j; break;}
		      }
		  }
		else if (strcmp(p[i].l,"be")==0)
		  {
		    for (j=0;j<control->nBE;j++)
		      {
			if (strcmp(p[i].v,control->be[j].name)==0)
			  {setBE=j; break;}
		      }
		  }
		else if (strcmp(p[i].l,"outlier")==0)
		  {
		    char *tok;
		    setOutlier=1;
		    tok = strtok(p[i].v,";");
		    printf("Trying: %s\n",tok+1);
		    strcpy(outlierAmp,tok+1);
		    tok = strtok(NULL,"]");
		    printf("now Trying %s\n",tok);
		    strcpy(outlierProb,tok);
		  }
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
;
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			//			if (setlabel==1)
			//			  printf("label: in here\n");
			no = control->sched[ns].nObsSched;
			//			strcpy(control->sched[ns].obs[no].psrName,control->psr[i].name);
			control->sched[ns].obs[no].psrNum = i;
			strcpy(control->sched[ns].obs[no].toaErr.inVal,toaErr);
			strcpy(control->sched[ns].obs[no].efac.inVal,efac);
			strcpy(control->sched[ns].obs[no].equad.inVal,equad);
			strcpy(control->sched[ns].obs[no].freq.inVal,freq);
			if (setOutlier==-1)
			  {
			    control->sched[ns].obs[no].outlierAmp.set=0;
			    control->sched[ns].obs[no].outlierProb.set=0;
			  }
			else
			  {
			    strcpy(control->sched[ns].obs[no].outlierAmp.inVal,outlierAmp);
			    strcpy(control->sched[ns].obs[no].outlierProb.inVal,outlierProb);
			    control->sched[ns].obs[no].outlierAmp.set=1;
			    control->sched[ns].obs[no].outlierProb.set=1;

			  }
			if (setTobs==1)
			    strcpy(control->sched[ns].obs[no].tobs.inVal,tobs);
			if (setHa==1)
			  {
			    strcpy(control->sched[ns].obs[no].ha.inVal,ha);
			    control->sched[ns].obs[no].ha.set=1;
			  }
			else
			  control->sched[ns].obs[no].ha.set=0;

			if (setStart==1)
			  {
			    strcpy(control->sched[ns].obs[no].start.inVal,start);
			    control->sched[ns].obs[no].start.set=1;
			  }
			else
			  control->sched[ns].obs[no].start.set=0;

			if (setFinish==1)
			  {
			    strcpy(control->sched[ns].obs[no].finish.inVal,finish);
			    control->sched[ns].obs[no].finish.set=1;
			  }
			else
			  control->sched[ns].obs[no].finish.set=0;

			control->sched[ns].obs[no].rcvrNum=-1;
			control->sched[ns].obs[no].beNum=-1;
			control->sched[ns].obs[no].obsSysNum=-1;


			if (setRcvr>-1)
			  control->sched[ns].obs[no].rcvrNum=setRcvr;
			if (setBE>-1)
			  control->sched[ns].obs[no].beNum=setBE;
			if (setObsSys>-1)			 
			  control->sched[ns].obs[no].obsSysNum=setObsSys;
			  
			(control->sched[ns].nObsSched)++;
		      }
		  }
	      }
	    else
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->sched[ns].obs[no].psrNum = i; break;}
		  }
		strcpy(control->sched[ns].obs[no].toaErr.inVal,toaErr);
		strcpy(control->sched[ns].obs[no].efac.inVal,efac);
		strcpy(control->sched[ns].obs[no].equad.inVal,equad);
		strcpy(control->sched[ns].obs[no].freq.inVal,freq);
		if (setOutlier==-1)
		  {
		    control->sched[ns].obs[no].outlierAmp.set=0;
		    control->sched[ns].obs[no].outlierProb.set=0;
		  }
		else
		  {
		    control->sched[ns].obs[no].outlierAmp.set=1;
		    control->sched[ns].obs[no].outlierProb.set=1;

		    strcpy(control->sched[ns].obs[no].outlierAmp.inVal,outlierAmp);
		    strcpy(control->sched[ns].obs[no].outlierProb.inVal,outlierProb);
		  }
		
		if (setTobs==1)
		  strcpy(control->sched[ns].obs[no].tobs.inVal,tobs);
		if (setHa==1)
		  {
		    strcpy(control->sched[ns].obs[no].ha.inVal,ha);
		    control->sched[ns].obs[no].ha.set=1;
		  }
		else
		  control->sched[ns].obs[no].ha.set=0;
		
		if (setStart==1)
		  {
		    strcpy(control->sched[ns].obs[no].start.inVal,start);
		    control->sched[ns].obs[no].start.set=1;
		  }
		else
		  control->sched[ns].obs[no].start.set=0;

		if (setFinish==1)
		  {
		    strcpy(control->sched[ns].obs[no].finish.inVal,finish);
		    control->sched[ns].obs[no].finish.set=1;
		  }
		else
		  control->sched[ns].obs[no].finish.set=0;

		control->sched[ns].obs[no].rcvrNum=-1;
		control->sched[ns].obs[no].beNum=-1;
		control->sched[ns].obs[no].obsSysNum=-1;
		
		if (setRcvr>-1)
		  control->sched[ns].obs[no].rcvrNum=setRcvr;
		if (setBE>-1)
		  control->sched[ns].obs[no].beNum=setBE;
		printf("SET OBSSYS %d\n",setObsSys);
		if (setObsSys>-1)
		  {

		    control->sched[ns].obs[no].obsSysNum=setObsSys;
		  }
		(control->sched[ns].nObsSched)++;
	      }
	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);

  // 
  //  printf("SCHED: %s %d %d\n",control->sched[ns].name,control->nSched,control->sched[ns].nObsSched);
  //  for (i=0;i<control->sched[ns].nObsSched;i++)
  //    {
  //      printf("SCHED: ..  %d %s\n",i,control->sched[ns].obs[i].psrName);
  //    }
  (control->nSched)++;
}

void readDefineFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int i,np;

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the 'define' information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</define>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"name:")==0)
	  strcpy(control->name,p[0].v);
	else if (strcmp(label,"nproc:")==0)
	  sscanf(p[0].v,"%d",&(control->nproc));
	else if (strcmp(label,"cut:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"label")==0)
		  strcpy(control->cutName[control->nCut],p[i].v);
		else if (strcmp(p[i].l,"mjd")==0)
		  sscanf(p[i].v,"%f",&(control->mjdCut[control->nCut]));
	      }
	    control->nCut++;
	  }	 
	else if (strcmp(label,"t2exe:")==0)
	  strcpy(control->t2exe,p[0].v);
	else if (strcmp(label,"shell:")==0)
	  strcpy(control->shell,p[0].v);
	else if (strcasecmp(label,"shellpth:")==0)
	  strcpy(control->shellPth,p[0].v);
	else if (strcmp(label,"nreal:")==0)
	  sscanf(p[0].v,"%d",&(control->nreal));
	else if (strcmp(label,"output:")==0)
	  {
	    int no = control->nOutput;
	    control->output[no].nAdd=0;
	    for (i=0;i<np;i++)
	      strcpy(control->output[no].label[i],p[i].v);
	    control->output[no].nAdd=np;
	    (control->nOutput)++;
	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
  
}

void readT2filesFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int i,np;

  printf("Reading T2 files from script\n");

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the 't2files' information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</t2files>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"clk:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  strcpy(control->useClock,p[i].v);
		else if (strcmp(p[i].l,"sim")==0)
		  strcpy(control->simClock,p[i].v);
	      }
	  }
	else if (strcmp(label,"eop:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  strcpy(control->useEOP,p[i].v);
		else if (strcmp(p[i].l,"sim")==0)
		  strcpy(control->simEOP,p[i].v);
	      }
	  }
	else if (strcmp(label,"ephem:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  {
		    strcpy(control->useEphem,p[i].v);
		    control->useTypeEphem=1;
		  }
		else if (strcmp(p[i].l,"sim")==0)
		  {
		    strcpy(control->simEphem,p[i].v);
		    control->simTypeEphem=1;
		  }
	      }
	  }
	else if (strcmp(label,"eph_file:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  {
		    strcpy(control->useEphem,p[i].v);
		    control->useTypeEphem=2;
		  }
		else if (strcmp(p[i].l,"sim")==0)
		  {
		    strcpy(control->simEphem,p[i].v);
		    control->simTypeEphem=2;
		  }
	      }
	  }
	else if (strcmp(label,"swm:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  strcpy(control->useSWM,p[i].v);
		else if (strcmp(p[i].l,"sim")==0)
		  strcpy(control->simSWM,p[i].v);
	      }
	  }
	else if (strcmp(label,"ne_sw:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"use")==0)
		  strcpy(control->useNE_SW,p[i].v);
		else if (strcmp(p[i].l,"sim")==0)
		  strcpy(control->simNE_SW,p[i].v);
	      }
	  }
	else
	  {
	    printf("line: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    printf("Exiting\n");
	    exit(1);
	  }
      }
  } while (endit==0);
  printf("Complete reading\n");
}


void readAdditionsFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int i,np;

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the 'add' information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</add>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	if (strcmp(label,"gwb:")==0)
	  {
	    control->gw[control->nGW].type=1;
	    strcpy(control->gw[control->nGW].alpha.inVal,"-0.666666");
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"amp")==0)
		  strcpy(control->gw[control->nGW].amp.inVal,p[i].v);
	      }	    
	    (control->nGW)++;
	  }
	else if (strcmp(label,"gwsingle:")==0)
	  {
	    control->gw[control->nGW].type=2;
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"ra")==0)
		  strcpy(control->gw[control->nGW].ra.inVal,p[i].v);
		else if (strcmp(p[i].l,"dec")==0)
		  strcpy(control->gw[control->nGW].dec.inVal,p[i].v);
		else if (strcmp(p[i].l,"ap")==0)
		  strcpy(control->gw[control->nGW].ap.inVal,p[i].v);
		else if (strcmp(p[i].l,"ac")==0)
		  strcpy(control->gw[control->nGW].ac.inVal,p[i].v);
		else if (strcmp(p[i].l,"gwm_amp")==0)
		  {strcpy(control->gw[control->nGW].gwmAmp.inVal,p[i].v); control->gw[control->nGW].type=3;}
		else if (strcmp(p[i].l,"gwm_epoch")==0)
		  {strcpy(control->gw[control->nGW].gwmEpoch.inVal,p[i].v); control->gw[control->nGW].type=3;}
		else if (strcmp(p[i].l,"gwm_phi")==0)
		  {strcpy(control->gw[control->nGW].gwmPhi.inVal,p[i].v); control->gw[control->nGW].type=3;}
		else if (strcmp(p[i].l,"gwcs_amp1")==0)
		  {strcpy(control->gw[control->nGW].gwcsAmp1.inVal,p[i].v); control->gw[control->nGW].type=6;}
		else if (strcmp(p[i].l,"gwcs_amp2")==0)
		  {strcpy(control->gw[control->nGW].gwcsAmp2.inVal,p[i].v); control->gw[control->nGW].type=6;}
		else if (strcmp(p[i].l,"gwcs_epoch")==0)
		  {strcpy(control->gw[control->nGW].gwcsEpoch.inVal,p[i].v); control->gw[control->nGW].type=6;}
		else if (strcmp(p[i].l,"gwcs_width")==0)
		  {strcpy(control->gw[control->nGW].gwcsWidth.inVal,p[i].v); control->gw[control->nGW].type=6;}
		else if (strcmp(p[i].l,"cgw_freq")==0)
		  {strcpy(control->gw[control->nGW].cgw_freq.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"cgw_h0")==0)
		  {strcpy(control->gw[control->nGW].cgw_h0.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"cgw_epoch")==0)
		  {strcpy(control->gw[control->nGW].cgw_epoch.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"cgw_cosinc")==0)
		  {strcpy(control->gw[control->nGW].cgw_cosinc.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"cgw_angpol")==0)
		  {strcpy(control->gw[control->nGW].cgw_angpol.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"cgw_mc")==0)
		  {strcpy(control->gw[control->nGW].cgw_mc.inVal,p[i].v); control->gw[control->nGW].type=4;}
		else if (strcmp(p[i].l,"file")==0) // Type = 5 - file listing of source parameters
		  {strcpy(control->gw[control->nGW].fname,p[i].v); control->gw[control->nGW].type=5;}
	      }	    
	    (control->nGW)++;
	  }
	else if (strcmp(label,"clknoise:")==0)
	  {
	    control->nClkNoise=1;
	    control->clkNoise.gwAmp.set=0;
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"alpha")==0)
		  strcpy(control->clkNoise.alpha.inVal,p[i].v);
		else if (strcmp(p[i].l,"p0")==0)
		  strcpy(control->clkNoise.p0.inVal,p[i].v);
		else if (strcmp(p[i].l,"fc")==0)
		  strcpy(control->clkNoise.fc.inVal,p[i].v);
		else if (strcmp(p[i].l,"gwamp")==0)
		  {
		    strcpy(control->clkNoise.gwAmp.inVal,p[i].v);
		    control->clkNoise.gwAmp.set=1;
		  }
	      }
	  }
	else if (strcmp(label,"ephemnoise:")==0)
	  {
	    control->nEphemNoise=1;
	    control->ephemNoise.gwAmp.set=0;
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"alpha")==0)
		  strcpy(control->ephemNoise.alpha.inVal,p[i].v);
		else if (strcmp(p[i].l,"p0")==0)
		  strcpy(control->ephemNoise.p0.inVal,p[i].v);
		else if (strcmp(p[i].l,"fc")==0)
		  strcpy(control->ephemNoise.fc.inVal,p[i].v);
		else if (strcmp(p[i].l,"gwamp")==0)
		  {
		    strcpy(control->ephemNoise.gwAmp.inVal,p[i].v);
		    control->ephemNoise.gwAmp.set=1;
		  }
	      }
	  }
	else if (strcmp(label,"tnoise:")==0)
	  {
	    int nt = control->nTnoise;
	    char alpha[1024],p0[1024],fc[1024],beta[1024];
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;
	    char idLabel[1024];
	    int setIDlabel=0;

	    strcpy(beta,"0");

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"label")==0)
		  {strcpy(idLabel,p[i].v); setIDlabel=1;}
		else if (strcmp(p[i].l,"gwamp")==0)
		  {
		    strcpy(alpha,"gwamp_auto");
		    strcpy(p0,p[i].v);
		    strcpy(fc,"gwamp_auto");
		  }
		else if (strcmp(p[i].l,"alpha")==0)
		  strcpy(alpha,p[i].v);
		else if (strcmp(p[i].l,"beta")==0)
		  strcpy(beta,p[i].v);
		else if (strcmp(p[i].l,"p0")==0)
		  {strcpy(p0,p[i].v);}
		else if (strcmp(p[i].l,"fc")==0)
		  strcpy(fc,p[i].v);
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			nt = control->nTnoise;
			if (setIDlabel==0)
			  strcpy(control->tnoise[nt].label,"UNSET");
			else
			  strcpy(control->tnoise[nt].label,idLabel);

			control->tnoise[nt].psrNum = i;
			strcpy(control->tnoise[nt].alpha.inVal,alpha);
			strcpy(control->tnoise[nt].beta.inVal,beta);
			strcpy(control->tnoise[nt].p0.inVal,p0);
			strcpy(control->tnoise[nt].fc.inVal,fc);
 			(control->nTnoise)++;
		      }
		  }
	      }
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->tnoise[nt].psrNum = i; break;}
		  }
		if (setIDlabel==0)
		  strcpy(control->tnoise[nt].label,"UNSET");
		else
		  strcpy(control->tnoise[nt].label,idLabel);

		strcpy(control->tnoise[nt].alpha.inVal,alpha);
		strcpy(control->tnoise[nt].beta.inVal,beta);
		strcpy(control->tnoise[nt].p0.inVal,p0);
		strcpy(control->tnoise[nt].fc.inVal,fc);
		(control->nTnoise)++;	
	      }
	  }
	else if (strcmp(label,"planet:")==0)
	  {
	    int npl = control->nPlanets;

	    char pb[1024],ecc[1024],a1[1024],t0[1024],om[1024];
	    char pname[1024]="unknown";
	    char label[1024];
	    char idLabel[1024];
	    int setIDlabel=0;
	    int setlabel=0;


	    for (i=0;i<np;i++)
	      {
		printf("planets: test %s %s\n",p[i].l,p[i].v);
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"label")==0)
		  {strcpy(idLabel,p[i].v); setIDlabel=1;}
		else if (strcmp(p[i].l,"pb")==0)
		  strcpy(pb,p[i].v);
		else if (strcmp(p[i].l,"ecc")==0)
		  strcpy(ecc,p[i].v);
		else if (strcmp(p[i].l,"a1")==0)
		  {strcpy(a1,p[i].v);}
		else if (strcmp(p[i].l,"t0")==0)
		  strcpy(t0,p[i].v);
		else if (strcmp(p[i].l,"om")==0)
		  strcpy(om,p[i].v);
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {

		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			npl = control->nPlanets;
			if (setIDlabel==0)
			  strcpy(control->planets[npl].label,"UNSET");
			else
			  strcpy(control->planets[npl].label,idLabel);
			control->planets[npl].psrNum = i;
			strcpy(control->planets[npl].pb.inVal,pb);
			strcpy(control->planets[npl].ecc.inVal,ecc);
			strcpy(control->planets[npl].a1.inVal,a1);
			strcpy(control->planets[npl].t0.inVal,t0);
			strcpy(control->planets[npl].om.inVal,om);
 			(control->nPlanets)++;
		      }
		  }
	      }
	    else 
	      {
		printf("planets: npsr = %d\n",control->npsr);
		for (i=0;i<control->npsr;i++)
		  {
		    printf("planets: checking >%s< >%s<\n",control->psr[i].name,pname);
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->planets[npl].psrNum = i; break;}
		  }
		control->planets[npl].psrNum = i;
		if (setIDlabel==0)
		  strcpy(control->planets[npl].label,"UNSET");
		else
		  strcpy(control->planets[npl].label,idLabel);

		strcpy(control->planets[npl].pb.inVal,pb);
		strcpy(control->planets[npl].ecc.inVal,ecc);
		strcpy(control->planets[npl].a1.inVal,a1);
		strcpy(control->planets[npl].t0.inVal,t0);
		strcpy(control->planets[npl].om.inVal,om);
		(control->nPlanets)++;	
	      }
	  }
	else if (strcmp(label,"glitch:")==0)
	  {
	    int ng = control->nGlitches;
	    char glep[1024],glph[1024],glf0[1024],glf1[1024],glf0d[1024],gltd[1024];
	    int setGlep=0;
	    int setGlph=0;
	    int setGlf0=0;
	    int setGlf1=0;
	    int setGlf0d=0;
	    int setGltd=0;
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"glep")==0)
		  {strcpy(glep,p[i].v); setGlep=1;}
		else if (strcmp(p[i].l,"glf0")==0)
		  {strcpy(glf0,p[i].v); setGlf0=1;}
		else if (strcmp(p[i].l,"glf1")==0)
		  {strcpy(glf1,p[i].v); setGlf1=1;}
		else if (strcmp(p[i].l,"glph")==0)
		  {strcpy(glph,p[i].v); setGlph=1;}
		else if (strcmp(p[i].l,"glf0d")==0)
		  {strcpy(glf0d,p[i].v); setGlf0d=1;}
		else if (strcmp(p[i].l,"gltd")==0)
		  {strcpy(gltd,p[i].v); setGltd=1;}
	      }
	    printf("ADDING GLITCH\n");
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			ng = control->nGlitches;
			control->glitches[ng].psrNum = i;
			if (setGlep==1)
			  {
			    strcpy(control->glitches[ng].glep.inVal,glep);
			    control->glitches[ng].glep.set = 1;
			  }
			else
			  control->glitches[ng].glep.set = 0;
			if (setGlf0==1)
			  {
			    strcpy(control->glitches[ng].glf0.inVal,glf0);
			    control->glitches[ng].glf0.set = 1;
			  }
			else
			  control->glitches[ng].glf0.set = 0;
			  
			if (setGlf1==1)
			  {
			    strcpy(control->glitches[ng].glf1.inVal,glf1);
			    control->glitches[ng].glf1.set = 1;
			  }
			else
			  control->glitches[ng].glf1.set = 0;

			if (setGlph==1)
			  {
			    strcpy(control->glitches[ng].glph.inVal,glph);
			    control->glitches[ng].glph.set = 1;
			  }
			else
			  control->glitches[ng].glph.set = 0;

			if (setGlf0d==1)
			  {
			    strcpy(control->glitches[ng].glf0d.inVal,glf0d);
			    control->glitches[ng].glf0d.set = 1;
			  }
			else
			  control->glitches[ng].glf0d.set = 0;

			if (setGltd==1)
			  {
			    strcpy(control->glitches[ng].gltd.inVal,glep);
			    control->glitches[ng].gltd.set = 1;
			  }
			else
			  control->glitches[ng].gltd.set = 0;


 			(control->nGlitches)++;
		      }
		  }
	      } 
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->glitches[ng].psrNum = i; break;}
		  }
		if (setGlep==1)
		  {
		    strcpy(control->glitches[ng].glep.inVal,glep);
		    control->glitches[ng].glep.set = 1;
		  }
		else
		  control->glitches[ng].glep.set = 0;
		if (setGlf0==1)
		  {
		    strcpy(control->glitches[ng].glf0.inVal,glf0);
		    control->glitches[ng].glf0.set = 1;
		  }
		else
		  control->glitches[ng].glf0.set = 0;
		
		if (setGlf1==1)
		  {
		    strcpy(control->glitches[ng].glf1.inVal,glf1);
		    control->glitches[ng].glf1.set = 1;
		  }
		else
		  control->glitches[ng].glf1.set = 0;
		
		if (setGlph==1)
		  {
		    strcpy(control->glitches[ng].glph.inVal,glph);
		    control->glitches[ng].glph.set = 1;
		  }
		else
		  control->glitches[ng].glph.set = 0;
		
		if (setGlf0d==1)
		  {
		    strcpy(control->glitches[ng].glf0d.inVal,glf0d);
		    control->glitches[ng].glf0d.set = 1;
		  }
		else
		  control->glitches[ng].glf0d.set = 0;
		
		if (setGltd==1)
		  {
		    strcpy(control->glitches[ng].gltd.inVal,glep);
		    control->glitches[ng].gltd.set = 1;
		  }
		else
		  control->glitches[ng].gltd.set = 0;
		
		
		(control->nGlitches)++;
		
		
	      }
	  }
	else if (strcmp(label,"dmvar:")==0)
	  {
	    int ndm = control->nDMvar;
	    char d_tscale[1024],dVal[1024],refFreq[1024];
	    int setD=0;
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"D")==0)
		  {
		    char str[MAX_STRLEN];
		    char *tok;

		    setD=1;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(d_tscale,tok);
		    tok = strtok(NULL,";]");
		    strcpy(refFreq,tok);
		    tok = strtok(NULL,";]");
		    strcpy(dVal,tok);
		  }
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			ndm = control->nDMvar;
			control->dmVar[ndm].psrNum = i;
			if (setD==1)
			  {
			    strcpy(control->dmVar[ndm].d_tscale.inVal,d_tscale);
			    strcpy(control->dmVar[ndm].dVal.inVal,dVal);
			    strcpy(control->dmVar[ndm].refFreq.inVal,refFreq);
			    control->dmVar[ndm].type=1;
			  }
 			(control->nDMvar)++;
		      }
		  }
	      } 
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->dmVar[ndm].psrNum = i; break;}
		  }
		strcpy(control->dmVar[ndm].d_tscale.inVal,d_tscale);
		strcpy(control->dmVar[ndm].dVal.inVal,dVal);
		strcpy(control->dmVar[ndm].refFreq.inVal,refFreq);
		control->dmVar[ndm].type=1;
		
		(control->nDMvar)++;	
	      }
	  }
	else if (strcmp(label,"dmCovar:")==0)
	  {
	    int ndm = control->nDMcovar;
	    char alpha[1024],a[1024],b[1024];
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"alpha")==0)
		  {strcpy(alpha,p[i].v);}
		else if (strcmp(p[i].l,"a")==0)
		  {strcpy(a,p[i].v);}
		else if (strcmp(p[i].l,"b")==0)
		  {strcpy(b,p[i].v);}
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			ndm = control->nDMcovar;
			control->dmCovar[ndm].psrNum = i;
			strcpy(control->dmCovar[ndm].alpha.inVal,alpha);
			strcpy(control->dmCovar[ndm].a.inVal,a);
			strcpy(control->dmCovar[ndm].b.inVal,b);
			control->dmCovar[ndm].type=1;

 			(control->nDMcovar)++;
		      }
		  }
	      }
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->dmCovar[ndm].psrNum = i; break;}
		  }
		printf("Setting dmCovar: %i %s %s %s \n",control->dmCovar[ndm].psrNum,alpha,a,b);
		strcpy(control->dmCovar[ndm].alpha.inVal,alpha);
		strcpy(control->dmCovar[ndm].a.inVal,a);
		strcpy(control->dmCovar[ndm].b.inVal,b);
		control->dmCovar[ndm].type=1;
		
		(control->nDMcovar)++;	
	      }

	  }
	else if (strcmp(label,"dmFunc:")==0)
	  {
	    int ndm = control->nDMfunc;
	    char ddm[1024];
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"ddm")==0)
		  {strcpy(ddm,p[i].v);}
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			ndm = control->nDMfunc;
			control->dmFunc[ndm].psrNum = i;
			strcpy(control->dmFunc[ndm].ddm.inVal,ddm);
 			(control->nDMfunc)++;
		      }
		  }
	      }
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->dmFunc[ndm].psrNum = i; break;}
		  }
		strcpy(control->dmFunc[ndm].ddm.inVal,ddm);
		(control->nDMfunc)++;	
	      }

	  }
	else if (strcmp(label,"jitter:")==0)
	  {
	    int nJitter = control->nJitter;
	    char d_tscale[1024],sigma_j[1024],refFreq[1024];
	    int setSJ=0;
	    char pname[1024];
	    char label[1024];
	    int setlabel=0;

	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"psr")==0)
		  strcpy(pname,p[i].v);
		else if (strcmp(p[i].l,"psrLabel")==0)
		  {strcpy(label,p[i].v); setlabel=1;}
		else if (strcmp(p[i].l,"SJ")==0)
		  {
		    char str[MAX_STRLEN];
		    char *tok;

		    setSJ=1;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(d_tscale,tok);
		    tok = strtok(NULL,";]");
		    strcpy(refFreq,tok);
		    tok = strtok(NULL,";]");
		    strcpy(sigma_j,tok);
		  }
	      }
	    if (strcmp(pname,"all")==0 || setlabel==1)
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (setlabel==0 || strcmp(control->psr[i].label,label)==0)
		      {
			nJitter = control->nJitter;
			control->jitter[nJitter].psrNum = i;
			if (setSJ==1)
			  {
			    strcpy(control->jitter[nJitter].t0.inVal,d_tscale);
			    strcpy(control->jitter[nJitter].sigma_j.inVal,sigma_j);
			    strcpy(control->jitter[nJitter].refFreq.inVal,refFreq);
			    control->jitter[nJitter].type=1;
			  }
 			(control->nJitter)++;
		      }
		  }
	      }
	    else 
	      {
		for (i=0;i<control->npsr;i++)
		  {
		    if (strcmp(control->psr[i].name,pname)==0)
		      {control->jitter[nJitter].psrNum = i; break;}
		  }
		if (setSJ==1)
		  {
		    strcpy(control->jitter[nJitter].t0.inVal,d_tscale);
		    strcpy(control->jitter[nJitter].sigma_j.inVal,sigma_j);
		    strcpy(control->jitter[nJitter].refFreq.inVal,refFreq);
		    control->jitter[nJitter].type=1;
		  }
	      
		(control->nJitter)++;	
	      }
	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	    exit(1);
	  }
      }
  } while (endit==0);
    
}

void readPulsarsFromScript(controlStruct *control,FILE *fin)
{
  int endit=0;
  char line[1024];
  char trimLine[1024];
  char label[1024];
  paramStruct p[MAX_LINE_PARAMS];
  int np,i,npsr;

  do {
    if (fgets(line,1024,fin)==NULL)
      {
	printf("Error: script ended whilst still reading the pulsar information\n");
	fclose(fin);
	finishOff(control);
      }
    trimwhitespace(trimLine,1024,line);
    if (trimLine[0]=='#' || strlen(trimLine)==0)
      {
	// Do nothing
      }
    else if (strcmp(trimLine,"</pulsars>")==0)
      endit=1;
    else
      {
	np = getParams(trimLine,label,p);
	printf("processing: %s\n",label);
	if (strcmp(label,"psr:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		printf("... %s\n",p[i].l);
		if (strcmp(p[i].l,"name")==0){
		  initialisePulsar(&(control->psr[control->npsr]));
		  strcpy(control->psr[control->npsr].name,p[i].v);
		  control->psr[control->npsr].setName=1;
		  control->psr[control->npsr].requireCatRead=1;
		}
		else if (strcmp(p[i].l,"label")==0)
		  strcpy(control->psr[control->npsr].label,p[i].v);
		else if (strcmp(p[i].l,"profileFile")==0)
		  {
		    int nf = control->psr[control->npsr].nProfileFile;
		    char str[MAX_STRLEN];
		    char *tok;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(control->psr[control->npsr].freqProfileFile[nf].inVal,tok);
		    tok = strtok(NULL,";]");
		    strcpy(control->psr[control->npsr].profileFile[nf],tok);
		    (control->psr[control->npsr].nProfileFile)++;
		  }
		else if (strcmp(p[i].l,"profileEqn")==0)
		  {
		    strcpy(control->psr[control->npsr].profileEqn,p[i].v);
		    control->psr[control->npsr].setProfileEqn=1;
		  }
		else if (strcmp(p[i].l,"flux")==0)
		  {
		    int nf = control->psr[control->npsr].nFlux;
		    char str[MAX_STRLEN];
		    char *tok;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(control->psr[control->npsr].freqFlux[nf].inVal,tok);
		    tok = strtok(NULL,";]");
		    strcpy(control->psr[control->npsr].flux[nf].inVal,tok);
		    (control->psr[control->npsr].nFlux)++;
		  }
		else if (strcmp(p[i].l,"tsky")==0)
		  {
		    int nt = control->psr[control->npsr].nTsky;
		    char str[MAX_STRLEN];
		    char *tok;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(control->psr[control->npsr].freqTsky[nt].inVal,tok);
		    tok = strtok(NULL,";]");
		    strcpy(control->psr[control->npsr].tsky[nt].inVal,tok);
		    (control->psr[control->npsr].nTsky)++;
		  }
		else if (strcmp(p[i].l,"diff_df")==0)
		  {
		    char str[MAX_STRLEN];
		    char *tok;
		    printf("a1\n");
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(control->psr[control->npsr].diff_dfFreq.inVal,tok);
		    tok = strtok(NULL,";]");
		    strcpy(control->psr[control->npsr].diff_df.inVal,tok);
		    control->psr[control->npsr].setDiff_df=1;
		    printf("a2\n");
		  }
		else if (strcmp(p[i].l,"diff_ts")==0)
		  {
		    char str[MAX_STRLEN];
		    char *tok;
		    strcpy(str,p[i].v+1);
		    tok = strtok(str,";]");
		    strcpy(control->psr[control->npsr].diff_tsFreq.inVal,tok);
		    tok = strtok(NULL,";]");
		    strcpy(control->psr[control->npsr].diff_ts.inVal,tok);
		    control->psr[control->npsr].setDiff_ts=1;
		  }
		else
		  {
		    int npv;
		    npsr = control->npsr;
		    printf("In here\n");
		    npv = control->psr[npsr].nSetParam;
		    strcpy(control->psr[npsr].setParamName[npv],p[i].l);
		    strcpy(control->psr[npsr].paramVal[npv].inVal,p[i].v);
		    control->psr[npsr].paramVal[npv].constant = 0;
		    (control->psr[npsr].nSetParam)++;
		  }
	      }

	    (control->npsr)++;
	  }
	else if (strcmp(label,"ephem:")==0)
	  {
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"name")==0){

		  initialisePulsar(&(control->psr[control->npsr]));
		  strcpy(control->psr[control->npsr].ephem,p[i].v);
		  control->psr[control->npsr].setEphem=1;
		  control->psr[control->npsr].requireEphemRead=1;
		}
		else if (strcmp(p[i].l,"label")==0)
		  strcpy(control->psr[control->npsr].label,p[i].v);
	      }
	    (control->npsr)++;

	  }
	else if (strcmp(label,"createpsr:")==0)
	  {
	    int npsr = control->npsr;
	    int npv;
	    int j;
	    int npsrReq=1;
	    // Search for multiple pulsar request
	    for (i=0;i<np;i++)
	      {
		if (strcmp(p[i].l,"n")==0)
		  sscanf(p[i].v,"%d",&npsrReq);
	      }
	    for (j=0;j<npsrReq;j++)
	      {
		npsr = control->npsr;
		for (i=0;i<np;i++)
		  {
		    if (strcmp(p[i].l,"name")==0)
		      {
			strcpy(control->psr[npsr].name,p[i].v);
			control->psr[npsr].setName=1;
		      }
		    else if (strcmp(p[i].l,"label")==0)
		      {
			strcpy(control->psr[npsr].label,p[i].v);
			control->psr[npsr].setLabel=1;
		      }
		    else if (strcmp(p[i].l,"pos")==0)
		      {
			printf("Setting position %s %d\n",p[i].v,npsr);
			//			exit(1);
			if (strstr(p[i].v,"isotropic(")!=0)
			  {
			    double theta,phi;
			    char raS[1024],decS[1024];
			    int inp=0;
			    char iparam[10][128]; // SHOULD SET PROPERLY
			    char *itok;
			    char t[1024];

			    int notGood,ii;
			    double cval;

			    strcpy(t,p[i].v);
			    // Get parameters
			    itok = strtok(t,"(");
			    while ((itok = strtok(NULL,";)"))!=NULL)
			      {
				strcpy(iparam[inp],itok);
				inp++;
			      }
			    //			    printf("Have %d parameters\n",inp);
			    do {
			      theta = (acos((TKranDev(&(control->seed))-0.5)*2)-M_PI/2.0)/(2.0*M_PI);
			      phi = TKranDev(&(control->seed));
			      notGood=0;
			      for (ii=0;ii<inp;ii++)
				{
				  if (strstr(iparam[ii],"decj")!=NULL)
				    {
				      if ((itok = strstr(iparam[ii],"<"))!=NULL)
					{
					  sscanf(itok+1,"%lf",&cval);
					  if (theta*360.0 >= cval)
					    notGood=1;
					}
				      if ((itok = strstr(iparam[ii],">"))!=NULL)
					{
					  sscanf(itok+1,"%lf",&cval);
					  if (theta*360.0 <= cval)
					    notGood=1;
					}
				    }
				  if (strstr(iparam[ii],"raj")!=NULL)
				    {
				      if ((itok = strstr(iparam[ii],"<"))!=NULL)
					{
					  sscanf(itok+1,"%lf",&cval);
					  if (phi*360.0 >= cval)
					    notGood=1;
					}
				      if ((itok = strstr(iparam[ii],">"))!=NULL)
					{
					  sscanf(itok+1,"%lf",&cval);
					  if (phi*360.0 <= cval)
					    notGood=1;
					}
				    }
				}
			    } while (notGood==1);
			    turn_hms(phi,raS);
			    turn_dms(theta,decS);
			    printf("Fixed isotropic distribution %g\n",theta);
			    printf("raS and decS = %s %s\n",raS,decS);
			    strcpy(control->psr[npsr].rajStr,raS);
			    npv = control->psr[npsr].nSetParam;
			    strcpy(control->psr[npsr].setParamName[npv],"RAJ");
			    strcpy(control->psr[npsr].paramVal[npv].inVal,raS);
			    (control->psr[npsr].nSetParam)++;

			    strcpy(control->psr[npsr].decjStr,decS);
			    npv = control->psr[npsr].nSetParam;
			    strcpy(control->psr[npsr].setParamName[npv],"DECJ");
			    strcpy(control->psr[npsr].paramVal[npv].inVal,decS);
			    (control->psr[npsr].nSetParam)++;
			  }
			else
			  {
			    printf("Unknown source position parameter: %s\n",p[i].v);
			    finishOff(control);
			  }
		      }
		    else if (strcmp(p[i].l,"n")==0)
		      {  // Do nothing
		      }
		    else if (strcmp(p[i].l,"fix")==0)
		      { 
			int ii,found=-1;
			printf("Fixing: %s\n",p[i].v);
			// Find parameter
			for (ii=0;ii<control->psr[npsr].nSetParam;ii++)
			  {
			    if (strcmp(control->psr[npsr].setParamName[ii],p[i].v)==0)
			      {found = ii; break;}
			  }
			if (found == -1)
			  {
			    printf("fixing a parameter that doesn't exist: %s\n",p[i].v);
			    finishOff(control);
			  }
			control->psr[npsr].paramVal[found].constant = 1;
		      }
		    else
		      {
			npv = control->psr[npsr].nSetParam;
			strcpy(control->psr[npsr].setParamName[npv],p[i].l);
			strcpy(control->psr[npsr].paramVal[npv].inVal,p[i].v);
			control->psr[npsr].paramVal[npv].constant = 0;
			(control->psr[npsr].nSetParam)++;
		      }
		  }
		(control->npsr)++;
	      }
	  }
	else
	  {
	    printf("psr: %s\n",trimLine);
	    printf("Unkown label: %s\n",label);
	    for (i=0;i<np;i++)
	      printf(" ... [%s] %d %s %s\n",label,p[i].type,p[i].l,p[i].v);
	  }
      }
  } while (endit==0);
}

void initialisePulsar(psrStruct *psr)
{
  psr->setName = 0;
  psr->setEphem = 0;
  psr->setLabel = 0;
  psr->nSetParam=0;
  psr->nToAs = 0;
  psr->nFlux = 0;
  psr->nTsky = 0;
  psr->requireCatRead=0;
  psr->requireEphemRead=0;
  psr->setProfileEqn=0;
  psr->nProfileFile=0;
  psr->setDiff_df=0;
  psr->setDiff_ts=0;

}

int getParams(char *line,char *label,paramStruct *p)
{
  char *tok;
  char *new;
  char *tok2,*new2;
  char trim[1024],trim2[1024];
  int n=0;

  strcpy(label,strtok_r(line," ",&tok));
  // Read label

  // Read parameter
  while ((new = strtok_r(NULL,",\n",&tok))!=NULL)
    {
      trimwhitespace(trim,1024,new);
      // Search for an equal sign
      if (strstr(trim,"=")==NULL)
	{
	  strcpy(p[n].v,trim);
	  p[n].type=1;
	  n++;
	}
      else
	{
	  new2 = strtok_r(trim,"=",&tok2);
	  trimwhitespace(trim2,1024,new2);
	  strcpy(p[n].l,trim2);
	  new2 = strtok_r(NULL,"",&tok2);
	  trimwhitespace(trim2,1024,new2);
	  strcpy(p[n].v,trim2);
	  p[n].type = 2;
	  n++;
	}
    } 
  return n;
}

void loadInputs(controlStruct *control,int argc,char *argv[])
{
  if (argc!=2)
    {
      printf("Usage: ptaSimulate scriptName\n");
      finishOff(control);
    }
  strcpy(control->inputScript,argv[1]);
}

void finishOff(controlStruct *control)
{
  free(control);
  exit(1);
}

void initialiseControl(controlStruct *control)
{
  int i;
  control->nCut=0;
  control->nOutput=1;
  strcpy(control->output[0].label[0],"DEFAULT");
  control->output[0].nAdd=1;
  control->npsr=0;
  control->nObsRun=0;
  control->nSched=0;
  control->nTnoise=0;
  control->nPlanets=0;
  control->nClkNoise=0;
  control->nEphemNoise=0;
  control->nObsSys=0;
  control->nGlitches=0;
  control->nGW=0;
  control->nBE=0;
  control->nDMvar=0;
  control->nDMcovar=0;
  control->nDMfunc=0;
  control->nJitter=0;
  control->seed = TKsetSeed();
  control->nreal = 1;
  control->minT = 99999;
  control->maxT = 0;
  for (i=0;i<MAX_GWS;i++)
    {
      strcpy(control->gw[i].ap.inVal,"0");
      strcpy(control->gw[i].ac.inVal,"0");
    }
  strcpy(control->simEphem,"DE421");
  strcpy(control->useEphem,"DE421");
  control->simTypeEphem = 1;
  control->useTypeEphem = 1;

  strcpy(control->simClock,"TT(TAI)");
  strcpy(control->useClock,"TT(TAI)");
  strcpy(control->simEOP,"/earth/eopc04_IAU2000.62-now");
  strcpy(control->useEOP,"/earth/eopc04_IAU2000.62-now");
  strcpy(control->t2exe,"tempo2");
  strcpy(control->shell,"tcsh");
  strcpy(control->shellPth,"/usr/bin/");
  strcpy(control->useSWM,"DEFAULT");
  strcpy(control->simSWM,"DEFAULT");
  strcpy(control->useNE_SW,"4.0");
  strcpy(control->simNE_SW,"4.0");
}

// Stores the trimmed input string into the given output buffer, which must be
// large enough to store the result.  If it is too small, the output is
// truncated.
size_t trimwhitespace(char *out, size_t len, const char *str)
{
  if(len == 0)
    return 0;

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    {
      *out = 0;
      return 1;
    }

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;
  // Copy trimmed string and add null terminator
  memcpy(out, str, out_size);
  out[out_size] = 0;

  return out_size;
}

int turn_hms(double turn, char *hms){
 
  /* Converts double turn to string " hh:mm:ss.ssss" */
  
  int hh, mm, isec;
  double sec;

  hh = turn*24.;
  mm = (turn*24.-hh)*60.;
  sec = ((turn*24.-hh)*60.-mm)*60.;
  isec = (sec*10000. +0.5)/10000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        hh=hh+1;
        if(hh==24){
          hh=0;
        }
      }
    }

    //  sprintf(hms," %02d:%02d:%010.7f",hh,mm,sec);
  sprintf(hms," %02d:%02d:%05.2f",hh,mm,sec);
 
}

int turn_dms(double turn, char *dms){
  
  /* Converts double turn to string "sddd:mm:ss.sss" */
  
  int dd, mm, isec;
  double trn, sec;
  char sign;
  
  sign=' ';
  if (turn < 0.){
    sign = '-';
    trn = -turn;
  }
  else{
    sign = '+';
    trn = turn;
  }
  dd = trn*360.;
  mm = (trn*360.-dd)*60.;
  sec = ((trn*360.-dd)*60.-mm)*60.;
  isec = (sec*1000. +0.5)/1000;
    if(isec==60){
      sec=0.;
      mm=mm+1;
      if(mm==60){
        mm=0;
        dd=dd+1;
      }
    }
    //  sprintf(dms,"%c%02d:%02d:%010.7f",sign,dd,mm,sec);
  sprintf(dms,"%c%02d:%02d:%05.2f",sign,dd,mm,sec);
 
}

void createDirectoryStructure(controlStruct *control)
{
  char resDir[MAX_STRLEN];
  char dir[MAX_STRLEN];
  int i,j,k;
  FILE *fout;

  // Make filenames for output directories
  for (j=0;j<control->nOutput;j++)
    {
      if (j==0)
	strcpy(control->output[j].fname,"");
      else
	{
	  strcpy(control->output[j].fname,"");
	  for (k=0;k<control->output[j].nAdd;k++)
	    {
	      strcat(control->output[j].fname,control->output[j].label[k]);
	      if (k<control->output[j].nAdd-1)
		strcat(control->output[j].fname,"_");
	    }
	  printf("Made filename: %s\n",control->output[j].fname);
	}
    }


  printf("Creating directory 1\n");
  sprintf(resDir,control->name);
  mkdir(resDir,0700);
  sprintf(dir,"%s/%s",resDir,"output");
  mkdir(dir,0700);
  sprintf(dir,"%s/%s",resDir,"scripts");
  mkdir(dir,0700);
  sprintf(dir,"%s/scripts/status",resDir);
  mkdir(dir,0700);
  printf("Creating directory 2\n");
  for (i=0;i<control->nreal;i++)
    {
      sprintf(dir,"%s/output/real_%d",resDir,i);
      mkdir(dir,0700);
      
      for (j=1;j<control->nOutput;j++)
	{
	  sprintf(dir,"%s/output/real_%d/%s",resDir,i,control->output[j].fname);
	  mkdir(dir,0700);
	}
    }
  printf("Creating directory 3\n");
  sprintf(dir,"%s/%s",resDir,"workFiles");
  mkdir(dir,0700);
  sprintf(dir,"%s/workFiles/common",resDir,i);
  mkdir(dir,0700);
  printf("Creating directory 4\n");
  for (i=0;i<control->nreal;i++)
    {
      sprintf(dir,"%s/workFiles/real_%d",resDir,i);
      mkdir(dir,0700);
    }
  printf("Creating directory 5\n");
  sprintf(dir,"%s/%s",resDir,"setup");
  mkdir(dir,0700);
  sprintf(dir,"%s/setup/useParams",resDir);

  //  printf("Creating directory61\n");  
  //  fout = fopen(dir,"w");
  //  fprintf(fout,"Parameters used by ptaSimulate\n");
  //  fclose(fout);
  printf("Creating directory 7\n");
}

void createDMvar(controlStruct *control,int r)
{
  int i,nit,j,p;
  char fname[MAX_STRLEN];
  double globalParameter;
  long double result;

  double secperyear=365*86400.0;
  // my parameters
  double alpha= -8.0/3.0;
  int npts=1024;
  float D_d=1; // us
  float ref_freq=1400; // MHz
  float d=1000;
  char writeTextFiles=0;
  double lastMJD=1e99;
  char name[1024];
  
  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double offsets[MAX_TOAS]; // Will change to doubles - should use malloc
  double dms[MAX_TOAS]; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  float beta=0;
  int dd;
  
  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  nit = 1;

  for (dd=0;dd<control->nDMvar;dd++)
    {
      p = control->dmVar[dd].psrNum;
      D_d = control->dmVar[dd].dVal.dval;
      ref_freq = control->dmVar[dd].refFreq.dval;
      d = control->dmVar[dd].d_tscale.dval;

      header = toasim_init_header();
      strcpy(header->short_desc,"addDmVar");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = nit;
      
      sprintf(fname,"%s/workFiles/real_%d/%s.dmvar.%d",control->name,r,control->psr[control->dmVar[dd].psrNum].name,dd);
      // First we write the header...
      file = toasim_write_header(header,fname);
      
      double mjd_start=1000000.0;
      double mjd_end=-10000000.0;
      for (j=0;j<control->psr[p].nToAs;j++){
	// find the start and end times
	if(control->psr[p].obs[j].sat < mjd_start)mjd_start=(double)control->psr[p].obs[j].sat;
	if(control->psr[p].obs[j].sat > mjd_end)mjd_end=(double)control->psr[p].obs[j].sat;
      }
      
      
      printf("start    = %f (mjd)\n",mjd_start);
      printf("end      = %f (mjd)\n",mjd_end  );
      printf("npts     = %d (days)\n",npts     );
      printf("D_d(%f)  = %f (us^2)\n",d,D_d    );
      printf("ref_freq = %f (MHz)\n",ref_freq    );     
      printf("seed     = %d\n",control->seed);
      
      
      D_d *=1e-12; // convert us to seconds
      d*=86400.0;  // convert days to seconds.
            
      double pism = 0.0112 * D_d * pow(d,(-5.0/3.0)) * pow(secperyear,-1.0/3.0);
      
      printf("pism(1yr)  = %g (yr^3) \n",pism );
                  
      double yr2dm = secperyear * DM_CONST*pow(ref_freq,2.0);
      printf("yr2dm   = %f \n",yr2dm    );
      
      pism *= pow(yr2dm,2); // convert yr^3 to yr.cm^-3.pc
      
      printf("pism(1yr)  = %g ((cm^-3pc)^2 yr) \n",pism );
      
      printf("\n");
      printf("Generating red noise...\n");
      
      rednoisemodel_t* model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,pism,alpha,beta);
      populateRedNoiseModel(model,&(control->seed));
            
      int itjmp=nit/50;
      if (itjmp<1)itjmp=1;
      int dots=0;
      printf("v");
      for (i=0;i<nit/itjmp;i++){
	printf("_");
      }
      printf("v\n");
      printf("[");
      fflush(stdout);
      for (i=0;i<nit;i++)
	{
	  if (i%itjmp==0){
	    int v = i/itjmp;
	    v-=dots;
	    while (v > 0){
	      printf(".");
	      fflush(stdout);
	      v--;
	      dots++;
	    }
	  }
	  
	  for (j=0;j<control->psr[p].nToAs;j++){
	    double t = (double)(control->psr[p].obs[j].sat);
	    if(t > lastMJD)t=lastMJD;
	    dms[j]=getRedNoiseValue(model,t,i);
	  }
	  FILE *log_ts;
	  double sum=0;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    sum+=dms[j];
	  }
	  sum/=control->psr[p].nToAs;
	  int mm=-1;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    dms[j]-=sum;
	    double ofreq=control->psr[p].obs[j].freq.dval*1e6;
	    offsets[j] = (double)(dms[j]/DM_CONST/ofreq/ofreq)*1e12;
	  }
	  toasim_write_corrections(corr,header,file);
	}
      fclose(file);
    }
  free(corr);
}

void createBEoffsets(controlStruct *control,int r)
{
  int i,nit,j,p;
  char fname[MAX_STRLEN];
  double globalParameter;
  long double result;

  char name[1024];
  
  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double *offsets; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

  int dd;
  int k;
  int beNum;

  offsets = (double *)malloc(sizeof(double)*MAX_TOAS);

  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  nit = 1;
  printf("ABC\n");
  for (p=0;p<control->npsr;p++)
    {
      printf("p = %d\n",p);
      header = toasim_init_header();
      strcpy(header->short_desc,"addBEoffsets");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = nit;
      printf("Opening file%s/workFiles/real_%d/%s.addBEoffsets \n",control->name,r,control->psr[p].name);
      sprintf(fname,"%s/workFiles/real_%d/%s.addBEoffsets",control->name,r,control->psr[p].name);
      // First we write the header...
      printf("writing file header\n");
      file = toasim_write_header(header,fname);
      printf("noffsets = %d\n",control->be[beNum].nOffset);
      for (j=0;j<control->psr[p].nToAs;j++){
	beNum = control->psr[p].obs[j].beNum;
	offsets[j]=0.0;
	for (k=0;k<control->be[beNum].nOffset;k++)
	  {
	    if (control->be[beNum].offsetMJD[k].dval < control->psr[p].obs[j].sat)
	      offsets[j] = control->be[beNum].offsetVal[k].dval;
	  }
      }
      printf("Writing corrections\n");
      toasim_write_corrections(corr,header,file);
    
      fclose(file);
    }
  printf("Complete create BE\n");
  free(offsets);
  free(corr);
}

void createDMcovar(controlStruct *control,int r)
{
  int i,nit,j,p;
  char fname[MAX_STRLEN];
  double globalParameter;
  long double result;

  double secperyear=365*86400.0;
  // my parameters
  int npts=1024;
  double alpha,a,b;
  char writeTextFiles=0;
  double lastMJD=1e99;
  char name[1024];
  
  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double offsets[MAX_TOAS]; // Will change to doubles - should use malloc
  double dms[MAX_TOAS]; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

  int dd;
  //  double *covar,x;
  double x;
  int ndays;
  fftwf_complex *covar;
  fftwf_complex *out;
  fftwf_complex *spectrum;
  fftwf_complex *data;
  fftw_plan plan;
  fftwf_plan planf;
  double scale;

  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  nit = 1;

  for (dd=0;dd<control->nDMcovar;dd++)
    {
      p = control->dmCovar[dd].psrNum;
      alpha = control->dmCovar[dd].alpha.dval;
      a = control->dmCovar[dd].a.dval;
      b = control->dmCovar[dd].b.dval;

      header = toasim_init_header();
      strcpy(header->short_desc,"addDmCovar");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = nit;
      
      sprintf(fname,"%s/workFiles/real_%d/%s.dmcovar.%d",control->name,r,control->psr[control->dmCovar[dd].psrNum].name,dd);
      // First we write the header...
      file = toasim_write_header(header,fname);
      
      double mjd_start=1000000.0;
      double mjd_end=-10000000.0;
      for (j=0;j<control->psr[p].nToAs;j++){
	// find the start and end times
	if(control->psr[p].obs[j].sat < mjd_start)mjd_start=(double)control->psr[p].obs[j].sat;
	if(control->psr[p].obs[j].sat > mjd_end)mjd_end=(double)control->psr[p].obs[j].sat;
      }
      
      ndays=ceil((mjd_end-mjd_start)+1e-10);
      //      covar=(double*)malloc(sizeof(double)*(ndays+1)*2);
      covar = (fftwf_complex*) fftwf_malloc((2*ndays+1)*sizeof(fftwf_complex));
      out = (fftwf_complex*) fftwf_malloc((2*ndays+1)*sizeof(fftwf_complex));
      spectrum = (fftwf_complex*) fftwf_malloc((2*ndays+1)*sizeof(fftwf_complex));
      data = (fftwf_complex*) fftwf_malloc((2*ndays+1)*sizeof(fftwf_complex));
      //     out=fftw_malloc(sizeof(fftw_complex)*ndays);
      //      printf("start    = %f (mjd)\n",mjd_start);
      //      printf("end      = %f (mjd)\n",mjd_end  );
      //      printf("npts     = %d (days)\n",npts     );
      //      printf("D_d(%f)  = %f (us^2)\n",d,D_d    );
      //      printf("ref_freq = %f (MHz)\n",ref_freq    );     
      //      printf("seed     = %d\n",control->seed);
      
      
      // Form the covariance function
      for (i=0; i <= ndays; i++){
	x = (i+1e-10);
	covar[i][0]=a*exp(-pow(x/b,alpha));
	covar[i][1]=0;
      }
      for (i=ndays+1;i<=2*ndays;i++)
	{
	  covar[i][0]=covar[2*ndays+1-i][0];
	  covar[i][1]=0;
	}
      //      plan = fftw_plan_dft_r2c_1d(ndays*2+1,covar,out,FFTW_ESTIMATE);
      planf = fftwf_plan_dft_1d(ndays*2+1,covar,out,FFTW_FORWARD,FFTW_ESTIMATE);
      fftwf_execute(planf);
      for (i=0;i<ndays*2;i++)
	{
	  printf("fft: %d %g %g %g\n",i,out[i][0],out[i][1],covar[i][0]);
	}
      fftwf_destroy_plan(planf);
      //      exit(1);
      // Checking with an inverse
      //      planf=fftwf_plan_dft_1d(ndays*2+1,out,covar,FFTW_BACKWARD,FFTW_ESTIMATE);
      //      fftwf_execute(planf);
      //      for (i=0;i<ndays*2+1;i++)
      //	printf("reverse: %d %g\n",i,covar[i][0]);

      //      fftwf_destroy_plan(planf);
      //      exit(1);
      spectrum[0][0]=0;
      spectrum[0][1]=0;

      for (i=0;i<2*ndays+1;i++)
	{
	  scale = sqrt(fabs(out[i][0]))/(double)(sqrt(2*ndays+1));
	  //	  scale=1;
	  spectrum[i][0] = (scale*TKgaussDev(&(control->seed)));
      	  spectrum[i][1] = (scale*TKgaussDev(&(control->seed)));
	  printf("spectrum %d %g %g\n",i,spectrum[i][0],spectrum[i][1]);
	}
      planf=fftwf_plan_dft_1d(2*ndays+1,spectrum,data,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftwf_execute(planf);

      for (i=0;i<2*ndays+1;i++)
	printf("tseries: %d %g %g\n",i,data[i][0],data[i][1]);

      printf("A\n");
      fftwf_destroy_plan(planf);
      printf("B\n");
      fftwf_free(spectrum);
      printf("C\n");
      //      free(data);
      printf("D\n");
      //fftwf_free(out);
      printf("E\n");
      /*      double pism = 0.0112 * D_d * pow(d,(-5.0/3.0)) * pow(secperyear,-1.0/3.0);
	      
	      printf("pism(1yr)  = %g (yr^3) \n",pism );
	      
	      double yr2dm = secperyear * DM_CONST*pow(ref_freq,2.0);
	      printf("yr2dm   = %f \n",yr2dm    );
	      
	      pism *= pow(yr2dm,2); // convert yr^3 to yr.cm^-3.pc
	      
	      printf("pism(1yr)  = %g ((cm^-3pc)^2 yr) \n",pism );
	      
	      printf("\n");
	      printf("Generating red noise...\n");
	      
	      rednoisemodel_t* model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,pism,alpha);
	      populateRedNoiseModel(model,&(control->seed));
	      
	      int itjmp=nit/50;
	      if (itjmp<1)itjmp=1;
	      int dots=0;
	      printf("v");
	      for (i=0;i<nit/itjmp;i++){
	      printf("_");
	      }
	      printf("v\n");
	      printf("[");
	      fflush(stdout);
      */
      for (i=0;i<nit;i++)
	{ 	     	  
	  for (j=0;j<control->psr[p].nToAs;j++){
	    double t = (double)(control->psr[p].obs[j].sat);
	    if(t > lastMJD)t=lastMJD;
	    //	    dms[j]=getRedNoiseValue(model,t,i);
	    dms[j]=data[(int)(t-mjd_start)][0];
	  }
	  FILE *log_ts;
	  double sum=0;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    sum+=dms[j];
	  }
	  sum/=control->psr[p].nToAs;
	  int mm=-1;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    dms[j]-=sum;
	    double ofreq=control->psr[p].obs[j].freq.dval*1e6;
	    offsets[j] = (double)(dms[j]/DM_CONST/ofreq/ofreq)*1e12;
	  }
	  toasim_write_corrections(corr,header,file);
	} 
      free(covar);
      
      fclose(file);
    }
  free(corr);
}

void createDMfunc(controlStruct *control,int r)
{
  int i,nit,j,p;
  char fname[MAX_STRLEN];
  double globalParameter;
  long double result;
  double res;
  char expression[1024];
  char expression1[1024];
  int errorFlag=0;

  double secperyear=365*86400.0;
  double ofreq;
  int dd;
  //
  // For the output file
  //
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  double offsets[MAX_TOAS]; // Will change to doubles - should use malloc
  double dms[MAX_TOAS]; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  char name[1024];

  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  nit = 1;

  for (dd=0;dd<control->nDMfunc;dd++)
    {
      p = control->dmFunc[dd].psrNum;

      header = toasim_init_header();
      strcpy(header->short_desc,"addDmFunc");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Desciprtion of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = nit;
      
      sprintf(fname,"%s/workFiles/real_%d/%s.dmfunc.%d",control->name,r,control->psr[control->dmFunc[dd].psrNum].name,dd);
      // First we write the header...
      file = toasim_write_header(header,fname);
      
      for (i=0;i<nit;i++)
	{ 	     	  
	  strcpy(expression1,control->dmFunc[dd].ddm.inVal);
	  changeRandomOnce(expression1,control);
	  for (j=0;j<control->psr[p].nToAs;j++){
	    sprintf(expression,"x=%g; v=%s;",(double)control->psr[p].obs[j].sat,expression1);
	    errorFlag = runEvaluateExpression(expression,control);      
	    res = variable[0].value;
	    printf("dmfunc: %g %g\n",(double)control->psr[p].obs[j].sat,res);
	    ofreq=control->psr[p].obs[j].freq.dval*1e6;
	    offsets[j] = (double)(res/DM_CONST/ofreq/ofreq)*1e12;
	  }
	  toasim_write_corrections(corr,header,file);
	} 
      fclose(file);
    }
  free(corr);
}

void createTnoise(controlStruct *control,int r)
{
  int p,i,j;
  FILE *file;
  int npts=1024;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  double offsets[MAX_TOAS]; // should use malloc
  double mjds[MAX_TOAS]; //  should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  char name[MAX_STRLEN];
  float cnr_cut=0;
  float cnr_flat=0;
  float old_fc=-1;
  int nit=1;
  double alpha,beta;
  float p_1yr=-1; // s^2 yr
  double secperyear = 86400.0*365.25;
  rednoisemodel_t* model;
  FILE *fout;
  int t;
  char fn[1024];
  sprintf(fname,"%s/setup/useParams",control->name);
  fout = fopen(fname,"a"); 

  for (t=0;t<control->nTnoise;t++)
    {
      p = control->tnoise[t].psrNum;
      alpha = control->tnoise[t].alpha.dval;
      beta = control->tnoise[t].beta.dval;
      p_1yr = control->tnoise[t].p0.dval*secperyear*secperyear;
      old_fc = control->tnoise[t].fc.dval;
      cnr_flat = old_fc;
      fprintf(fout,"rednoise: [%d] %s %g %g %g\n",r,control->psr[p].name,control->tnoise[t].p0.dval,control->tnoise[t].alpha.dval,control->tnoise[t].fc.dval);
      corr->offsets=offsets;
      corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
      // Same length string in every iteration - defined in r_param_length see below
      corr->a0=0; // constant
      corr->a1=0; // a1*x
      corr->a2=0; // a2*x*X
      
      header = toasim_init_header();
      strcpy(header->short_desc,"tNoise");
      strcpy(header->invocation,"ptaSimulate");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Description of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = 1;
      
      // First we write the header...
      sprintf(fname,"%s/workFiles/real_%d/%s.tnoise.%d",control->name,r,control->psr[control->tnoise[t].psrNum].name,t);
      file = toasim_write_header(header,fname);
      
      double mjd_start=(double)control->minT;
      double mjd_end=(double)control->maxT;
      
      printf("Setting up the rednoise model %g %g\n",mjd_start,mjd_end);
      model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha,beta);
      model->cutoff=cnr_cut;
      model->flatten=cnr_flat;
      if(old_fc>0)
	model->mode=MODE_T2CHOL;
      populateRedNoiseModel(model,&(control->seed));
      
      int itjmp=nit/50;
      if (itjmp<1)itjmp=1;
      int dots=0;
      printf("v");
      for (i=0;i<nit/itjmp;i++){
	printf("_");
      }
      //      printf("v\n");
      //      printf("[");
      fflush(stdout);
      for (i=0;i<nit;i++)
	{
	  if (i%itjmp==0){
	    int v = i/itjmp;
	    v-=dots;
	    while (v > 0){
	      //	      printf(".");
	      fflush(stdout);
	      v--;
	      dots++;
	    }
	  }
	  
	  for (j=0;j<control->psr[p].nToAs;j++){
	    offsets[j]=getRedNoiseValue(model,control->psr[p].obs[j].sat,i);
	    printf("rednoise1 offsets = %g %g\n",offsets[j],(double)control->psr[p].obs[j].sat);
	  }
	  //	  exit(1);
	  FILE *log_ts;
	  double sum=0;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    sum+=offsets[j];
	  }
	  sum/=control->psr[p].nToAs;
	  for (j=0;j<control->psr[p].nToAs;j++){
	    mjds[j]=(double)control->psr[p].obs[j].sat;
	    offsets[j]-=sum;
	  }
	  TKremovePoly_d(mjds,offsets,control->psr[p].nToAs,2); // remove a quadratic to reduce the chances of phase wraps
	  // The above is ok because it's linear with F0/F1
	  //	  for (j=0;j<control->psr[p].ntoas;j++){
	  //	    	    printf("offsets: %g\n",offsets[j]);
	  //	  }
	  toasim_write_corrections(corr,header,file);
	}
      int v = i/itjmp;
      v-=dots;
      while (v > 0){
	//	printf(".");
	v--;
	dots++;
      }
      
      //      printf("]\n");
      
      printf("Close file\n");
      fclose(file);
    }
  fclose(fout);
  free(corr);
}

void createPlanets(controlStruct *control,int r)
{
  int p,i,j,nit=1;
  FILE *file;
  int npts=1024;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  double offsets[MAX_TOAS]; // should use malloc
  double mjds[MAX_TOAS]; //  should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  char name[MAX_STRLEN];
  double pb,ecc,a1,t0,om;
  double phase;
  FILE *fout;
  int t;
  char fn[1024];
  sprintf(fname,"%s/setup/useParams",control->name);
  fout = fopen(fname,"a"); 

  printf("planets: in here with %d\n",control->nPlanets);

  for (t=0;t<control->nPlanets;t++)
    {
      p = control->planets[t].psrNum;
      pb = control->planets[t].pb.dval;
      ecc = control->planets[t].ecc.dval;
      t0 = control->planets[t].t0.dval;
      om = control->planets[t].om.dval;     
      a1 = control->planets[t].a1.dval;       

      corr->offsets=offsets;
      corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
      // Same length string in every iteration - defined in r_param_length see below
      corr->a0=0; // constant
      corr->a1=0; // a1*x
      corr->a2=0; // a2*x*X
      
      header = toasim_init_header();
      strcpy(header->short_desc,"planets");
      strcpy(header->invocation,"ptaSimulate");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Description of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed;
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = 1;
      
      // First we write the header...
      sprintf(fname,"%s/workFiles/real_%d/%s.planets.%d",control->name,r,control->psr[control->planets[t].psrNum].name,t);
      file = toasim_write_header(header,fname);
      
      
      int itjmp=nit/50;
      if (itjmp<1)itjmp=1;
      int dots=0;
      printf("v");
      for (i=0;i<nit/itjmp;i++){
	printf("_");
      }
      //      printf("v\n");
      //      printf("[");
      fflush(stdout);
      for (i=0;i<nit;i++)
	{
	  if (i%itjmp==0){
	    int v = i/itjmp;
	    v-=dots;
	    while (v > 0){
	      //	      printf(".");
	      fflush(stdout);
	      v--;
	      dots++;
	    }
	  }
	  printf("planets: number of toas = %d, psr = %d\n",control->psr[p].nToAs,p);
	  for (j=0;j<control->psr[p].nToAs;j++){
	    offsets[j]=-BTmodel(pb,ecc,a1,t0,om,control->psr[p].obs[j].sat);
	    printf("planets: offsets = %g\n",offsets[j]);
	  }
	  //	  exit(1);
	  toasim_write_corrections(corr,header,file);
	}
      int v = i/itjmp;
      v-=dots;
      while (v > 0){
	//	printf(".");
	v--;
	dots++;
      }
      
      //      printf("]\n");
      
      printf("Close file\n");
      fclose(file);
    }
  fclose(fout);
  free(corr);
}


void createClkNoise(controlStruct *control,int r)
{
  int p,i,j;
  FILE *file;
  int npts=1024;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  double offsets[MAX_TOAS]; // should use malloc
  double mjds[MAX_TOAS]; //  should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  char name[MAX_STRLEN];
  float cnr_cut=0;
  float cnr_flat=0;
  float old_fc=-1;
  int nit=1;
  double alpha,beta;
  float p_1yr=-1; // s^2 yr
  double secperyear = 86400.0*365.25;
  rednoisemodel_t* model;
  FILE *fout;
  int t;
  char fn[1024];
  sprintf(fname,"%s/setup/useParams",control->name);
  fout = fopen(fname,"a"); 

  for (t=0;t<control->nClkNoise;t++)
    {
      for (p=0;p<control->npsr;p++)
	{
	  alpha = control->clkNoise.alpha.dval;
	  beta = 0;
	  p_1yr = control->clkNoise.p0.dval*secperyear*secperyear;
	  old_fc = control->clkNoise.fc.dval;
	  cnr_flat = old_fc;
	  corr->offsets=offsets;
	  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
	  // Same length string in every iteration - defined in r_param_length see below
	  corr->a0=0; // constant
	  corr->a1=0; // a1*x
	  corr->a2=0; // a2*x*X
	  
	  header = toasim_init_header();
	  strcpy(header->short_desc,"clkNoise");
	  strcpy(header->invocation,"ptaSimulate");
	  sprintf(name,"%s.sim",control->psr[p].name);
	  strcpy(header->timfile_name,name);
	  
	  strcpy(header->parfile_name,"Unknown");
	  header->idealised_toas="NotSet"; // What should this be
	  header->orig_parfile="NA";
	  header->gparam_desc=""; // Global parameters
	  header->gparam_vals="";
	  header->rparam_desc=""; // Description of the parameters
	  header->rparam_len=0; // Size of the string
	  header->seed = control->seed;
	  
	  header->ntoa = control->psr[p].nToAs;
	  header->nrealisations = 1;
	  
	  // First we write the header...
	  sprintf(fname,"%s/workFiles/real_%d/%s.clknoise.%d",control->name,r,control->psr[p].name,t);
	  file = toasim_write_header(header,fname);

	  if (p==0)
	    {
	      double mjd_start=(double)control->minT;
	      double mjd_end=(double)control->maxT;
	      
	      
	      printf("Setting up the noise model %g %g\n",mjd_start,mjd_end);
	      model = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha,beta);
	      model->cutoff=cnr_cut;
	      model->flatten=cnr_flat;
	      if(old_fc>0)
		model->mode=MODE_T2CHOL;
	      populateRedNoiseModel(model,&(control->seed));
	    }	      
	  int itjmp=nit/50;
	  if (itjmp<1)itjmp=1;
	  int dots=0;
	  printf("v");
	  for (i=0;i<nit/itjmp;i++){
	    printf("_");
	  }
	  //      printf("v\n");
	  //      printf("[");
	  fflush(stdout);
	  for (i=0;i<nit;i++)
	    {
	      if (i%itjmp==0){
		int v = i/itjmp;
		v-=dots;
		while (v > 0){
		  //	      printf(".");
		  fflush(stdout);
		  v--;
		  dots++;
		}
	      }
	      
	      for (j=0;j<control->psr[p].nToAs;j++){
		offsets[j]=getRedNoiseValue(model,control->psr[p].obs[j].sat,i);
		//	    printf("offsets = %g\n",offsets[j]);
	      }
	      //	  exit(1);
	      FILE *log_ts;
	      double sum=0;
	      for (j=0;j<control->psr[p].nToAs;j++){
		sum+=offsets[j];
	      }
	      sum/=control->psr[p].nToAs;
	      for (j=0;j<control->psr[p].nToAs;j++){
		mjds[j]=(double)control->psr[p].obs[j].sat;
		offsets[j]-=sum;
	      }
	      TKremovePoly_d(mjds,offsets,control->psr[p].nToAs,2); // remove a quadratic to reduce the chances of phase wraps
	      // The above is ok because it's linear with F0/F1
	      //	  for (j=0;j<control->psr[p].ntoas;j++){
	      //	    	    printf("offsets: %g\n",offsets[j]);
	      //	  }
	      toasim_write_corrections(corr,header,file);
	    }
	  int v = i/itjmp;
	  v-=dots;
	  while (v > 0){
	    //	printf(".");
	    v--;
	    dots++;
	  }
	  
	  //      printf("]\n");
	  
	  printf("Close file\n");
	  fclose(file);
	}
    }
  fclose(fout);
  free(corr);
}




void createRadiometerNoise(controlStruct *control, int r)
{
  int i,p,j;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  char name[MAX_STRLEN];
  double offsets[MAX_TOAS];
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  double err,efac,equad;
  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X


  for (p=0;p<control->npsr;p++)
    {
      printf("Processing pulsar: %s\n",control->psr[p].name);
      header = toasim_init_header();
      printf("... Created header\n");
      strcpy(header->short_desc,"addGaussian");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Description of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed; // SHOULD SET THIS
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = 1;

      sprintf(fname,"%s/workFiles/real_%d/%s.addGauss",control->name,r,control->psr[p].name);
      printf("... Opening file\n");
      file = toasim_write_header(header,fname);
      printf("... Creating offsets: %d\n",control->psr[p].nToAs);

      // ADD IN EFAC/EQUAD
      // MUST DO
      for (j=0;j<control->psr[p].nToAs;j++)
	{
	  err = control->psr[p].obs[j].toaErr.dval;
	  efac = control->psr[p].obs[j].efac.dval;
	  equad = control->psr[p].obs[j].equad.dval;
	  printf("have efac %g %g %g >%s<\n",err,efac,equad,control->psr[p].obs[j].efac.inVal);

	  err = sqrt(pow(err,2)+pow(equad,2))*efac;
	  offsets[j] = err*TKgaussDev(&(control->seed));
	}
      printf(" ... Outputing file\n");
      toasim_write_corrections(corr,header,file);
      printf("... Closing file\n");
      fclose(file);
    }
  free(corr);
}

void createJitter(controlStruct *control, int r)
{
  int i,p,j;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  char name[MAX_STRLEN];
  double offsets[MAX_TOAS];
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));
  int dd;
  double jLevel,tobs;

  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X

  for (dd=0;dd<control->nJitter;dd++)
    {
      p = control->jitter[dd].psrNum;

      header = toasim_init_header();
      strcpy(header->short_desc,"addJitter");
      strcpy(header->invocation,"");
      sprintf(name,"%s.sim",control->psr[p].name);
      strcpy(header->timfile_name,name);
      strcpy(header->parfile_name,"Unknown");
      header->idealised_toas="NotSet"; // What should this be
      header->orig_parfile="NA";
      header->gparam_desc=""; // Global parameters
      header->gparam_vals="";
      header->rparam_desc=""; // Description of the parameters
      header->rparam_len=0; // Size of the string
      header->seed = control->seed; // SHOULD SET THIS
      
      header->ntoa = control->psr[p].nToAs;
      header->nrealisations = 1;

      sprintf(fname,"%s/workFiles/real_%d/%s.jitter.%d",control->name,r,control->psr[control->jitter[dd].psrNum].name,dd);
      // First we write the header...
      file = toasim_write_header(header,fname);

      for (j=0;j<control->psr[p].nToAs;j++)
	{
	  tobs = control->psr[p].obs[j].tobs.dval;
	  jLevel = control->jitter[dd].sigma_j.dval*sqrt(control->jitter[dd].t0.dval/tobs);
	  offsets[j] = jLevel*TKgaussDev(&(control->seed));
	  printf("Adding jitter: %g\n",offsets[j]);
	}
      toasim_write_corrections(corr,header,file);

      fclose(file);
    }
  free(corr);
}


double hms_turn(char *line){

  /* Converts string " hh:mm:ss.ss" or " hh mm ss.ss" to double turn */
  
  int i;int turn_hms(double turn, char *hms);
  double hr, min, sec, turn=0;
  char hold[MAX_STRLEN];

  strcpy(hold,line);

  /* Get rid of ":" */
  for(i=0; *(line+i) != '\0'; i++)if(*(line+i) == ':')*(line+i) = ' ';

  i = sscanf(line,"%lf %lf %lf", &hr, &min, &sec);
  if(i > 0){
    turn = hr/24.;
    if(i > 1)turn += min/1440.;
    if(i > 2)turn += sec/86400.;
  }
  if(i == 0 || i > 3)turn = 1.0;


  strcpy(line,hold);

  return turn;
}

double dms_turn(char *line){

  /* Converts string "-dd:mm:ss.ss" or " -dd mm ss.ss" to double turn */
  
  int i;
  char *ic, ln[40];
  double deg, min, sec, sign, turn=0;

  /* Copy line to internal string */
  strcpy(ln,line);

  /* Get rid of ":" */
  for(i=0; *(ln+i) != '\0'; i++)if(*(ln+i) == ':')*(ln+i) = ' ';

  /* Get sign */
  if((ic = strchr(ln,'-')) == NULL)
     sign = 1.;
  else {
     *ic = ' ';
     sign = -1.;
  }

  /* Get value */
  i = sscanf(ln,"%lf %lf %lf", &deg, &min, &sec);
  if(i > 0){
    turn = deg/360.;
    if(i > 1)turn += min/21600.;
    if(i > 2)turn += sec/1296000.;
    if(turn >= 1.0)turn = turn - 1.0;
    turn *= sign;
  }
  if(i == 0 || i > 3)turn =1.0;

  return turn;
}

void createGW(controlStruct *control, int r)
{
  int i,p,j,k;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  FILE* file;
  char expression[MAX_STRLEN];
  int errorFlag=0;
  char name[MAX_STRLEN];
  double offsets[MAX_TOAS];
  double epochs[MAX_TOAS]; // Will change to doubles - should use malloc
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

  gwSrc *gw;
  long double timeOffset=0.0; 
  long double scale;
  long double alpha;
  long double gwAmp;
  long double ra_p,dec_p;
  long double flo=0.0,fhi=0.0;
  long double kp[3];            /* Vector pointing to pulsar           */
  long double tspan = (control->maxT - control->minT)*86400.0L;
  long double time;
  long double gwRes[MAX_TOAS];
  double gwRes_p[MAX_TOAS];
  double gwRes_c[MAX_TOAS];
  long double dist[control->npsr];
  long double mean;
  int distNum=0;
  int logspacing=1;
  int ngw=1000;
  char readGW=0;
  char writeGW=0;
  char gwFileName[MAX_STRLEN];
  FILE *gwFile;
  int kk;
  int nCW=0;

  double resp,resc;
  double lambda_p,beta_p,lambda,beta;
  double *h0_g,*omega_g,*angpol_g,*lambda_g,*beta_g,*skyInc_g;
  long sSeed = TKsetSeed();

  
  for (kk=0;kk<control->nGW;kk++)
    {
      alpha = control->gw[kk].alpha.dval;
      gwAmp = control->gw[kk].amp.dval;

      if (control->gw[kk].type==1)
	{
	  scale = pow(86400.0*365.25,alpha);
	  gwAmp *= scale;
	  if ((gw = (gwSrc *)malloc(sizeof(gwSrc)*ngw))==NULL)
	    {
	      printf("Unable to allocate memory for %d GW sources\n",ngw);
	      exit(1);
	    }
	  
	  if (flo==0)
	    flo=0.01/tspan;
	  if (fhi==0)
	    fhi = 1.0/(long double)86400.0L;
	  
	  timeOffset = 0.5*(control->maxT + control->minT);
	  
	  
	  GWbackground(gw,ngw,&(control->seed),flo,fhi,gwAmp,alpha,logspacing);
	  for (i=0;i<ngw;i++)
	    setupGW(&gw[i]);
	}

      if (control->gw[kk].type==5)
	{
	  double res_r,res_i;
	  double SPEED_LIGHT = 299792458.0; /*!< Speed of light (m/s)                       */
	  double GM = 1.3271243999e20;      /*!< Gravitational constant * mass sun          */
	  double n1,n2,n3;
	  double B1,B2,deltaPhi,res;
	  double cosTheta;
	  double redshift,skyPhi,skyTheta,logfgw,polaris_ang,cm_dist_mpc,cw_phase,log_ch_mass,cosinc,mc,angpol;
	  double massRatio,logObsFreq,temp,distance;
	  double H0=70;
	  FILE *fin;


	  h0_g = (double *)malloc(sizeof(double)*MAX_CWS);
	  omega_g = (double *)malloc(sizeof(double)*MAX_CWS);
	  if (!(angpol_g = (double *)malloc(sizeof(double)*MAX_CWS)))
	    {printf("Error in allocating memory\n"); exit(1);}
       
	  lambda_g = (double *)malloc(sizeof(double)*MAX_CWS);
	  beta_g = (double *)malloc(sizeof(double)*MAX_CWS);
	  if (!(skyInc_g = (double *)malloc(sizeof(double)*MAX_CWS)))
	    {printf("Error in allocating memory\n"); exit(1);}
	
	  if (!(fin = fopen(control->gw[kk].fname,"r")))
	    {
	      printf("Unable to open file >%s< to read the GW soure parameters\n",control->gw[kk].fname);
	      exit(1);
	    }		 
	  
	  // Read all the sources for this pulsar	  
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&log_ch_mass,&massRatio,&redshift,&distance,&logObsFreq,
			 &temp,&temp,&temp,&temp,&(skyInc_g[nCW]))==10)
		//		  if (fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&redshift,&skyPhi,&skyTheta,
		//			     &skyInc,&logfgw,&polaris_ang,&cm_dist_mpc,&cw_phase,&log_ch_mass)==9)
		
		{
		  cosinc = cos(skyInc_g[nCW]);
		  mc = pow(10,log_ch_mass);
		  cm_dist_mpc = distance*1000.0; // Now read this from the input file		 		  
		  h0_g[nCW] = pow(GM*mc*(1+redshift),(5.0/3.0))/pow(SPEED_LIGHT,4)/(cm_dist_mpc*1e6*3.08568025e16)*pow((M_PI*pow(10,logObsFreq)),(2.0/3.0));		 
		  //		      omega_g = 2*M_PI*pow(10,logfgw);
		  omega_g[nCW] = 2*M_PI*pow(10,logObsFreq);		  
		  angpol_g[nCW] = 2*M_PI*TKranDev(&sSeed);
		  lambda_g[nCW] = acos((TKranDev(&sSeed)-0.5)*2);  
		  beta_g[nCW]   = TKranDev(&sSeed)*2*M_PI;  
		  nCW++;
		  if (nCW == MAX_CWS)
		    {
		      printf("ERROR: Must increase number of CW sources in ptaSimulate.h\n");
		      exit(1);
		    }		  
		}
	    }
	  fclose(fin);
	}
		  
      for (p=0;p<control->npsr;p++)
	{
	  printf("Creating GWB noise file for pulsar: %d\n",p);
	  corr->offsets=offsets;
	  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
	  // Same length string in every iteration - defined in r_param_length see below
	  corr->a0=0; // constant
	  corr->a1=0; // a1*x
	  corr->a2=0; // a2*x*X
	  //  printf("Here with %g %g %g %g\n",(double)flo,(double)fhi,(double)gwAmp,(double)alpha);

	  header = toasim_init_header();
	  strcpy(header->short_desc,"addGW");
	  strcpy(header->invocation,"");
	  sprintf(name,"%s.sim",control->psr[p].name);
	  strcpy(header->timfile_name,name);
	  strcpy(header->parfile_name,"Unknown");
	  header->idealised_toas="NotSet"; // What should this be
	  header->orig_parfile="NA";
	  header->gparam_desc=""; // Global parameters
	  header->gparam_vals="";
	  header->rparam_desc=""; // Description of the parameters
	  header->rparam_len=0; // Size of the string
	  header->seed = 0; // SHOULD SET THIS
	  
	  header->ntoa = control->psr[p].nToAs;
	  header->nrealisations = 1;
	  
	  sprintf(fname,"%s/workFiles/real_%d/%s.addGW.%d",control->name,r,control->psr[p].name,kk);
	  file = toasim_write_header(header,fname);
	  
	  dist[p] = control->psr[p].dist;
	  printf("dist = %Lg\n",dist[p]);
	  ra_p = control->psr[p].rajd*M_PI/180.0;
	  dec_p = control->psr[p].decjd*M_PI/180.0;
	  if (control->gw[kk].type==1)
	    setupPulsar_GWsim(ra_p,dec_p,kp);
	  else if (control->gw[kk].type==2)
	    {
	      double n1,n2,n3;
	      double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
	      double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
	      double cosTheta;

	      lambda_p = (double)ra_p;
	      beta_p   = (double)dec_p;
	      lambda   = control->gw[kk].ra.dval;
	      beta   = control->gw[kk].dec.dval;
	      // Pulsar vector
	      n1 = cosl(lambda_p)*cosl(beta_p);
	      n2 = sinl(lambda_p)*cosl(beta_p);
	      n3 = sinl(beta_p);
	      
	      cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
		sinl(beta)*sinl(beta_p);
	      
	      // From KJ's paper
	      // Gravitational wave matrix
	      
	      // NOTE: This is for the plus terms.  For cross should use different terms
	      e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
	      e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
	      e31p = cosl(lambda)*sinl(beta)*cosl(beta);
	      
	      e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
	      e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
	      e32p = sinl(lambda)*sinl(beta)*cosl(beta);
	      
	      e13p = cosl(lambda)*sinl(beta)*cosl(beta);
	      e23p = sinl(lambda)*sinl(beta)*cosl(beta);
	      e33p = -powl(cosl(beta),2);
	      
	      resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
		      n2*(n1*e21p+n2*e22p+n3*e23p)+
		      n3*(n1*e31p+n2*e32p+n3*e33p));
	      
	      if ((1-cosTheta)==0.0)
		resp = 0.0;  // Check if this is sensible
	      else
		resp = 1.0L/(2.0L*(1.0L-cosTheta))*(resp); 

	      e11c = sin(2*lambda)*sin(beta);
	      e21c = -cos(2*lambda)*sin(beta);
	      e31c = -sin(lambda)*cos(beta);
	      
	      e12c = -cos(2*lambda)*sin(beta);
	      e22c = -sin(2*lambda)*sin(beta);
	      e32c = cos(lambda)*cos(beta);
	      
	      e13c = -sin(lambda)*cos(beta);
	      e23c = cos(lambda)*cos(beta);
	      e33c  = 0;
	      
	      resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
		      n2*(n1*e21c+n2*e22c+n3*e23c)+
		      n3*(n1*e31c+n2*e32c+n3*e33c));
	      
	      if ((1-cosTheta)==0.0)
		resc = 0.0;  // Check if this is sensible
	      else
		resc = 1.0L/(2.0L*(1.0L-cosTheta))*(resc); 
	    }
	    
	  if (control->gw[kk].type==5)
	    {
	      double res_r,res_i;
	      double SPEED_LIGHT = 299792458.0; /*!< Speed of light (m/s)                       */
	      double GM = 1.3271243999e20;      /*!< Gravitational constant * mass sun          */
	      double n1,n2,n3;
	      double B1,B2,deltaPhi,res;
	      double cosTheta;
	      double redshift,skyPhi,skyTheta,logfgw,polaris_ang,cm_dist_mpc,cw_phase,log_ch_mass,h0,cosinc,mc,angpol;
	      double massRatio,logObsFreq,temp,distance;
	      double H0=70;

	      mean=0;
	      for (i=0;i<control->psr[p].nToAs;i++)
		gwRes[i]=0.0;

	      // Read all the sources for this pulsar
	    	      
	      for (j=0;j<nCW;j++)
		{
		  lambda_p = (double)ra_p;
		  beta_p   = (double)(dec_p); 
		  
		  // Lee et al. (2011) equation 12
		  B1 = (1+pow(sin(beta_g[j]),2))*pow(cos(beta_p),2)*cos(2*(lambda_g[j]-lambda_p))-
		    sin(2*beta_g[j])*sin(2*beta_p)*cos(lambda_g[j]-lambda_p)+(2.0-3.0*pow(cos(beta_p),2))*pow(cos(beta_g[j]),2);
		  // Lee et al. (2011) equation 13
		  B2 = 2*cos(beta_g[j])*sin(beta_p)*sin(lambda_g[j]-lambda_p)-2*sin(beta_g[j])*pow(cos(beta_p),2)*sin(2*(lambda_g[j]-lambda_p));
		  // Lee et al. (2011) equation 14
		  cosTheta = -(cosl(beta_g[j])*cosl(beta_p)*cosl(lambda_g[j]-lambda_p) + sinl(beta_g[j])*sinl(beta_p));
		  deltaPhi = omega_g[j]*dist[p]*(1.0-cosTheta);
		  
		  // NOTE THAT h0 IS NOT THE SAME AS THE GW AMPLITUDE GIVEN IN LEE ET AL. UNDER EQUATION 14
		  
		  // Now loop through all the residuals for this pulsar
		  mean = 0.0;
		  for (i=0;i<control->psr[p].nToAs;i++)
		    {
		      time = (control->psr[p].obs[i].sat - timeOffset)*86400.0L;
		      // Equation 11 in Lee et al. (2011)
		      res = h0_g[j]/2.0/omega_g[j]*sin(deltaPhi/2.0)/(1.0-cosTheta)*
			((B1*cos(2*angpol_g[j])+B2*sin(2*angpol_g[j]))*cos(omega_g[j]*time-deltaPhi/2.0)*(1+pow(cos(skyInc_g[j]),2))
			 + 2*(B2*cos(2*angpol_g[j])-B1*sin(2*angpol_g[j]))*sin(omega_g[j]*time-deltaPhi/2.0)*cos(skyInc_g[j]));
		      
		      gwRes[i] += (res);
		      mean+=gwRes[i];
		    }
		}
	    }
	  else
	    {
	      mean = 0.0;
	      for (i=0;i<control->psr[p].nToAs;i++)
		{
		  time = (control->psr[p].obs[i].sat - timeOffset)*86400.0L;
		  
		  if (control->gw[kk].type==1)
		    {
		      gwRes[i] = 0.0;	      
		      for (k=0;k<ngw;k++)
			{
			  gwRes[i]+=calculateResidualGW(kp,&gw[k],time,dist[p]);
			  //	      printf("Have: %d %g %g %g %g %g\n",k,(double)gwRes[i],(double)kp[0],(double)kp[1],(double)kp[2],(double)dist[p]);
			}
		      mean += gwRes[i];
		    }
		  else if (control->gw[kk].type==2)
		    {
		      // Must process Ap and Ac at this time
		      sprintf(expression,"x=%g; v=%s;",(double)control->psr[p].obs[i].sat,control->gw[kk].ap.inVal);
		      errorFlag = runEvaluateExpression(expression,control);      
		      gwRes_p[i] = variable[0].value;
		      printf("Results = %g %g %g\n",(double)control->psr[p].obs[i].sat, variable[0].value,variable[1].value);
		      
		      sprintf(expression,"x=%g; v=%s;",(double)control->psr[p].obs[i].sat,control->gw[kk].ac.inVal);
		      errorFlag = runEvaluateExpression(expression,control);      
		      gwRes_c[i] = variable[0].value;
		      
		      gwRes[i] = gwRes_p[i]*resp+gwRes_c[i]*resc;
		      mean+=gwRes[i];
		    }
		  else if (control->gw[kk].type==3) // GWM
		    {
		      // Equations here from Wang et al.
		      long double dt,scale;
		      double cos2Phi;
		      double cosPhi;
		      double l1,l2,l3,m1,m2,m3;
		      //double beta_m;
		      double d1,d2,d3,md;
		      double a1,a2,a3,ma;
		      long double time;
		      double g1,g2,g3;
		      double n1,n2,n3;
		      double cosTheta;
		      
		      /* define the GW coordinate system see Hobbs,G. (2009)*/
		      /* the d vector point to north pole or south pole*/ 
		      
		      time    = (control->psr[p].obs[i].sat - control->gw[kk].gwmEpoch.dval)*86400.0L;
		      if (time > 0)
			{
			  lambda_p = (double)ra_p;
			  beta_p   = (double)dec_p;
			  lambda   = control->gw[kk].ra.dval;
			  beta   = control->gw[kk].dec.dval;
			  
			  // GW vector
			  g1 = -cosl(lambda)*cosl(beta);
			  g2 = -sinl(lambda)*cosl(beta);
			  g3 = -sinl(beta);
			  
			  // Pulsar vector
			  n1 = cosl(lambda_p)*cosl(beta_p);
			  n2 = sinl(lambda_p)*cosl(beta_p);
			  n3 = sinl(beta_p);
			  //	       printf("n = %g %g %g\n",n1,n2,n3);
			  cosTheta = -(cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p) + sinl(beta)*sinl(beta_p));
			  
			  if (beta == 0.0 )
			    {
			      d1 = 0.0;
			      d2 = 0.0;
			      d3 = 1.0;
			    }
			  
			  if ( beta > 0)
			    {
			      d1 = g1*cosl(0.5*M_PI - beta);
			      d2 = g2*cosl(0.5*M_PI - beta);
			      d3 = 1.0 + g3*cos(0.5*M_PI - beta);
			      md = sqrt(d1*d1 + d2*d2 + d3*d3);
			      d1 = d1/md;
			      d2 = d2/md;
			      d3 = d3/md;
			      /*covert d to unit vector */
			    } 
			  else if (beta < 0)
			    {
			      d1 = g1*cosl(-0.5*M_PI - beta);
			      d2 = g2*cosl(-0.5*M_PI - beta);
			      d3 = -1.0 + g3*cos(-0.5*M_PI - beta);
			      md = sqrt(d1*d1 + d2*d2 + d3*d3);
			      d1 = d1/md;
			      d2 = d2/md;
			      d3 = d3/md;
			    } 
		      
			  a1 =  (d2*g3-d3*g2);
			  a2 =  (d3*g1-d1*g3);
			  a3 =  (d1*g2-d2*g1);
			  /* conver it to unit vector */
			  ma = sqrt(a1*a1 +a2*a2 + a3*a3);
			  a1 = a1/ma;
			  a2 = a2/ma;
			  a3 = a3/ma;
			  
			  /* polarisation vector of GW source */
			  m1 = d1*cosl(control->gw[kk].gwmPhi.dval)	+ a1*sinl(control->gw[kk].gwmPhi.dval);   
			  m2 = d2*cosl(control->gw[kk].gwmPhi.dval)	+ a2*sinl(control->gw[kk].gwmPhi.dval);
			  m3 = d3*cosl(control->gw[kk].gwmPhi.dval)	+ a3*sinl(control->gw[kk].gwmPhi.dval);
			  
			  if  (cosTheta != 1.0 && cosTheta != -1.0)
			    {g1 = g1*cosTheta; 
			      g2 = g2*cosTheta;
			      g3 = g3*cosTheta;
			      
			      /*l is  the projection of pulsar vector on the plane which pependicular to the GW source direction */
			      l1 = n1 - g1;
			      l2 = n2 - g2;
			      l3 = n3 - g3;
			      cosPhi = (l1*m1 + l2*m2 + l3*m3)/sqrt(l1*l1 + l2*l2 + l3*l3);
			      //		        if  (cosPhi >= 1.0/sqrt(2.0))
			      cos2Phi = 2*cosPhi*cosPhi - 1.0;
			      //		       else
			      //		      	   cos2Phi = 2*sqrt(1.0 - cosPhi*cosPhi)*sqrt(1.0 - cosPhi*cosPhi) - 1.0;
			    }
			  else 
			    {cos2Phi = 0;}
			  
			  scale = -0.5*cos2Phi*(1-cosTheta);
			  //		   scale=1.0;
			  gwRes[i] = scale*time*control->gw[kk].gwmAmp.dval;
			  mean+=gwRes[i];
			} 
		      else 
			{
			  gwRes[i] =0.0;
			  mean+=gwRes[i];
			}
		    }
		  else if (control->gw[kk].type==4) // CW source
		    {
		      double res_r,res_i;
		      double SPEED_LIGHT = 299792458.0; /*!< Speed of light (m/s)                       */
		      double GM = 1.3271243999e20;      /*!< Gravitational constant * mass sun          */
		      double n1,n2,n3;
		      double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
		      double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
		      double cosTheta;
		      double omega_g;
		      
		      omega_g = 2*M_PI*control->gw[kk].cgw_freq.dval;
		      lambda_p = (double)ra_p;
		      beta_p   = (double)dec_p;
		      lambda   = control->gw[kk].ra.dval;
		      beta   = control->gw[kk].dec.dval;
		      // Pulsar vector
		      n1 = cosl(lambda_p)*cosl(beta_p);
		      n2 = sinl(lambda_p)*cosl(beta_p);
		      n3 = sinl(beta_p);
		      
		      cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
			sinl(beta)*sinl(beta_p);
		      
		      // From KJ's paper
		      // Gravitational wave matrix
		      
		      // NOTE: This is for the plus terms.  For cross should use different terms
		      e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
		      e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		      e31p = cosl(lambda)*sinl(beta)*cosl(beta);
		      
		      e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		      e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
		      e32p = sinl(lambda)*sinl(beta)*cosl(beta);
		      
		      e13p = cosl(lambda)*sinl(beta)*cosl(beta);
		      e23p = sinl(lambda)*sinl(beta)*cosl(beta);
		      e33p = -powl(cosl(beta),2);
		      
		      resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
			      n2*(n1*e21p+n2*e22p+n3*e23p)+
			      n3*(n1*e31p+n2*e32p+n3*e33p));
		      
		      if ((1-cosTheta)==0.0)
			resp = 0.0;  // Check if this is sensible
		      else
			resp = 1.0L/(2.0L*(1.0L-cosTheta))*(resp); 
		      
		      e11c = sin(2*lambda)*sin(beta);
		      e21c = -cos(2*lambda)*sin(beta);
		      e31c = -sin(lambda)*cos(beta);
		      
		      e12c = -cos(2*lambda)*sin(beta);
		      e22c = -sin(2*lambda)*sin(beta);
		      e32c = cos(lambda)*cos(beta);
		      
		      e13c = -sin(lambda)*cos(beta);
		      e23c = cos(lambda)*cos(beta);
		      e33c  = 0;
		      
		      resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
			      n2*(n1*e21c+n2*e22c+n3*e23c)+
			      n3*(n1*e31c+n2*e32c+n3*e33c));
		      
		      if ((1-cosTheta)==0.0)
			resc = 0.0;  // Check if this is sensible
		      else
			resc = 1.0L/(2.0L*(1.0L-cosTheta))*(resc); 
		      
		      res_r = (control->gw[kk].cgw_h0.dval/omega_g*((1+pow(control->gw[kk].cgw_cosinc.dval,2))*cos(2*control->gw[kk].cgw_angpol.dval)*sin(omega_g*time)+2*control->gw[kk].cgw_cosinc.dval*sin(2*control->gw[kk].cgw_angpol.dval)*cos(omega_g*time)))*resp + (control->gw[kk].cgw_h0.dval/omega_g*((1+pow(control->gw[kk].cgw_cosinc.dval,2))*sin(2*control->gw[kk].cgw_angpol.dval)*sin(omega_g*time)-2*control->gw[kk].cgw_cosinc.dval*cos(2*control->gw[kk].cgw_angpol.dval)*cos(omega_g*time)))*resc; 
		      res_i = 0.0;
		      if (dist[p]>0) // Add in the pulsar term  (NOTE: using subtraction here)
			{
			  double omega_prime_g;
			  double h0_prime;
			  
			  if (control->gw[kk].cgw_mc.dval == 0) {omega_prime_g = omega_g; h0_prime = control->gw[kk].cgw_h0.dval;}
			  else {
			    //			 omega_prime_g = omega_g - 2*M_PI*2.77e-8*pow(psr[p].cgw_mc/1e8,5.0/3.0)*pow(omega_g/2.0/M_PI/1e-7,11.0/3.0)*(psr[p].gwsrc_psrdist/PCM/1000.0)*(1-cosTheta);
			    
			    omega_prime_g = 2*M_PI*pow((1-cosTheta)*dist[p]/SPEED_LIGHT*256.0/5.0/pow(SPEED_LIGHT,5)*pow(M_PI,8.0/3.0)*pow(GM*control->gw[kk].cgw_mc.dval,5.0/3.0)+pow(omega_g/2.0/M_PI,-8.0/3.0),-3.0/8.0);
			    
			    h0_prime = control->gw[kk].cgw_h0.dval*pow(omega_prime_g/omega_g,2.0/3.0);
			    //			  printf("Using: omega_prime_g = %g, omega_g = %g, h0_prime = %g, h0 = %g\n",omega_prime_g,omega_g,h0_prime,psr[p].cgw_h0);
			  }
			  res_r -= ((h0_prime/omega_prime_g*((1+pow(control->gw[kk].cgw_cosinc.dval,2))*cos(2*control->gw[kk].cgw_angpol.dval)*sin(omega_prime_g*time-(1-cosTheta)*dist[p]/SPEED_LIGHT*omega_prime_g)+2*control->gw[kk].cgw_cosinc.dval*sin(2*control->gw[kk].cgw_angpol.dval)*cos(omega_prime_g*time-(1-cosTheta)*dist[p]/SPEED_LIGHT*omega_prime_g)))*resp 
				    + (h0_prime/omega_prime_g*((1+pow(control->gw[kk].cgw_cosinc.dval,2))*sin(2*control->gw[kk].cgw_angpol.dval)*sin(omega_prime_g*time-(1-cosTheta)*dist[p]/SPEED_LIGHT*omega_prime_g)-2*control->gw[kk].cgw_cosinc.dval*cos(2*control->gw[kk].cgw_angpol.dval)*cos(omega_prime_g*time-(1-cosTheta)*dist[p]/SPEED_LIGHT*omega_prime_g)))*resc); 
			}
		      if ((1-cosTheta)==0.0)
			{
			  res_r = 0.0;
			  res_i = 0.0;
			}
		      else
			{
			  res_r = -(1.0)/((2.0)*((1.0)-cosTheta))*(res_r); 
			  res_i = -(1.0)/((2.0)*((1.0)-cosTheta))*(res_i); 
			}
		      gwRes[i] = (res_r+res_i);
		      mean+=gwRes[i];
		    }
		  else if (control->gw[kk].type==6)  // Cosmic string burst
		    {
		      long double dt,width,width_day;
		      double res_r,res_i;
		      double SPEED_LIGHT = 299792458.0; /*!< Speed of light (m/s)                       */
		      double GM = 1.3271243999e20;      /*!< Gravitational constant * mass sun          */
		      double n1,n2,n3;
		      double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
		      double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
		      double cosTheta;
		      double omega_g;
		      
		      omega_g = 2*M_PI*control->gw[kk].cgw_freq.dval;
		      lambda_p = (double)ra_p;
		      beta_p   = (double)dec_p;
		      lambda   = control->gw[kk].ra.dval;
		      beta   = control->gw[kk].dec.dval;
		      // Pulsar vector
		      n1 = cosl(lambda_p)*cosl(beta_p);
		      n2 = sinl(lambda_p)*cosl(beta_p);
		      n3 = sinl(beta_p);
		      
		      cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
			sinl(beta)*sinl(beta_p);
		      
		      // From KJ's paper
		      // Gravitational wave matrix
		      
		      // NOTE: This is for the plus terms.  For cross should use different terms
		      e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
		      e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		      e31p = cosl(lambda)*sinl(beta)*cosl(beta);
		      
		      e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
		      e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
		      e32p = sinl(lambda)*sinl(beta)*cosl(beta);
		      
		      e13p = cosl(lambda)*sinl(beta)*cosl(beta);
		      e23p = sinl(lambda)*sinl(beta)*cosl(beta);
		      e33p = -powl(cosl(beta),2);
		      
		      resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
			      n2*(n1*e21p+n2*e22p+n3*e23p)+
			      n3*(n1*e31p+n2*e32p+n3*e33p));
		      
		      if ((1-cosTheta)==0.0)
			resp = 0.0;  // Check if this is sensible
		      else
			resp = 1.0L/(2.0L*(1.0L-cosTheta))*(resp); 
		      
		      e11c = sin(2*lambda)*sin(beta);
		      e21c = -cos(2*lambda)*sin(beta);
		      e31c = -sin(lambda)*cos(beta);
		      
		      e12c = -cos(2*lambda)*sin(beta);
		      e22c = -sin(2*lambda)*sin(beta);
		      e32c = cos(lambda)*cos(beta);
		      
		      e13c = -sin(lambda)*cos(beta);
		      e23c = cos(lambda)*cos(beta);
		      e33c  = 0;
		      
		      resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
			      n2*(n1*e21c+n2*e22c+n3*e23c)+
			      n3*(n1*e31c+n2*e32c+n3*e33c));
		      
		      if ((1-cosTheta)==0.0)
			resc = 0.0;  // Check if this is sensible
		      else
			resc = 1.0L/(2.0L*(1.0L-cosTheta))*(resc); 

		      
		      dt = (control->psr[p].obs[i].sat - control->gw[kk].gwcsEpoch.dval)*86400.0L;
		      width = control->gw[kk].gwcsWidth.dval*86400.0;
		      width_day = control->gw[kk].gwcsWidth.dval;

		      if (control->psr[p].obs[i].sat < control->gw[kk].gwcsEpoch.dval-width_day/2.0)
			gwRes_p[i] = gwRes_c[i] = 0;
		      else if (control->psr[p].obs[i].sat <= control->gw[kk].gwcsEpoch.dval)
			{
			  gwRes_p[i] = control->gw[kk].gwcsAmp1.dval*(3.0/4.0*(pow(0.5*width,4.0/3.0)-pow(fabs(dt),4.0/3.0))-pow(0.5*width,1.0/3.0)*(dt+0.5*width));
			  gwRes_c[i] = control->gw[kk].gwcsAmp2.dval*(3.0/4.0*(pow(0.5*width,4.0/3.0)-pow(fabs(dt),4.0/3.0))-pow(0.5*width,1.0/3.0)*(dt+0.5*width));
			}
		      else if (control->psr[p].obs[i].sat <= control->gw[kk].gwcsEpoch.dval+width_day/2.0)
			{
			  gwRes_p[i] = control->gw[kk].gwcsAmp1.dval*(3.0/4.0*(pow(0.5*width,4.0/3.0)+pow(fabs(dt),4.0/3.0))-pow(0.5*width,1.0/3.0)*(dt+0.5*width));
			  gwRes_c[i] = control->gw[kk].gwcsAmp2.dval*(3.0/4.0*(pow(0.5*width,4.0/3.0)+pow(fabs(dt),4.0/3.0))-pow(0.5*width,1.0/3.0)*(dt+0.5*width));
			}
		      else
			{
			  gwRes_p[i]=-0.25*(pow(0.5,1.0/3.0)*control->gw[kk].gwcsAmp1.dval*pow(width,4.0/3.0));
			  gwRes_c[i]=-0.25*(pow(0.5,1.0/3.0)*control->gw[kk].gwcsAmp2.dval*pow(width,4.0/3.0));
			}
		      printf("Have gwcs %g %g %g %g %g %g\n",gwRes_p[i],resp,gwRes_c[i],resc,(double)control->psr[p].obs[i].sat,(double)control->gw[kk].gwcsEpoch.dval);
		      gwRes[i] = gwRes_p[i]*resp+gwRes_c[i]*resc;
		      mean+=gwRes[i];
		    }

		}
		
	    }
	
	  mean /= (double)control->psr[p].nToAs;
	  for (i=0;i<control->psr[p].nToAs;i++)
	    {
	      epochs[i]=(double)control->psr[p].obs[i].sat;
	      offsets[i] = (double)((gwRes[i]-mean));
	      printf("offsetsGW = %g\n",offsets[i]);
	    }
	  //      exit(1);
	  // remove quadratic to make the total variation smaller.
	  TKremovePoly_d(epochs,offsets,control->psr[p].nToAs,2);
	  printf("Writing corr\n");
	  toasim_write_corrections(corr,header,file);
	  printf("Done\n");
	  fclose(file);
	}
    }
  if (control->gw[kk].type==5)
    {
      free(h0_g);
      free(omega_g);
      free(angpol_g);
      free(lambda_g);
      free(beta_g);
      free(skyInc_g);
    }
}

double calcDiffractiveScint(controlStruct *control,int s0,int j,int sys)
{
  double nScintD;
  int nScint;
  double fFactor;
  double expon;
  int i,jj;
  double f0;
  double df,dfp;
  double tau,taup;
  double freq;
  double fillFactor=0.2; // From Shannon paper on jitter (eq 14)
  double tobs;
  double bw;
  double diff;
  int p;
  int r,b;
  int s;
  double flux;
  int fClosest;

  if (s0==-1)
    {
      p=sys;
      s=-1;
      freq=control->psr[p].obs[j].freq.dval;
      r=control->psr[p].obs[j].rcvrNum;
      b=control->psr[p].obs[j].beNum;      
    }
  else
    {
      p=control->sched[s0].obs[j].psrNum;
      s=control->sched[s0].obs[j].obsSysNum;
      if (s != -1)
	{
	  freq=control->obsSys[s].freq[sys].dval*1e6;
	  r=control->obsSys[s].rcvrNum[sys];
	  b=control->obsSys[s].beNum[sys];
	}
      else
	{
	  freq=control->sched[s0].obs[j].freq.dval*1e6;
	  r=control->sched[s0].obs[j].rcvrNum;
	  b=control->sched[s0].obs[j].beNum;
	}

    }

  tobs = control->psr[p].obs[control->psr[p].nToAs].tobs.dval;
  printf("Calc scint: tobs = %g ntoas = %d p=%d j=%d freq=%g\n",tobs,control->psr[p].nToAs,p,j,freq);

  bw = control->be[b].bw.dval*1e6;

  // Finding closest flux density value;

  fClosest=0;
  diff=fabs(control->psr[p].freqFlux[0].dval-freq);
  for (i=1;i<control->psr[p].nFlux;i++)
    {
      if (diff > fabs(control->psr[p].freqFlux[i].dval-freq))
	{
	  diff=fabs(control->psr[p].freqFlux[i].dval-freq);
	  fClosest=i;
	}
    }
  flux = control->psr[p].flux[fClosest].dval*1e-3;
  printf("Flux = %g\n",flux);


  df = control->psr[p].diff_df.dval;
  tau = control->psr[p].diff_ts.dval;
  printf("df = %g, tau = %g\n",df,tau);
   
  f0 =  control->psr[p].diff_dfFreq.dval*1e6;
  dfp = df*pow(freq/f0,4.4); // CHECK THESE NUMBERS  SHOULD THIS BE POSITIVE??

  f0 =  control->psr[p].diff_tsFreq.dval*1e6;
  taup = tau*pow(freq/f0,1.2);
  nScintD = (1+fillFactor*tobs/taup)*(1+fillFactor*bw/dfp);

  printf("Diffractive scintillation %g\n",nScintD);
  nScint = (int)(nScintD+0.5);
  if (nScint < 1) nScint = 1; // CHECK THIS

  expon=0;
  for (jj=0;jj<nScint;jj++)
    expon += -flux*log(TKranDev(&control->seed));  
  expon/=(double)nScint;
  printf("Expon = %g %g\n",expon,expon/flux);
  return expon/flux;
}

double calculateToaErrRadiometer(controlStruct *control,int s0,int j,int sys,double scale,int realisation)
{
  static int counter=0;
  double flux;
  double freq;
  double tsys=20; // Recevier parameter
  double gain=1.1; // Receiver parameter
  double tobs; // Observation parameter  THIS SHOULD BE UPDATED EACH OBSERVATION
  // AND RECORDED IN THE PSR OBSERVATION VALUE!
  double deltaf; // Observation parameter
  double tsky=0; // External information
  int nbin=4096; // Backend parameter
  double prof[nbin];
  double templ[nbin];
  double phase;
  int i;
  char expression[MAX_STRLEN];
  double radNoise;
  int errorFlag=0;
  double sum=0.0;
  double toaErr;
  int p;
  int r,b;
  int s=-1;
  int fClosest;
  FILE *fout;
  char fileOut[128];

  double diff;
  int closestProf=-1;
  int tskyClosest;

  if (s0!=-1)
    {
      p = control->sched[s0].obs[j].psrNum;
      s =control->sched[s0].obs[j].obsSysNum;
    }
  else
    {
      p = sys;
    }
  tobs = control->psr[p].obs[control->psr[p].nToAs].tobs.dval;
  if (s0 == -1)
    {
      freq = control->psr[p].obs[j].freq.dval;
      r=control->psr[p].obs[j].rcvrNum;
      b=control->psr[p].obs[j].beNum;
      printf("Using: %g %d %d\n",freq,r,b);
    }
  else
    {
      if (s != -1)
	{
	  freq=control->obsSys[s].freq[sys].dval;
	  r=control->obsSys[s].rcvrNum[sys];
	  b=control->obsSys[s].beNum[sys];
	}
      else
	{
	  freq=control->sched[s0].obs[j].freq.dval;
	  r=control->sched[s0].obs[j].rcvrNum;
	  b=control->sched[s0].obs[j].beNum;
	}
    }

  // Finding closest Tsky value;

  tskyClosest=0;
  diff=fabs(control->psr[p].freqTsky[0].dval-freq);
  for (i=1;i<control->psr[p].nTsky;i++)
    {
      if (diff > fabs(control->psr[p].freqTsky[i].dval-freq))
	{
	  diff=fabs(control->psr[p].freqTsky[i].dval-freq);
	  tskyClosest=i;
	}
    }
  tsky = control->psr[p].tsky[tskyClosest].dval;
  printf("Using tsky = %g\n",tsky);

  // Finding closest profile file
  //  printf("profile = %d\n",control->psr[p].nProfileFile);
  if (control->psr[p].nProfileFile > 0)
    {
      closestProf=0;
      diff=fabs(control->psr[p].freqProfileFile[0].dval-freq);
      for (i=1;i<control->psr[p].nProfileFile;i++)
	{
	  if (diff > fabs(control->psr[p].freqProfileFile[i].dval-freq))
	    {
	      diff=fabs(control->psr[p].freqProfileFile[i].dval-freq);
	      closestProf=i;
	    }
	}
      //      printf("closest profile = %d\n",closestProf);

    }
  

  // Finding closest flux density value;
  fClosest=0;
  diff=fabs(control->psr[p].freqFlux[0].dval-freq);
  for (i=1;i<control->psr[p].nFlux;i++)
    {
      if (diff > fabs(control->psr[p].freqFlux[i].dval-freq))
	{
	  diff=fabs(control->psr[p].freqFlux[i].dval-freq);
	  fClosest=i;
	}
    }
  flux = control->psr[p].flux[fClosest].dval*1e-3*scale;
  //  printf("Using flux = %g %g %d %d\n",freq,control->psr[p].flux[fClosest].dval,fClosest,control->psr[p].nFlux);

  tsys = control->rcvr[r].tsys.dval;
  gain = control->rcvr[r].gain.dval;
  nbin = control->be[b].nbin.dval;
  deltaf = control->be[b].bw.dval*1e6;
  printf("Trying: %d %d %d %d %d %g %g %d %g %g\n",p,r,b,s0,j,tsys,gain,nbin,deltaf,flux);
  radNoise = (tsys+tsky)/(gain)/sqrt(2.0*(tobs/nbin)*deltaf);
  printf("radNoise = %s %g %g %d %d\n",control->psr[p].name,radNoise,tobs,p,j);
  
  if (control->psr[p].nProfileFile==0)
    {      
      for (i=0;i<nbin;i++)
	{
	  sprintf(expression,"x=%g; v=%s;",i/(double)nbin,control->psr[p].profileEqn);
	  errorFlag = runEvaluateExpression(expression,control);      
	  prof[i] = variable[0].value;
	}
    }
  else
    loadProfileFromFile(prof,templ,nbin,control->psr[p].profileFile[closestProf]);

  for (i=0;i<nbin;i++)
    {
      sum+=prof[i];
    }
  for (i=0;i<nbin;i++)
    {
      //      templ[i]=prof[i]; //*(flux/sum)*nbin;

      prof[i] = prof[i]*(flux/sum)*nbin + TKgaussDev(&(control->seed))*radNoise;
      //       printf("prof = %d %g %g\n",i,variable[0].value,prof[i]);
    }
  sprintf(fileOut,"%s.%g.%g.%d.%d.%d.%d.%d.prof",control->psr[p].name,freq,tobs,s0,j,sys,realisation,counter);
  fout = fopen(fileOut,"w");
  for (i=0;i<nbin;i++)
    fprintf(fout,"%d %g\n",i,prof[i]);
  fclose(fout);


  //
  //   exit(1);
  // Now cross correlate the template with itself to get the error bar size (use a fast fft method)
  // Either do this for every observation and every scintillation state and receivers or 
  // just do it once for each profile and then scale for different noise levels
  toaErr = getToaErr(prof,templ,nbin)*control->psr[p].p0;
  //  printf("toaErr = %g\n",toaErr);
  //  exit(1);
  counter++;
  return toaErr;
}

double getToaErr(double *prof,double *templ,int nbin)
{
  double toaErr;//=0.1e-6;
  double shift,eshift,snr,esnr,b,errb;
  int ngood;

  fftfit(prof,templ,nbin,&shift,&eshift,&snr,&esnr,&b,&errb,&ngood);
  toaErr = eshift/nbin;
  printf("toaErr = %g in phase\n",toaErr);

  return toaErr;
}




void fftfit(double *prof,double *standard,int nmax,double *shift,double *eshift,
	    double *snr,double *esnr,double *b,double *errb,int *ngood)
{
  int nh,i,j,k;
  double sum,ave,errtau;
  double s[10000],phi[10000];   /* Transform of template */
  double p[10000],theta[10000]; /* Transform of profile  */
  double tmp[10000],r[10000],fac;
  double tau,s1,s2,s3,cosfac,sq,rms;
  int isum;
  int nsum,ntries;
  double edtau,ftau,a,fa,dtau,fb;
  int low,high;

  /* Obtain Fourier transform of template */
  fft(standard,nmax,s,phi);

  nh = nmax/2;
  sum=0.0;
  for (i=nh/2+1;i<=nh;i++)
    sum+=s[i];
  ave=2.0*sum/nh;

  for (i=1;i<=nh;i++)
    if (s[i] < ave) break;
  *ngood=i-1;
  //  printf("Test 1: %g %d\n",ave,*ngood);

  /*  *ngood = 89; */ /* SET IT EQUAL TO FORTRAN VERSION __ WHAT IS WRONG ? */

  //  printf("Ngood = %d %d %f %f\n",*ngood,nmax,ave,sum);

  /* Obtain Fourier transform of profile */
  fft(prof,nmax,p,theta);

  for (k=0;k<nh;k++)
    {
      tmp[k]=p[k]*s[k];
      r[k]=theta[k]-phi[k];
    }

  fac = nmax/(2.0*M_PI);
  //  printf("Test 2: %g\n",fac);
  //  printf("fac = %f\n",fac);
  fccf(tmp,r,shift,nmax);
  //  printf("From fccf, shift = %f\n",*shift);
  //  printf("Test 3: %g\n",*shift);
  tau = *shift;
  for (isum=5;isum<99;isum++)
    {
      nsum=(int)pow(2,isum);
      if (nsum>nh) break;
      dtau = (2.0*M_PI)/(nsum*5);
      edtau = 1.0/(2.0*nsum+1.0);
      if (nsum > (nh/2.+0.5)) edtau = 1.0e-4;

      ntries = 0;
      low = -1;
      high = -1;
      do {
	ftau = dchisqr(tau,tmp,r,nsum);
	ntries++;
	if (ftau < 0)
	  {
	    a=tau;
	    fa=ftau;
	    tau+=dtau;
	    low=1;
	  }
	else
	  {
	    *b=tau;
	    fb=ftau;
	    tau=tau-dtau;
	    high = 1;	  
	  }
	if (ntries>10)
	  {
	    *shift=0.0;
	    *eshift=999.0;
	    *snr=0.0;
	    *esnr=0.0;
	    return;
	  }
      }	while (low!=high);
      tau = zbrent(a,*b,fa,fb,edtau,tmp,r,nsum);
    }
  s1=0.0;
  s2=0.0;
  s3=0.0;
  for (k=1;k<=nh;k++)  /* SHOULD THIS START FROM 0 */
    {
      cosfac = cos(-r[k]+k*tau);
      s1=s1+tmp[k]*cosfac;
      s2=s2+s[k]*s[k];
      s3=s3+k*k*tmp[k]*cosfac;
    }

  *b=s1/s2;
  //  printf("b= %g\n",*b);
  s1=0.0;
  
  for (k=1;k<=nh;k++)
    {
      sq = p[k]*p[k]-2.0*(*b)*p[k]*s[k]*cos(r[k]-k*tau)+pow((*b)*s[k],2);
      s1+=sq;
    }
  rms =sqrt(s1/nh);
  *errb = rms/sqrt(2.0*s2);
  errtau = rms/sqrt(2.0*(*b)*s3);
  *snr = 2.0*sqrt(2.0*nh)*(*b)/rms;
  *shift = fac*tau;
  *eshift = fac*(errtau);
  *esnr = *snr*(*errb)/(*b);
  
  return;
}

double dchisqr(double tau,double *tmp,double *r,int nsum)
{
  int k;
  double s;

  s=0.0;
  for (k=1;k<=nsum;k++)
    s+=k*tmp[k]*sin(-r[k]+k*tau);

  return s;

}


void fccf(double *amp,double *pha,double *shift,int nprof)
{
  /* nprof = 64 in the fccf.f program ??? */
  int nh,i;
  double ccf[5000];
  double cmax,rc,fb,fa,fc;
  double cReal[5000],cImag[5000];
  int imax,ia,ic;

  nprof = 64;

  nh = nprof/2;
  ccf[0]=0;
  ccf[1]=0;
  for (i=1;i<=nh/2;i++)
    {
      ccf[2*i] = amp[i]*cos(pha[i]);
      ccf[2*i+1] = amp[i]*sin(pha[i]);
      ccf[2*nprof-(2*i)] = amp[i]*cos(pha[i]);
      ccf[2*nprof-(2*i-1)] = -amp[i]*sin(pha[i]);
    }
  for (i=nh/2+1;i<=nh;i++)
    {
      ccf[2*i]=0.0;
      ccf[2*i+1]=0.0;
      ccf[2*nprof-(2*i)]=0.0;
      ccf[2*nprof-(2*i-1)]=0.0;
    }

  four1(ccf-1,nprof,-1);

  for (i=0;i<=nprof;i++)
    {
      cReal[i]=ccf[2*i];
      cImag[i]=ccf[2*i+1];
    }
  
  cmax = -1.0e30;
  for (i=0;i<nprof;i++)
    {
      rc = ccf[2*i];
      if (rc>cmax)
	{
	  cmax= rc;
	  imax = i;
	}
    }
  fb = cmax;
  //  printf("fb = %f\n",fb);
  ia = imax-1;
  if (ia==-1) ia=nprof-1;
  fa=ccf[ia*2];
  ic=imax+1;
  if (ic==nprof) ic=0;
  fc=ccf[ic*2];
  if ((2*fb-fc-fa)!=0.0)
    *shift=imax+0.5*(fa-fc)/(2*fb-fc-fa);
  else
    *shift=imax;
  if (*shift>nh) *shift-=nprof;
  *shift=*shift*(2.0*M_PI)/nprof;

}

void fft(double *y,int nmax,double *amp,double *pha)
{
  double temp[5000];
  int i;

  for (i=0;i<nmax;i++)
    {
      temp[2*i]  =y[i];
      temp[2*i+1]=0.0;  /* profile is real */
    }
  four1(temp-1,nmax,1);
  for (i=0;i<nmax/2;i++)
    {
      amp[i] = 2*sqrt(pow(temp[2*i],2)+pow(temp[2*i+1],2));
      pha[i] = atan2(temp[2*i+1],temp[2*i]);
      /*      printf("%d %f\n",i,temp[2*i]);*/
    }
  /*  exit(1); */
}


void four1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

double zbrent(double x1,double x2,double f1,double f2,double tol,double *tmp,
	      double *pha,int nsum)
{
  double a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r;
  int iter;
  int itmax=100;
  double eps=6.0e-8,s;

  a=x1;
  b=x2;
  fa=f1;
  fb=f2;
  fc=fb;

  for (iter=1;iter<=itmax;iter++)
    {
      if (fb*fc>0)
	{
	  c=a;
	  fc=fa;
	  d=b-a;
	  e=d;
	}
      if (fabs(fc)<fabs(fb))
	{
	  a=b;
	  b=c;
	  c=a;
	  fa=fb;
	  fb=fc;
	  fc=fa;
	}
      tol1 = 2.0*eps*fabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if (fabs(xm)<=tol1 || fb==0)
	return b;
      
      if (fabs(e)>=tol1 && fabs(fa)>fabs(fb))
	{
	  s=fb/fa;
	  if (a==c)
	    {
	      p=2.0*xm*s;
	      q=1.0-s;
	    }
	  else
	    {
	      q=fa/fc;
	      r=fb/fc;
	      p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	      q=(q-1.0)*(r-1.0)*(s-1.0);     
	    }
	  if (p>0.0)q=-q;
	  p=fabs(p);
	  if (2.0*p < min(3.0*xm*q-fabs(tol1*q),fabs(e*q)))
	    {
	      e=d;
	      d=p/q;
	    }
	  else
	    {
	      d=xm;
	      e=d;
	    }
	}
      else
	{
	  d=xm;
	  e=d;
	}
      a=b;
      fa=fb;
      if (fabs(d) > tol1)
	b=b+d;
      else
	b=b+sign(tol1,xm);

      fb=dchisqr(b,tmp,pha,nsum);
    }
  return b;
}

double sign(double a,double b)
{
  if (b>0)
    return fabs(a);
  else
    return -fabs(a);
}

double min(double a,double b)
{
  if (a<b)
    return a;
  else
    return b;
}

void loadProfileFromFile(double *prof,double *templ,int nbin,char *fname)
{
  tmplStruct template;
  int i;

  initialiseTemplate(&template);
  readTemplate(fname,&template);
  
  for (i=0;i<nbin;i++)
    prof[i] = evaluateTemplateChannel(&template,(double)i/(double)nbin,0,0,0);

  //  template.channel[0].pol[0].comp[0].concentration = 500;
  for (i=0;i<nbin;i++)
    templ[i] = evaluateTemplateChannel(&template,(double)i/(double)nbin,0,0,0.2);


}

void readT2TimFile(controlStruct *control,int or,int t2Num,int r)
{
  FILE *fin;
  int p=-1;
  int i,j;
  char line[2048];
  char first[1024];

  // Find correct pulsar


  printf("Looking for pulsar\n");
  for (i=0;i<control->npsr;i++)
    {
      if (strcmp(control->psr[i].name,control->obsRun[or].T2Tim[t2Num].psrName)==0)
	{p=i; break;}
    }
  printf("Found pulsar: %d\n",p);

  if (!(fin = fopen(control->obsRun[or].T2Tim[t2Num].fileName,"r")))
    {
      printf("Unable to open file: %s\n",control->obsRun[or].T2Tim[t2Num].fileName);
      finishOff(control);
    }
  while (!feof(fin))
    {
      if (fgets(line,2048,fin)!=NULL)
	{
	  if (line[0]!='#')
	    {
	      sscanf(line,"%s",first);
	      if (strcasecmp(first,"FORMAT")==0 || (strcasecmp(first,"MODE")==0) ||
		  strcasecmp(first,"EFAC")==0 || strcasecmp(first,"EQUAD")==0)
		printf("Ignoring: %s\n",line);
	      else
		{
		  long double sat;
		  double freq,toaErr;
		  char telCode[1024];
		  char temp[1024];
		  int nobs;
		  nobs = control->psr[p].nToAs;

		  if (sscanf(line,"%s %lf %Lf %lf %s",temp,&freq,&sat,&toaErr,telCode)==5)
		    {
		      if (sat > control->maxT)
			control->maxT = (double)sat;
		      if (sat < control->minT)
			control->minT = (double)sat;

		      // Now look for flags
		      if (strstr(line,"-tobs")!=NULL)
			{
			  char *tok;
			  char tt[1024];
			  strcpy(tt,strstr(line,"-tobs"));
			  tok = strtok(tt," \n");
			  tok = strtok(NULL," \n");
			  sscanf(tok,"%lf",&(control->psr[p].obs[nobs].tobs.dval));
			  control->psr[p].obs[nobs].tobs.set = 1;
			}
		      if (strstr(line,"-rcvr")!=NULL)
			{
			  char *tok;
			  char tt[1024];
			  int k,r0=-1;
			  strcpy(tt,strstr(line,"-rcvr"));
			  tok = strtok(tt," \n");
			  tok = strtok(NULL," \n");
			  printf("Got recevier : %s\n",tok);
			  for (k=0;k<control->nRCVR;k++)
			    {
			      if (strcmp(control->rcvr[k].name,tok)==0)
				{r0 = k; break;}
			    }
			  if (r0 == -1)
			    {
			      printf("Unable to find receiver with name %s\n",tok);
			      finishOff(control);
			    }
			  control->psr[p].obs[nobs].rcvrNum = r0;
			}
		      if (strstr(line,"-ptaSimbe")!=NULL)
			{
			  char *tok;
			  char tt[1024];
			  int k,b0=-1;
			  strcpy(tt,strstr(line,"-ptaSimbe"));
			  tok = strtok(tt," \n");
			  tok = strtok(NULL," \n");
			  printf("Got recevier : %s\n",tok);
			  for (k=0;k<control->nBE;k++)
			    {
			      if (strcmp(control->be[k].name,tok)==0)
				{b0 = k; break;}
			    }
			  if (b0 == -1)
			    {
			      printf("Unable to find backend with name %s\n",tok);
			      finishOff(control);
			    }
			  control->psr[p].obs[nobs].beNum = b0;
			}

		      strcpy(control->psr[p].obs[nobs].tel,telCode);
		      control->psr[p].obs[nobs].freq.dval = freq;
		      control->psr[p].obs[nobs].freq.set = 1;
		      control->psr[p].obs[nobs].sat = sat;
		      control->psr[p].obs[nobs].satSet = 1;
		      fillDval(&(control->obsRun[or].T2Tim[t2Num].efac),control);
		      fillDval(&(control->obsRun[or].T2Tim[t2Num].equad),control);
		      control->psr[p].obs[nobs].efac.dval = control->obsRun[or].T2Tim[t2Num].efac.dval;
		      control->psr[p].obs[nobs].equad.dval = control->obsRun[or].T2Tim[t2Num].equad.dval;

		      if (control->obsRun[or].T2Tim[t2Num].toaErr.set==0)
			{
			  control->psr[p].obs[nobs].toaErr.dval = toaErr/1e6;
			  control->psr[p].obs[nobs].toaErr.set = 1;
			}
		      else
			{
			  double scale=1;

			  if (strcmp(control->obsRun[or].T2Tim[t2Num].toaErr.inVal,"radiometer")==0)
			    {
			      printf("Trying radiometer through t2file:\n");
			      printf("need s0, sys\n");
			      printf("Have observation %d, freq = %g\n",nobs,freq);

			      if (control->psr[p].setDiff_df==1
				  && control->psr[p].setDiff_ts==1)
				//				scale = calcDiffractiveScint(control,s0,j,sys);
				scale = calcDiffractiveScint(control,-1,nobs,p);
			      else
				scale=1;

			      
			      control->psr[p].obs[nobs].toaErr.dval = calculateToaErrRadiometer(control,-1,nobs,p,scale,r);
			    }
			  else
			    {
			      fillDval(&(control->obsRun[or].T2Tim[t2Num].toaErr),control);
			      control->psr[p].obs[nobs].toaErr.dval = control->obsRun[or].T2Tim[t2Num].toaErr.dval;
			    }
			  control->psr[p].obs[nobs].toaErr.set = 1;
			}
		      (control->psr[p].nToAs)++;
		    }
		}
	    }
	}
    }
  fclose(fin);

  //	  printf("Processing psr %s\n",control->obsRun[p].T2Tim[k].psrName);
  //	  readT2TimFile(control,p,k);
}

void readObservatoryPositions(controlStruct *control)
{
  FILE *fin;
  int n=0,i;
  char fname[1024];
  char line[1024];

  sprintf(fname,"%s/observatory/observatories.dat",getenv("TEMPO2"));
  printf("Reading observatory positions from %s\n",fname);
  if (!(fin = fopen(fname,"r")))
    {
      printf("Unable to read observatories.dat file from %s. Perhaps you need to install tempo2\n",fname);
      finishOff(control);
    }
  while (!feof(fin))
    {
      if (fgets(line,1024,fin)!=NULL)
	{
	  if (line[0]!='#')
	    {
	      if (sscanf(line,"%lf %lf %lf %s %s",&(control->observatory[n].posX),
			 &(control->observatory[n].posY),&(control->observatory[n].posZ),
			 control->observatory[n].name1,control->observatory[n].name2)==5)
		n++;
	    }
	}
    }
  fclose(fin);
  control->nObservatory=n;
  printf("Completed reading observatory positions\n");
}

long double getTimeHA(controlStruct *control,double ha,long double sat,int telID,long double ra)
{
  long double solsid  = 1.002737909;
  long double obslong,almst;
  long double ha0;
  double x,y,z;
  long double sat2;

  x = control->observatory[telID].posX;
  y = control->observatory[telID].posY;
  z = control->observatory[telID].posZ;
  obslong = atan2(y,x)*180.0/M_PI; // In degrees
  almst = fortran_mod((sat-47892.0)*solsid+0.276105324+obslong/360.0,(long double)1.0);
  printf("Coords %g %g %g %g\n",x,y,z,(double)obslong);
  printf("RA = %g\n",(double)ra);
  ra = ra/360.0;
  ha0 = almst-ra;
  if (ha0 < 0.0) ha0 +=1.0;
  printf("almst = %g, ha0 = %g, ha = %g\n",(double)almst,(double)ha0,(double)ha);  
  ha0 -= ha/24.0; // This is the request ha
  sat2 = sat + 1.0-ha0;
  //  exit(1);

  return sat2;
}

long double fortran_mod(long double a,long double p)
{
  long double ret;
  ret = a - (int)(a/p)*p;
  return ret;
}

double BTmodel(double pb,double ecc,double a1,double t0,double om,long double t)
{
  long double tt0;
  double edot=0;
  double pbdot=0;
  double xpbdot=0;
  double xdot=0.0;
  double omdot=0;
  double asini;
  double omega;
  double gamma=0;
  double torb;
  double orbits,phase;
  double ep,dep,bige,som,com,alpha,beta,sbe,cbe,q,r,s,tt;
  int norbits;

  tt0 = (t - t0)*SECDAY;

  pb     = pb*SECDAY;
  ecc    = ecc + edot*tt0;

  asini  = (double)(a1 + xdot*tt0);
  printf("asini = %g %g %g\n",a1,xdot,(double)tt0);
  omega  = (om + omdot*tt0/(SECDAY * 365.25))/(180.0/M_PI);

  torb = 0.0;
  orbits = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2); 
  norbits = (int)orbits;
  if (orbits < 0.0) norbits--;
  
  phase = 2.0*M_PI * (orbits-norbits);

  /* Using Pat Wallace's method of solving Kepler's equation -- code based on bnrybt.f */
  ep = phase + ecc*sin(phase)*(1.0+ecc*cos(phase));

  /* This line is wrong in the original tempo: should be inside the do loop */
  /*  denom = 1.0 - ecc*cos(ep);*/
  
  do {
    dep = (phase - (ep-ecc*sin(ep)))/(1.0 - ecc*cos(ep));
    ep += dep;
  } while (fabs(dep) > 1.0e-12);
  bige = ep;

  tt = 1.0-ecc*ecc;
  som = sin(omega);
  com = cos(omega);

  alpha = asini*som;
  beta = asini*com*sqrt(tt);
  sbe = sin(bige);
  cbe = cos(bige);
  q = alpha * (cbe-ecc) + (beta+gamma)*sbe;
  r = -alpha*sbe + beta*cbe;
  s = 1.0/(1.0-ecc*cbe);

  torb = -q+(2*M_PI/pb)*q*r*s + torb;
  printf("torb = %g %g %g %g\n",torb,ecc,com,tt);
  return torb;
}
