#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ptaSimulate.h"
#include "toasim.h"
#include "makeRedNoise.h"

void processEphemNoise(controlStruct *control,int r)
{
  int i,p;
  char str[1024];
 
  // Evaluate expressions
  if (control->ephemNoise.gwAmp.set==1)
    {
      double fc,alpha,p0;
      double secperyear = 86400.0*365.25;
      
      alpha = -13.0/3.0;
      fc = 1.0/(((control->maxT-control->minT)/365.25)*2);
      fillDval(&(control->ephemNoise.gwAmp),control);
      p0 = control->ephemNoise.gwAmp.dval;

      control->ephemNoise.alpha.dval = alpha;
      control->ephemNoise.fc.dval = fc;
      control->ephemNoise.p0.dval = pow(p0,2)/12.0/M_PI/M_PI*pow(fc,alpha);//*secperyear*secperyear;
      
    }
  else
    {
      fillDval(&(control->ephemNoise.alpha),control);
      fillDval(&(control->ephemNoise.p0),control);
      fillDval(&(control->ephemNoise.fc),control);
    }
}

void createEphemNoise(controlStruct *control,int r)
{
  int p,i,j;
  FILE *file;
  int npts=1024;
  char fname[MAX_STRLEN];
  toasim_header_t* header;
  toasim_header_t* read_header;
  double offsets[MAX_TOAS]; // should use malloc
  double dx,dy,dz;
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
  rednoisemodel_t* modelx,*modely,*modelz;
  FILE *fout;
  int t;
  char fn[1024];
  long double ra_p,dec_p;
  long double kp[3];            /* Vector pointing to pulsar           */

  sprintf(fname,"%s/setup/useParams",control->name);
  fout = fopen(fname,"a"); 

  for (t=0;t<control->nEphemNoise;t++)
    {
      for (p=0;p<control->npsr;p++)
	{
	  alpha = control->ephemNoise.alpha.dval;
	  beta = 0;
	  p_1yr = control->ephemNoise.p0.dval*secperyear*secperyear;
	  old_fc = control->ephemNoise.fc.dval;
	  cnr_flat = old_fc;
	  corr->offsets=offsets;
	  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
	  // Same length string in every iteration - defined in r_param_length see below
	  corr->a0=0; // constant
	  corr->a1=0; // a1*x
	  corr->a2=0; // a2*x*X
	  
	  header = toasim_init_header();
	  strcpy(header->short_desc,"ephemNoise");
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
	  sprintf(fname,"%s/workFiles/real_%d/%s.ephemnoise.%d",control->name,r,control->psr[p].name,t);
	  file = toasim_write_header(header,fname);

	  if (p==0)
	    {
	      double mjd_start=(double)control->minT;
	      double mjd_end=(double)control->maxT;
	      
	      
	      printf("Setting up the noise model %g %g\n",mjd_start,mjd_end);
	      modelx = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha,beta);
	      modely = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha,beta);
	      modelz = setupRedNoiseModel(mjd_start,mjd_end,npts,nit,p_1yr,alpha,beta);

	      modelx->cutoff=cnr_cut;
	      modelx->flatten=cnr_flat;
	      if(old_fc>0)
		modelx->mode=MODE_T2CHOL;
	      populateRedNoiseModel(modelx,&(control->seed));

	      modely->cutoff=cnr_cut;
	      modely->flatten=cnr_flat;
	      if(old_fc>0)
		modely->mode=MODE_T2CHOL;
	      populateRedNoiseModel(modely,&(control->seed));

	      modelz->cutoff=cnr_cut;
	      modelz->flatten=cnr_flat;
	      if(old_fc>0)
		modelz->mode=MODE_T2CHOL;
	      populateRedNoiseModel(modelz,&(control->seed));
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
	      
	      ra_p = control->psr[p].rajd*M_PI/180.0;
	      dec_p = control->psr[p].decjd*M_PI/180.0;
	      setupPulsar_GWsim(ra_p,dec_p,kp);

	      for (j=0;j<control->psr[p].nToAs;j++){
		dx = getRedNoiseValue(modelx,control->psr[p].obs[j].sat,i);
		dy = getRedNoiseValue(modely,control->psr[p].obs[j].sat,i);
		dz = getRedNoiseValue(modelz,control->psr[p].obs[j].sat,i);

		offsets[j]=dx*kp[0] + dy*kp[1] + dz*kp[2]; // Assume equatorial coordinates
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
}
