#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ptaSimulate.h"
#include "toasim.h"

void createOutliers(controlStruct *control,int r)
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
  double offsets[MAX_TOAS]; // Will change to doubles - should use malloc
  // Create a set of corrections.
  toasim_corrections_t* corr = (toasim_corrections_t*)malloc(sizeof(toasim_corrections_t));

  int dd;
  int k;
  int fail;
  corr->offsets=offsets;
  corr->params=""; // Normally leave as NULL. Can store this along with each realisation. 
  // Same length string in every iteration - defined in r_param_length see below
  corr->a0=0; // constant
  corr->a1=0; // a1*x
  corr->a2=0; // a2*x*X
  
  nit = 1;

  for (p=0;p<control->npsr;p++)
    {
      header = toasim_init_header();
      strcpy(header->short_desc,"addOutliers");
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
      
      sprintf(fname,"%s/workFiles/real_%d/%s.addOutliers",control->name,r,control->psr[p].name);
      // First we write the header...
      file = toasim_write_header(header,fname);
      
      for (j=0;j<control->psr[p].nToAs;j++){
	offsets[j]=0.0;
	//	printf("outliers %d %d\n",control->psr[p].obs[j].outlierAmp.set,control->psr[p].nToAs);
	if (control->psr[p].obs[j].outlierAmp.set==1)
	  {
	    fillDval(&(control->psr[p].obs[j].outlierAmp),control);
	    fillDval(&(control->psr[p].obs[j].outlierProb),control);
	    fail =0;
	    fail = checkProbability(control->psr[p].obs[j].outlierProb,control);
	    if (fail==1)
	      offsets[j] = control->psr[p].obs[j].outlierAmp.dval;
	  }
      }
      printf("Writing corrections\n");
      toasim_write_corrections(corr,header,file);
      printf("Complete writing\n");
      fclose(file);
    }
}
