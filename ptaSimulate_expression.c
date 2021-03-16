#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ptaSimulate.h"
#include "evaldefs.h"
#include "T2toolkit.h"

int runEvaluateExpression(char *expression,controlStruct *control)
{
  int result;
  char express2[1024];
  char temp[1024],temp2[1024];
  char *tok,*tok2,*e;
  int add;
  double change;
  char changeStr[1024];
  printf("Evaluating expression: %s\n",expression);
  if (strstr(expression,"ran(")!=NULL)
    {
      //      printf("Must process\n");
      strcpy(express2,expression);
      while (strstr(express2,"ran(")!=NULL)
	{
	  strcpy(temp,express2);
	  tok = strstr(temp,"ran(");
	  express2[tok-temp]='\0';
	  // Now split between brackets to get the information
	  tok2 = strtok(tok,"(");
	  tok2 = strtok(NULL,";)");
	  // Should check to see the type of random variable
	  if (strcmp(tok2,"gaussian")==0)
	    {
	      double v;
	      char param[1024];
	      strcpy(param,tok2);
	      change = TKgaussDev(&(control->seed));
	      sprintf(changeStr,"%20.20g",change);  // THIS IS A BIG PROBLEM!! HOW TO SET THIS CORRECTLY???
	    }
	  else if (strcmp(tok2,"linear")==0)
	    {
	      char param[1024];
	      strcpy(param,tok2);
	      change = TKranDev(&(control->seed));
	      sprintf(changeStr,"%20.20g",change);  // THIS IS A BIG PROBLEM!! HOW TO SET THIS CORRECTLY???
	    }
	  else
	    {
	      printf("Unknown random parameter: %s\n",tok2);
	      finishOff(control);
	    }

	  add = tok2+strlen(tok2)+1-temp;
	  strcpy(temp2,temp + add);

	  strcat(express2,changeStr);
	  strcat(express2,temp2);
	}
      strcpy(expression,express2);
    }

  if (strstr(expression,"fdist(")!=NULL)
    {
      FILE *fin;
      int n,r;
      double v[MAX_FDIST_VALS];
      printf("Must process fdist\n");
      strcpy(express2,expression);
      while (strstr(express2,"fdist(")!=NULL)
	{
	  strcpy(temp,express2);
	  tok = strstr(temp,"fdist(");
	  express2[tok-temp]='\0';
	  // Now split between brackets to get the information
	  tok2 = strtok(tok,"(");
	  tok2 = strtok(NULL,";)");
	  // Should check to see the type of random variable
	  if (!(fin = fopen(tok2,"r")))
	    {
	      printf("Unable to open file: %s in fdist\n",tok2);
	      finishOff(control);
	    }
	  n=0;
	  while (!feof(fin))
	    {
	      if (fscanf(fin,"%lf",&v[n])==1)
		{
		  n++;
		  if (n > MAX_FDIST_VALS)
		    {
		      printf("Too many fdist values in %s\n",tok2);
		      finishOff(control);
		    }
		}
	    }
	  fclose(fin);
	  r = TKranDev(&(control->seed))*n;
	  printf("Random value from %s = %d %g\n",tok2,r,v[r]);
	  sprintf(changeStr,"%20.20g",v[r]);  // THIS IS A BIG PROBLEM!! HOW TO SET 
	  add = tok2+strlen(tok2)+1-temp;
	  strcpy(temp2,temp + add);

	  strcat(express2,changeStr);
	  strcat(express2,temp2);
	}
      strcpy(expression,express2);
    }



  strcpy(express2,expression);


  printf("Evaluated to: %s\n",express2);
  //      printf("Evalulating expression: %s\n",express2);
  result = evaluateExpression(express2);
  return result;
}

//
// Routine to remove random numbers that should 
// only be set once
//
int changeRandomOnce(char *expression,controlStruct *control)
{
  int result;
  char express2[1024];
  char temp[1024],temp2[1024];
  char *tok,*tok2,*e;
  int add;
  double change;
  char changeStr[1024];

  printf("Expression at top = %s\n",expression);

  if (strstr(expression,"ranOnce(")!=NULL)
    {
      //      printf("Must process\n");
      strcpy(express2,expression);
      while (strstr(express2,"ranOnce(")!=NULL)
	{
	  strcpy(temp,express2);
	  tok = strstr(temp,"ranOnce(");
	  express2[tok-temp]='\0';
	  // Now split between brackets to get the information
	  tok2 = strtok(tok,"(");
	  tok2 = strtok(NULL,";)");
	  // Should check to see the type of random variable
	  if (strcmp(tok2,"gaussian")==0)
	    {
	      double v;
	      char param[1024];
	      strcpy(param,tok2);
	      change = TKgaussDev(&(control->seed));
	      sprintf(changeStr,"%20.20g",change);  // THIS IS A BIG PROBLEM!! HOW TO SET THIS CORRECTLY???
	    }
	  else if (strcmp(tok2,"linear")==0)
	    {
	      char param[1024];
	      strcpy(param,tok2);
	      change = TKranDev(&(control->seed));
	      sprintf(changeStr,"%20.20g",change);  // THIS IS A BIG PROBLEM!! HOW TO SET THIS CORRECTLY???
	    }
	  else
	    {
	      printf("Unknown random parameter: %s\n",tok2);
	      finishOff(control);
	    }

	  add = tok2+strlen(tok2)+1-temp;
	  strcpy(temp2,temp + add);

	  strcat(express2,changeStr);
	  strcat(express2,temp2);
	}
      strcpy(expression,express2);
    }

  //  strcpy(express2,expression);
  printf("Expression from ranOnce = %s\n",expression);
  //  printf("express2 = %s\n",express2);


}

int checkProbability(valStruct in,controlStruct *control)
{
  int res;
  char process[1024];
  strcpy(process,"v = ");
  strcat(process,in.inVal);
  res = runEvaluateExpression(process,control);
  if (variable[0].value > 0.1)
    return 1;
  else
    return 0;
}
