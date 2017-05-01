#pragma once
#include "psmcreader.h"

//from psmc
typedef struct{
  int n; // $n$ in psmc.tex. number of intervals equals to $n+1$
  int n_free; // number of free lambdas
  int *par_map; // parameter groups
  char *pattern;
  double *times;
  double *params;
}psmc_par;


typedef struct {
  char *chooseChr;
  int start;
  int stop;
  size_t nSites;
  int maxIter;
  double tole;
  perpsmc * perc;
  char *fname;
  int onlyOnce;
  long seed;//<-seed=-1 old version;seed=0 means time; othervise it will be used as seed
  int block;
  psmc_par *par;
}args;
args * getArgs(int argc,char **argv);
void destroy_args(args *p);
int main_psmc(int argc,char **argv);
