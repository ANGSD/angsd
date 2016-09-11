#pragma once
#include "psmcreader.h"

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
}args;
args * getArgs(int argc,char **argv);
void destroy_args(args *p);
int main_psmc(int argc,char **argv);
