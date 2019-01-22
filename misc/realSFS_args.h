#pragma once
#include "../prep_sites.h"
#include "safreader.h"
#include <vector>
typedef struct {
  char *chooseChr;
  int start;
  int stop;
  size_t nSites;
  int maxIter;
  double tole;
  int nThreads;
  std::vector<char *> sfsfname;
  std::vector<persaf *> saf;
  //  char *fname;
  int onlyOnce;
  int emAccl;
  char *outname;
  int win;
  int step;
  int type;
  int oldout;
  long seed;//<-seed=-1 old version;seed=0 means time; othervise it will be used as seed
  int bootstrap;
  int resample_chr;//nspope; 1=bootstrap sites by resampling chromosomes, 0=regular bootstrap
  filt *fl;
  int whichFst;
  int fold;
  char *ref;
  char *anc;
}args;
args * getArgs(int argc,char **argv);
void destroy_args(args *p);
