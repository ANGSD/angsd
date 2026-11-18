#pragma once
#include "abc.h"
#include "prep_sites.h"

class abcFilter : public abc{
  int strict;
  int doMajorMinor;
  char *fname;
  int minInd;
  int setMinIndDepth;
  int capDepth;
public:
  filt *fl;//<-DRAGON this is being used in abcFreq.cpp for pop freq posts
  void readSites(int refId);
  abcFilter(argStruct *arguments);
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  ~abcFilter();
};
