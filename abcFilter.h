#include "prep_sites.h"
class abcFilter : public abc{

  int doMajorMinor;
  char *fname;
  FILE *fp;

  char *keepsChr;//<- a char array of length=reflength, used for .keeps
  int minInd;
  //  int nInd;
  int curChr;
  int capDepth;
public:
  filt *fl;//<-DRAGON this is being used in abcFreq.cpp for pop freq posts
  void readSites(int refId);
  //none optional stuff
  FILE *outfile;
  abcFilter(argStruct *arguments);
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  ~abcFilter();
};
