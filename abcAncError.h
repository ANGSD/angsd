#include "abc.h"

class abcAncError:public abc{
private:
  char *refName;
  char *ancName;
  int doAncError;
  int nInd;
  int sample;
  FILE *outfile;
  FILE *outfile2;
  int currentChr;

  size_t **alleleCounts; //[ind][125]; 
  size_t **alleleCountsChr; //[ind][125]; 

  void model1(funkyPars *pars);
  char *tsk_outname;
public:
  abcAncError(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcAncError();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};
