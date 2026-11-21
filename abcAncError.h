#include "abc.h"

//anc, ref and the observed base, each 0-4 (A,C,G,T,N): anc*25+ref*5+base
#define ANCERR_NCAT 125

class abcAncError:public abc{
private:
  char *refName;
  char *ancName;
  int doAncError;
  int nInd;
  FILE *outfile;
  FILE *outfile2;
  int currentChr;

  size_t **alleleCounts; //[ind][ANCERR_NCAT], whole run
  size_t **alleleCountsChr; //[ind][ANCERR_NCAT], current chromosome

  void model1(funkyPars *pars);
  void count(int ind,funkyPars *pars,int s,char base);
  void printChr();
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
