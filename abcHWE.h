
#include "abc.h"


typedef struct {
  double *freq;
  double *freqHWE;
  double *like0;
  double *likeF;
  double *pval;
  
  double *F;
}funkyHWE;


class abcHWE:public abc{
private:
  Chisqdist *chisq;
  //  int doSnpStat;
  BGZF* outfileZ;
  kstring_t bufstr;
  void estHWE(double *x,double *loglike,int nInd);
  double HWE_like(double *x,double *loglike,int nInd);
  void HWE_EM(double *x,double *loglike,int nInd);
  //  double HWE_pval;
  int doHWE;
  double minHWEpval;
  double maxHWEpval;
  int testMe;
  double tolStop;
  double differ(double *x,double *y);

public:
  abcHWE(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcHWE();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};

