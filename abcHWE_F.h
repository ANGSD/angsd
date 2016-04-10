
#include "abc.h"


typedef struct {
  double *freq;
  double *like0;
  double *likeF;
  double *F;
}funkyHWE_F;


class abcHWE_F:public abc{
private:
  Chisqdist *chisq;
  int doHWE_F;
  BGZF* outfileZ;
  kstring_t bufstr;
  void estHWE(double *x,double *loglike,int nInd);
  double HWE_like(double *x,double *loglike,int nInd);
  void HWE_EM(double *x,double *loglike,int nInd);
  double HWE_pval_F;
  double LRT_thres;
  int testMe;
  double tolStop;
  double differ(double *x,double *y);
public:
  abcHWE_F(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcHWE_F();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};

