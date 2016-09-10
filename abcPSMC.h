#include "abc.h"

typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  double *pLikes;
}psmcRes;


class abcPSMC:public abc{
private:
  int noTrans;
  int dopsmc;
  void writeAll();
  BGZF *outfileSAF;
  FILE *outfileSAFIDX;
  BGZF *outfileSAFPOS;
  char *tmpChr;
  int nnnSites;
  int64_t offs[2];
  void algoJoint(double **liks,int nsites,int *keepSites,psmcRes *r,int noTrans);
 public:
  abcPSMC(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcPSMC();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void changeChr(int refId);
};
