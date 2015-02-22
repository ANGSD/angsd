
#include "abc.h"
#include <htslib/bgzf.h>

class abcSmartCounts:public abc{
private:
  int doSmartCounts;
  FILE *fidx;
  BGZF *fbin;
  int curChr;
  int len;
public:
  abcSmartCounts(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcSmartCounts();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void changeChr(int refId);
};

