#include "abc.h"
class abcPSMC:public abc{
private:
  int dopsmc;
public:
  abcPSMC(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcPSMC();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
