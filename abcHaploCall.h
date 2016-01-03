#include "abc.h"
typedef struct{
  int **dat;
  int *major;
}haploCalls;

class abcHaploCall:public abc{
private:

  //non optional arguments
  int doHaploCall;
  int doCount;

  //optional arguments
  int maxMis;
  int minMinor;

  //out file
  BGZF* outfileZ;

  //functions
  void printHaplo(funkyPars *pars);
  void getHaplo(funkyPars *pars);

  //print buffer
  kstring_t bufstr;

public:
  abcHaploCall(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcHaploCall();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};
