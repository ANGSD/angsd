#include "abc.h"
typedef struct{
  int **dat;
  int *major;
}IBSstruct;

class abcIBS:public abc{
private:

  //non optional arguments
  int doIBS;
  int majorminor;
  int doCount;
  int output01;
  int intToMajorMinor[5];
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
  abcIBS(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcIBS();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};
