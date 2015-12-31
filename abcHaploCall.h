#include "abc.h"
typedef struct{
  int **dat;
  int *major;
}haploCalls;

class abcHaploCall:public abc{
private:
  char *refName;
  char *ancName;
  int doHaploCall;
  int nInd;
  int sample;
  int doCount;

  int maxMis;
  int minMinor;

  FILE *outfile;
  FILE *outfile2;
  BGZF* outfileZ;

  size_t **alleleCounts; //[ind][125]; 
  size_t **alleleCountsChr; //[ind][125]; 

  void model1(funkyPars *pars);
  void printHaplo(funkyPars *pars);
  void getHaplo(funkyPars *pars);
  char *tsk_outname;

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
