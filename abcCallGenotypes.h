#pragma once
#include "abc.h"

/*
  contains -1,0,1,2.
  -1 no call
  0 AA
  1 Aa
  2 aa
 */
typedef struct{
  int **dat;
}genoCalls;


class abcCallGenotypes:public abc{
private:
  int doGeno;
  float postCutoff;
  BGZF* outfileZ;
  double geno_minMM;
  int geno_minDepth;
  int geno_maxDepth;
  int minInd;
  kstring_t bufstr;
public:
  abcCallGenotypes(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcCallGenotypes();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void getGeno(funkyPars *pars);
  void printGeno(funkyPars *pars);
  
};
