#pragma once
#include "abc.h"

class abcWritePlink:public abc{
private:
  int doPlink;
  FILE *fp;//tped
  kstring_t kstr;
public:
  abcWritePlink(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWritePlink();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
