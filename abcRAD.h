#pragma once
#include "abc.h"

class abcRAD:public abc{
private:
  int doRad;
  char *refName;
  FILE *outfile;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcRAD(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcRAD();
};


