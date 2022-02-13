#pragma once
#include "abc.h"

class abcAMOVA:public abc{
private:
  int doAMOVA;
  char *refName;
  FILE *outfile;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcAMOVA(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcAMOVA();
};


