#pragma once
#include "abc.h"

class abcAncestry:public abc{
private:
  int doAncestry;
  FILE *outfile;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcAncestry(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcAncestry();
};


