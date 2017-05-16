#pragma once
#include "abc.h"

class abcScounts:public abc{
private:
  int doScounts;
  BGZF *outfile;
char *vcfname;
aMap am;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcScounts(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcScounts();
};


