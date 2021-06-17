#pragma once
#include "abc.h"

class abcTemplate2:public abc{
private:
  int doTemplate2;
  char *refName;
  FILE *outfile;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcTemplate2(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcTemplate2();
};


