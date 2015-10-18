#pragma once
#include "abc.h"

class abcTemplate:public abc{
private:
  int doTemplate;
  char *refName;
  FILE *outfile;
public:
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcTemplate(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcTemplate();
};


