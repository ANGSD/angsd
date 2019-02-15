#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "analysisFunction.h"

class vcfReader2{
private:
  int nInd;
  gzFile gz;
  int len;
  char *buf,*saveptr,*original;
  const aMap *revMap;
  int parseline(double **lk,double **gp,char &major,char &minor);
  int curChr;
public:
  funkyPars *fetch(int chunkSize);
  vcfReader2(int &nInd_a,gzFile gz_a, const aMap *revMap);
  ~vcfReader2();
};
