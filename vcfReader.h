#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "analysisFunction.h"

class vcfReader{
private:
  int nInd;
  gzFile gz;
  int len;
  char *buf,*saveptr;
  const aMap *revMap;
  int parseline(double **lk,double **gp,char &major,char &minor);
  int curChr;
public:
  funkyPars *fetch(int chunkSize);
  vcfReader(int &nInd_a,gzFile gz_a, const aMap *revMap);
  ~vcfReader();
};
