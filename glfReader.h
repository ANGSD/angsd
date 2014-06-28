#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "analysisFunction.h"

class glfReader{
private:
  int nInd;
  int nInd2;//<-used when slicing out samples, totally addhoc, should be generealized
  int from;
  int to;
  gzFile gz;
  int isSim;
public:
  funkyPars *fetch(int chunkSize);
  glfReader(int &nInd_a,int nInd2_a,int from_a, int to_a,gzFile gz_a,int isSim_a);
  ~glfReader();
};
