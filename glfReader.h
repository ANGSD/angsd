#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "analysisFunction.h"

//binary glf: per site nInd*nGL doubles, no positions (sites are numbered from 0 on the first reference)
class glfReader{
private:
  int nInd;
  gzFile gz;
  int nGL;//10 for -glf, 3 for -glf3
  int pos;//position of the next site
public:
  funkyPars *fetch(int chunkSize);
  glfReader(int &nInd_a,gzFile gz_a,int nGL_a);
  ~glfReader();
};
