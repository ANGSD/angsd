#pragma once
#include "zlib.h"
#include "argStruct.h"
#include "analysisFunction.h"

class glfReader{
private:
  int nInd;
  gzFile gz;
  int nGL;
public:
  funkyPars *fetch(int chunkSize);
  glfReader(int &nInd_a,gzFile gz_a,int nGL_a);
  ~glfReader();
};
