#pragma once
#include <vector>
#include <zlib.h>
#include "shared.h"

//text glf (as written by -doGlf 4): chr pos and 10 log likelihoods per individual.
//chr and pos can also be one column chr_pos.
class glfReader_text{
private:
  const aMap *revMap;
  int l;//capacity of buffer
  char *buffer;//current line
  int bufSize;//capacity of buf
  char *buf;//copy of buffer that is tokenized
  std::vector<char *> fields;
  gzFile gz;
  int nInd;
  int lineNo;
  int changed;//buffer holds the first line of the next chromosome
  int parseLine(int *refId,int *posi,double *likes);
public:
  glfReader_text(int nInd,gzFile gz_a,const aMap *revMap_a);
  ~glfReader_text();

  funkyPars *fetch(int chunkSize);
};
