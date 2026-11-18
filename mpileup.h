#pragma once
#include <vector>
#include <zlib.h>
#include "shared.h"

//reads samtools mpileup text: chr pos ref and depth bases quals for each individual
class mpileup{
  int nInd;
  gzFile gz;
  int l;//capacity of buffer
  char *buffer;//current line
  int bufSize;//capacity of buf
  char *buf;//copy of buffer that is split into fields
  std::vector<char *> fields;
  const aMap *revMap;
  int minQ;
  int lineNo;
  int changed;//buffer holds the first line of the next chromosome
  int nIndelWarn;
  int splitLine();
  tNode **parseNd(char ref);
public:
  funkyPars *fetch(int chunkSize);
  mpileup(int nInd_a,gzFile gz,const aMap *revMap,int minQ_a);
  ~mpileup();
};
