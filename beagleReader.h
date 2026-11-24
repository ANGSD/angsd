#pragma once
#include <zlib.h>
#include "shared.h"
class beagle_reader{
private:
  const aMap *revMap;
  int l;//capacity of buffer
  char *buffer;//current line
  int bufSize;//capacity of buf
  char *buf;//copy of buffer that is tokenized
  gzFile gz;
  int intName;//marker is chr_pos
  int nInd;
  int lineNo;
  int positions;//position counter when not intName
  int changed;//buffer holds the first line of the next chromosome
  int isGL;//the file has genotype likelihoods (e.g. from -doGlf 2), stored as log likelihoods in likes
  int parseLine(int *refId,int *posi,char *major,char *minor,double *post);
public:
  beagle_reader(gzFile gz_a,const aMap *revMap_a,int,int&,int isGL_a);
  ~beagle_reader();

  funkyPars *fetch(int chunkSize);
};
