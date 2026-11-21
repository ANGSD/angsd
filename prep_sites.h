#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <htslib/bgzf.h>
#include "argStruct.h"//ltstr
#include "analysisFunction.h"

//one chromosome in the index (.idx) of a -sites file
typedef struct{
  int64_t offs;//offset of the chromosome in the .bin file
  size_t len;//last position (0-based) +1, the number of entries stored for the chromosome
}filt_chr;

typedef std::map<char*,filt_chr,ltstr> pMap;

typedef struct{
  BGZF *bg;
  FILE *fp;//the .idx file, closed after it has been read in filt_read
  pMap offs;
  char *keeps;
  char *major;
  char *minor;
  float *af;
  int *an;
  int *ac;
  int hasExtra;// 0=only sites,1=majorminor,2=also an ac af
  size_t curLen;//<-length of the char arrays;NOT the last position filtersite
  char *curNam;//not a copy, points to the name given to filt_readSites
}filt;


filt *filt_read(const char *fname);
void filt_dalloc(filt *f);
void filt_readSites(filt*fl,char *chr,size_t hint);
