#include <cstdio>
#include <cstdlib>
#include "analysisFunction.h"

typedef struct{
  int64_t offs;
  size_t len;
}asdf_dats;

typedef std::map<char*,asdf_dats,ltstr> pMap;

typedef struct{
  BGZF *bg;
  FILE *fp;
  pMap offs;
  char *keeps;
  char *major;
  char *minor;
  int hasMajMin;
  size_t curLen;//<-length of the char arrays;
  char *curNam;
}filt;


filt *filt_read(const char *fname);
void dalloc(filt *f);
void filt_readSites(filt*fl,char *chr,size_t hint) ;
