#include <cstdio>
#include <cstdlib>
#include "razf.h"
#include "analysisFunction.h"

typedef struct{
  int64_t offs;
  int len;
}asdf_dats;



typedef struct{
  BGZF *bg;
  FILE *fp;
  std::map<int,asdf_dats> offs;
  char *keeps;
  char *major;
  char *minor;
  int hasMajMin;
  int curLen;//<-length of the char arrays;
}filt;


filt *filt_init(const char *fname,const aMap* revMap,const aHead *hd,int isBed);
