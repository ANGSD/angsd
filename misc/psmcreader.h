#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "header.h"

typedef struct{
  size_t nSites;
  myMap mm;
  BGZF *pos;
  BGZF *saf;
  int version;
  int *ppos;
  char *fname;
  size_t first;
  size_t last;
  double *gls;
}perpsmc;

perpsmc* perpsmc_init(char *fname);
void writesaf_header(FILE *fp,perpsmc *pp);
void perpsmc_destroy(perpsmc *pp);
myMap::iterator iter_init(perpsmc *,char *,int,int);
