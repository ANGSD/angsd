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
  BGZF *bgzf_pos;
  BGZF *bgzf_gls;
  int version;//is version1, now
  int *pos;//contains the positions
  char *fname;//input.saf.idx?
  size_t first;//if we have specified a region, then this is the first index to use
  size_t last;//if we have specified a region, then this is the last index to use
  double *gls;//2*nSites long homogl_1 hetgl_1 homogl_2 hetgl_2
}perpsmc;

perpsmc* perpsmc_init(char *fname);
void writepsmc_header(FILE *fp,perpsmc *pp);
void perpsmc_destroy(perpsmc *pp);
myMap::iterator iter_init(perpsmc *,char *,int,int);
