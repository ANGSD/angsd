#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "keep.hpp"
#include "header.h"

typedef struct{
  size_t nSites;//sum of all sites
  size_t sumBand;
  size_t nChr;//nr of chromosomes(bins)
  myMap mm;
  BGZF *pos;
  BGZF *saf;
  int version;
  keep<char> *toKeep;
  int at;
  int *ppos;
  int kind;//0 = only llh; 1 = only pos ; 2 = both
  char *fname;
  int dontRead;
}persaf;

template <typename T> size_t iter_read(persaf *saf, T *&data, T *buffer, int *pos);
template <typename T> persaf* persaf_init(char *fname, int verbose);
void writesaf_header(FILE *fp, persaf *pp);
void persaf_destroy(persaf *pp);
myMap::iterator iter_init(persaf*, char*, int, int);
