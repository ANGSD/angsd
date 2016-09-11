#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>
#include "keep.hpp"
#include "header.h"

typedef struct{
  size_t nSites;
  size_t nChr;
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

size_t iter_read(persaf *saf, void *data, size_t length,int *pos);
template <typename T>
persaf* persaf_init(char *fname);
void writesaf_header(FILE *fp,persaf *pp);
void persaf_destroy(persaf *pp);
myMap::iterator iter_init(persaf *,char *,int,int);
