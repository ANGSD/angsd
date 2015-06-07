#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <vector>
#include <htslib/bgzf.h>
#include "header.h"

typedef struct{
  size_t nSites;
  int64_t off;
}dat;

typedef std::map<char*,dat,ltstr> myFstMap;

typedef struct{
  size_t nSites;
  myFstMap mm;
  BGZF *fp;
  std::vector<char*> names;
  int version;
}perfst;

perfst* perfst_init(char *fname);
void perfst_destroy(perfst *pp);
void writefst_header(FILE *fp,perfst *pp);
