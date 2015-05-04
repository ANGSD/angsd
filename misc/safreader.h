#pragma once
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <map>
#include <htslib/bgzf.h>

#include "keep.hpp"

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef struct{
  size_t nSites;
  int64_t pos;
  int64_t saf;
}datum;

typedef std::map<char*,datum,ltstr> myMap;

typedef struct{
  size_t nSites;
  size_t nChr;
  myMap mm;
  BGZF *pos;
  BGZF *saf;
  size_t fsize;//contains an estimate of the uncompressed fsize
  int version;
  keep<char> *toKeep;
  int at;
  int *ppos;
  int with_iter;
  int start;
  int stop;
}persaf;

size_t iter_read(persaf *saf, void *data, size_t length);
template <typename T>
persaf* readsaf(char *fname);
void writesaf_header(FILE *fp,persaf *pp);
void persaf_destroy(persaf *pp);
myMap::iterator iter_init(persaf *,char *,int,int);
size_t iter_read(persaf *saf, void *data, size_t length);
