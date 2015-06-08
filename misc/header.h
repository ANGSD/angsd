#pragma once
#include <cstdlib>
#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <cstring>


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};


void normalize(double *tmp,int len);
size_t fsize(const char* fname);
BGZF *openFileBG(const char* a,const char* b);
FILE *openFile(const char* a,const char* b);
int fexists(const char* str);
