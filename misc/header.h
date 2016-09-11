#pragma once
#include "../prep_sites.h"
#include <cstdlib>
#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <cstring>

void normalize(double *tmp,size_t len);
size_t fsize(const char* fname);
BGZF *openFileBG(const char* a,const char* b);
FILE *openFile(const char* a,const char* b);
int fexists(const char* str);

typedef struct{
  size_t nSites;
  int64_t pos;
  int64_t saf;
}datum;

typedef std::map<char*,datum,ltstr> myMap;


void my_bgzf_write(BGZF *fp, const void *data, size_t length);
void my_bgzf_seek(BGZF *fp, int64_t pos, int whence);
void my_bgzf_read(BGZF *fp, void *data, size_t length);
size_t getTotalSystemMemory();
char * get_region(char *extra,int &start,int &stop) ;
