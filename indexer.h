#include "argStruct.h"
int parse_region(char *extra,const aHead *hd,int &ref,int &start,int &stop,const aMap *revMap);
void getOffsets(char *bamname,const aHead *hd,iter_t &ITER,int ref,int start,int stop);
void dalloc(tindex idx);
