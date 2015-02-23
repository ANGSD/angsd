#include <htslib/hts.h>
#include "argStruct.h"
int parse_region(char *extra,const aHead *hd,int &ref,int &start,int &stop,const aMap *revMap);
void getOffsets(htsFile *fp, char *fn,const aHead *hd,iter_t &ITER,int ref,int start,int stop,bam_hdr_t *hdr);
void dalloc(tindex idx);
