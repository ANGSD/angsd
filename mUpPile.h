#ifndef muppile_
#define muppile_

#include "argStruct.h"
#include "bambi_interface.h"

#include <htslib/kstring.h>
#include <htslib/hts.h>
#include "bams.h"



/*
  node for uppile for a single individual for a single site
  this is used for textoutput emulating samtools mpileup
 */

typedef struct{
  int len;//length of seq,qs,pos
  int maxLen;//maxlength of seq,qs.pos, maybe do realloc
  kstring_t seq;
  kstring_t qs;
  kstring_t pos;
  int depth;
  int refPos;
}node;

typedef struct{
  char *refName;
  int regStart;
  int regStop;
  int nSites;
  int nSamples;
  node **nd;//nd[site][ind]
  int *refPos;//length is nSites
  int first;
  int last;
  int start;
  int length;
  int refId;
}chunky;

typedef struct{
  htsFile *fp;
  char *fn;
  bam_hdr_t *hd;
  int isEOF;
  int regionDone;
  iter_t it;
  regs regions;
}bufReader;



int uppile(int show,int nThreads,bufReader *rd,int NLINES,int nFiles,std::vector<regs> regions);
void dalloc_node(node &n);
tNode initNodeT(int l);


#endif
