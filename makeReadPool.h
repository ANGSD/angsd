#pragma once

#include <cassert>
#include <cstdlib>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "pop1_read.h"
#include "argStruct.h"

typedef struct{
  int l;//<- number of elements
  int m;//<-maxnumber of elements

  int *first;//first postition of read, relative to reference
  int *last;//last position of read,relative to reference
  bam1_t **reads;
  int readIDstop;
  int lowestStart;//the lowestStartacroos an array of sglPool
  bam1_t *bufferedRead;//this is used for buffering a read, in the case of chromosome change
}readPool;


typedef struct{
  htsFile *fp;
  char *fn;
  bam_hdr_t *hdr;
  int isEOF;
  int regionDone;
  hts_idx_t *idx;
  hts_itr_t *itr;
  regs regions;
}bufReader;

readPool makePoolb(int l);
void dalloc (readPool *ret);

int collect_reads(bufReader *rd,int nFiles,int &notDone,readPool *ret,int &readNlines,int ref,int &pickStop);
int collect_reads2(bufReader *rd,int nFiles,int &notDone,readPool *ret,int &readNlines,int ref,int &pickStop);
