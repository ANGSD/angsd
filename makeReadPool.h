#pragma once

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
  int readIDstop;//an intervalue represeing how many reads we want to process
  int lowestStart;//the lowestStartacroos an array of sglPool
  bam1_t *bufferedRead;//this is used for buffering a read, in the case of chromosome change
}sglPoolb;


typedef struct{
  htsFile *fp;
  char *fn;
  bam_hdr_t *hdr;
  int isEOF;
  int regionDone;
  iter_t it;
  regs regions;
}bufReader;

sglPoolb makePoolb(int l);
void dalloc (sglPoolb *ret);

int collect_reads(bufReader *rd,int nFiles,int &notDone,sglPoolb *ret,int &readNlines,int ref,int &pickStop);
