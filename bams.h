//based on SAMtools code

#ifndef BAMS_H
#define BAMS_H
 
#include <stdint.h>
#include <cstdio>
#include <htslib/hts.h>
#include <htslib/sam.h>
#define MAX_SEQ_LEN 200 //this is used for getting some secure value for when a site is securely covered for sample
//#define NTHREADS 10
#define UPPILE_LEN 8



//uint32_t bam_calend(const aRead& rd, const uint32_t *cigar);
htsFile *openBAM(const char *fname);

typedef struct{
  int nReads;
  int *first;//first postition of read, relative to reference
  int *last;//last position of read,relative to reference
  int l;
  int m;
  bam1_t **reads;
  
  int readIDstop;//an intervalue represeing how many reads we want to process
  int lowestStart;//the lowestStartacroos an array of sglPool
  
  //int regionDone;//can be either chromosomeDone or regionDone (whens supplying regions)
  bam1_t *bufferedRead;//this is used for buffering a read, in the case of chromosome change
  //int isEOF;
}sglPoolb;




typedef struct{
  hts_idx_t *hts_idx;
  hts_itr_t *hts_itr;
}iter_t;

int bam_iter_read2(htsFile *fp, iter_t *iter,bam1_t *b,bam_hdr_t *hdr);
int restuff(bam1_t *b);
#endif
