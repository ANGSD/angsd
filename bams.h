//based on SAMtools code

#ifndef BAMS_H
#define BAMS_H
 
#include <stdint.h>
#include <map>
#include <vector>
#include <htslib/sam.h>
#include <cstdio>
#include <htslib/hts.h>

#define MAX_SEQ_LEN 200 //this is used for getting some secure value for when a site is securely covered for sample
//#define NTHREADS 10
#define UPPILE_LEN 8

// macros taken from samtools
#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam1_strand(b) (((b)->flag&BAM_FREVERSE) != 0)


/*
  strute to represent a single read, names taken from bam spec
 */
typedef struct {
  int refID;//4
  int pos;// 4
  //  unsigned int flag_nc;//4
  int l_seq;//4 <-  length of sequence
  int next_refID;//4
  int next_pos;//4
  int tlen;//4 sum to here=32bytes
  uint8_t *vDat;//this size varies..
  uint32_t block_size;//32 has been subtracted from this value
 
  //below are the bitparsed attributes (to make life easy)
  uint32_t nBin:16, mapQ:8, l_qname:8;
  uint32_t flag:16, nCig:16;
}aRead;


#define getSeq(b) ((b)->vDat + (b)->nCig*sizeof(uint32_t) + (b)->l_qname)
#define getCig(b) ((uint32_t*)((b)->vDat + (b)->l_qname))
#define getQuals(b) ((b)->vDat + (b)->nCig*sizeof(uint32_t) + (b)->l_qname + (((b)->l_seq + 1)>>1))

//uint32_t bam_calend(const aRead& rd, const uint32_t *cigar);
int is_overlap(uint32_t beg, uint32_t end, const aRead &b);
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
