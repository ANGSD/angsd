//based on SAMtools code

#ifndef BAMS_H
#define BAMS_H
 
#include <stdint.h>
#include <map>
#include <vector>
#include "bgzf.h"
#include "kstring.h"


#define MAX_SEQ_LEN 200 //this is used for getting some secure value for when a site is securely covered for sample
//#define NTHREADS 10
#define UPPILE_LEN 8

// macros taken from samtools
#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)


#define RLEN 512
/*
  CIGAR operations.
 */
/*! @abstract CIGAR: M = match or mismatch*/
#define BAM_CMATCH      0
/*! @abstract CIGAR: I = insertion to the reference */
#define BAM_CINS        1
/*! @abstract CIGAR: D = deletion from the reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: S = clip on the read with clipped sequence
  present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: P = padding */
#define BAM_CPAD        6
/*! @abstract CIGAR: equals = match */
#define BAM_CEQUAL        7
/*! @abstract CIGAR: X = mismatch */
#define BAM_CDIFF        8


/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024



#define bam1_strand(b) (((b)->flag&BAM_FREVERSE) != 0)


/*
  structure to represent the header of a read, names taken from bam spec
 */
typedef struct {
  int l_text;
  char *text;//len is l_text
  int n_ref;//number of chromosomes
  int *l_name;//length of chromoname name
  char **name;//name[n_ref][l_name]
  int *l_ref;//length of the chromosome
}aHead;

/*
  strute to represent a single read, names taken from bam spec
 */
typedef struct {
  int refID;//4
  int pos;// 4
  unsigned int bin_mq_nl;//4
  unsigned int flag_nc;//4
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






/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)





#define getSeq(b) ((b)->vDat + (b)->nCig*sizeof(uint32_t) + (b)->l_qname)
#define getCig(b) ((uint32_t*)((b)->vDat + (b)->l_qname))
#define getQuals(b) ((b)->vDat + (b)->nCig*sizeof(uint32_t) + (b)->l_qname + (((b)->l_seq + 1)>>1))
#define getAuxStart(b) ((b)->vDat + (b)->nCig*sizeof(uint32_t) + (b)->l_qname + (b)->l_seq + ((b)->l_seq + 1)/2)

//utility forprinting out problematic stuff
#define printErr(b) (({fprintf(stderr,"\t-> Problems: in func:%s at line:%d in file:%s \n",__FUNCTION__,__LINE__,__FILE__),  exit(0);}))



uint32_t bam_calend(const aRead& rd, const uint32_t *cigar);
int is_overlap(uint32_t beg, uint32_t end, const aRead &b);
int bam_validate1(const aHead *header, const aRead b);


int bam_read1(BGZF *fp,aRead & b);
aHead *getHd(BGZF *gz);
int getAlign(BGZF *gz,int block_size,aRead &st);
void dalloc(const aHead *hd);
BGZF *openBAM(const char *fname);




typedef struct{
  int nReads;

  int *first;//first postition of read, relative to reference
  int *last;//last position of read,relative to reference
  int l;
  int m;
  aRead *reads;
  
  int readIDstop;//an intervalue represeing how many reads we want to process
  int lowestStart;//the lowestStartacroos an array of sglPool
  
  //int regionDone;//can be either chromosomeDone or regionDone (whens supplying regions)
  aRead bufferedRead;//this is used for buffering a read, in the case of chromosome change
  //int isEOF;
}sglPool;




sglPool getPool(BGZF *fp,int nReads,aRead &bufRead,int &keepGoing);
void printSglPool(const sglPool& sgl,FILE *fp,aHead *hd);
void dalloc (sglPool &ret);

typedef struct{
  uint64_t chunk_beg;
  uint64_t chunk_end;
}pair64_t;


typedef struct __tindex *tindex;


typedef struct{
  int from_first;
  int tid,beg,end,n_off,i,finished;
  uint64_t curr_off;
  pair64_t *off;//offsets for alignements to loop through
  tindex dasIndex;
}iter_t;

void printIter(const iter_t& it,FILE *fp);

int bam_iter_read(BGZF *fp, iter_t *iter, aRead &b);

aHead *getHd_andClose(const char *fname);
void printHd(const aHead *hd,FILE *fp);
int restuff(aRead &b);
int compHeader(aHead *hd1,aHead *hd2);
#endif
