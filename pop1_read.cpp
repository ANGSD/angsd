/*
  angsd thorfinn@binf.ku.dk
  27 feb 2015
  
  functions to read one single read using htslib
*/


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctype.h>
#include <cmath>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include "pop1_read.h"
#include "abcGetFasta.h"
#include "cigstat.h"
#include "from_samtools.h"
#include "aio.h"

//from parseArgs_bambi.cpp
extern int uniqueOnly;
extern int only_proper_pairs;
extern int remove_bads;
extern int minMapQ;
extern float downSample;
extern int adjustMapQ;
extern int baq;
extern void *rghash;

//number of reads read and removed by each filter, printed at the end by printReadCounts().
//pop1_read() and restuff() are only called from the main thread, so no lock is needed
static size_t nReadsRead = 0;
static size_t nReadsUsed = 0;
static size_t nReadsNoData = 0;//unmapped, duplicate (flag 1024), no aligned bases or no sequence
static size_t nReadsRG = 0;
static size_t nReadsDownSample = 0;
static size_t nReadsBad = 0;
static size_t nReadsNotProperPair = 0;
static size_t nReadsNotUnique = 0;
static size_t nReadsMapQ = 0;

void printReadCounts(FILE *fp){
  if(nReadsRead==0)
    return;
  const char *names[] = {"unmapped, duplicates or without aligned bases","+RG/+LB","-downSample","-remove_bads","-only_proper_pairs","-uniqueOnly","-minMapQ/-C"};
  size_t counts[] = {nReadsNoData,nReadsRG,nReadsDownSample,nReadsBad,nReadsNotProperPair,nReadsNotUnique,nReadsMapQ};
  //the numbers are right aligned below the total
  int w = snprintf(NULL,0,"%zu",nReadsRead);
  fprintf(fp,"\t-> Total number of reads read: %zu\n",nReadsRead);
  fprintf(fp,"\t-> Number of reads used after the read filters: %zu\n",nReadsUsed);
  for(int i=0;i<7;i++)
    if(counts[i]>0)
      fprintf(fp,"\t\t%*zu %s%s\n",w,counts[i],i==0?"":"removed by ",names[i]);
  if(nReadsNotProperPair>0&&nReadsUsed==0)
    fprintf(fp,"\t-> No reads were used. Paired reads without the proper pair flag (2) are removed by -only_proper_pairs 1, use -only_proper_pairs 0 to keep them\n");
}
//number of best hits from the X0 tag (made by bwa aln), used by -uniqueOnly
//reads without X0 (eg. mate rescued reads) are kept
static int getNumBest(bam1_t *b) {
  static int warned = 0;
  uint8_t *x0 = bam_aux_get(b,"X0");
  if(x0==NULL){
    if(warned==0){
      fprintf(stderr,"\t-> Warning: -uniqueOnly uses the X0 tag (made by bwa aln), reads without it are kept (first one: \'%s\')\n",bam_get_qname(b));
      warned = 1;
    }
    return 1;
  }
  return (int) bam_aux2i(x0);
}

//simple function to perform additional analysis, only used when we have to redo a single read when we change chr

int restuff(bam1_t *b){
  //  fprintf(stderr,"ohh yes: b.refID=%d b.pos=%d\n",b.refID,b.pos);
 extern abcGetFasta *gf;
 aio::doAssert(gf->ref->curChr==b->core.tid,1,AT,"");

 if(baq){
   aio::doAssert(gf->ref!=NULL,1,AT,"");
   aio::doAssert(gf->ref->curChr==b->core.tid,1,AT,"");
   sam_prob_realn(b,gf->ref->seqs,gf->ref->chrLen,baq);
 }

 if(adjustMapQ!=0){
   aio::doAssert(gf->ref!=NULL,1,AT,"");
   int q = sam_cap_mapq(b,gf->ref->seqs,gf->ref->chrLen,adjustMapQ);
   if(q<0){
     nReadsUsed--;nReadsMapQ++;
     return 0;
   }
   if(q<b->core.qual)
     b->core.qual = q;
 }
 if(b->core.qual<minMapQ){
   nReadsUsed--;nReadsMapQ++;
   return 0;
 }
  
  return 1;

}
extern int cigstat;

int pop1_read(htsFile *fp, hts_itr_t *itr,bam1_t *b,bam_hdr_t *hdr) {
  int r;
 bam_iter_reread:
  
  if(itr==NULL)
    r= sam_read1(fp,hdr,b);
  else
    r = sam_itr_next(fp, itr, b);

  //-1 is the end of the file or region, below -1 is an error (eg. a truncated file)
  if(r<-1){
    fprintf(stderr,"\t-> Problem reading: \'%s\' (error %d), the file looks truncated or damaged\n",fp->fn,r);
    exit(1);
  }
  
  if(r!=-1) {
    nReadsRead++;

    //pathologial case with no data in CIGAAR. I dont see how this can happen, but we have observed this is some data.
    //previously it was only one insertion, but in new data it has been a combination of softclip and insertion
    int hasdata = 0;
    for(int i=0;i<b->core.n_cigar;i++){
      uint32_t *cigs = bam_get_cigar(b);
      int opCode = cigs[i]&BAM_CIGAR_MASK; //what to do
      if(opCode==BAM_CMATCH||opCode==BAM_CEQUAL||opCode==BAM_CDIFF){
	hasdata = 1;
	break;
      }
    }
    if(hasdata==0){
      nReadsNoData++;
      goto bam_iter_reread;
    }

    //extract reads from that has correct RG
    if(rghash){
      uint8_t *rg = bam_aux_get(b, "RG");
      int keep = (rg && khash_str2int_get(rghash, (const char*)(rg+1), NULL)==0);
      if (keep==0){
	nReadsRG++;
	goto bam_iter_reread;
      }

    }
    //    uint8_t *rg = bam_aux_get(b, "RG"); fprintf(stderr,"%s\n",rg);
    if((downSample>0 )&& (drand48()>downSample)){
      nReadsDownSample++;
      goto bam_iter_reread;
    }
  
  extern abcGetFasta *gf;
  if(b->core.flag&4||b->core.n_cigar==0||b->core.l_qseq==0||b->core.flag&1024){//discard unmapped, nocigar,no seq and dups
      nReadsNoData++;
      goto bam_iter_reread;
  }
    
    
    if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
      //fprintf(stderr,"r:%d\n",r);
      if(cigstat)
	aio::doAssert(cigstat_calc(b)==0,1,AT,"");
      nReadsUsed++;
      return r;
    }
    
    //check if read was bad
    if(remove_bads==1 && b->core.flag &512){
      nReadsBad++;
      goto bam_iter_reread;
    }
    
    //check if orphan read, 
    if(only_proper_pairs&&(b->core.flag%2) ){//only check if pairend
      if(!(b->core.flag&BAM_FPROPER_PAIR)){
	nReadsNotProperPair++;
	goto bam_iter_reread;
      }
    }
    if(uniqueOnly&& getNumBest(b)!=1){
      nReadsNotUnique++;
      goto bam_iter_reread;
    }
    
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(baq){
	sam_prob_realn(b,gf->ref->seqs,gf->ref->chrLen,baq);
      }
    }
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(adjustMapQ!=0){
	aio::doAssert(gf->ref!=NULL,1,AT,"");
	int q = sam_cap_mapq(b,gf->ref->seqs,gf->ref->chrLen,adjustMapQ);
	if(q<0){
	  nReadsMapQ++;
	  goto bam_iter_reread;
	}
	if(q<b->core.qual)
	  b->core.qual = q;
      }
    }

    if(b->core.qual<minMapQ){
      nReadsMapQ++;
      goto bam_iter_reread;
    }
    nReadsUsed++;
  }
  // fprintf(stderr,"r:%d\n",r);
  if(cigstat&&r>=0)//at the end of the file b is not a new read
    aio::doAssert(cigstat_calc(b)==0,1,AT,"");
  return r; 
}

