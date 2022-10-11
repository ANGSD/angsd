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

//three externs below are from
extern int uniqueOnly;
extern int only_proper_pairs;
extern int remove_bads;
extern int minMapQ;
extern float downSample;
extern int adjustMapQ;
extern int baq;
extern void *rghash;
int getNumBest(bam1_t *b) {
  uint8_t *s= bam_get_aux(b);
  uint8_t *sStop = s+bam_get_l_aux(b);
  int retVal =1;
  
  while (s < sStop) {
    uint8_t type;
    int doRet =0;
    if(s[0]=='X'&&s[1]=='0')
      doRet = 1;

    s += 2; type = *s; ++s;
    if (type == 'A') {if(doRet) return (int) *s ;else  ++s; }
    else if (type == 'C') {if(doRet) return *s;else ++s; }
    else if (type == 'c') {if(doRet) return *(int8_t*)s; else  ++s; }
    else if (type == 'S') {if(doRet) return *(uint16_t*)s; else s += 2; }
    else if (type == 's') {if(doRet) return *(int16_t*)s; else  s += 2; }
    else if (type == 'I') {if(doRet) return *(uint32_t*)s;else  s += 4; }
    else if (type == 'i') {if(doRet) return *(int32_t*)s; else  s += 4; }
    else if (type == 'f') { s += 4; continue; }
    else if (type == 'd') {s += 8;continue; }
    else if (type == 'Z' || type == 'H') return retVal;
  }
  fprintf(stderr,"\t-> Need XO tag to calculate the number of best hits\n");
  exit(0);
  return retVal; //<- will never happen
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
   if(q<0)
     return 0;
   if(q<b->core.qual)
     b->core.qual = q;
 }
 if(b->core.qual<minMapQ)
   return 0;
  
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
  
  if(r!=-1) {

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
    if(hasdata==0)
      	goto bam_iter_reread;
    if(b->core.n_cigar==1){

      uint32_t *cigs = bam_get_cigar(b);
      int opCode = cigs[0]&BAM_CIGAR_MASK; //what to do
      if(opCode==BAM_CINS)
	goto bam_iter_reread;
    }
    
    //extract reads from that has correct RG
    if(rghash){
      uint8_t *rg = bam_aux_get(b, "RG");
      int keep = (rg && khash_str2int_get(rghash, (const char*)(rg+1), NULL)==0);
      if(0&&keep)
	fprintf(stderr,"rg:%s keep:%d\n",rg,keep);
      if (keep==0) goto bam_iter_reread;

    }
    //    uint8_t *rg = bam_aux_get(b, "RG"); fprintf(stderr,"%s\n",rg);
    if((downSample>0 )&& (drand48()>downSample))
    goto bam_iter_reread;
  
  extern abcGetFasta *gf;
  if(b->core.flag&4||b->core.n_cigar==0||b->core.l_qseq==0||b->core.flag&1024)//discard unmapped, nocigar,no seq and dups
      goto bam_iter_reread;
    
    
    if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
      //fprintf(stderr,"r:%d\n",r);
      if(cigstat)
	aio::doAssert(cigstat_calc(b)==0,1,AT,"");
      return r;
    }
    
    //check if read was bad
    if(remove_bads==1 && b->core.flag &512)
	goto bam_iter_reread;
    
    //check if orphan read, 
    if(only_proper_pairs&&(b->core.flag%2) ){//only check if pairend
      if(!(b->core.flag&BAM_FPROPER_PAIR))
	goto bam_iter_reread;
    }if(uniqueOnly&& getNumBest(b)!=1)
       goto bam_iter_reread;
    
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(baq){
	sam_prob_realn(b,gf->ref->seqs,gf->ref->chrLen,baq);
      }
    }
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(adjustMapQ!=0){
	aio::doAssert(gf->ref!=NULL,1,AT,"");
	int q = sam_cap_mapq(b,gf->ref->seqs,gf->ref->chrLen,adjustMapQ);
	if(q<0)
	  goto bam_iter_reread;
	if(q<b->core.qual)
	  b->core.qual = q;
      }
    }

    if(b->core.qual<minMapQ)
      goto bam_iter_reread;
  }
  // fprintf(stderr,"r:%d\n",r);
  if(cigstat)
    aio::doAssert(cigstat_calc(b)==0,1,AT,"");
  return r; 
}

