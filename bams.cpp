//based on SAMtools code
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctype.h>
#include <cmath>
#include <htslib/hts.h>
#include "bams.h"
#include "bam_md.h"
#include "abcGetFasta.h"

//three externs below are from
extern int uniqueOnly;
extern int only_proper_pairs;
extern int remove_bads;
extern int minMapQ;
extern int adjustMapQ;
extern int baq;

int getNumBest2(bam1_t *b) {
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
 assert(gf->ref->curChr==b->core.tid);

 if(baq){
   assert(gf->ref!=NULL);
   assert(gf->ref->curChr==b->core.tid);
   bam_prob_realn_core(b,gf->ref->seqs,baq);
 }

 if(adjustMapQ!=0){
   assert(gf->ref!=NULL);
   int q = bam_cap_mapQ(b,gf->ref->seqs,adjustMapQ);
   if(q<0)
     return 0;
   if(q<b->core.qual)
     b->core.qual = q;
 }
 if(b->core.qual<minMapQ)
   return 0;
  
  return 1;

}

htsFile *openBAM(const char *fname){
  htsFile *fp =NULL;
  if((fp=sam_open(fname,"r"))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }
  const char *str = strrchr(fname,'.');
  if(str&&strcasecmp(str,".bam")!=0&&str&&strcasecmp(str,".cram")!=0){
    fprintf(stderr,"\t-> file:\"%s\" should be suffixed with \".bam\" or \".cram\"\n",fname);
    exit(0);
  }
  return fp;
}

int bam_iter_read2(htsFile *fp, iter_t *iter,bam1_t *b,bam_hdr_t *hdr) {
  int r;
 bam_iter_reread:
  
  if(iter->hts_itr==NULL)
    r= sam_read1(fp,hdr,b);
  else
    r = sam_itr_next(fp, iter->hts_itr, b);
  if(r!=-1) {
    extern abcGetFasta *gf;
    if(b->core.flag&4||b->core.n_cigar==0)
      goto bam_iter_reread;
    
    
    if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
      return 1;
    }
    
    //check if read was bad
    if(remove_bads ){
      if(b->core.flag>=256)
	goto bam_iter_reread;
    }
    
    
    //check if orphan read, 
    if(only_proper_pairs&&(b->core.flag%2) ){//only check if pairend
      if(!(b->core.flag&BAM_FPROPER_PAIR))
	goto bam_iter_reread;
    }if(uniqueOnly&& getNumBest2(b)!=1)
       goto bam_iter_reread;
    
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(baq){
	bam_prob_realn_core(b,gf->ref->seqs,baq);
      }
    }
    
    if(gf->ref!=NULL&&gf->ref->curChr==b->core.tid){
      if(adjustMapQ!=0){
	assert(gf->ref!=NULL);
	int q = bam_cap_mapQ(b,gf->ref->seqs,adjustMapQ);
	if(q<0)
	  goto bam_iter_reread;
	if(q<b->core.qual)
	  b->core.qual = q;
      }
    }
    if(b->core.qual<minMapQ)
      goto bam_iter_reread;
  }
  return r; 
}

