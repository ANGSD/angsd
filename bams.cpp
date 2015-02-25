//based on SAMtools code
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctype.h>
#include <cmath>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <cram/cram.h>

#include <htslib/hts.h>
#include "bams.h"
#include "baq_adjustMapQ.h"
#include "abcGetFasta.h"

//three externs below are from
extern int uniqueOnly;
extern int only_proper_pairs;
extern int remove_bads;
extern int minMapQ;
extern int adjustMapQ;
extern int baq;

uint32_t bam_calend_old(const aRead& rd, const uint32_t *cigar){
  uint32_t end = rd.pos;
  for (int k = 0; k < rd.nCig; ++k) {
    int op = cigar[k] & BAM_CIGAR_MASK;
    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
      end += cigar[k] >> BAM_CIGAR_SHIFT;
  }
  return end;
}
uint32_t bam_calend(const aRead& rd, const uint32_t *cigar){
  uint32_t end = rd.pos;
  for (int k = 0; k < rd.nCig; ++k) {
    int op = cigar[k] & BAM_CIGAR_MASK;
    if (op==BAM_CEQUAL||op==BAM_CDIFF||op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
      end += cigar[k] >> BAM_CIGAR_SHIFT;
  }
  //  fprintf(stderr,"end=%d\n",end);
  return end;
}



int is_overlap(uint32_t beg, uint32_t end, const aRead &b) {
  uint32_t rbeg = b.pos;
  uint32_t rend = b.nCig ? bam_calend(b, getCig(&b)) : b.pos + 1;
  int ret = (rend > beg && rbeg < end);
  return ret;
}



int bam_validate1(const bam_hdr_t *header, const aRead b)
{
  
  if (b.refID < -1 || b.next_refID < -1){
    fprintf(stderr,"error first\n");
    return 0;
  }
  if (header && (b.refID >= header->n_targets || b.next_refID >= header->n_targets)){
    fprintf(stderr,"error second\n");
    return 0;
  }

  if (b.block_size < b.l_qname){
      fprintf(stderr,"error third\n");
      return 0;

    }
    char *s =(char *) memchr(b.vDat, '\0', b.l_qname);
 
    //    fprintf(stderr,"%d\n",(int )s-(b.vDat+1));
    if (s[0] != b.vDat[b.l_qname-1]) {
      fprintf(stderr,"error forth\n");
    return 0;
    }
    return 1;
}

/*
  given the current offset read an alignment and put the data in st
 */
int getAlign(BGZF *gz,int block_size,aRead &st) {

  uint32_t tmp[8];
  bgzf_read(gz,tmp,32);
  //copy vals to st
  st.refID = tmp[0];
  st.pos = tmp[1];
  unsigned int bin_mq_nl = tmp[2];
  unsigned int flag_nc = tmp[3];
  st.l_seq = tmp[4];
  st.next_refID = tmp[5];
  st.next_pos = tmp[6];
  st.tlen = tmp[7];
  
   //begin vars
  // fprintf(stderr,"reading\t%d\n",block_size-32);
  if(block_size-32>=RLEN){//fix if blocks to big
    delete [] st.vDat;
    int tmpSize = block_size-32;
    kroundup32(tmpSize);
    st.vDat = new uint8_t[tmpSize];
  }
    
  bgzf_read(gz,st.vDat,block_size-32);
  st.block_size=block_size-32;
   
  //plugin bitparsed attributes maybe we don't want to do this...
  st.nBin = bin_mq_nl >>16;
  st.mapQ = bin_mq_nl >>8&0xff;
  st.l_qname = bin_mq_nl &0xff;
  st.flag = flag_nc >>16;
  st.nCig = flag_nc &0xffff;
  return 0;
}

int getNumBest(uint8_t *s, uint8_t *sStop) {
 
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
int restuff(aRead &b){
  //  fprintf(stderr,"ohh yes: b.refID=%d b.pos=%d\n",b.refID,b.pos);
 extern abcGetFasta *gf;
 assert(gf->ref->curChr==b.refID);

 if(baq){
   assert(gf->ref!=NULL);
   assert(gf->ref->curChr==b.refID);
   bam_prob_realn_core(b,gf->ref->seqs,baq);
 }

 if(adjustMapQ!=0){
   assert(gf->ref!=NULL);
   int q = bam_cap_mapQ(b,gf->ref->seqs,adjustMapQ);
   if(q<0)
     return 0;
   if(q<b.mapQ)
     b.mapQ = q;
 }
 if(b.mapQ<minMapQ)
   return 0;
  
  return 1;

}


int bam_t_to_aRead (bam1_t *from,aRead &to){
 to.refID = from->core.tid;
 to.pos = from->core.pos;
 to.nBin = from->core.bin;
 to.mapQ = from->core.qual;
 to.l_qname = from->core.l_qname;
 to.flag = from->core.flag;
 to.nCig = from->core.n_cigar;
 to.l_seq = from->core.l_qseq;
 to.next_refID = from->core.mtid;
 to.next_pos= from->core.mpos;
 to.tlen = from->core.isize;
 // fprintf(stderr,"vDat:Rlen:%d l_data:%d m_data:%d\n",RLEN,from->l_data,from->m_data);
 if(from->l_data+32>RLEN){
   delete [] to.vDat;
   to.vDat = new uint8_t[from->m_data];
 }
 memcpy(to.vDat,from->data,sizeof(uint8_t)*from->l_data);
 // fprintf(stderr,"refId:%d pos:%d flag:%d\n",to.refID,to.pos,to.flag);
 to.block_size= from->l_data;
 return 1;
}

int bam_read2(htsFile *fp,aRead & b,bam_hdr_t *hdr){
  bam1_t *bb = bam_init1();
 bam_read2_reread :
  
  int r = sam_read1(fp,hdr, bb);
  //fprintf(stderr,"r: %d\n",r);
  if(r==-1){
    bam_destroy1(bb);
    return r;
  }
  bam_t_to_aRead(bb,b);
  {
    extern abcGetFasta *gf;
    if(b.flag&4)
      goto bam_read2_reread;
    
    
    if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
      return 1;
    }
    
    //check if read was bad
    if(remove_bads ){
      if(b.flag>=256)
	goto bam_read2_reread;
    }
      
      
    //check if orphan read, 
    if(only_proper_pairs&&(b.flag%2) ){//only check if pairend
      if(!(b.flag&BAM_FPROPER_PAIR))
	goto bam_read2_reread;
    }if(uniqueOnly&& getNumBest(getAuxStart(&b),b.vDat+b.block_size)!=1)
       goto bam_read2_reread;
      
    
    if(gf->ref!=NULL&&gf->ref->curChr==b.refID){
      if(baq){
	bam_prob_realn_core(b,gf->ref->seqs,baq);
      }
    }
    
    if(gf->ref!=NULL&&gf->ref->curChr==b.refID){
      if(adjustMapQ!=0){
	assert(gf->ref!=NULL);
	int q = bam_cap_mapQ(b,gf->ref->seqs,adjustMapQ);
	if(q<0)
	  goto bam_read2_reread;
	if(q<b.mapQ)
	  b.mapQ = q;
      }
    }
    if(b.mapQ<minMapQ)
      goto bam_read2_reread;
  }

  
  bam_destroy1(bb);
  return r;
  
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

int bam_iter_read1(htsFile *fp, iter_t *iter, aRead &b,bam_hdr_t *hdr) {
  int r;
  bam1_t *bb=bam_init1();
  
  if(iter->hts_itr==NULL){
    r= bam_read2(fp, b,hdr);
  }else {
  bam_iter_reread :
    r = sam_itr_next(fp, iter->hts_itr, bb);
    fprintf(stderr,"rit:%d\n");
    if(r!=-1) {
      bam_t_to_aRead(bb,b);
      {
	extern abcGetFasta *gf;
	if(b.flag&4)
	  goto bam_iter_reread;
	
	
	if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
	  return 1;
	}
	
	//check if read was bad
	if(remove_bads ){
	  if(b.flag>=256)
	    goto bam_iter_reread;
	}
	
	
	//check if orphan read, 
	if(only_proper_pairs&&(b.flag%2) ){//only check if pairend
	  if(!(b.flag&BAM_FPROPER_PAIR))
	    goto bam_iter_reread;
	}if(uniqueOnly&& getNumBest(getAuxStart(&b),b.vDat+b.block_size)!=1)
	   goto bam_iter_reread;
	
	
	if(gf->ref!=NULL&&gf->ref->curChr==b.refID){
	  if(baq){
	    bam_prob_realn_core(b,gf->ref->seqs,baq);
	  }
	}
	
	if(gf->ref!=NULL&&gf->ref->curChr==b.refID){
	  if(adjustMapQ!=0){
	    assert(gf->ref!=NULL);
	    int q = bam_cap_mapQ(b,gf->ref->seqs,adjustMapQ);
	    if(q<0)
	      goto bam_iter_reread;
	    if(q<b.mapQ)
	      b.mapQ = q;
	  }
	}
	if(b.mapQ<minMapQ)
	  goto bam_iter_reread;
      }
      
      
    }
    
  }
  bam_destroy1(bb); 
  
  return r; 
}


void dalloc (sglPool &ret){
  for(int i=0;i<ret.l;i++)
    delete [] ret.reads[i].vDat;
  if(0){
    free( ret.reads);
    free(ret.first);
    free(ret.last);
  }else{
    delete [] ret.reads;
    delete [] ret.first;
    delete [] ret.last;


  }
}
