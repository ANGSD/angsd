//based on SAMtools code
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctype.h>
#include <cmath>

#include "bams.h"
#include "baq_adjustMapQ.h"
#include "bgzf.h"
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



int bam_validate1(const aHead *header, const aRead b)
{
  
  if (b.refID < -1 || b.next_refID < -1){
    fprintf(stderr,"error first\n");
    return 0;
  }
  if (header && (b.refID >= header->n_ref || b.next_refID >= header->n_ref)){
    fprintf(stderr,"error second\n");
    return 0;
  }
  //  int l_qname = b.bin_mq_nl &0xff;

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


 void printIter(const iter_t& it,FILE *fp){
   fprintf(stderr,"--------------\n");
   fprintf(stderr,"from_first=%d\n",it.from_first);
   fprintf(stderr,"tid=%d\tbeg=%d\tend=%d\tn_off=%d\ti=%d\tfinished=%d\n",it.tid,it.beg,it.end,it.n_off,it.i,it.finished);
   fprintf(stderr,"i=%d\tcurr_off:%zu\n",it.i,(size_t)it.curr_off);
   for(int i=0;i<it.n_off;i++)
     fprintf(stderr,"%d: (%lu-%lu)\n",i,(size_t)it.off[i].chunk_beg,(size_t)it.off[i].chunk_end);
   fprintf(stderr,"--------------\n");
 }




aHead *getHd(BGZF *gz){
  aHead *h = new aHead;
  
  const char * magic ="BAM\1"; 
  char head[4];
  if(4!=bgzf_read(gz,head,4)){
    fprintf(stderr,"[%s] Error reading\n",__FUNCTION__);
  }
  //  print(stderr,head,4);
  if(strncmp(head,magic,4)!=0) {
    fprintf(stderr,"[%s]\t-> Problem with magic number in header from BAM file\n",__FUNCTION__);
    exit(0);
  }
  bgzf_read(gz,&h->l_text,sizeof(int));
  h->text = new char[h->l_text];
  bgzf_read(gz,h->text,h->l_text);
  #if 0
  fprintf(stdout,"[%s] hd:%s\n",__FILE__,h->text);
  exit(0);
  #endif
  bgzf_read(gz,&h->n_ref,sizeof(int));
  h->l_name=new int[h->n_ref];
  h->name = new char*[h->n_ref];
  h->l_ref = new int[h->n_ref];
  
  for(int i=0;i<h->n_ref;i++){
    bgzf_read(gz,&h->l_name[i],sizeof(int));
    h->name[i] =(char*) malloc(h->l_name[i]);
    bgzf_read(gz,h->name[i],h->l_name[i]);
    bgzf_read(gz,&h->l_ref[i],sizeof(int));
  } 
 
  return h;
}

void printHd(const aHead *hd,FILE *fp){
  fprintf(fp,"htext=%s\n",hd->text);
  fprintf(fp,"n_ref=%d\n",hd->n_ref);
  for(int i=0;i<hd->n_ref;i++)
    fprintf(fp,"i=%d name=%s length=%d\n",i,hd->name[i],hd->l_ref[i]);

}



aHead *getHd_andClose(const char *fname){
  BGZF *fp = openBAM(fname);
  aHead *hd = getHd(fp);
  bgzf_close(fp);
  return hd;
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
  st.bin_mq_nl = tmp[2];
  st.flag_nc = tmp[3];
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
  st.nBin = st.bin_mq_nl >>16;
  st.mapQ = st.bin_mq_nl >>8&0xff;
  st.l_qname = st.bin_mq_nl &0xff;
  st.flag = st.flag_nc >>16;
  st.nCig = st.flag_nc &0xffff;
  return 0;
}


int bam_read1(BGZF *fp,aRead & b){

  int block_size;
  if(4!=bgzf_read(fp,&block_size,sizeof(int))){
    fprintf(stderr,"[%s] error reading some file stuff in assuming EOF\n",__FUNCTION__);
    //    return 0;
    return -1;
  }
  getAlign(fp,block_size,b);
  return 1;
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
  assert(1==0);
  return retVal;
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


#if 1
//this is the old
//read a single 'read'
int bam_read2(BGZF *fp,aRead & b){
  extern abcGetFasta *gf;
 
 reread: //oh yes baby a goto statement

  int block_size;
  if(4!=bgzf_read(fp,&block_size,sizeof(int))){
    return -2;
  }
  getAlign(fp,block_size,b);
  if(b.flag&4)
    goto reread;


  if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
    return 1;
  }

  //check if read was bad
  if(remove_bads ){
    if(b.flag>=256)
      goto reread;
  }


  //check if orphan read, 
  if(only_proper_pairs&&(b.flag%2) ){//only check if pairend
    if(!(b.flag&BAM_FPROPER_PAIR))
      goto reread;
  }if(uniqueOnly&& getNumBest(getAuxStart(&b),b.vDat+b.block_size)!=1)
    goto reread;
  

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
	goto reread;
      if(q<b.mapQ)
	b.mapQ = q;
    }
  }
  if(b.mapQ<minMapQ)
    goto reread;
  
  return 1;
}
#else
//this is the new, that I haven't fixed yet
int bam_read2(BGZF *fp,aRead & b){
  extern getFasta *gf;
  extern unsigned int includeflags,discardflags;
 reread: //oh yes baby a goto statement

  int block_size;
  if(4!=bgzf_read(fp,&block_size,sizeof(int))){
    return -2;
  }
  getAlign(fp,block_size,b);
  if(b.flag&discardflags||((b.flag&includeflags)!=includeflags))
    goto reread;


  if(uniqueOnly==0&&only_proper_pairs==0 &&remove_bads==0&&minMapQ==0&&adjustMapQ==0&&baq==0){
    return 1;
  }

  //check if read was bad
  if(remove_bads ){
    if(b.flag>=256)
      goto reread;
  }


  //check if orphan read, 
  if(only_proper_pairs&&(b.flag%2) ){//only check if pairend
    if(!(b.flag&BAM_FPROPER_PAIR))
      goto reread;
  }if(uniqueOnly&& getNumBest(getAuxStart(&b),b.vDat+b.block_size)!=1)
    goto reread;
  

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
	goto reread;
      if(q<b.mapQ)
	b.mapQ = q;
    }
  }
  if(b.mapQ<minMapQ)
    goto reread;
  
  return 1;
}
#endif

//DRAGON, make sure that all hd->names are allocated with malloc or new. Dont' remember what is used.
void dalloc(const aHead *hd){
  delete [] hd->text;
  delete [] hd->l_name;
  delete [] hd->l_ref;
  for(int i=0;i<hd->n_ref;i++)
    free(hd->name[i]);
  delete [] hd->name;
  delete hd;
  hd=NULL;
}


BGZF *openBAM(const char *fname){
  BGZF* fp =NULL;
  if((fp=bgzf_open(fname,"r"))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }
  const char *str = strrchr(fname,'.');
  if(str&&strcasecmp(str,".bam")!=0){
    fprintf(stderr,"file:\"%s\" should be suffixed with \".bam\"\n",fname);
    exit(0);
  }
  return fp;
}



int bam_iter_read(BGZF *fp, iter_t *iter, aRead &b) {
  int ret;
  if (iter && iter->finished) return -1;
  if (iter == 0 || iter->from_first) {
    ret = bam_read1(fp, b);
    if (ret < 0 && iter) iter->finished = 1;
    return ret;
  }
  if (iter->off == 0) return -1;
  for (;;) {
    if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].chunk_end) { // then jump to the next chunk
      if (iter->i == iter->n_off - 1) {
	ret = -1; break; 
      } // no more chunks
      if (iter->i >= 0) assert(iter->curr_off == iter->off[iter->i].chunk_end); // otherwise bug
      if (iter->i < 0 || iter->off[iter->i].chunk_end != iter->off[iter->i+1].chunk_beg) { // not adjacent chunks; then seek
	bgzf_seek(fp, iter->off[iter->i+1].chunk_beg, SEEK_SET);
	iter->curr_off = bgzf_tell(fp);
      }
      ++iter->i;
    }
    if ((ret = bam_read1(fp, b)) >= 0) {
      iter->curr_off = bgzf_tell(fp);
      if (b.refID != iter->tid || b.pos >= iter->end) { // no need to proceed
	ret = bam_validate1(NULL, b)? -1 : -5; // determine whether end of region or error
	break;
      }
      else if (is_overlap(iter->beg, iter->end, b)) return ret;
    } else break; // end of file or error
  }
  iter->finished = 1;
  return ret;
}


int bam_iter_read1(BGZF *fp, iter_t *iter, aRead &b) {
  int ret;
  if (iter && iter->finished) return -2;
  if (iter == 0 || iter->from_first) {
    //fprintf(stderr,"[%s] reading\n",__FUNCTION__);
    ret = bam_read2(fp, b);
    if (ret < 0 && iter) iter->finished = 1;
    return ret;
  }
  if (iter->off == 0) return -1;
  for (;;) {
    if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].chunk_end) { // then jump to the next chunk
      if (iter->i == iter->n_off - 1) {
	ret = -1; break; 
      } // no more chunks
      //      if (iter->i >= 0) assert(iter->curr_off == iter->off[iter->i].chunk_end); // otherwise bug
      if (iter->i < 0 || iter->off[iter->i].chunk_end != iter->off[iter->i+1].chunk_beg) { // not adjacent chunks; then seek
	bgzf_seek(fp, iter->off[iter->i+1].chunk_beg, SEEK_SET);
	iter->curr_off = bgzf_tell(fp);
      }
      ++(iter->i);
    }
    //    fprintf(stderr,"in the for looop2 curr_off=%lu\n",iter->curr_off);
    if ((ret = bam_read2(fp, b)) >= 0) {
      iter->curr_off = bgzf_tell(fp);
      //      fprintf(stderr,"in the for looop3 curr_off=%lu\n",iter->curr_off);
      if (b.refID != iter->tid || b.pos >= iter->end) { // no need to proceed
	//	fprintf(stderr,"in the for looop4 no ned to procedd:%lu\n",iter->curr_off);
	//	fprintf(stderr,"bamret=%d\n", bam_validate1(NULL, b));
	ret = bam_validate1(NULL, b)? -1 : -5; // determine whether end of region or error

	break;
      }
      else if (is_overlap(iter->beg, iter->end, b)) return ret;
    } else break; // end of file or error
  }
  //  fprintf(stderr,"asdf;laksdf;alsdkfj ret=%d\n",ret);
  iter->finished = 1;
  return ret;
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


/*
  compare all entries in the 2 headers, if difference return 0;
*/
int compHeader(aHead *hd1,aHead *hd2){
  if(0){
    if(hd1->l_text!=hd2->l_text)
      fprintf(stderr,"problem with l_text in header\n");
    if(memcmp(hd1->text,hd2->text,hd1->l_text)!=0)
      fprintf(stderr,"problem with text in header\n");
  }
  if(hd1->n_ref!=hd2->n_ref){
    fprintf(stderr,"Difference in BAM headers: Problem with number of chromosomes in header\n");
    exit(0);
  }
  for(int i=0;i<hd1->n_ref;i++){
    if(strcasecmp(hd1->name[i],hd2->name[i])!=0){
      fprintf(stderr,"Difference in BAM headers: Problem with chromosome ordering");
      exit(0);
    }
    if(hd1->l_ref[i]!=hd2->l_ref[i]){
      fprintf(stderr,"Difference in BAM headers: Problem with length of chromosomes");
      exit(0);
    }
  }
  return 0;
}
