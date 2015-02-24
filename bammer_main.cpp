/*
  2 leaks
  1) when choosing region
  2) when doing strdup in indexing
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <csignal>
#include <vector>
#include <stdint.h>
#include <algorithm>
#include <htslib/hts.h>
#include "bams.h"
#include "mUpPile.h"
#include "parseArgs_bambi.h"
#include "indexer.h"
#include "abc.h"
#include "abcGetFasta.h"
#include "analysisFunction.h"


extern abc **allMethods;
abcGetFasta *gf=NULL;


extern int SIG_COND;

static const char *bam_nt16_rev_table2 = "=ACMGRSVTWYHKDBN";



void printHd(const bam_hdr_t *hd,FILE *fp){
  fprintf(fp,"htext=%s\n",hd->text);
  fprintf(fp,"n_ref=%d\n",hd->n_targets);
  for(int i=0;i<hd->n_targets;i++)
    fprintf(fp,"i=%d name=%s length=%d\n",i,hd->target_name[i],hd->target_len[i]);

}




/*
  compare all entries in the 2 headers, if difference return 1;
*/

int compHeader(bam_hdr_t *hd1,bam_hdr_t *hd2){
  if(0){
    if(hd1->l_text!=hd2->l_text)
      fprintf(stderr,"problem with l_text in header\n");
    if(memcmp(hd1->text,hd2->text,hd1->l_text)!=0)
      fprintf(stderr,"problem with text in header\n");
  }
  if(hd1->n_targets!=hd2->n_targets){
    fprintf(stderr,"Difference in BAM headers: Problem with number of chromosomes in header\n");
    return 1;
  }
  for(int i=0;i<hd1->n_targets;i++){
    if(strcasecmp(hd1->target_name[i],hd2->target_name[i])!=0){
      fprintf(stderr,"Difference in BAM headers: Problem with chromosome ordering");
      return 1;
    }
    if(hd1->target_len[i]!=hd2->target_len[i]){
      fprintf(stderr,"Difference in BAM headers: Problem with length of chromosomes");
      return 1;
    }
  }
  return 0;
}



void dalloc_bufReader(bufReader &ret){
  if(ret.hd)
    bam_hdr_destroy(ret.hd);
  if(ret.it.hts_itr)
    hts_itr_destroy(ret.it.hts_itr);
  if(ret.it.hts_idx)
    hts_idx_destroy(ret.it.hts_idx);
  free(ret.fn);
  free(ret.it.off);//cleanup offsets if used
  hts_close(ret.fp);
}



bufReader initBufReader2(const char*fname){
  bufReader ret;
  ret.fn = strdup(fname);
  int newlen=strlen(fname);//<-just to avoid valgrind -O3 uninitialized warning
  ret.fp = openBAM(ret.fn);
  //  fprintf(stderr,"[%s] Not implemented fromat.format:%d cram:%d\n",__FUNCTION__,ret.fp->format.format,cram);  

  ret.isEOF =0;
  ret.it.from_first=1;//iterator, used when supplying regions
  ret.it.finished=0;//iterator, used when supplying regions
  ret.it.off = NULL;
  ret.it.dasIndex = NULL;

  ret.it.hts_itr=NULL;
  ret.it.hts_idx=NULL;
 
  ret.hd = sam_hdr_read(ret.fp);
  if(ret.hd==NULL) {
    fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", ret.fn);
    exit(0);
  }
  
  return ret;
}

 typedef struct{
    size_t val;
    int key;
  }pair;

bool operator < (const pair& v1, const pair& v2)
  {
    return v1.val > v2.val;
  }

int *bamSortedIds = NULL;
bufReader *initializeBufReaders2(const std::vector<char *> &vec,int exitOnError){
  bufReader *ret = new bufReader[vec.size()];

  for(size_t i =0;i<vec.size();i++)
    ret[i] = initBufReader2(vec[i]);
  //now all readers are inialized, lets validate the header is the same
  for(size_t i=1;i<vec.size();i++)
    if(compHeader(ret[0].hd,ret[i].hd)){
      fprintf(stderr,"Difference in BAM headers for \'%s\' and \'%s\'\n",vec[0],vec[i]);
      fprintf(stderr,"HEADER BAM1\n");
      printHd(ret[0].hd,stderr);
      printHd(ret[i].hd,stderr);
      if(exitOnError)
	exit(0);
    }
 
  
  pair *val_keys = new pair[vec.size()];
  size_t fsize(const char* fname);
  for(size_t i=0;i<vec.size();i++){
    pair p;
    if (strstr(vec[i], "ftp://") == vec[i] || strstr(vec[i], "http://") == vec[i]){
      p.val = 100+i;//doesntmatter
      p.key = i;

    }else{
      p.val = aio::fsize(vec[i]);
      p.key = i;
      val_keys[i]=p;
    }
  }
  std::sort(val_keys,val_keys+vec.size());

  #if 0
  for(int i=0;i<vec.size();i++)
    fprintf(stderr,"i:%d key:%zu val:%d\n",i,val_keys[i].val,val_keys[i].key);
  #endif

  bamSortedIds = new int[vec.size()];
  for(int i=0;i<vec.size();i++)
    bamSortedIds[i] = val_keys[i].key;
  delete [] val_keys;

  return ret;
}



void printAuxBuffered(uint8_t *s, uint8_t *sStop,kstring_t &str ) {
  //  fprintf(stderr,"\ncomp:%p vs %p\n",s,sStop);
  
  while (s < sStop) {
    uint8_t type;
    kputc('\t', &str);kputc(s[0], &str);kputc(s[1], &str); kputc(':', &str); 
    //    fprintf(stderr,"\t%c%c:",s[0],s[1]);
    s += 2; type = *s; ++s;
    //    fprintf(stderr,"\ntype=%c\n",type);//,(char)*s);
    //    kputc('\t', &str); kputsn((char*)key, 2, &str); kputc(':', &str);
    if (type == 'A') { kputsn("A:", 2, &str); kputc(*s, &str); ++s; }
    else if (type == 'C') { kputsn("i:", 2, &str); kputw(*s, &str); ++s; }
    else if (type == 'c') { kputsn("i:", 2, &str); kputw(*(int8_t*)s, &str); ++s; }
    else if (type == 'S') { kputsn("i:", 2, &str); kputw(*(uint16_t*)s, &str); s += 2; }
    else if (type == 's') { kputsn("i:", 2, &str); kputw(*(int16_t*)s, &str); s += 2; }
    else if (type == 'I') { kputsn("i:", 2, &str); kputuw(*(uint32_t*)s, &str); s += 4; }
    else if (type == 'i') { kputsn("i:", 2, &str); kputw(*(int32_t*)s, &str); s += 4; }
    else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
    else if (type == 'd') { ksprintf(&str, "d:%lg", *(double*)s); s += 8; }
    else if (type == 'Z' || type == 'H') { kputc(type, &str); kputc(':', &str); while (*s) kputc(*s++, &str); ++s; }
    else if (type == 'B') {
      uint8_t sub_type = *(s++);
      int32_t n;
      memcpy(&n, s, 4);
      s += 4; // no point to the start of the array
      kputc(type, &str); kputc(':', &str); kputc(sub_type, &str); // write the typing
      for (int i = 0; i < n; ++i) {
	kputc(',', &str);
	if ('c' == sub_type || 'c' == sub_type) { kputw(*(int8_t*)s, &str); ++s; }
	else if ('C' == sub_type) { kputw(*(uint8_t*)s, &str); ++s; }
	else if ('s' == sub_type) { kputw(*(int16_t*)s, &str); s += 2; }
	else if ('S' == sub_type) { kputw(*(uint16_t*)s, &str); s += 2; }
	else if ('i' == sub_type) { kputw(*(int32_t*)s, &str); s += 4; }
	else if ('I' == sub_type) { kputuw(*(uint32_t*)s, &str); s += 4; }
	else if ('f' == sub_type) { ksprintf(&str, "%g", *(float*)s); s += 4; }
      }
    }
  }
  //  fprintf(stderr,"done\n");
}

void printChunky2(const chunky* chk,FILE *fp,char *refStr) {
  //  fprintf(stderr,"[%s] nsites=%d region=(%d,%d) itrReg=(%d,%d)\n",__FUNCTION__,chk->nSites,chk->refPos[0],chk->refPos[chk->nSites-1],chk->regStart,chk->regStop);
  if(chk->refPos[0]>chk->regStop){
    fprintf(stderr,"\t->Problems with stuff\n");
    exit(0);
  }
  int refId = chk->refId;
  for(int s=0;s<chk->nSites;s++) {
    if(chk->refPos[s]<chk->regStart || chk->refPos[s]>chk->regStop-1 ){
      for(int i=0;i<chk->nSamples;i++)
	dalloc_node(chk->nd[s][i]);
      delete [] chk->nd[s];
      continue;
    }
    fprintf(fp,"%s\t%d",refStr,chk->refPos[s]+1);     
    if(gf->ref!=NULL){
      if(refId!=gf->ref->curChr)
	gf->loadChr(gf->ref,refStr,refId);
      if(gf->ref->seqs!=NULL){
	fprintf(fp,"\t%c",gf->ref->seqs[chk->refPos[s]]);
      }
    }
    for(int i=0;i<chk->nSamples;i++) {

      //      fprintf(stderr,"seqlen[%d,%d]=%lu\t",s,i,chk->nd[s][i].seq->l);
      if(chk->nd[s][i].seq.l!=0){
	fprintf(fp,"\t%d\t",chk->nd[s][i].depth);
	for(size_t l=0;l<chk->nd[s][i].seq.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].seq.s[l]);
	fprintf(fp,"\t");
		
	for(size_t l=0;l<chk->nd[s][i].qs.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].qs.s[l]);
	//	fprintf(fp,"\t");

      }else
	fprintf(fp,"\t0\t*\t*");
      dalloc_node(chk->nd[s][i]);
    }
    //    fprintf(stderr,"\n");
    fprintf(fp,"\n");
    delete [] chk->nd[s];
  }
  delete [] chk->nd;
  delete [] chk->refPos;
  delete chk;
} 

void (*func)(void *) = NULL;

void printReg(FILE *fp,std::vector<regs> &regions){
  fprintf(fp,"-------------\n");
  fprintf(fp,"regions.size()=%lu\n",regions.size());
  for(size_t i=0;i<regions.size();i++)
    fprintf(fp,"reg[%zu]= %d %d %d\n",i,regions[i].refID,regions[i].start,regions[i].stop);
  fprintf(fp,"-------------\n");
}



int bammer_main(argStruct *args){

  gf=(abcGetFasta *) allMethods[1];

  //read bamfiles
  extern int checkBamHeaders;
  bufReader *rd = initializeBufReaders2(args->nams,checkBamHeaders);

  extern int maxThreads;
  
  uppile(args->show,maxThreads,rd,args->nLines,args->nams.size(),args->regions);

  //cleanup stuff
  for(unsigned i=0;i<args->nams.size();i++)
    dalloc_bufReader(rd[i]);
  
  delete [] rd;
  return 0;
}
