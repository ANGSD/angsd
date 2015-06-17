#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <csignal>
#include <vector>
#include <stdint.h>
#include <algorithm>
#include <htslib/hts.h>
#include "mUpPile.h"
#include "parseArgs_bambi.h"
#include "abc.h"
#include "abcGetFasta.h"
#include "analysisFunction.h"



extern int SIG_COND;
#define bam_nt16_rev_table seq_nt16_str



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
  if(ret.hdr)
    bam_hdr_destroy(ret.hdr);
  if(ret.itr)
    hts_itr_destroy(ret.itr);
  //  fprintf(stderr,"idx:%p\n",ret.idx);
  //  exit(0);
  if(ret.idx)
    hts_idx_destroy(ret.idx);
  free(ret.fn);
  hts_close(ret.fp);
}




void printHd(const bam_hdr_t *hd,FILE *fp){
  fprintf(fp,"htext=%s\n",hd->text);
  fprintf(fp,"n_ref=%d\n",hd->n_targets);
  for(int i=0;i<hd->n_targets;i++)
    fprintf(fp,"i=%d name=%s length=%d\n",i,hd->target_name[i],hd->target_len[i]);

}


int checkIfSorted(char *str){
  //check if proper header exists
  if(strncmp(str,"@HD",3)!=0){
    fprintf(stderr,"\t-> We require a proper header starting with @HD for ANGSD\n");
    fprintf(stderr,"\t-> We observed: \'%.10s\' will exit\n",str);
    return 1;
  }
  //check if SO:coordinate exists
  char *so = strstr(str,"SO:coordinate");
  if(so==NULL){
    fprintf(stderr,"\t-> We require files to be sorted by coordinate\n");
    return 2;
  }
  if(strchr(str,'\n')<so){
    fprintf(stderr,"\t-> We require a SO:coordinate tag in the firs line of header\n");
    return 3;
  }
  return 0;
}



bufReader initBufReader2(const char*fname){
  bufReader ret;
  ret.fn = strdup(fname);
  int newlen=strlen(fname);//<-just to avoid valgrind -O3 uninitialized warning
  ret.fp = openBAM(ret.fn);
  ret.isEOF =0;
  ret.itr=NULL;
  ret.idx=NULL;
 
  ret.hdr = sam_hdr_read(ret.fp);
  if(checkIfSorted(ret.hdr->text))
    exit(0);
  if(ret.hdr==NULL) {
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
    if(compHeader(ret[0].hdr,ret[i].hdr)){
      fprintf(stderr,"Difference in BAM headers for \'%s\' and \'%s\'\n",vec[0],vec[i]);
      fprintf(stderr,"HEADER BAM1\n");
      printHd(ret[0].hdr,stderr);
      printHd(ret[i].hdr,stderr);
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

void printReg(FILE *fp,std::vector<regs> &regions){
  fprintf(fp,"-------------\n");
  fprintf(fp,"regions.size()=%lu\n",regions.size());
  for(size_t i=0;i<regions.size();i++)
    fprintf(fp,"reg[%zu]= %d %d %d\n",i,regions[i].refID,regions[i].start,regions[i].stop);
  fprintf(fp,"-------------\n");
}


extern abc **allMethods;
abcGetFasta *gf=NULL;

int bammer_main(argStruct *args){

  gf=(abcGetFasta *) allMethods[1];

  //read bamfiles
  extern int checkBamHeaders;
  bufReader *rd = initializeBufReaders2(args->nams,checkBamHeaders);

  extern int maxThreads;
  
  uppile(args->show,maxThreads,rd,args->nLines,args->nams.size(),args->regions,gf);

  //cleanup stuff
  for(unsigned i=0;i<args->nams.size();i++)
    dalloc_bufReader(rd[i]);
  
  delete [] rd;
  return 0;
}
