#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <csignal>
#include <vector>
#include <stdint.h>
#include <algorithm>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include "mUpPile.h"
#include "parseArgs_bambi.h"
#include "abc.h"
#include "abcGetFasta.h"
#include "analysisFunction.h"
#include "sample.h"
#include "aio.h"

extern int SIG_COND;
#define bam_nt16_rev_table seq_nt16_str

//map from libid to vector of associated readgroups
std::map<char*,std::vector<char *>,ltstr > get_lb_rg(bam_hdr_t *h){
  std::map<char*,std::vector<char *>,ltstr > ret;
  int nrg = sam_hdr_count_lines(h, "RG");
  //  fprintf(stderr,"nrg: %d\n",nrg);
  for(int i=0;i<nrg;i++){
    const char *rgid = sam_hdr_line_name(h,"RG",i);
    kstring_t lib = { 0, 0, NULL };
    if (sam_hdr_find_tag_id(h, "RG", "ID", rgid, "LB", &lib)  < 0){
      fprintf(stderr,"Header problem\n");
      continue;
    }
    //fprintf(stderr,"rgid: %s lib.s:%s\n",rgid,lib.s);
    std::map<char*,std::vector<char *>,ltstr >::iterator it=ret.find(lib.s);
    if(it==ret.end()){
      std::vector<char *> vec;vec.push_back(strdup(rgid));
      ret[strdup(lib.s)] = vec;
    }else
      it->second.push_back(strdup(rgid));

  }
  return ret;
}

htsFile *openBAM(const char *fname,int doCheck){

  htsFile *fp =NULL;
  extern htsFormat *dingding2;//<-externed from abcGetFasta. This is very bad style. Should Fix, dragon
  if((fp=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }
  const char *str = strrchr(fname,'.');
  if(doCheck==1){
    if(str&&strcasecmp(str,".bam")!=0&&str&&strcasecmp(str,".cram")!=0&&strcasecmp(str,".sam")!=0){
      fprintf(stderr,"\t-> file:\"%s\" should be suffixed with \".bam\" or \".cram\"\n",fname);
      fprintf(stderr,"\t-> If you know what you are doing you can disable with -doCheck 0\"\n");
      exit(0);
    }
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
    fprintf(stderr,"\t-> We require a SO:coordinate tag in the first line of header\n");
    return 3;
  }
  return 0;
}



bufReader initBufReader2(const char*fname,int doCheck,char *fai_fname,bam_sample_t *sm){
  bufReader ret;
  ret.fn = strdup(fname);
  int newlen=strlen(fname);//<-just to avoid valgrind -O3 uninitialized warning
  ret.fp = openBAM(ret.fn,doCheck);
  if (fai_fname && hts_set_fai_filename(ret.fp, fai_fname) != 0) {
    fprintf(stderr, "[%s] failed to process %s\n",
	    __func__, fai_fname);
    exit(EXIT_FAILURE);
  }
  ret.isEOF =0;
  ret.itr=NULL;
  ret.idx=NULL;
 
  ret.hdr = sam_hdr_read(ret.fp);
  if(strlen(ret.hdr->text)==0){
    fprintf(stderr,"\t-> No header information could be found for BAM/CRAM file: \'%s\' will exit\n",fname);
    exit(1);
  }

  checkIfSorted(ret.hdr->text);
  extern int MPLP_IGNORE_RG;
  bam_smpl_add(sm, fname, MPLP_IGNORE_RG? 0 : ret.hdr->text);
  //  fprintf(stderr,"sm->smpl[0]:%s\n",sm->smpl[0]);
  if(ret.hdr==NULL) {
    fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", ret.fn);
    exit(0);
  }
  //OK this is really bad coding style.
  //But i am here pullin in the char *libname. This contains either filename for libs to include or the name of the lib to include
  extern std::map<char *,int,ltstr> LBS;//<- contains the libs to include.
  extern void *rghash;//this is the readgroups that should have the libraries from LBS
  //  fprintf(stderr,"LBS.size():%lu rghash:%p\n",LBS.size(),rghash);
  if(LBS.size()>0){
    if(!rghash)
      rghash = khash_str2int_init();//initialize if its not defined
    //first build a lib1->{rd1,rd2},lib2->{rd3,rd4} map;
    std::map<char*,std::vector<char *> ,ltstr> lb_rg = get_lb_rg(ret.hdr);
    fprintf(stderr,"\t-> Number of different libs in BAM/CRAM file: %lu\n",lb_rg.size());

    //loop over the librarie we want to have included
    for(std::map<char *,int,ltstr>::iterator it=LBS.begin();it!=LBS.end();it++){
      //then seek in the lib->rd map
      std::map<char*,std::vector<char *> ,ltstr>::iterator it2 = lb_rg.find(it->first);
      if(it2!=lb_rg.end()){
	//	fprintf(stderr,"lib: %s has %lu readgroups\n",it2->first,it2->second.size());
	//all the read groups in the vector below should be added
	std::vector <char *> rgs_to_add = it2->second;
	for(int i=0;i<rgs_to_add.size();i++){
	  int ret = khash_str2int_inc(rghash,strdup(rgs_to_add[i]));
	  if(ret==-1)
	    fprintf(stderr,"\t-> [READGROUP info] Problems adding readgroup: %s for lib: %s\n",rgs_to_add[i],it2->first);
	  else{
	    int tmp;
	    assert(khash_str2int_get(rghash,rgs_to_add[i],&tmp)==0);
	    fprintf(stderr,"\t-> [READGROUP info] Added readgroup %s for lib: %s (check: %s) -> %d\n",rgs_to_add[i],it2->first,it->first,tmp);
	    if(tmp>254){
	      fprintf(stderr,"\t-> Readgroup representation is only valid for fewer than 255 different readgroups\n");
	      exit(0);
	    }
	  }
	}
      }
    }
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
bufReader *initializeBufReaders2(const std::vector<char *> &vec,int exitOnError,int doCheck,char *fai_fname,bam_sample_t *sm){
  bufReader *ret = new bufReader[vec.size()];

  for(size_t i =0;i<vec.size();i++)
    ret[i] = initBufReader2(vec[i],doCheck,fai_fname,sm);
  //now all readers are inialized, lets validate the header is the same
  for(size_t i=1;i<vec.size();i++)
    if(compHeader(ret[0].hdr,ret[i].hdr)){
      fprintf(stderr,"Difference in BAM headers for \'%s\' and \'%s\'\n",vec[0],vec[i]);
      fprintf(stderr,"HEADER BAM1\n");
      printHd(ret[0].hdr,stderr);
      fprintf(stderr,"HEADER BAM2\n");
      printHd(ret[i].hdr,stderr);
      fprintf(stderr,"Difference in BAM headers for \'%s\' and \'%s\'\n",vec[0],vec[i]);
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

void init_bamreaders(argStruct *args){


}


int bammer_main(argStruct *args){
  gf=(abcGetFasta *) allMethods[1];
  extern int maxThreads;
  uppile(args->show,maxThreads,args->rd,args->nReads,args->nams.size(),args->regions,gf);

  //cleanup stuff
  for(unsigned i=0;i<args->nams.size();i++)
    dalloc_bufReader(args->rd[i]);
  
  delete [] args->rd;
  return 0;
}
