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
#include <htslib/thread_pool.h>
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
  kstring_t lib = { 0, 0, NULL };
  for(int i=0;i<nrg;i++){
    const char *rgid = sam_hdr_line_name(h,"RG",i);
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
  free(lib.s);
  return ret;
}

htsFile *openBAM(const char *fname,int doCheck){

  htsFile *fp =NULL;
  extern htsFormat *dingding2;//<-externed from abcGetFasta. This is very bad style. Should Fix, dragon
  if((fp=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(1);
  }
  const char *str = strrchr(fname,'.');
  if(doCheck==1){
    if(str&&strcasecmp(str,".bam")!=0&&str&&strcasecmp(str,".cram")!=0&&strcasecmp(str,".sam")!=0){
      fprintf(stderr,"\t-> file:\"%s\" should be suffixed with \".bam\", \".cram\" or \".sam\"\n",fname);
      fprintf(stderr,"\t-> If you know what you are doing you can disable with -doCheck 0\n");
      exit(1);
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



/*
  The name of the sample for a bam file, used in the bcf header. angsd always has one sample
  per input file. With --ignore-RG 1 (default) it is the filename, with --ignore-RG 0 the SM
  tag of the @RG lines if the file has exactly one SM, otherwise the filename.
  A name that is already used gets _2, _3 ... added, so the number of samples is always the
  number of files. Returns a new string.
*/
static char *sampleName(const char *fname,const char *hdrtext,int ignoreRG,const bam_sample_t *sm){
  kstring_t name = {0,0,NULL};
  if(!ignoreRG){
    char *sm1 = NULL;//the first SM
    int nDistinct = 0;
    for(const char *p=hdrtext;p&&*p;){
      const char *eol = strchr(p,'\n');
      size_t len = eol?(size_t)(eol-p):strlen(p);
      if(len>3&&strncmp(p,"@RG\t",4)==0){
	//find \tSM: within this line
	for(const char *q=p;q+4<=p+len;q++){
	  if(strncmp(q,"\tSM:",4)==0){
	    const char *v = q+4;
	    size_t vl = 0;
	    while(v+vl<p+len&&v[vl]!='\t')
	      vl++;
	    if(sm1==NULL){
	      sm1 = strndup(v,vl);
	      nDistinct = 1;
	    }else if(strlen(sm1)!=vl||strncmp(sm1,v,vl)!=0)
	      nDistinct = 2;
	    break;
	  }
	}
      }
      p = eol?eol+1:NULL;
    }
    if(nDistinct==1&&sm1[0]!='\0')
      kputs(sm1,&name);
    else if(nDistinct>1)
      fprintf(stderr,"\t-> \'%s\' has more than one SM tag, angsd analyses it as one sample called \'%s\'\n",fname,fname);
    else
      fprintf(stderr,"\t-> \'%s\' has no SM tag, the sample is called \'%s\'\n",fname,fname);
    free(sm1);
  }
  if(name.l==0)
    kputs(fname,&name);

  //make the name unique
  size_t baseLen = name.l;
  for(int k=2;;k++){
    int used = 0;
    for(int i=0;i<sm->n;i++)
      if(strcmp(sm->smpl[i],name.s)==0)
	used = 1;
    if(!used)
      break;
    name.l = baseLen;
    name.s[baseLen] = '\0';
    ksprintf(&name,"_%d",k);
  }
  if(name.l!=baseLen)
    fprintf(stderr,"\t-> The sample name for \'%s\' is already used, it is called \'%s\'\n",fname,name.s);
  return name.s;
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
  if(ret.hdr==NULL){
    fprintf(stderr,"\t-> Could not read the header of: \'%s\', is it a SAM/BAM/CRAM file?\n",fname);
    exit(1);
  }
  if(ret.hdr->text==NULL||strlen(ret.hdr->text)==0){
    fprintf(stderr,"\t-> No header information could be found for BAM/CRAM file: \'%s\' will exit\n",fname);
    exit(1);
  }

  checkIfSorted(ret.hdr->text);
  extern int MPLP_IGNORE_RG;
  char *smName = sampleName(fname,ret.hdr->text,MPLP_IGNORE_RG,sm);
  bam_smpl_add(sm, smName, 0);//one sample per file
  free(smName);
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
	    aio::doAssert(khash_str2int_get(rghash,rgs_to_add[i],&tmp)==0,1,AT,"");
	    fprintf(stderr,"\t-> [READGROUP info] Added readgroup %s for lib: %s (check: %s) -> %d\n",rgs_to_add[i],it2->first,it->first,tmp);
	    if(tmp>254){
	      fprintf(stderr,"\t-> Readgroup representation is only valid for fewer than 255 different readgroups\n");
	      exit(0);
	    }
	  }
	}
      }
    }
    for(std::map<char*,std::vector<char *> ,ltstr>::iterator it2 = lb_rg.begin();it2!=lb_rg.end();it2++){
      std::vector <char *> rgs_to_add = it2->second;
      for(int i=0;i<rgs_to_add.size();i++)
	free(rgs_to_add[i]);
      free(it2->first);
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
	exit(1);
    }
 
  
  pair *val_keys = new pair[vec.size()];
  for(size_t i=0;i<vec.size();i++){
    //the largest file is tried first when picking how far to read. fsize is 0 for urls (ftp, http, https, s3, ...),
    //they were not put in val_keys before, so the order had uninitialized file ids
    val_keys[i].val = aio::fsize(vec[i]);
    val_keys[i].key = i;
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

  //with -P >1 the files are decompressed by a htslib thread pool shared by all files,
  //decompression was done by the reading thread and was a large part of its time.
  //Each file gets a queue of 2 blocks, the htslib default (2*threads) used ~3.5MB per file with -P 8, this ~0.6MB
  htsThreadPool tpool = {NULL, 2};
  if(maxThreads>1){
    tpool.pool = hts_tpool_init(maxThreads);
    if(tpool.pool==NULL){
      fprintf(stderr,"\t-> Could not make a thread pool for reading the files, will exit\n");
      exit(1);
    }
    for(unsigned i=0;i<args->nams.size();i++)
      if(hts_set_opt(args->rd[i].fp,HTS_OPT_THREAD_POOL,&tpool)!=0){
	fprintf(stderr,"\t-> Could not use the thread pool for reading \'%s\', will exit\n",args->rd[i].fn);
	exit(1);
      }
  }

  uppile(args->show,maxThreads,args->rd,args->nReads,args->nams.size(),args->regions,gf);

  //cleanup stuff
  for(unsigned i=0;i<args->nams.size();i++)
    dalloc_bufReader(args->rd[i]);
  if(tpool.pool)
    hts_tpool_destroy(tpool.pool);//after the files are closed

  delete [] args->rd;
  return 0;
}
