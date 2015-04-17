/*
  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

  april 13, safv3 added, safv2 removed for know. Will be reintroduced later.
  
 */

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>
#include <htslib/bgzf.h>
#include <map>
#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#define LENS 4096

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};

int fexists(const char* str){
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.                             
}

int SIG_COND =1;
pthread_t *thd=NULL;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef struct{
  size_t nSites;
  int64_t pos;
  int64_t saf;
}datum;

typedef std::map<char*,datum,ltstr> myMap;

typedef struct{
  size_t nSites;
  size_t nChr;
  myMap mm;
  BGZF *pos;
  BGZF *saf;
  size_t fsize;//contains an estimate of the uncompressed fsize
}perpop;

void dalloc(myMap &mm){
  for(myMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}

void dalloc(perpop &pp){
  bgzf_close(pp.pos);
  bgzf_close(pp.saf);
  dalloc(pp.mm);
}



void writePerPop(FILE *fp,perpop &pp){
  fprintf(fp,"\t\tInformation from index file: nChr:%lu nSites:%lu\n",pp.nChr,pp.nSites);
  
  int i=0;
  for(myMap::const_iterator it=pp.mm.begin();it!=pp.mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%ld\n",i++,it->first,d.nSites,(long int)d.pos,(long int)d.saf);
  }

}


perpop getMap(const char *fname){
  perpop ret;
  ret.pos=ret.saf=NULL;

  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }
  FILE *fp = fopen(fname,"r");
  char buf[8];
  assert(fread(buf,1,8,fp)==8);
 
  if(1!=fread(&ret.nChr,sizeof(size_t),1,fp)){
    fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
    exit(0);
  }
  ret.nSites =0;
  while(fread(&clen,sizeof(size_t),1,fp)){
    char *chr = (char*)malloc(clen+1);
    assert(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';
    
    datum d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    ret.nSites += d.nSites;
    if(1!=fread(&d.pos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    if(1!=fread(&d.saf,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
  
    myMap::iterator it = ret.mm.find(chr);
    if(it==ret.mm.end())
      ret.mm[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists\n",chr);
      exit(0);
    }
  }
  fclose(fp);
  char *tmp =(char*)calloc(strlen(fname)+100,1);//that should do it
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  //  fprintf(stderr,"tmp:%s\n",tmp);
  
  char *tmp2 = (char*)calloc(strlen(fname)+100,1);//that should do it
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  fprintf(stderr,"\t-> Assuming .saf.gz file: %s\n",tmp2);
  ret.saf = bgzf_open(tmp2,"r");bgzf_seek(ret.saf,8,SEEK_SET);
  snprintf(tmp2,strlen(fname)+100,"%spos.gz",tmp);
  fprintf(stderr,"\t-> Assuming .saf.pos.gz: %s\n",tmp2);
  ret.pos = bgzf_open(tmp2,"r");bgzf_seek(ret.pos,8,SEEK_SET);
  assert(ret.pos!=NULL&&ret.saf!=NULL);
  free(tmp);free(tmp2);
  
  ret.fsize = sizeof(float)*ret.nSites*(ret.nChr+1)+sizeof(double *)*ret.nSites;
  

  return ret;
}




/*
  loads first 8bytes and checks magic value

  //return 0 if it is old format
  //returns 1 if safv2
  //returns 2 if safv3

 */

int isNewFormat(const char *fname){
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",fname);
    exit(0);
  }
  char buf[8];
  gzread(gz,buf,8*sizeof(char));
  //if(0==strcmp(buf,"safv2"))
  //  fprintf(stderr,"File is new format... Opening string is: %s\n",buf);
  gzclose(gz);
  if(0==strcmp(buf,"safv2"))
    return 1;
  else if(0==strcmp(buf,"safv3"))
    return 2;
  else 
    return 0;
}





template <typename T>
Matrix<float> alloc(size_t x,size_t y){
  //fprintf(stderr,"def=%f\n",def);
  Matrix<float> ret;
  ret.x=x;
  ret.y=y;
  ret.mat= new float*[x];
  for(size_t i=0;i<x;i++)
    ret.mat[i]=new float[y];
  return ret;
};

void dalloc(Matrix<float> &ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret.mat[i];
  delete [] ret.mat;
}



typedef struct emPars_t{
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
  Matrix<float> *GL1;
  Matrix<float> *GL2;
  int from;
  int to;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  int dim;
} emPars;

emPars *emp = NULL;



void normalize(double *tmp,int len){
  double s=0;
  for(int i=0;i<len;i++)
    s += tmp[i];
  for(int i=0;i<len;i++)
    tmp[i] /=s;
}



#ifdef __APPLE__
size_t getTotalSystemMemory(){
  uint64_t mem;
  size_t len = sizeof(mem);
  sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
  return mem;
}
#else
size_t getTotalSystemMemory(){
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
#endif





char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

template<typename T>
void readGL(BGZF *fp,size_t nSites,int dim,Matrix<T> &ret){
  ret.x=nSites;
  ret.y=dim;
  size_t i;
  for(i=0;SIG_COND&&i<nSites;i++){
    //    fprintf(stderr,"i:%lu\n",i);
    int bytes_read = bgzf_read(fp,ret.mat[i],sizeof(float)*dim);

    if(bytes_read!=0 && bytes_read<sizeof(float)*dim){
      fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0)
      break;
    
    for(size_t j=0;j<dim;j++)
      ret.mat[i][j] = exp(ret.mat[i][j]);
    
  }
  
  ret.x=i;
  if(SIG_COND==0)
    exit(0);
  
}


size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

void readSFS(const char*fname,int hint,double *ret){
  fprintf(stderr,"reading: %s\n",fname);
  FILE *fp = NULL;
  if(((fp=fopen(fname,"r")))){
    fprintf(stderr,"problems opening file:%s\n",fname);
    exit(0);
  }
  char buf[fsize(fname)+1];
  if(fsize(fname)!=fread(buf,sizeof(char),fsize(fname),fp)){
    fprintf(stderr,"Problems reading file: %s\n will exit\n",fname);
    exit(0);
  }
  buf[fsize(fname)]='\0';
  std::vector<double> res;
  char *tok=NULL;
  tok = strtok(buf,"\t\n ");
  if(!tok){
    fprintf(stderr,"File:%s looks empty\n",fname);
    exit(0);
  }
  res.push_back(atof(tok));

  while((tok=strtok(NULL,"\t\n "))) {  
    //fprintf(stderr,"%s\n",tok);
    res.push_back(atof(tok));

  }
  //  fprintf(stderr,"size of prior=%lu\n",res.size());
  if(hint!=res.size()){
    fprintf(stderr,"problem with size of dimension of prior %d vs %lu\n",hint,res.size());
    for(size_t i=0;0&&i<res.size();i++)
      fprintf(stderr,"%zu=%f\n",i,res[i]);
    exit(0);
  }
  for(size_t i=0;i<res.size();i++){
    
    ret[i] = exp(res[i]);
    // fprintf(stderr,"i=%d %f\n",i,ret[i]);
  }
  fclose(fp);
}

typedef struct {
  char *chooseChr;
  int calcLike;
  int nSites;
  int maxIter;
  double tole;
  int nThreads;
  char *sfsfname;
}args;

void dalloc(args *p){
  free(p->chooseChr);
}

args * getArgs(int argc,char **argv){
  args *p =(args*) calloc(1,sizeof(args));
  p->sfsfname=p->chooseChr=NULL;
  p->maxIter=1e2;
  p->tole=1e-6;
  p->nThreads=4;
  p->nSites =0;
  if(argc==0)
    return p;
  if(argc%2){
    fprintf(stderr,"Extra args must be given as -par VALUE\n");
    exit(0);
  }
  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-P"))
      p->nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-calcLike"))
      p->calcLike = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = strdup(*(++argv));
    }
    else  if(!strcasecmp(*argv,"-start")){
      p->sfsfname = *(++argv);
    }else{
      fprintf(stderr,"Unknown arg: %s\n",*argv);
      exit(0);
    }
    argv++;
  }
  fprintf(stderr,"args: tole:%f nthreads:%d maxiter:%d nsites:%lu calcLike:%d chooseChr:%s start:%s\n",p->tole,p->nThreads,p->maxIter,p->nSites,p->calcLike,p->chooseChr,p->sfsfname);
  return p;
}

template <typename T>
double lik1(double *sfs,Matrix<T> *ret,int from,int to){
  double r =0;
  for(int s=from;s<to;s++){
    double tmp =0;
    for(int i=0;i<ret->y;i++)
      tmp += sfs[i]* ret->mat[s][i];
    r += log(tmp);
  }
  return r;
}



void *lik1_slave(void *p){
  emPars &pars = emp[(size_t) p];

  pars.lik = lik1(pars.sfs,pars.GL1,pars.from,pars.to);
  //fprintf(stderr," thdid=%d lik=%f\n",pars.threadId,pars.lik);
  return NULL;
}



double lik1_master(int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,lik1_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  double res=0;
  for(int i=0;i<nThreads;i++){
    //    fprintf(stderr,"lik=%f\n",emp[i].lik);
    res += emp[i].lik;
  
  }
  return res;
}


void emStep1(double *pre,Matrix<float> *GL1,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    for(int x=0;x<dim;x++)
      inner[x] = pre[x]*GL1->mat[s][x];
  
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}



void *emStep1_slave(void *p){
  emPars &pars = emp[(size_t) p];

  emStep1(pars.sfs,pars.GL1,pars.post,pars.from,pars.to,pars.dim);

  return NULL;
}


void emStep1_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep1_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  memcpy(post,emp[0].post,emp[0].dim*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<emp[0].dim;j++)
      post[j] += emp[i].post[j];
  }
  
  normalize(post,emp[0].dim);

#if 0
  for(int i=0;i<nThreads;i++){
    for(int j=0;j<dim;j++)
      fprintf(stdout,"%f ",emp[i].post[j]);
    fprintf(stdout,"\n");
  }
#endif
  
}






void em1(double *sfs,Matrix<float> *GL1,double tole,int maxIter,int nThreads,int dim){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik1_master(nThreads);
  else
    oldLik = lik1(sfs,GL1,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);
  fflush(stderr);

  double tmp[dim];
  
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep1_master(tmp,nThreads);
    else
      emStep1(sfs,GL1,tmp,0,GL1->x,dim);
    
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];

    if(nThreads>1)
      lik = lik1_master(nThreads);
    else
      lik = lik1(sfs,GL1,0,GL1->x);

    fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  
}


double lik2(double *sfs,Matrix<float> *GL1,Matrix<float> *GL2,size_t start,size_t stop){
  double res =0;
  for(int s=start;SIG_COND &&s<stop;s++){
    //    fprintf(stderr,"s=%d\n",s);
    double tmp =0;
    int inc =0;
    for(int x=0;x<GL1->y;x++)
      for(int y=0;y<GL2->y;y++)
	tmp += sfs[inc++]* GL1->mat[s][x]*GL2->mat[s][y];
    res +=log(tmp);
  }
  return res;
}


void setThreadPars(Matrix<float> *GL1,Matrix<float> *GL2,double *sfs,int nThreads,int dim){
  emp = new emPars[nThreads];
  int blockSize = GL1->x/nThreads;
  for(int i=0;i<nThreads;i++){
    emp[i].threadId = i;
    emp[i].GL1=GL1;
    emp[i].GL2=GL2;
    emp[i].from =0;
    emp[i].to=blockSize;
    emp[i].sfs = sfs;
    emp[i].post=new double[dim];
  }
  //redo the from,to
  for(int i=1;i<nThreads;i++){
    emp[i].from = emp[i-1].to;
    emp[i].to = emp[i].from+blockSize;
  }
  //fix last end point
  emp[nThreads-1].to=GL1->x;
#if 0
  for(int i=0;i<nThreads;i++)
    fprintf(stderr,"%d:(%d,%d)=%d ",emp[i].threadId,emp[i].from,emp[i].to,emp[i].to-emp[i].from);
  fprintf(stderr,"\n");
#endif 
  
  thd= new pthread_t[nThreads];
}



void *lik2_slave(void *p){
  emPars &pars = emp[(size_t) p];

  pars.lik = lik2(pars.sfs,pars.GL1,pars.GL2,pars.from,pars.to);
  //fprintf(stderr," thdid=%d lik=%f\n",pars.threadId,pars.lik);
  return NULL;
}



double lik2_master(int nThreads){
  
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,lik2_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  double res=0;
  for(int i=0;i<nThreads;i++){
    //    fprintf(stderr,"lik=%f\n",emp[i].lik);
    res += emp[i].lik;
  
  }
  return res;
}




void emStep2(double *pre,Matrix<float> *GL1,Matrix<float> *GL2,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(int x=0;x<GL1->y;x++)
      for(int y=0;y<GL2->y;y++){
	inner[inc] = pre[inc]*GL1->mat[s][x]*GL2->mat[s][y];
	inc++;
      }
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}

void *emStep2_slave(void *p){
  emPars &pars = emp[(size_t) p];

  emStep2(pars.sfs,pars.GL1,pars.GL2,pars.post,pars.from,pars.to,pars.dim);

  return NULL;
}


void emStep2_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep2_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  memcpy(post,emp[0].post,emp[0].dim*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<emp[0].dim;j++)
      post[j] += emp[i].post[j];
  }
  
  normalize(post,emp[0].dim);

#if 0
  for(int i=0;i<nThreads;i++){
    for(int j=0;j<dim;j++)
      fprintf(stdout,"%f ",emp[i].post[j]);
    fprintf(stdout,"\n");
  }
#endif
  
}


int really_kill =3;
int VERBOSE = 1;
void handler(int s) {

  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want \'realSFS\' to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;

}





void em2(double *sfs,Matrix<float> *GL1,Matrix<float> *GL2,double tole,int maxIter,int nThreads,int dim){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik2_master(nThreads); 
  else
    oldLik = lik2(sfs,GL1,GL2,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);

  double tmp[dim];//<- wont be cleaned up, but is only allocated once
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep2_master(tmp,nThreads);
    else
      emStep2(sfs,GL1,GL2,tmp,0,GL1->x,dim);

    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    if(nThreads>1)
      lik = lik2_master(nThreads);
    else
      lik = lik2(sfs,GL1,GL2,0,GL1->x);
    
    fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  
  
}

template <typename T>
int main_2dsfs(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./emOptim2 2dsfs pop1.saf.idx pop2.saf.idx [-start FNAME -P nThreds -tole tole -maxIter -r chrnam ] (only works if the two saf files covers the same region)\n");
    return 0;
  }
  perpop p1 = getMap(argv[1]);
  perpop p2 = getMap(argv[2]);
  assert(p1.fsize==p2.fsize);
  argc -=2;
  argv += 2;
  args *p = getArgs(argc,argv);
  if(p->nSites==0){
    if(p1.fsize+p2.fsize > getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n");
    //this doesnt make sense if ppl supply a filelist containing safs
    p->nSites=p1.fsize;
  }
  fprintf(stderr,"chr1:%lu chr2:%lu startsfs:%s nThreads=%d tole=%f maxIter=%d nSites:%lu\n",p1.nChr,p2.nChr,p->sfsfname,p->nThreads,p->tole,p->maxIter,p->nSites);
  float bytes_req_megs = p->nSites*(sizeof(T)*(p1.nChr+1) + sizeof(double)*(p2.nChr+1)+2*sizeof(T*))/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);
  
  int dim=(p1.nChr+1)*(p2.nChr+1);
  
  Matrix<T> GL1=alloc<T>(p->nSites,p1.nChr+1);
  Matrix<T> GL2=alloc<T>(p->nSites,p2.nChr+1);
  dim=GL1.y*GL2.y;
  
  double *sfs = new double[dim];
  while(1){
    readGL(p1.saf,p->nSites,p1.nChr,GL1);
    readGL(p2.saf,p->nSites,p2.nChr,GL2);
          
    assert(GL1.x==GL2.x);
    if(GL1.x==0)
      break;
    
    if(p->sfsfname!=NULL){
      readSFS(p->sfsfname,dim,sfs);
    }else{
      for(int i=0;i<dim;i++)
	sfs[i] = (i+1)/((double)dim);
      normalize(sfs,dim);
    }
    
    setThreadPars(&GL1,&GL2,sfs,p->nThreads,dim);
    if(p->calcLike==0){
      if(SIG_COND) 
	em2(sfs,&GL1,&GL2,p->tole,p->maxIter,p->nThreads,dim);
    }
    double lik;
    if(p->nThreads>1)
      lik = lik1_master(p->nThreads);
    else
      lik = lik1(sfs,&GL1,0,GL1.x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    int inc=0;
    for(int x=0;x<=p1.nChr;x++){
      for(int y=0;y<=p2.nChr+1;y++)
	fprintf(stdout,"%f ",log(sfs[inc++]));
      fprintf(stdout,"\n");
    }
#endif
   
  }
  dalloc(GL1,p->nSites);
  dalloc(GL2,p->nSites);
  return 0;
}



template <typename T>
int main_1dsfs(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return 0;
  }


  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming .saf.idx:%s\n",bname);

  perpop p1 = getMap(bname);
  writePerPop(stderr,p1);

  args *ar = getArgs(--argc,++argv);
  
  int nSites = 0;
  if(nSites == 0){//if no -nSites is specified
    if(p1.fsize>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    //this doesnt make sense if ppl supply a filelist containing safs
    nSites=p1.nSites;
  }
  fprintf(stderr,"nChr:%lu startsfs:%s nThreads:%d ",p1.nChr,ar->sfsfname,ar->nThreads);
fprintf(stderr," tole=%f maxIter=%d nSites=%lu\n",ar->tole,ar->maxIter,nSites);
  float bytes_req_megs = p1.fsize/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  

  Matrix<T> GL1=alloc<T>(nSites,p1.nChr+1);
  double *sfs=new double[p1.nChr+1];
  
  while(1) {
    readGL(p1.saf,nSites,p1.nChr+1,GL1);
    
    if(GL1.x==0)
      break;
    fprintf(stderr,"dim(GL1)=%zu,%zu\n",GL1.x,GL1.y);
   
    
  
    if(ar->sfsfname!=NULL){
      readSFS(ar->sfsfname,p1.nChr+1,sfs);
    }else{
      
      for(int i=0;i<p1.nChr+1;i++)
	sfs[i] = (i+1)/((double)(p1.nChr+1));
      normalize(sfs,p1.nChr+1);
    }
    //  em2_smart(sfs2,pops,1e-6,1e3);
    setThreadPars(&GL1,NULL,sfs,ar->nThreads,p1.nChr+1);
    if(ar->calcLike==0)
      em1(sfs,&GL1,ar->tole,ar->maxIter,ar->nThreads,p1.nChr+1);

    double lik;
    if(ar->nThreads>1)
      lik = lik1_master(ar->nThreads);
    else
      lik = lik1(sfs,&GL1,0,GL1.x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    for(int x=0;x<=p1.nChr;x++)
      fprintf(stdout,"%f ",log(sfs[x]));
    fprintf(stdout,"\n");
    fflush(stdout);
#endif
    
  }
  dalloc(GL1,ar->nSites);
  dalloc(p1);
  delete [] sfs;
  return 0;
}




int print(int argc,char **argv){
  //  fprintf(stderr,"chooseChr:%s argc:%d argv:%s\n",chooseChr,argc,*argv);
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [-r chr]\n");
    return 0;
  }

  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
  args *p = getArgs(--argc,++argv);
  perpop pp = getMap(bname);
  writePerPop(stderr,pp);
  
  float *flt = new float[pp.nChr+1];
  for(myMap::iterator it=pp.mm.begin();it!=pp.mm.end();++it){
    if(p->chooseChr!=NULL){
      it = pp.mm.find(p->chooseChr);
      if(it==pp.mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",p->chooseChr);
	break;
      }
    }
    bgzf_seek(pp.pos,it->second.pos,SEEK_SET);
    bgzf_seek(pp.saf,it->second.saf,SEEK_SET);
    int *ppos = new int[it->second.nSites];
    bgzf_read(pp.pos,ppos,sizeof(int)*it->second.nSites);
    for(int s=0;s<it->second.nSites;s++){
      bgzf_read(pp.saf,flt,sizeof(float)*(pp.nChr+1));
      fprintf(stdout,"%s\t%d",it->first,ppos[s]);
      for(int is=0;is<pp.nChr+1;is++)
	fprintf(stdout,"\t%f",flt[is]);
      fprintf(stdout,"\n");
    }
    delete [] ppos;
    if(p->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  dalloc(pp);
  return 0;
}





int main(int argc,char **argv){
  //set of signal handling

  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  if(argc==1){
    fprintf(stderr,"\n./realSFS afile.saf nChr [-start FNAME -P nThreads -tole tole -maxIter  -nSites ]\n");
    fprintf(stderr,"OR\n ./realSFS 2dsfs pop1.saf pop2.saf nchrPop1 nChrPop2 [-start FNAME -P nThreads -tole tole -maxIter  -nSites ]\n");
    fprintf(stderr,"\nnChr is the number of chromosomes. (twice the number of diploid invididuals)\n");    
    return 0;
  }
  ++argv;
  --argc;

  if(isatty(fileno(stdout))){
    fprintf(stderr,"\t-> You are printing the optimized SFS to the terminal consider dumping into a file\n");
    fprintf(stderr,"\t-> E.g.: \'./realSFS");
    for(int i=0;i<argc;i++)
      fprintf(stderr," %s",argv[i]);
    fprintf(stderr," >sfs.ml.txt\'\n");   

  }

  if(!strcasecmp(*argv,"2dsfs"))
    main_2dsfs<float>(argc,argv);
  if(!strcasecmp(*argv,"print"))
    print(--argc,++argv);
  else
    main_1dsfs<float>(argc,argv);
  

  return 0;
}
