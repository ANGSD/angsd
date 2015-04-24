/*
  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

  april 13, safv3 added, safv2 removed for know. Will be reintroduced later.
  april 20, removed 2dsfs as special scenario
  april 20, split out the safreader into seperate cpp/h
*/

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#include "safreader.h"

typedef struct {
  char *chooseChr;
  int nSites;
  int maxIter;
  double tole;
  int nThreads;
  char *sfsfname;
  std::vector<persaf *> saf;
}args;

void destroy(args *p){
  free(p->chooseChr);
  delete p;
}

args * getArgs(int argc,char **argv){
  args *p = new args;
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
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = strdup(*(++argv));
    }
    else  if(!strcasecmp(*argv,"-start")){
      p->sfsfname = *(++argv);
    }else{
      p->saf.push_back(readsaf(*(++argv)));
    }
    argv++;
  }
  fprintf(stderr,"args: tole:%f nthreads:%d maxiter:%d nsites:%lu chooseChr:%s start:%s\n",p->tole,p->nThreads,p->maxIter,p->nSites,p->chooseChr,p->sfsfname);
  for(int i =0;i<p->saf.size();i++)
    fprintf(stderr,"input safs[%d]: %s\n",i,p->saf[i]);
  return p;
}

void print(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [chrname]\n");
    return; 
  }

  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
  persaf *saf = readsaf(bname);
  writesaf_header(stderr,saf);
  
  args *pars = getArgs(--argc,++argv);
  if(argc>0)
    pars->chooseChr = argv[1];
  float *flt = new float[saf->nChr+1];
  for(myMap::iterator it=saf->mm.begin();it!=saf->mm.end();++it){
    if(pars->chooseChr!=NULL){
      it = saf->mm.find(pars->chooseChr);
      if(it==saf->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",pars->chooseChr);
	break;
      }
    }
    bgzf_seek(saf->pos,it->second.pos,SEEK_SET);
    bgzf_seek(saf->saf,it->second.saf,SEEK_SET);
    int *ppos = new int[it->second.nSites];
    bgzf_read(saf->pos,ppos,sizeof(int)*it->second.nSites);
    for(int s=0;s<it->second.nSites;s++){
      bgzf_read(saf->saf,flt,sizeof(float)*(saf->nChr+1));
      fprintf(stdout,"%s\t%d",it->first,ppos[s]);
      for(int is=0;is<saf->nChr+1;is++)
	fprintf(stdout,"\t%f",flt[is]);
      fprintf(stdout,"\n");
    }
    delete [] ppos;
    if(pars->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  destroy(pars);
  destroy(saf);
}

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};

int SIG_COND =1;
pthread_t *thd=NULL;

template <typename T>
void destroy(Matrix<T> *ret,int x){
  for(size_t i=0;i<x;i++)
    delete [] ret->mat[i];
  delete [] ret->mat;
  delete ret;
}




template <typename T>
Matrix<T> *alloc(size_t x,size_t y){
  //fprintf(stderr,"def=%f\n",def);
  Matrix<T> *ret = new Matrix<T>;
  ret->x=x;
  ret->y=y;
  ret->mat= new T*[x];
  for(size_t i=0;i<ret->x;i++)
    ret->mat[i]=new T[y];
  return ret;
};



template<typename T>
struct emPars{
  int threadId; //size_t is the largest primitive datatype.
  Matrix<T> *GL;//<-this will be new approach
  int from;
  int to;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  int dim;
};



emPars<float> *emp = NULL;




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

template<typename T>
void readGL(BGZF *fp,size_t nSites,int dim,Matrix<T> *ret){
  ret->x=nSites;
  ret->y=dim;
  size_t i;
  for(i=0;SIG_COND&&i<nSites;i++){
    //    fprintf(stderr,"i:%lu\n",i);
    int bytes_read = bgzf_read(fp,ret->mat[i],sizeof(T)*dim);

    if(bytes_read!=0 && bytes_read<sizeof(float)*dim){
      fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0)
      break;
    
    for(size_t j=0;j<dim;j++)
      ret->mat[i][j] = exp(ret->mat[i][j]);
    
  }
  
  ret->x=i;
  if(SIG_COND==0)
    exit(0);
  
}

template<typename T>
void readGL(std::vector<persaf> &adolf,size_t nSites,std::vector< Matrix<T> > &ret){

  for(int f=0;f<adolf.size();f++){
    ret[f].y=adolf[f].nChr+1;
    size_t i;  
    for(i=0;SIG_COND&&i<nSites;i++){
      //    fprintf(stderr,"i:%lu\n",i);
      int bytes_read = bgzf_read(adolf[f].saf,ret[f].mat[i],sizeof(T)*(adolf[f].nChr+1));
      
      if(bytes_read!=0 && bytes_read<sizeof(T)*(adolf[f].nChr+1)){
	fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
	exit(0);
      }
      if(bytes_read==0)
	break;
      
      for(size_t j=0;j<adolf[f].nChr+1;j++)
	ret[f].mat[i][j] = exp(ret[f].mat[i][j]);
      
    }
    ret[f].x=i;
  }
 
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



template <typename T>
void *lik1_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];

  pars.lik = lik1(pars.sfs,pars.GL,pars.from,pars.to);
  //fprintf(stderr," thdid=%d lik=%f\n",pars.threadId,pars.lik);
  return NULL;
}


template <typename T>
double lik1_master(int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,lik1_slave<T>,(void*) i);
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


template <typename T>
void emStep1(double *pre,Matrix<T> *GL1,double *post,int start,int stop,int dim){
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


template <typename T>
void *emStep1_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];

  emStep1(pars.sfs,pars.GL,pars.post,pars.from,pars.to,pars.dim);

  return NULL;
}

template <typename T>
void emStep1_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep1_slave<T>,(void*) i);
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





template <typename T>
void em1(double *sfs,Matrix<float> *GL1,double tole,int maxIter,int nThreads,int dim){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik1_master<T>(nThreads);
  else
    oldLik = lik1(sfs,GL1,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);
  fflush(stderr);

  double tmp[dim];
  
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep1_master<T>(tmp,nThreads);
    else
      emStep1<T>(sfs,GL1,tmp,0,GL1->x,dim);
    
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];

    if(nThreads>1)
      lik = lik1_master<T>(nThreads);
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

template<typename T>
emPars<T> *setThreadPars(Matrix<T> *GL1,double *sfs,int nThreads,int dim){
  emPars<T> *temp = new emPars<T>[nThreads];
  int blockSize = GL1->x/nThreads;
  for(int i=0;i<nThreads;i++){
    temp[i].threadId = i;
    temp[i].GL=GL1;
    temp[i].from =0;
    temp[i].to=blockSize;
    temp[i].sfs = sfs;
    temp[i].post=new double[dim];
    temp[i].dim = dim;
  }
  //redo the from,to
  for(int i=1;i<nThreads;i++){
    temp[i].from = temp[i-1].to;
    temp[i].to = temp[i].from+blockSize;
  }
  //fix last end point
  temp[nThreads-1].to=GL1->x;
#if 0
  for(int i=0;i<nThreads;i++)
    fprintf(stderr,"%d:(%d,%d)=%d ",temp[i].threadId,temp[i].from,temp[i].to,temp[i].to-temp[i].from);
  fprintf(stderr,"\n");
#endif 
  
  thd= new pthread_t[nThreads];
  return temp;
}

template<typename T>
void destroy(emPars<T> *a,int nThreads ){
  for(int i=0;i<nThreads;i++)
    delete [] a[i].post;
  delete [] a;
}


template<typename T>
void setThreadPars(std::vector<Matrix < T > >  &gls,double *sfs,int nThreads,int tdim){
  emp = new emPars<T>[nThreads];
  int blockSize = gls[0]->x/nThreads;
  for(int i=0;i<nThreads;i++){
    emp[i].threadId = i;
    emp[i].GL=&gls[i];
    emp[i].from =0;
    emp[i].to=blockSize;
    emp[i].sfs = sfs;
    emp[i].post=new double[tdim];
  }
  //redo the from,to
  for(int i=1;i<nThreads;i++){
    emp[i].from = emp[i-1].to;
    emp[i].to = emp[i].from+blockSize;
  }
  //fix last end point
  emp[nThreads-1].to=gls[0].x;
#if 0
  for(int i=0;i<nThreads;i++)
    fprintf(stderr,"%d:(%d,%d)=%d ",emp[i].threadId,emp[i].from,emp[i].to,emp[i].to-emp[i].from);
  fprintf(stderr,"\n");
#endif 
  
  thd= new pthread_t[nThreads];
}

template<typename T>
double lik2(double *sfs,Matrix<T> *GL1,Matrix<T> *GL2,size_t start,size_t stop){
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



template<typename T>
void *lik2_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];

  pars.lik = lik2(pars.sfs,pars.GL1,pars.GL2,pars.from,pars.to);
  //fprintf(stderr," thdid=%d lik=%f\n",pars.threadId,pars.lik);
  return NULL;
}


template<typename T>
double lik2_master(int nThreads){
  
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,lik2_slave<T>,(void*) i);
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



template <typename T>
void emStep2(double *pre,Matrix<T> *GL1,Matrix<T> *GL2,double *post,int start,int stop,int dim){
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
template<typename T>
void *emStep2_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];

  emStep2(pars.sfs,pars.GL1,pars.GL2,pars.post,pars.from,pars.to,pars.dim);

  return NULL;
}

template<typename T>
void emStep2_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep2_slave<T>,(void*) i);
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



template <typename T>
int main_1dsfs(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return 0;
  }


  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming .saf.idx:%s\n",bname);

  persaf *saf = readsaf(bname);
  //  writesaf_header(stderr,p1);

  args *ar = getArgs(--argc,++argv);
  
  int nSites = 0;
  if(nSites == 0){//if no -nSites is specified
    if(saf->fsize>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    //this doesnt make sense if ppl supply a filelist containing safs
    nSites=saf->nSites;
  }
  fprintf(stderr,"nChr:%lu startsfs:%s nThreads:%d ",saf->nChr,ar->sfsfname,ar->nThreads);
fprintf(stderr," tole=%f maxIter=%d nSites=%lu\n",ar->tole,ar->maxIter,nSites);
  float bytes_req_megs = saf->fsize/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  

  Matrix<T> *GL1=alloc<T>(nSites,saf->nChr+1);
  double *sfs=new double[saf->nChr+1];
  
  while(1) {
    readGL(saf->saf,nSites,saf->nChr+1,GL1);
    
    if(GL1->x==0)
      break;
    fprintf(stderr,"dim(GL1)=%zu,%zu\n",GL1->x,GL1->y);
   
    
  
    if(ar->sfsfname!=NULL){
      readSFS(ar->sfsfname,saf->nChr+1,sfs);
    }else{
      
      for(int i=0;i<saf->nChr+1;i++)
	sfs[i] = (i+1)/((double)(saf->nChr+1));
      normalize(sfs,saf->nChr+1);
    }
    emp = setThreadPars<T>(GL1,sfs,ar->nThreads,saf->nChr+1);
    em1<T>(sfs,GL1,ar->tole,ar->maxIter,ar->nThreads,saf->nChr+1);

    double lik;
    if(ar->nThreads>1)
      lik = lik1_master<T>(ar->nThreads);
    else
      lik = lik1(sfs,GL1,0,GL1->x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    for(int x=0;x<=saf->nChr;x++)
      fprintf(stdout,"%f ",log(sfs[x]));
    fprintf(stdout,"\n");
    fflush(stdout);
#endif
    
  }
  destroy<T>(emp,ar->nThreads);
  destroy(GL1,nSites);
  destroy(saf);

  destroy(ar);
  delete [] sfs;
  delete [] thd;

  return 0;
}

size_t fsizes(std::vector<persaf> &pp){
  size_t res = 0;
  for(int i=0;i<pp.size();i++)
    res += pp[i].fsize;
  return res;
}



int main(int argc,char **argv){
  //start of signal handling
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  if(argc==1){
    fprintf(stderr,"Make better output. All info below is potential wrong should be upated\n");
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
  
  
  if(!strcasecmp(*argv,"print"))
    print(--argc,++argv);
  else {
    args *arg = getArgs(argc,argv);
    if(arg->saf.size()==1)
      main_1dsfs<float>(argc,argv);
    //    else
    //  main_2dsfs<float>(argc,argv);
  }

  return 0;
}
