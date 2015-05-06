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
  may 5, seems to work well now
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
#include <zlib.h>
#include <htslib/bgzf.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "safreader.h"
#include "keep.hpp"
 
typedef struct {
  char *chooseChr;
  int start;
  int stop;
  int nSites;
  int maxIter;
  double tole;
  int nThreads;
  char *sfsfname;
  std::vector<persaf *> saf;
  int posOnly;
}args;


int SIG_COND =1;
pthread_t *thd=NULL;
int howOften =5e6;//how often should we print out (just to make sure something is happening)
void destroy_safvec(std::vector<persaf *> &saf){
  //fprintf(stderr,"destroy &saf\n");
  for(int i=0;i<saf.size();i++)
    persaf_destroy(saf[i]);
}


void destroy_args(args *p){
  //  fprintf(stderr,"destroy args\n");
  destroy_safvec(p->saf);
  delete p;
}

//just approximate
template <typename T>
size_t fsizes(std::vector<persaf *> &pp, int nSites){
  size_t res = 0;
  for(int i=0;i<pp.size();i++){
    res += nSites*(pp[i]->nChr+1)*sizeof(T)+nSites*sizeof( T*);
  }
  return res;
}

size_t helper(persaf * pp,char *chr){
  if(chr==NULL)
    return pp->nSites;
  myMap::iterator it=pp->mm.find(chr);
  if(it==pp->mm.end()){
    fprintf(stderr,"\t-> Problem finding chromosome: %s\n",chr);
    exit(0);
  }
  return it->second.nSites;
}

size_t nsites(std::vector<persaf *> &pp,args *ar){
  size_t res = helper(pp[0],ar->chooseChr);
  for(int i=1;i<pp.size();i++)
    if(helper(pp[i],ar->chooseChr) > res)
      res = helper(pp[i],ar->chooseChr);
  return res;
}


char * get_region(char *extra,int &start,int &stop) {
  if(!extra){
    fprintf(stderr,"Must supply parameter for -r option\n");
    return NULL;
  }
  if(strrchr(extra,':')==NULL){//only chromosomename
    char *ref = extra;
    start = stop = -1;;
    return ref;
  }
  char *tok=NULL;
  tok = strtok(extra,":");

  char *ref = tok;

  start =stop=-1;

  tok = extra+strlen(tok)+1;//tok now contains the rest of the string
 
  if(strlen(tok)==0)//not start and/or stop ex: chr21:
    return ref;
  

  if(tok[0]=='-'){//only contains stop ex: chr21:-stop
    tok =strtok(tok,"-");
    stop = atoi(tok);
  }else{
    //catch single point
    int isProper =0;
    for(size_t i=0;i<strlen(tok);i++)
      if(tok[i]=='-'){
	isProper=1;
	 break;
      }
    //fprintf(stderr,"isProper=%d\n",isProper);
    if(isProper){
      tok =strtok(tok,"-");
      start = atoi(tok)-1;//this is important for the zero offset
      tok = strtok(NULL,"-");
      if(tok!=NULL)
	stop = atoi(tok);
    }else{
      //single point
      stop = atoi(tok);
      start =stop -1;
      
    }
    
  }
  if(stop!=-1&&stop<start){
    fprintf(stderr,"endpoint:%d is larger than startpoint:%d\n",start,stop);
    exit(0);
    
  }
  if(0){
    fprintf(stderr,"[%s] ref=%s,start=%d,stop=%d\n",__FUNCTION__,ref,start,stop);
    exit(0);
  }
  return ref;
}




args * getArgs(int argc,char **argv){
  args *p = new args;
  p->sfsfname=p->chooseChr=NULL;
  p->start=p->stop=-1;
  p->maxIter=1e2;
  p->tole=1e-6;
  p->nThreads=4;
  p->nSites =0;
  p->posOnly = 0;

  if(argc==0)
    return p;

  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-P"))
      p->nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-posOnly"))
      p->posOnly = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = get_region(*(++argv),p->start,p->stop);
      if(!p->chooseChr)
	return NULL;
    }
    else  if(!strcasecmp(*argv,"-start")){
      p->sfsfname = *(++argv);
    }else{
      p->saf.push_back(readsaf<float>(*argv));
      //   fprintf(stderr,"toKeep:%p\n",p->saf[p->saf.size()-1]->toKeep);
    }
    argv++;
  }
  fprintf(stderr,"\t-> args: tole:%f nthreads:%d maxiter:%d nsites:%d chooseChr:%s start:%s chr:%s start:%d stop:%d\n",p->tole,p->nThreads,p->maxIter,p->nSites,p->chooseChr,p->sfsfname,p->chooseChr,p->start,p->stop);
  return p;
}

int print(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [chrname]\n");
    return 0; 
  }
  
  args *pars = getArgs(argc,argv);
  if(!pars)
    return 1;
  if(pars->saf.size()!=1){
    fprintf(stderr,"Print only implemeted for single safs\n");
    exit(0);
  }
  writesaf_header(stderr,pars->saf[0]);
  
  float *flt = new float[pars->saf[0]->nChr+1];
  for(myMap::iterator it=pars->saf[0]->mm.begin();it!=pars->saf[0]->mm.end();++it){
    if(pars->chooseChr!=NULL){
      it = pars->saf[0]->mm.find(pars->chooseChr);
      if(it==pars->saf[0]->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",pars->chooseChr);
	break;
      }
    }
    bgzf_seek(pars->saf[0]->pos,it->second.pos,SEEK_SET);
    bgzf_seek(pars->saf[0]->saf,it->second.saf,SEEK_SET);
    int *ppos = new int[it->second.nSites];
    bgzf_read(pars->saf[0]->pos,ppos,sizeof(int)*it->second.nSites);
    
    int first=0;
    if(pars->start!=-1)
      while(ppos[first]<pars->start) 
	first++;
    
    int last=it->second.nSites;
    //    fprintf(stderr,"pars-.stop:%d ppos:%d\n",pars->stop,ppos[last-1]);
    if(pars->stop!=-1&&pars->stop<=ppos[last-1]){
      last=first;
      while(ppos[last]<pars->stop) 
	last++;
    }
    //fprintf(stderr,"first:%d last:%d\n",first,last);
    int at=0;
    for(int s=0;SIG_COND&&s<it->second.nSites;s++) {
      bgzf_read(pars->saf[0]->saf,flt,sizeof(float)*(pars->saf[0]->nChr+1));
      if(at>=first&&at<last){
	if(pars->posOnly==0){
	  fprintf(stdout,"%s\t%d",it->first,ppos[s]+1);
	  for(int is=0;is<pars->saf[0]->nChr+1;is++)
	    fprintf(stdout,"\t%f",flt[is]);
	}else
	  fprintf(stdout,"%d",ppos[s]+1);
	  fprintf(stdout,"\n");
      }
      at++;
    }
    delete [] ppos;
    if(pars->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  destroy_args(pars);
  return 0;
}

#if 1
void print2(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return; 
  }
  
  args *pars = getArgs(argc,argv);
  
  if(pars->saf.size()!=1){
    fprintf(stderr,"Print only implemeted for single safs\n");
    exit(0);
  }

  writesaf_header(stderr,pars->saf[0]);
  
  float *flt = new float[pars->saf[0]->nChr+1];
  for(myMap::iterator it=pars->saf[0]->mm.begin();it!=pars->saf[0]->mm.end();++it){
    
    if(pars->chooseChr!=NULL)
      it = iter_init(pars->saf[0],pars->chooseChr,pars->start,pars->stop);

    int *ppos = new int [it->second.nSites];
    bgzf_seek(pars->saf[0]->pos,it->second.pos,SEEK_SET);
    bgzf_read(pars->saf[0]->pos,ppos,sizeof(int)*it->second.nSites);
    int ret;
    fprintf(stderr,"in print2 first:%lu last:%lu\n",pars->saf[0]->toKeep->first,pars->saf[0]->toKeep->last);
    while((ret=iter_read(pars->saf[0],flt,sizeof(float)*(pars->saf[0]->nChr+1)))){
   
      //      fprintf(stderr,"[%s] pars->saf[0]->at:%d nSites: %lu ret:%d\n",__FUNCTION__,pars->saf[0]->at,it->second.nSites,ret);
      fprintf(stdout,"%s\t%d",it->first,ppos[pars->saf[0]->at]+1);
      for(int is=0;is<pars->saf[0]->nChr+1;is++)
	fprintf(stdout,"\t%f",flt[is]);
      fprintf(stdout,"\n");
    }
    delete [] ppos;
    //fprintf(stderr,"[%s] after while:%d\n",__FUNCTION__,ret);
    if(pars->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  destroy_args(pars);
}
#endif

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};


template <typename T>
void destroy(Matrix<T> *ret,int x){
  for(size_t i=0;i<x;i++)
    delete [] ret->mat[i];
  delete [] ret->mat;
  delete ret;
}

template <typename T>
void destroy(std::vector< Matrix<T> * > &gls,int x){
  for(size_t i=0;i<gls.size();i++)
    destroy(gls[i],x);
}



template <typename T>
void matrix_print( Matrix<T> *gls){
  for(size_t s=0;s<gls->x;s++){
    for(size_t i=0;i<gls->y;i++)
      fprintf(stderr,"\t%f",gls->mat[s][i]);
  fprintf(stderr,"\n");
  }
}




template <typename T>
Matrix<T> *alloc(size_t x,size_t y){
  Matrix<T> *ret = new Matrix<T>;
  ret->x=x;
  ret->y=y;
  ret->mat= new T*[x];
  for(size_t i=0;i<ret->x;i++)
    ret->mat[i]=new T[y];
  ret->x=0;
  return ret;
};



template<typename T>
struct emPars{
  int threadId;
  std::vector <Matrix<T> *> gls;
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
void readGL(persaf *fp,size_t nSites,int dim,Matrix<T> *ret){

  // ret->x=nSites;
  ret->y=dim;
  size_t i;
  for(i=ret->x;SIG_COND&&i<nSites;i++){
    if(i>0 &&(i% howOften)==0  )
      fprintf(stderr,"\r\t-> Has read 5mio sites now at: %lu      ",i);
    //
    
    int bytes_read= iter_read(fp,ret->mat[i],sizeof(T)*dim);//bgzf_read(fp,ret->mat[i],sizeof(T)*dim);
    if(bytes_read!=0 && bytes_read<sizeof(T)*dim){
      fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0){
      //fprintf(stderr,"[%s] bytes_read==0 at i:%lu\n",__FUNCTION__,i);
      break;
    }
    for(size_t j=0;j<dim;j++)
      ret->mat[i][j] = exp(ret->mat[i][j]);
  }
  //fprintf(stderr,"[%s] i:%lu\n",__FUNCTION__,i);
  ret->x=i;
  if(SIG_COND==0)
    exit(0);
  //  matrix_print<T>(ret);
  //  exit(0);
  fprintf(stderr,"\r");
}


//returns the number of sites read
template<typename T>
int readGLS(std::vector<persaf *> &adolf,size_t nSites,std::vector< Matrix<T> *> &ret){
  int pre=ret[0]->x;
  for(int i=0;i<adolf.size();i++){
    readGL(adolf[i],nSites,adolf[i]->nChr+1,ret[i]);
    //    fprintf(stderr,"adolf:%d\t%lu\n",i,ret[i]->x);
  }
  
  return ret[0]->x-pre;
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
double lik1(double *sfs,std::vector< Matrix<T> *> &gls,int from,int to){
  //  fprintf(stderr,"[%s] from:%d to:%d\n",__FUNCTION__,from,to);
  double r =0;
  for(int s=from;s<to;s++){
    double tmp =0;
    for(int i=0;i<gls[0]->y;i++)
      tmp += sfs[i]* gls[0]->mat[s][i];
    r += log(tmp);
  }
  return r;
}
template <typename T>
double lik2(double *sfs,std::vector< Matrix<T> *> &gls,int from,int to){
  double r =0;
  for(int s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(int i=0;i<gls[0]->y;i++)
      for(int j=0;j<gls[1]->y;j++)
	tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j];
    r += log(tmp);
  }
  return r;
}
template <typename T>
double lik3(double *sfs,std::vector< Matrix<T> *> &gls,int from,int to){
  double r =0;
  for(int s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(int i=0;i<gls[0]->y;i++)
      for(int j=0;j<gls[1]->y;j++)
	for(int k=0;k<gls[2]->y;k++)
	tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j]*gls[2]->mat[s][k];
    r += log(tmp);
  }
  return r;
}
template <typename T>
double lik4(double *sfs,std::vector< Matrix<T> *> &gls,int from,int to){
  double r =0;
  for(int s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(int i=0;i<gls[0]->y;i++)
      for(int j=0;j<gls[1]->y;j++)
	for(int k=0;k<gls[2]->y;k++)
	  for(int m=0;m<gls[3]->y;m++)
	tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j]*gls[2]->mat[s][k]*gls[3]->mat[s][m];
    r += log(tmp);
  }
  return r;
}

template <typename T>
void *like_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];
  if(pars.gls.size()==1)
    pars.lik = lik1(pars.sfs,pars.gls,pars.from,pars.to);
  else if(pars.gls.size()==2)
    pars.lik = lik2(pars.sfs,pars.gls,pars.from,pars.to);
  else if(pars.gls.size()==3)
    pars.lik = lik3(pars.sfs,pars.gls,pars.from,pars.to);
  else if(pars.gls.size()==4)
    pars.lik = lik4(pars.sfs,pars.gls,pars.from,pars.to);
  
  pthread_exit(NULL);
}


template <typename T>
double like_master(int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,like_slave<T>,(void*) i);
    if(rc)
      fprintf(stderr,"Error creating thread\n");
    
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
void emStep1(double *pre,std::vector< Matrix<T> * > &gls,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    for(int x=0;x<dim;x++)
      inner[x] = pre[x]*gls[0]->mat[s][x];
  
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}


template <typename T>
void emStep2(double *pre,std::vector<Matrix<T> *> &gls,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(int x=0;x<gls[0]->y;x++)
      for(int y=0;y<gls[1]->y;y++){
	inner[inc] = pre[inc]*gls[0]->mat[s][x]*gls[1]->mat[s][y];
	inc++;
      }
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}

template <typename T>
void emStep3(double *pre,std::vector<Matrix<T> *> &gls,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(int x=0;x<gls[0]->y;x++)
      for(int y=0;y<gls[1]->y;y++)
	for(int i=0;i<gls[2]->y;i++){
	  inner[inc] = pre[inc]*gls[0]->mat[s][x] * gls[1]->mat[s][y] * gls[2]->mat[s][i];
	  inc++;
	}
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
   
}

template <typename T>
void emStep4(double *pre,std::vector<Matrix<T> *> &gls,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(int x=0;x<gls[0]->y;x++)
      for(int y=0;y<gls[1]->y;y++)
	for(int i=0;i<gls[2]->y;i++)
	  for(int j=0;j<gls[3]->y;j++){
	    inner[inc] = pre[inc]*gls[0]->mat[s][x] * gls[1]->mat[s][y] * gls[2]->mat[s][i]* gls[3]->mat[s][j];
	    inc++;
	  }

  }
  normalize(inner,dim);
  for(int x=0;x<dim;x++)
    post[x] += inner[x];
  
  normalize(post,dim);
   
}

template <typename T>
void *emStep_slave(void *p){
  emPars<T> &pars = emp[(size_t) p];
  if(pars.gls.size()==1)
    emStep1<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim);
  else if(pars.gls.size()==2)
    emStep2<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim);
  else if(pars.gls.size()==3)
    emStep3<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim);
  else if(pars.gls.size()==4)
    emStep4<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim);
  pthread_exit(NULL);
}


template<typename T>
void emStep_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep_slave<T>,(void*) i);
    if(rc)
      fprintf(stderr,"Error creating thread\n");
    
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
double em(double *sfs,double tole,int maxIter,int nThreads,int dim){
  double oldLik,lik;
  oldLik = like_master<T>(nThreads);
  
  fprintf(stderr,"startlik=%f\n",oldLik);fflush(stderr);
  double tmp[dim];
  
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    emStep_master<T>(tmp,nThreads);
    
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];

    lik = like_master<T>(nThreads);

    fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return oldLik;
}

template<typename T>
emPars<T> *setThreadPars(std::vector<Matrix<T> * > &gls,double *sfs,int nThreads,int dim,int nSites){
  //  fprintf(stderr,"nSites:%d\n",nSites);
  emPars<T> *temp = new emPars<T>[nThreads];
  int blockSize = nSites/nThreads;
  for(int i=0;i<nThreads;i++){
    temp[i].threadId = i;
    temp[i].gls=gls;
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
  temp[nThreads-1].to=nSites;
#if 0
  fprintf(stderr,"--------------\n");
  for(int i=0;i<nThreads;i++)
    fprintf(stderr,"threadinfo %d)=(%d,%d)=%d \n",temp[i].threadId,temp[i].from,temp[i].to,temp[i].to-temp[i].from); //
  fprintf(stderr,"--------------\n");
#endif 
  
  thd= new pthread_t[nThreads];
  return temp;
}

template<typename T>
void destroy(emPars<T> *a,int nThreads ){
  for(int i=0;i<nThreads;i++)
    delete [] a[i].post;
  delete [] a;
  delete [] thd;
}

int really_kill =3;
int VERBOSE = 1;
void handler(int s) {
  if(s==13)//this is sigpipe
    exit(0);
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


size_t parspace(std::vector<persaf *> &saf){
  size_t ndim = 1;
  for(int i=0;i<saf.size();i++)
    ndim *= saf[i]->nChr+1;
  fprintf(stderr,"\t-> Dimension of parameter space: %lu\n",ndim);
  return ndim;
}

//unthreaded
//this will populate the keep vector by
// 1) set the chooseChr and populate toKeep
// 2) find over lap between different positions
// this is run once for each chromsome
int set_intersect_pos(std::vector<persaf *> &saf,char *chooseChr,int start,int stop){
  //  fprintf(stderr,"[%s] chooseChr:%s, start:%d stop:%d\n",__FUNCTION__,chooseChr,start,stop );

  if(saf.size()==1&&chooseChr==NULL){//use entire genome, then don't do any strange filtering
    //fprintf(stderr,"herer\n");
    return 0 ;
  }
  /*
    What happens here? very good question

   */
  static int firstTime =1;
  static myMap::iterator it_outer=saf[0]->mm.begin();
  if(chooseChr==NULL){
    if(it_outer==saf[0]->mm.end())//if we are done
      return -2;
    else if(firstTime==0){
      it_outer++;
      if(it_outer==saf[0]->mm.end()){
	//	fprintf(stderr,"done reading will exit \n");
	return -3;
      }
    }else
      firstTime =0;
    chooseChr = it_outer->first;
    fprintf(stderr,"\t-> Is in multi sfs, will now read data from chr:%s\n",chooseChr);
  }
  
  fprintf(stderr,"\t-> hello Im the master merge part of realSFS. and I'll now do a tripple bypass to find intersect \n");
  fprintf(stderr,"\t-> 1) Will set iter according to chooseChr and start and stop\n");
  assert(chooseChr!=NULL);

 //hit will contain the depth for the different populations
  keep<char> *hit =NULL;
  //  assert(saf.size()>1);

  if(saf.size()>1)
    hit =keep_alloc<char>();//
  
  
  //this loop will populate a 'hit' array containing the effective (differnt pops) depth
  //if we only have one population, then just return after iter_init
  for(int i=0;i<saf.size();i++){
    myMap::iterator it = iter_init(saf[i],chooseChr,start,stop);
    assert(it!=saf[i]->mm.end());  
    if(saf.size()==1)
      return 0;
    
    bgzf_seek(saf[i]->pos,it->second.pos,SEEK_SET);
    saf[i]->ppos = new int[it->second.nSites];
    bgzf_read(saf[i]->pos,saf[i]->ppos,it->second.nSites*sizeof(int));
    if(saf[i]->ppos[it->second.nSites-1] > hit->m)
      realloc(hit,saf[i]->ppos[it->second.nSites-1]+1);
    assert(hit->m>0);
    for(int j=saf[i]->toKeep->first;j<saf[i]->toKeep->last;j++)
      if(saf[i]->toKeep->d[j])
	hit->d[j]++;
  }
  
  // keep_info(hit,stderr,0,saf.size());
  
  //hit now contains the genomic position (that is the index).
  
  //let us now modify the the persaf toKeep char vector
  int tsk[saf.size()];
  for(int i=0;i<saf.size();i++){
    tsk[i] =0;
    for(int j=0;j<saf[i]->toKeep->last;j++)
      if(hit->d[j]!=saf.size())
	saf[i]->toKeep->d[j] =0;
      else
	tsk[i]++;
    fprintf(stderr,"\t-> Sites to keep from pop%d:\t%d\n",i,tsk[i]);
    if(i>0)
      assert(tsk[i]==tsk[i-1]);
#if 0
    keep_info(saf[i]->toKeep,stderr,0,1);
    //print out overlapping posiitons for all pops
    //    fprintf(stderr,"saf:%d\thit->m:%lu\td->m:%luppo")
    for(int j=0;0&&j<saf[i]->toKeep->m;j++){
      if(hit->d[saf[i]->ppos[j]]==saf.size())
	fprintf(stdout,"saf%d\t%d\n",i,saf[i]->ppos[j]);
    }
#endif
  }
  keep_destroy(hit);
}
/*
  return value 
  -3 indicates that we are doing multi sfs and that we are totally and should flush

 */



template <typename T>
int readdata(std::vector<persaf *> &saf,std::vector<Matrix<T> *> &gls,int nSites,char *chooseChr,int start,int stop){

  static int lastread=0;
  //  fprintf(stderr,"[%s] nSites:%d lastread:%d\n",__FUNCTION__,nSites,lastread);
  if(lastread==0 ){
    //    fprintf(stderr,"\t-> Done reading data from chromosome will prepare next chromosome\n");
    int ret = set_intersect_pos(saf,chooseChr,start,stop); 
    //fprintf(stderr,"ret:%d\n",ret);
    if(ret==-3)
      return -3;
  }
  lastread=readGLS(saf,nSites,gls);
  if(lastread==0)
    fprintf(stderr,"\t-> Only read nSites: %lu will therefore prepare next chromosome (or exit)\n",gls[0]->x);
  //fprintf(stderr,"readdata lastread:%d\n\n",lastread);
  // exit(0);
  if(chooseChr!=NULL&&lastread==0){
    //fprintf(stderr,"return -2\n");
    return -2;
  }
  else if(chooseChr==NULL &&lastread==0 )
    return -2;
  else
    return 1;
}


template <typename T>
int main_opt(args *arg){



  std::vector<persaf *> &saf =arg->saf;
  int nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=nsites(saf,arg);
  }
  if(fsizes<T>(saf,nSites)>getTotalSystemMemory())
    fprintf(stderr,"\t-> Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    
  fprintf(stderr,"\t-> nSites: %d\n",nSites);
  float bytes_req_megs = fsizes<T>(saf,nSites)/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"\t-> The choice of -nSites will require atleast: %f megabyte memory, that is atleast: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim= parspace(saf);
  double *sfs=new double[ndim];
  
  
  while(1) {
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop);//read nsites from data
    //    fprintf(stderr,"\t\tRET:%d\n",ret);
    if(ret==-2&gls[0]->x==0)//no more data in files or in chr, eith way we break;
      break;
    
    if(saf.size()==1){
      if(ret!=-2){
	if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	  //	  fprintf(stderr,"continue continue\n");
	  continue;
	}
      }
    }else{
      if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	//fprintf(stderr,"continue continue\n");
	continue;
      }

    }
  
      
    fprintf(stderr,"\t-> Will run optimization on nSites: %lu\n",gls[0]->x);
    
    if(arg->sfsfname!=NULL)
      readSFS(arg->sfsfname,ndim,sfs);
    else
      for(int i=0;i<ndim;i++)
	sfs[i] = (i+1)/((double)(ndim));

    normalize(sfs,ndim);
    emp = setThreadPars<T>(gls,sfs,arg->nThreads,ndim,gls[0]->x);
    fprintf(stderr,"------------\n");
    double lik = em<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim);
    fprintf(stderr,"likelihood: %f\n",lik);
    fprintf(stderr,"------------\n");
#if 1
    for(int x=0;x<ndim;x++)
      fprintf(stdout,"%f ",log(sfs[x]));
    fprintf(stdout,"\n");
    fflush(stdout);
#endif
    destroy<T>(emp,arg->nThreads);
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
  }

  destroy(gls,nSites);
  destroy_args(arg);
  delete [] sfs;
  

  return 0;
}

int main(int argc,char **argv){
#if 0
  char **reg = new char*[6];
  reg[0] = strdup("avsdf:");
  reg[1] = strdup("avsdf");
  reg[2] = strdup("avsdf:200");
  reg[3] = strdup("avsdf:200-300");
  reg[4] = strdup("avsdf:-300");
  reg[5] = strdup("avsdf:300-");

  int a,b;
  char *ref;
  for(int i=0;i<6;i++){
  fprintf(stderr,"%d) string=%s ",i,reg[i]);
  get_region(reg[i],&ref,a,b);
  fprintf(stderr," parsed as: ref:\'%s\' a:\'%d\' b:\'%d\'\n",ref,a,b);
}
  return 0;
#endif


  
  //start of signal handling
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  
\
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
  else if(!strcasecmp(*argv,"print2"))
    print2(--argc,++argv);
  else {
    args *arg = getArgs(argc,argv);
    if(!arg)
      return 0;
    main_opt<float>(arg);
    
  }

  return 0;
}
