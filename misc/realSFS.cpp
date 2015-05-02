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

void destroy(std::vector<persaf *> &saf){
  for(int i=0;i<saf.size();i++)
    persaf_destroy(saf[i]);
}


void destroy(args *p){
  destroy(p->saf);
  delete p;
}


size_t fsizes(std::vector<persaf *> &pp){
  size_t res = 0;
  for(int i=0;i<pp.size();i++)
    res += pp[i]->fsize;
  return res;
}

size_t nsites(std::vector<persaf *> &pp){
  size_t res = pp[0]->nSites;
  for(int i=1;i<pp.size();i++)
    if(pp[i]->nSites > res)
      res = pp[i]->nSites;
  return res;
}


char * get_region(char *extra,int &start,int &stop) {

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
    }
    else  if(!strcasecmp(*argv,"-start")){
      p->sfsfname = *(++argv);
    }else{
      p->saf.push_back(readsaf<float>(*argv));
      //   fprintf(stderr,"toKeep:%p\n",p->saf[p->saf.size()-1]->toKeep);
    }
    argv++;
  }
  fprintf(stderr,"args: tole:%f nthreads:%d maxiter:%d nsites:%d chooseChr:%s start:%s chr:%s start:%d stop:%d\n",p->tole,p->nThreads,p->maxIter,p->nSites,p->chooseChr,p->sfsfname,p->chooseChr,p->start,p->stop);
  return p;
}

void print(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [chrname]\n");
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
  destroy(pars);
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
    
    if(pars->chooseChr!=NULL){
      iter_init(pars->saf[0],pars->chooseChr,pars->start,pars->stop);
      it = pars->saf[0]->mm.find(pars->chooseChr);
    }
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
  destroy(pars);
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
Matrix<T> *alloc(size_t x,size_t y){
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
double lik1(double *sfs,std::vector< Matrix<T> *> &gls,int from,int to){
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
  
  fprintf(stderr,"startlik=%f\n",oldLik);
  fflush(stderr);

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
  fprintf(stderr,"nSites:%d\n",nSites);
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
#if 1
  for(int i=0;i<nThreads;i++)
    fprintf(stderr,"%d:(%d,%d)=%d ",temp[i].threadId,temp[i].from,temp[i].to,temp[i].to-temp[i].from); //
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

template <typename T>
void readdata(std::vector<persaf *> &saf,std::vector<Matrix<T> *> &gls,int nSites,char *chooseChr,int start,int stop){
  assert(saf.size()==1);
  
  if(chooseChr!=NULL){
    for(int f=0;f<saf.size();f++){
      myMap::iterator it = saf[f]->mm.find(chooseChr);
      if(it==saf[f]->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",chooseChr);
	break;
      }
      bgzf_seek(saf[f]->pos,it->second.pos,SEEK_SET);
      bgzf_seek(saf[f]->saf,it->second.saf,SEEK_SET);
    }
  }
  
  readGL(saf[0]->saf,nSites,saf[0]->nChr+1,gls[0]);
}


template <typename T>
int main_opt(args *arg){



  std::vector<persaf *> &saf =arg->saf;
  int nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    if(fsizes(saf)>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    //this doesnt make sense if ppl supply a filelist containing safs
    nSites=nsites(saf);
  }
  //  float bytes_req_megs = saf->fsize/1024/1024;
  //float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  //fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim= parspace(saf);
  double *sfs=new double[ndim];
  emp = setThreadPars<T>(gls,sfs,arg->nThreads,ndim,nSites);
  
  while(1) {
    readdata<T>(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop);//read nsites from data
    
    if(gls[0]->x==0)
      break;
    fprintf(stderr,"dim(GL1)=%zu,%zu\n",gls[0]->x,gls[0]->y);
     
    if(arg->sfsfname!=NULL)
      readSFS(arg->sfsfname,ndim,sfs);
    else
      for(int i=0;i<ndim;i++)
	sfs[i] = (i+1)/((double)(ndim));

    normalize(sfs,ndim);

    double lik = em<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    for(int x=0;x<ndim;x++)
      fprintf(stdout,"%f ",log(sfs[x]));
    fprintf(stdout,"\n");
    fflush(stdout);
#endif
    
  }
  destroy<T>(emp,arg->nThreads);
  destroy(gls,nSites);
  destroy(saf);

  destroy(arg);
  delete [] sfs;
  delete [] thd;

  return 0;
}

//unthreaded
keep<char> *merge(std::vector<persaf *> &saf,char *chooseChr){
  fprintf(stderr,"hello Im the master merge part of realSFS. and I'll now do a tripple bypass\n");
  assert(chooseChr!=NULL);
  
  //  static keep<int> *dat = alloc_keep<int>();//positions 
  static keep<char> *hit =alloc_keep<char>();//
   
  for(int i=0;i<saf.size();i++){
    myMap::iterator it = saf[i]->mm.find(chooseChr);
    assert(it!=saf[i]->mm.end());

    bgzf_seek(saf[i]->pos,it->second.pos,SEEK_SET);
    saf[i]->ppos = new int[it->second.nSites];
    bgzf_read(saf[i]->pos,saf[i]->ppos,it->second.nSites*sizeof(int));
    if(saf[i]->ppos[it->second.nSites-1] > hit->m)
      realloc(hit,saf[i]->ppos[it->second.nSites-1]+1);
   
    assert(hit->m>0);

    for(int j=0;j<it->second.nSites;j++)
      hit->d[saf[i]->ppos[j]]++;
  }
  
  
  return hit;
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
    merge(arg->saf,"18");
    return 0;
    if(arg->saf.size()==1)
      main_opt<float>(arg);
    
  }

  return 0;
}
