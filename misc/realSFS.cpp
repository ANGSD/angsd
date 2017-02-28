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
#include <zlib.h>
#include <htslib/tbx.h>
#include "Matrix.hpp"
#include "safstat.h"
#include <libgen.h>
#include <algorithm>
#include "realSFS_args.h"
#include "safreader.h"
#include "keep.hpp"
#include "header.h"
#include "safcat.h"

int SIG_COND =1;
int howOften =5e6;//how often should we print out (just to make sure something is happening)
#include "multisafreader.hpp"

double ttol = 1e-16; 
pthread_t *thd=NULL;

int **posiG  = NULL;
size_t *bootstrap = NULL;
int really_kill =3;
int VERBOSE = 1;

extern std::vector <char *> dumpedFiles;

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

/*
  over elaborate function to read a sfs. Assumption is that input file contains the expected values.
  output is plugged into ret
 */

void readSFS(const char*fname,size_t hint,double *ret){
  fprintf(stderr,"\t-> Reading: %s assuming counts (will normalize to probs internally)\n",fname);
  FILE *fp = NULL;
  if(((fp=fopen(fname,"r")))==NULL){
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
    fprintf(stderr,"problem with size of dimension of prior %lu vs %lu\n",hint,res.size());
    for(size_t i=0;0&&i<res.size();i++)
      fprintf(stderr,"%zu=%f\n",i,res[i]);
    exit(0);
  }
  for(size_t i=0;i<res.size();i++){
      ret[i] = res[i];
      //      fprintf(stderr,"i=%lu %f\n",i,ret[i]);
  }
  normalize(ret,(int)res.size());
  for(int i=0;0&&i<res.size();i++)
    ret[i] = log(ret[i]);
  fclose(fp);
}


size_t parspace(std::vector<persaf *> &saf){
  size_t ndim = 1;
  for(int i=0;i<saf.size();i++){
    ndim *= saf[i]->nChr+1;
    fprintf(stderr,"\t-> dim(%s):%lu\n",saf[i]->fname,saf[i]->nChr+1);
  }
  fprintf(stderr,"\t-> Dimension of parameter space: %lu\n",ndim);
  return ndim;
}

//just approximate
template <typename T>
size_t fsizes(std::vector<persaf *> &pp, size_t nSites){
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
/*
  returns the maxnumber of sites across all samples
 */
size_t nsites(std::vector<persaf *> &pp,args *ar){
  if(ar->start!=-1 &&ar->stop!=-1)
    return ar->stop-ar->start;
  size_t res = helper(pp[0],ar->chooseChr);
  for(int i=1;i<pp.size();i++)
    if(helper(pp[i],ar->chooseChr) > res)
      res = helper(pp[i],ar->chooseChr);
  return res;
}



int print_header(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return 0; 
  }
  
  args *pars = getArgs(argc,argv);
  if(!pars)
    return 1;
  if(pars->saf.size()!=1){
    fprintf(stderr,"print_header only implemeted for single safs\n");
    exit(0);
  }
  writesaf_header(stdout,pars->saf[0]);
  
  destroy_args(pars);
  return 0;
}

void setGloc(std::vector<persaf *> &saf,size_t nSites){
#if 1
  posiG = new int*[saf.size()];
  for(int i=0;i<saf.size();i++)
    posiG[i] = new int[nSites];
#endif
}

void delGloc(std::vector<persaf *> &saf,size_t nSites){
  for(int i=0;i<saf.size();i++)
    delete [] posiG[i];
  delete [] posiG;
}


template <typename T>
int printMulti(args *arg){
  fprintf(stderr,"[%s]\n",__FUNCTION__);
  std::vector<persaf *> &saf =arg->saf;
  for(int i=0;i<saf.size();i++)
    assert(saf[i]->pos!=NULL&&saf[i]->saf!=NULL);

  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=nsites(saf,arg);
  }
  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim=(int) parspace(saf);
  double *sfs=new double[ndim];
  
  //temp used for checking pos are in sync
  setGloc(saf,nSites);
  int *posiToPrint = new int[nSites];
  //used for printout old format
  FILE **oldfp = NULL;
  gzFile oldpos = Z_NULL;
  if(arg->oldout){
    oldfp = new FILE*[saf.size()];
    for(int i=0;i<saf.size();i++){
      size_t newlen = strlen(saf[i]->fname)+100;
      char *tmp =(char*) calloc(newlen,sizeof(char));
      tmp = strncpy(tmp,saf[i]->fname,strlen(saf[i]->fname)-4);
      fprintf(stderr,"\t-> Generating outputfile: %s\n",tmp);
      oldfp[i] = fopen(tmp,"wb");
      free(tmp);
    }
    size_t newlen = strlen(dirname(saf[0]->fname))+100;
    char *tmp = (char*) calloc(newlen,sizeof(char));
    snprintf(tmp,newlen,"%s/shared.pos.gz",dirname(saf[0]->fname));
    fprintf(stderr,"\t-> Generating outputfile: %s\n",tmp);
    oldpos = gzopen(tmp,"wb");
    free(tmp);
  }
  while(1) {
    static char *curChr=NULL;
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop,posiToPrint,&curChr,arg->fl,1);//read nsites from data
    if(arg->oldout==0){
      for(int s=0;s<gls[0]->x;s++){
	if(arg->chooseChr==NULL)
	  fprintf(stdout,"%s\t%d",curChr,posiToPrint[s]+1);
	else
	  fprintf(stdout,"%s\t%d",arg->chooseChr,posiToPrint[s]+1);
	for(int i=0;i<saf.size();i++)
	  for(int ii=0;ii<gls[i]->y;ii++)
	    fprintf(stdout,"\t%f",log(gls[i]->mat[s][ii]));
	fprintf(stdout,"\n");
      }
    }else{
      for(int s=0;s<gls[0]->x;s++){
	if(arg->chooseChr==NULL)
	  gzprintf(oldpos,"%s\t%d\n",curChr,posiToPrint[s]+1);
	else
	  gzprintf(oldpos,"%s\t%d\n",arg->chooseChr,posiToPrint[s]+1);
	for(int i=0;i<saf.size();i++){
	  double mytmp[gls[i]->y];
	  for(int ii=0;ii<gls[i]->y;ii++)
	    mytmp[ii] = log(gls[i]->mat[s][ii]);
	  fwrite(mytmp,sizeof(double),gls[i]->y,oldfp[i]);
	}
      }
    }
    if(ret==-3&&gls[0]->x==0){//no more data in files or in chr, eith way we break;g
      //fprintf(stderr,"breaking\n");
      break;
    }
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
#if 0
    if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
      fprintf(stderr,"continue continue\n");
      continue;
    }
#endif
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
    if(arg->onlyOnce)
      break;
  }
  delGloc(saf,nSites);
  destroy(gls,nSites);

  delete [] sfs;
  delete [] posiToPrint;

  if(arg->oldout==1){
    for(int i=0;i<saf.size();i++)
      fclose(oldfp[i]);
    delete [] oldfp;
    gzclose(oldpos);
  }
  destroy_args(arg);
  fprintf(stderr,"\t-> Run completed\n");
  return 0;
}


template<typename T>
void print(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"\t-> Must supply afile.saf.idx files \n");
    fprintf(stderr,"\t-> Examples \n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx \n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx -r chr1:10000000-12000000\n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx pop2.saf.idx -r chr2:10000000-12000000\n");
    fprintf(stderr,"\t-> You can generate the oldformat by appending the -oldout 1 to the print command like\n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx pop2.saf.idx -oldout 1\n");
    return; 
  }
  
  args *pars = getArgs(argc,argv);
  for(int i=0;i<pars->saf.size();i++)
    pars->saf[0]->kind = 2;
  if(1||pars->saf.size()!=1){
    fprintf(stderr,"\t-> Will jump to multisaf printer and will only print intersecting sites between populations\n");
    printMulti<T>(pars);
    return;
  }

  writesaf_header(stderr,pars->saf[0]);
  
  float *flt = new float[pars->saf[0]->nChr+1];
  for(myMap::iterator it=pars->saf[0]->mm.begin();it!=pars->saf[0]->mm.end();++it){

    if(pars->chooseChr!=NULL)
      it = iter_init(pars->saf[0],pars->chooseChr,pars->start,pars->stop);
    else
      it = iter_init(pars->saf[0],it->first,pars->start,pars->stop);
 
    size_t ret;
    int pos;

    while((ret=iter_read(pars->saf[0],flt,sizeof(float)*(pars->saf[0]->nChr+1),&pos))){
      fprintf(stdout,"%s\t%d",it->first,pos+1);
      for(int is=0;is<pars->saf[0]->nChr+1;is++)
	fprintf(stdout,"\t%f",flt[is]);
      fprintf(stdout,"\n");
    }
 
    if(pars->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  destroy_args(pars);
}

template<typename T>
struct emPars{
  int threadId;
  std::vector <Matrix<T> *> gls;
  size_t from;
  size_t to;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  double *inner;//allocated for every thread
  int dim;
};



emPars<float> *emp = NULL;

template <typename T>
double lik1(double *sfs,std::vector< Matrix<T> *> &gls,size_t from,size_t to){
  assert(from >=0 && to >=0);
  //  fprintf(stderr,"[%s] from:%d to:%d\n",__FUNCTION__,from,to);
  double r =0;
  if(!bootstrap)
    for(size_t s=from;s<to;s++){
      double tmp =0;
      for(size_t i=0;i<gls[0]->y;i++)
	tmp += sfs[i]* gls[0]->mat[s][i]; 
      r += log(tmp);
    }
  else
    for(size_t s=from;s<to;s++){
      double tmp =0;
      for(size_t i=0;i<gls[0]->y;i++)
	tmp += sfs[i]* gls[0]->mat[bootstrap[s]][i]; 
      r += log(tmp);
    }
  return r;
}
template <typename T>
double lik2(double *sfs,std::vector< Matrix<T> *> &gls,size_t from,size_t to){
  double r =0;
  if(bootstrap==NULL){
    for(size_t s=from;s<to;s++){
      double tmp =0;
      int inc =0;
      for(size_t i=0;i<gls[0]->y;i++)
	for(size_t j=0;j<gls[1]->y;j++)
	  tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j];
      r += log(tmp);
    }
  }
  else
    for(size_t s=from;s<to;s++){
      double tmp =0;
      int inc =0;
      for(size_t i=0;i<gls[0]->y;i++)
	for(size_t j=0;j<gls[1]->y;j++)
	  tmp += sfs[inc++]* gls[0]->mat[bootstrap[s]][i] *gls[1]->mat[bootstrap[s]][j];
      r += log(tmp);
    }
  return r;
}

template <typename T>
double lik3(double *sfs,std::vector< Matrix<T> *> &gls,size_t from,size_t to){
  double r =0;
  if(bootstrap==NULL)
  for(size_t s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(size_t i=0;i<gls[0]->y;i++)
      for(size_t j=0;j<gls[1]->y;j++)
	for(size_t k=0;k<gls[2]->y;k++)
	tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j]*gls[2]->mat[s][k];
    r += log(tmp);
  }
  else
  for(size_t s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(size_t i=0;i<gls[0]->y;i++)
      for(size_t j=0;j<gls[1]->y;j++)
	for(size_t k=0;k<gls[2]->y;k++)
	tmp += sfs[inc++]* gls[0]->mat[bootstrap[s]][i] *gls[1]->mat[bootstrap[s]][j]*gls[2]->mat[bootstrap[s]][k];
    r += log(tmp);
  }
  return r;
}
template <typename T>
double lik4(double *sfs,std::vector< Matrix<T> *> &gls,size_t from,size_t to){
  double r =0;
  if(bootstrap==NULL)
   for(size_t s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(size_t i=0;i<gls[0]->y;i++)
      for(size_t j=0;j<gls[1]->y;j++)
	for(size_t k=0;k<gls[2]->y;k++)
	  for(size_t m=0;m<gls[3]->y;m++)
	tmp += sfs[inc++]* gls[0]->mat[s][i] *gls[1]->mat[s][j]*gls[2]->mat[s][k]*gls[3]->mat[s][m];
    r += log(tmp);
  }
  else
  for(size_t s=from;s<to;s++){
    double tmp =0;
    int inc =0;
    for(size_t i=0;i<gls[0]->y;i++)
      for(size_t j=0;j<gls[1]->y;j++)
	for(size_t k=0;k<gls[2]->y;k++)
	  for(size_t m=0;m<gls[3]->y;m++)
	tmp += sfs[inc++]* gls[0]->mat[bootstrap[s]][i] *gls[1]->mat[bootstrap[s]][j]*gls[2]->mat[bootstrap[s]][k]*gls[3]->mat[bootstrap[s]][m];
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
void emStep1(double *pre,std::vector< Matrix<T> * > &gls,double *post,size_t start,size_t stop,int dim,double *inner){
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(size_t s=start;SIG_COND&&s<stop;s++){
    if(bootstrap==NULL)
      for(int x=0;x<dim;x++)
	inner[x] = pre[x]*gls[0]->mat[s][x];
    else
      for(int x=0;x<dim;x++)
	inner[x] = pre[x]*gls[0]->mat[bootstrap[s]][x];
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}


template <typename T>
void emStep2(double *pre,std::vector<Matrix<T> *> &gls,double *post,size_t start,size_t stop,int dim,double *inner){
  for(int x=0;x<dim;x++)
    post[x] =0.0;
   
  if(bootstrap==NULL)
    for(size_t s=start;SIG_COND&&s<stop;s++){
      int inc=0;
      for(size_t x=0;x<gls[0]->y;x++)
	for(size_t y=0;y<gls[1]->y;y++){
	  inner[inc] = pre[inc]*gls[0]->mat[s][x]*gls[1]->mat[s][y];
	  inc++;
	}
      normalize(inner,dim);
      for(int x=0;x<dim;x++)
	post[x] += inner[x];
    }
  else
    for(size_t s=start;SIG_COND&&s<stop;s++){
      int inc=0;
      for(size_t x=0;x<gls[0]->y;x++)
	for(size_t y=0;y<gls[1]->y;y++){
	  inner[inc] = pre[inc]*gls[0]->mat[s][x]*gls[1]->mat[bootstrap[s]][y];
	  inc++;
	}
      normalize(inner,dim);
      for(int x=0;x<dim;x++)
	post[x] += inner[x];
    }

  normalize(post,dim);
 
}

template <typename T>
void emStep3(double *pre,std::vector<Matrix<T> *> &gls,double *post,size_t start,size_t stop,int dim,double *inner){
  for(int x=0;x<dim;x++)
    post[x] =0.0;
  if(bootstrap==NULL)
  for(size_t s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(size_t x=0;x<gls[0]->y;x++)
      for(size_t y=0;y<gls[1]->y;y++)
	for(size_t i=0;i<gls[2]->y;i++){
	  inner[inc] = pre[inc]*gls[0]->mat[s][x] * gls[1]->mat[s][y] * gls[2]->mat[s][i];
	  inc++;
	}
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  else
  for(size_t s=start;SIG_COND&&s<stop;s++){
    int inc=0;
    for(size_t x=0;x<gls[0]->y;x++)
      for(size_t y=0;y<gls[1]->y;y++)
	for(size_t i=0;i<gls[2]->y;i++){
	  inner[inc] = pre[inc]*gls[0]->mat[bootstrap[s]][x] * gls[1]->mat[bootstrap[s]][y] * gls[2]->mat[bootstrap[s]][i];
	  inc++;
	}
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
   
}

template <typename T>
void emStep4(double *pre,std::vector<Matrix<T> *> &gls,double *post,size_t start,size_t stop,int dim,double *inner){

  for(int x=0;x<dim;x++)
    post[x] =0.0;
  if(bootstrap==NULL){
    for(size_t s=start;SIG_COND&&s<stop;s++){
      int inc=0;
      for(size_t x=0;x<gls[0]->y;x++)
	for(size_t y=0;y<gls[1]->y;y++)
	  for(size_t i=0;i<gls[2]->y;i++)
	    for(size_t j=0;j<gls[3]->y;j++){
	      inner[inc] = pre[inc]*gls[0]->mat[s][x] * gls[1]->mat[s][y] * gls[2]->mat[s][i]* gls[3]->mat[s][j];
	      inc++;
	    }
    }
  }
  else{
    for(size_t s=start;SIG_COND&&s<stop;s++){
      int inc=0;
      for(size_t x=0;x<gls[0]->y;x++)
	for(size_t y=0;y<gls[1]->y;y++)
	  for(size_t i=0;i<gls[2]->y;i++)
	    for(size_t j=0;j<gls[3]->y;j++){
	      inner[inc] = pre[inc]*gls[0]->mat[bootstrap[s]][x] * gls[1]->mat[bootstrap[s]][y] * gls[2]->mat[bootstrap[s]][i]* gls[3]->mat[bootstrap[s]][j];
	      inc++;
	    }
      
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
    emStep1<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim,pars.inner);
  else if(pars.gls.size()==2)
    emStep2<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim,pars.inner);
  else if(pars.gls.size()==3)
    emStep3<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim,pars.inner);
  else if(pars.gls.size()==4)
    emStep4<T>(pars.sfs,pars.gls,pars.post,pars.from,pars.to,pars.dim,pars.inner);
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


template<typename T>
emPars<T> *setThreadPars(std::vector<Matrix<T> * > &gls,double *sfs,int nThreads,int dim,size_t nSites){
  //  fprintf(stderr,"nSites:%d\n",nSites);
  emPars<T> *temp = new emPars<T>[nThreads];
  size_t blockSize = nSites/nThreads;
  for(int i=0;i<nThreads;i++){
    temp[i].threadId = i;
    temp[i].gls=gls;
    temp[i].from =0;
    temp[i].to=blockSize;
    temp[i].sfs = sfs;
    temp[i].post=new double[dim];
    temp[i].inner=new double[dim];
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
    fprintf(stderr,"threadinfo %d)=(%lu,%lu)=%d \n",temp[i].threadId,temp[i].from,temp[i].to,temp[i].to-temp[i].from); //
  fprintf(stderr,"--------------\n");
#endif 
  
  thd= new pthread_t[nThreads];
  return temp;
}

template<typename T>
void destroy(emPars<T> *a,int nThreads ){
  for(int i=0;i<nThreads;i++){
    delete [] a[i].inner;
    delete [] a[i].post;
  }
  delete [] a;
  delete [] thd;
}



template <typename T>
double em(double *sfs,double tole,int maxIter,int nThreads,int dim,std::vector<Matrix<T> *> &gls){
  emp = setThreadPars<T>(gls,sfs,nThreads,dim,gls[0]->x);
  fprintf(stderr,"------------\n");
  double oldLik,lik;
  oldLik = like_master<T>(nThreads);
  
  fprintf(stderr,"startlik=%f\n",oldLik);fflush(stderr);

  double tmp[dim];

  for(int it=0;SIG_COND&&it<maxIter;it++) {
    emStep_master<T>(tmp,nThreads);
    double sr2 = 0;
    for(int i=0;i<dim;i++){
      sr2 += (sfs[i]-tmp[i])*(sfs[i]-tmp[i]);
      sfs[i]= tmp[i];
      
    }

    lik = like_master<T>(nThreads);

    fprintf(stderr,"[%d] lik=%f diff=%e sr:%e\n",it,lik,fabs(lik-oldLik),sr2);

    if((fabs(lik-oldLik)<tole)||(0&&(sqrt(sr2)<tole))){//should update simfiles...
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
destroy<T>(emp,nThreads);
  return oldLik;
}


template <typename T>
double emAccl(double *p,double tole,int maxIter,int nThreads,int dim,std::vector<Matrix<T> *> &gls,int useSq){
  emp = setThreadPars<T>(gls,p,nThreads,dim,gls[0]->x);
  fprintf(stderr,"------------\n");
  double oldLik,lik;
  oldLik = like_master<T>(nThreads);
  
  fprintf(stderr,"startlik=%f\n",oldLik);fflush(stderr);

  double *p1 = new double[dim];
  double *p2 = new double[dim];
  double *q1 = new double[dim];
  double *q2 = new double[dim];
  double *pnew = new double[dim];
  double *oldp = new double[dim];
  double *tmp = new double[dim];
  int stepMax = 1;
  int mstep = 4;
  int stepMin = 1;
  
  int iter =0;

  while(SIG_COND&&iter<maxIter){
    emStep_master<T>(p1,nThreads);
    iter++;
    double sr2 =0;
    for(size_t i=0;i<dim;i++){
      q1[i] = p1[i]-p[i];
      sr2 += q1[i]*q1[i];
    }

    memcpy(oldp,p,sizeof(double)*dim);
    memcpy(p,p1,sizeof(double)*dim);  
    if(sqrt(sr2)<tole){
      fprintf(stderr,"\t-> Breaking EM(sr2) at iter:%d, sqrt(sr2):%e\n",iter,sqrt(sr2));
      break;
    }
    emStep_master<T>(p2,nThreads);
    iter++;


    double sq2 =0;
    for(size_t i=0;i<dim;i++){
      q2[i] = p2[i]-p1[i];
      sq2 += q2[i]*q2[i];
    }

    if(sqrt(sq2)<tole){
      fprintf(stderr,"\t-> Breaking EM(sq2) at iter:%d, sqrt(sq2):%e\n",iter,sqrt(sq2));
      break;
    }
    double sv2=0;
    for(size_t i=0;i<dim;i++)
      sv2 += (q2[i]-q1[i])*(q2[i]-q1[i]);

    double alpha = sqrt(sr2/sv2);
    alpha =std::max(stepMin*1.0,std::min(1.0*stepMax,alpha));
    
    //the magical step
    for(size_t j=0;j<dim;j++)
      pnew[j] = oldp[j] + 2.0 * alpha * q1[j] + alpha*alpha * (q2[j] - q1[j]);
#if 1 //fix for going out of bound
    for(size_t j=0;j<dim;j++){
      if(pnew[j]<ttol)
	pnew[j]=ttol;
      if(pnew[j]>1-ttol)
	pnew[j]=1-ttol;
    }
    normalize(pnew,dim);
#endif
    if(fabs(alpha-1) >0.01){
      //this is clearly to stabilize
      //      double tmp[dim];
      memcpy(p,pnew,sizeof(double)*dim);
      emStep_master<T>(tmp,nThreads);
      memcpy(pnew,tmp,sizeof(double)*dim);
      iter++;
    }
    memcpy(p,pnew,sizeof(double)*dim);
    
    //like_master is using sfs[] to calculate like
    lik = like_master<T>(nThreads);
    fprintf(stderr,"lik[%d]=%f diff=%e alpha:%f sr2:%e\n",iter,lik,fabs(lik-oldLik),alpha,sr2);
    if(std::isnan(lik)) {
      fprintf(stderr,"\t-> Observed NaN in accelerated EM, will use last reliable value. Consider using as input for ordinary em\n");
      fprintf(stderr,"\t-> E.g ./realSFS -sfs current.output -m 0 >new.output\n");//thanks morten rasmussen
      memcpy(p,p2,sizeof(double)*dim);
      break;
    }
    
    if(0&&lik<oldLik)//this should at some point be investigated further //
      fprintf(stderr,"\t-> New like is worse?\n");
#if 1
    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
#endif
    if (alpha == stepMax) 
      stepMax = mstep * stepMax;
    if(stepMin<0 &&alpha==stepMin)
      stepMin = mstep*stepMin;

  }
  destroy<T>(emp,nThreads);
  delete [] p1;
  delete [] p2;
  delete [] q1;
  delete [] q2;
  delete [] pnew;
  delete [] oldp;
  delete [] tmp;
  return oldLik;
}



/*
  return value 
  -3 indicates that we are doing multi sfs and that we are totally and should flush

 */

int fst_index(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [chrname, write more info]\n");
    return 0; 
  }
  args *arg = getArgs(argc,argv);
  if(!arg->outname){
    fprintf(stderr,"\t-> Must supply -fstout for doing fstindex\n");
    return 0;
  }

  std::vector<persaf *> &saf =arg->saf;
  //assert(saf.size()==2);
  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites = 100000;//<- set default to 100k sites, no need to load everything...
    // nSites=nsites(saf,arg);
  }
  fprintf(stderr,"\t-> nSites: %lu\n",nSites);
  std::vector<Matrix<float> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<float>(nSites,saf[i]->nChr+1));

  //  int ndim= parspace(saf);
  if(arg->sfsfname.size()!=choose(saf.size(),2)){
    fprintf(stderr,"\t-> You have supplied: %lu populations, that is %d pairs\n",saf.size(),choose(saf.size(),2));
    fprintf(stderr,"\t-> You therefore need to supply %d 2dsfs priors instead of:%lu\n",choose(saf.size(),2),arg->sfsfname.size());
    exit(0);
  }
  std::vector<double *> sfs;
  int inc =0;
  for(int i=0;i<saf.size();i++)
    for(int j=i+1;j<saf.size();j++){
      size_t pairdim = (saf[i]->nChr+1)*(saf[j]->nChr+1);
      double *ddd=new double[pairdim];
      readSFS(arg->sfsfname[inc],pairdim,ddd);
      normalize(ddd,pairdim);
      sfs.push_back(ddd);
      inc++;
    }

  
  double **a1,**b1;
  a1=new double*[choose(saf.size(),2)];
  b1=new double*[choose(saf.size(),2)];
  inc=0;
  for(int i=0;i<saf.size();i++)
    for(int j=i+1;j<saf.size();j++){
      calcCoef((int)saf[i]->nChr,(int)saf[j]->nChr,&a1[inc],&b1[inc],arg->whichFst);
      //      fprintf(stderr,"a1[%d]:%p b1[%d]:%p\n",inc,&a1[inc][0],inc,&b1[inc][0]);
      inc++;
    }

  BGZF *fstbg = openFileBG(arg->outname,".fst.gz");
  FILE *fstfp = openFile(arg->outname,".fst.idx");
  char buf[8]="fstv1";
  if(bgzf_write(fstbg,buf,8)!=8){
    fprintf(stderr,"\t-> Problem writing 8bytes into output file\n");
    exit(0);
  }
  fwrite(buf,1,8,fstfp);
#if 0
  for(int i=0;i<ndim;i++)
    fprintf(stdout,"%f %f\n",a1[i],b1[i]);
  exit(0);
#endif
#if 1
  size_t nsafs=saf.size();
  fwrite(&nsafs,sizeof(size_t),1,fstfp);
  for(int i=0;i<nsafs;i++){
    size_t clen= strlen(saf[i]->fname);
    fwrite(&clen,sizeof(size_t),1,fstfp);
    fwrite(saf[i]->fname,1,clen,fstfp);
  }
#endif
  int asdf = choose(saf.size(),2);
  std::vector<double> *ares = new std::vector<double> [choose(saf.size(),2)];
  std::vector<double> *bres = new std::vector<double> [choose(saf.size(),2)];
  //  for(int i=0;i<3;i++)
    //    fprintf(stderr,"ares.size():%lu bres.size():%lu sfs:%p\n",ares[i].size(),bres[i].size(),&sfs[i][0]);
  std::vector<int> posi;
  setGloc(saf,nSites);
  int *posiToPrint = new int[nSites];
  for(myMap::iterator it = saf[0]->mm.begin();it!=saf[0]->mm.end();++it) {
    //    fprintf(stderr,"doing chr:%s\n",it->first);
    if(arg->chooseChr!=NULL){
      it = saf[0]->mm.find(arg->chooseChr);
      if(it==saf[0]->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",arg->chooseChr);
	break;
      }
    }
    for(int i=0;i<choose(saf.size(),2);i++){
      ares[i].clear();
      bres[i].clear();
    }
    posi.clear();
    while(1) {
      int ret=readdata(saf,gls,nSites,it->first,arg->start,arg->stop,posiToPrint,NULL,arg->fl,1);//read nsites from data
      //  fprintf(stderr,"ret:%d glsx:%lu\n",ret,gls[0]->x);
      //if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	//fprintf(stderr,"continue continue\n");
      //	continue;
      //}
      
      fprintf(stderr,"\t-> Will now do fst temp dump using a chunk of %lu\n",gls[0]->x);
      int inc=0;
      for(int i=0;i<saf.size();i++)
	for(int j=i+1;j<saf.size();j++){
	  //	  fprintf(stderr,"i:%d j:%d inc:%d gls[i]:%p gls[j]:%p sfs:%p a1:%p b1:%p\n",i,j,inc,gls[i],gls[j],sfs[i],&a1[inc][0],&a1[inc][0]);
	  block_coef(gls[i],gls[j],sfs[inc],a1[inc],b1[inc],ares[inc],bres[inc]);
	  inc++;
	}
      for(int i=0;i<gls[0]->x;i++)
	posi.push_back(posiToPrint[i]);

      for(int i=0;i<gls.size();i++)
	gls[i]->x =0;
      if(ret==-2)//no more data in files or in chr, eith way we break;
	break;
    }
    size_t clen = strlen(it->first);
    fwrite(&clen,sizeof(size_t),1,fstfp);
    fwrite(it->first,1,clen,fstfp);
    size_t nit=posi.size();

    assert(1==fwrite(&nit,sizeof(size_t),1,fstfp));
    int64_t tell = bgzf_tell(fstbg);
    fwrite(&tell,sizeof(int64_t),1,fstfp);
    my_bgzf_write(fstbg,&posi[0],posi.size()*sizeof(int));
    int inc =0;
    for(int i=0;i<saf.size();i++)
      for(int j=i+1;j<saf.size();j++){
	my_bgzf_write(fstbg,&(ares[inc][0]),ares[inc].size()*sizeof(double));
	my_bgzf_write(fstbg,&(bres[inc][0]),bres[inc].size()*sizeof(double));
	inc++;
      }
    if(arg->chooseChr!=NULL)
      break;
  }
  delGloc(saf,nSites);
  destroy(gls,nSites);
  destroy_args(arg);
  for(int i=0;i<sfs.size();i++)
    delete [] sfs[i];
#if 0
  fprintf(stderr,"\n\t-> NB NB output is no longer log probs of the frequency spectrum!\n");
  fprintf(stderr,"\t-> Output is now simply the expected values! \n");
  fprintf(stderr,"\t-> You can convert to the old format simply with log(norm(x))\n");
#endif
  bgzf_close(fstbg);
  fclose(fstfp);
  fprintf(stderr,"\t-> fst index finished with no errors!\n");
  return 0;
}

template <typename T>
int main_opt(args *arg){

  std::vector<persaf *> &saf =arg->saf;
  for(int i=0;i<saf.size();i++)
    assert(saf[i]->pos!=NULL&&saf[i]->saf!=NULL);
  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=nsites(saf,arg);
  }
  if(fsizes<T>(saf,nSites)>getTotalSystemMemory())
    fprintf(stderr,"\t-> Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    
  fprintf(stderr,"\t-> nSites: %lu\n",nSites);
  float bytes_req_megs =(float) fsizes<T>(saf,nSites)/1024/1024;
  float mem_avail_megs =(float) getTotalSystemMemory()/1024/1024;//in percentile
  //fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"\t-> The choice of -nSites will require atleast: %f megabyte memory, that is at least: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim=(int) parspace(saf);
  double *sfs=new double[ndim];

  //temp used for checking pos are in sync
  setGloc(saf,nSites);
  while(1) {
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop,NULL,NULL,arg->fl,1);//read nsites from data
    int b=0;  
    //fprintf(stderr,"\t\tRET:%d gls->x:%lu\n",ret,gls[0]->x);
    if(ret==-2&&gls[0]->x==0)//no more data in files or in chr, eith way we break;
      break;
#if 0
    if(saf.size()==1){
      if(ret!=-2){
	if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	  //	  fprintf(stderr,"continue continue\n");
	  continue;
	}
      }
    }else
#endif
      {
      if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	//fprintf(stderr,"continue continue\n");
	continue;
      }

    }
    if(gls[0]->x==0)
      continue;
    
    fprintf(stderr,"\t-> Will run optimization on nSites: %lu\n",gls[0]->x);

  neverusegoto:
    if(arg->bootstrap)
      fprintf(stderr,"Will do bootstrap replicate %d/%d\n",b+1,arg->bootstrap);
    if(arg->sfsfname.size()!=0)
	readSFS(arg->sfsfname[0],ndim,sfs);
      else{
	if(arg->seed==-1){
	  for(int i=0;i<ndim;i++)
	    sfs[i] = (i+1)/((double)(ndim));
	}else{
	  for(int i=0;i<ndim;i++){
	    double r=drand48();
	    while(r==0.0)
	      r = drand48();
	    sfs[i] = r;
	  }
	}
	
      }
      normalize(sfs,ndim);
      
      if(bootstrap==NULL &&arg->bootstrap)
	bootstrap = new size_t[gls[0]->x];
      
      if(bootstrap){
	for(size_t i=0;i<gls[0]->x;i++)
	  bootstrap[i] = lrand48() % gls[0]->x;
	std::sort(bootstrap,bootstrap+gls[0]->x);
      }
      double lik;

      if(arg->emAccl==0)
	lik = em<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim,gls);
      else
	lik = emAccl<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim,gls,arg->emAccl);

      fprintf(stderr,"likelihood: %f\n",lik);
      fprintf(stderr,"------------\n");
#if 1
      //    fprintf(stdout,"#### Estimate of the sfs ####\n");
      //all gls have the same ->x. That means the same numbe of sites.
      for(int x=0;x<ndim;x++)
	fprintf(stdout,"%f ",((double)gls[0]->x)*sfs[x]);
      fprintf(stdout,"\n");
      fflush(stdout);
#endif
      if(++b<arg->bootstrap)
	goto neverusegoto;
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
    if(arg->onlyOnce)
      break;
  }
  delGloc(saf,nSites);
  destroy(gls,nSites);
  destroy_args(arg);
  delete [] sfs;
  
  fprintf(stderr,"\n\t-> NB NB output is no longer log probs of the frequency spectrum!\n");
  fprintf(stderr,"\t-> Output is now simply the expected values! \n");
  fprintf(stderr,"\t-> You can convert to the old format simply with log(norm(x))\n");
  return 0;
}

int fst(int argc,char**argv){
  if(argc==0){
    fprintf(stderr,"\t-> Possible options: index print\n");
    return 0;
  }
  if(!strcasecmp(*argv,"index"))  
    fst_index(--argc,++argv);
  else  if(!strcasecmp(*argv,"print"))  
    fst_print(--argc,++argv);
  else if(!strcasecmp(*argv,"stats"))  
    fst_stat(--argc,++argv);
  else if(!strcasecmp(*argv,"stats2"))  
    fst_stat2(--argc,++argv);
  else{
    fprintf(stderr,"unknown option: \'%s\'\n",*argv);
  }
  return 0;
}


void writeAllThetas(BGZF *dat,FILE *idx,char *tmpChr,int64_t &offs,std::vector<int> &p,std::vector<float> *res,int nChr){
  assert(dat!=NULL);
  assert(idx!=NULL);
  assert(p.size()==res[0].size());
  fprintf(stderr,"\t-> Writing %lu sites for chr:%s\n",p.size(),tmpChr);
  for(int i=1;i<5;i++)
    assert(p.size()==res[i].size());//DRAGON, might be discarded during compilation
      
  if(p.size()!=0&&tmpChr!=NULL){
    //write clen and chromoname for both index and bgzf
    size_t clen = strlen(tmpChr);
    fwrite(&clen,sizeof(size_t),1,idx);
    fwrite(tmpChr,1,clen,idx);
    if(sizeof(size_t)!=bgzf_write(dat,&clen,sizeof(size_t))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    if(clen!=bgzf_write(dat,tmpChr,clen)){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }

    //write number of sites for both index and bgzf
    size_t tt = p.size();
    fwrite(&tt,sizeof(size_t),1,idx);
    if(sizeof(size_t)!=bgzf_write(dat,&tt,sizeof(size_t))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    //write nChr for both index and bgzf
    fwrite(&nChr,sizeof(int),1,idx);
    if(sizeof(int)!=bgzf_write(dat,&nChr,sizeof(int))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    //write bgzf offset into idx
    fwrite(&offs,sizeof(int64_t),1,idx);
    for(int i=0;i<p.size();i++)
      aio::bgzf_write(dat,&p[i],sizeof(int));

    for(int i=0;i<5;i++)
      for(int j=0;j<p.size();j++)
	aio::bgzf_write(dat,&res[i][j],sizeof(float));
  }

  //reset
  offs = bgzf_tell(dat);
  p.clear();
  for(int i=0;i<5;i++)
    res[i].clear();
}


int saf2theta(int argc,char**argv){
  int fold =0;
  const char *THETAS =".thetas.gz";
  const char *THETASIDX =".thetas.idx";
  BGZF *theta_dat;
  FILE *theta_idx;
  std::vector<float> theta_res[5];// = new std::vector<float>[5];
  std::vector<int> theta_pos;

  if(argc==0){
    fprintf(stderr,"\t-> Possible options: addoptions\n");
    return 0;
  }
  args *arg = getArgs(argc,argv);
  if(arg->outname==NULL){
    fprintf(stderr,"\t-> Must supply -outname for generating outputfiles\n");
    return 0;
  }
  if(arg->saf.size()!=1){
    fprintf(stderr,"\t-> Must supply one, and only one saf.idx file\n");
    return 0;
  }
  if(arg->sfsfname.size()!=1){
    fprintf(stderr,"\t-> Must supply one, and only one -sfs argument which should contain the prior\n");
    return 0;
  }
  if(arg->nSites==0){
    int block = 4096;
    fprintf(stderr,"\t-> Will read chunks of size: %d\n",block);
    arg->nSites = block;
  }
  arg->saf[0]->kind=2; //<-important orhterwise we dont read positions from saffles\n
  double *prior = new double[arg->saf[0]->nChr+1];
  readSFS(arg->sfsfname[0],arg->saf[0]->nChr+1,prior);
  for(int i=0;i<arg->saf[0]->nChr+1;i++)
    prior[i] = log(prior[i]);
  char buf[8] = "thetav2";
  theta_dat = aio::openFileBG(arg->outname,THETAS);
  theta_idx = aio::openFile(arg->outname,THETASIDX);
  aio::bgzf_write(theta_dat,buf,8);
  fwrite(buf,1,8,theta_idx);
  //  theta_res = new std::vector<float>[5];
  int64_t offs_thetas = bgzf_tell(theta_dat);
  
  
  double aConst=0;
  int nChr = arg->saf[0]->nChr;
  fprintf(stderr,"\t-> nChr:%d\n",nChr);
  for(int i=1;i<nChr;i++)
    aConst += 1.0/i;
  aConst = log(aConst);//this is a1
  
  
  double aConst2 = log((nChr*(nChr-1))/2.0);//choose(nChr,2)
  double aConst3 = log((1.0*nChr-1.0));
  
  double *scalings = new double [nChr+1];
  for(int i=0;i<nChr+1;i++)
    scalings[i] = log(i)+log(nChr-i);

  std::vector<Matrix<float> *> gls;
  for(int i=0;i<arg->saf.size();i++)
    gls.push_back(alloc<float>(arg->nSites,arg->saf[i]->nChr+1));
  
  setGloc(arg->saf,arg->nSites);
  int *posiToPrint = new int[arg->nSites];

  char *tmpChr =NULL;
  static char *curChr=NULL;//why static?

  while(1) {
    int ret=readdata(arg->saf,gls,arg->nSites,arg->chooseChr,arg->start,arg->stop,posiToPrint,&curChr,arg->fl,0);//read nsites from data
    if(arg->chooseChr!=NULL){
      if(curChr==NULL)
	curChr=strdup(arg->chooseChr);
    }
    if(tmpChr==NULL)
      tmpChr = strdup(curChr);
    if(strcmp(tmpChr,curChr)!=0){
      writeAllThetas(theta_dat,theta_idx,tmpChr,offs_thetas,theta_pos,theta_res,nChr);
      free(tmpChr);
      tmpChr=strdup(curChr);
    }
    //calc thetas

    for(int s=0;s<gls[0]->x;s++){
      double workarray[nChr+1];
      for(int i=0;i<nChr+1;i++)//gls->mat is float lets pluginto double
	workarray[i] = gls[0]->mat[s][i];

      //calculate post probs
      double tsum =exp(workarray[0] + prior[0]);
      for(int i=1;i<nChr+1;i++)
	tsum += exp(workarray[i]+prior[i]);
      tsum = log(tsum);
      
      for(int i=0;i<nChr+1;i++){
	workarray[i] = workarray[i]+prior[i]-tsum;
	//      fprintf(stderr,"[%d]:%f\n",i,(workarray[i]));
      }
      //exit(0);
      //First find thetaW: nSeg/a1
      double pv,seq;
      if(fold)
	pv = 1-exp(workarray[0]);
      else
	pv = 1-exp(workarray[0])-exp(workarray[nChr]);
      //      fprintf(stderr,"pv:%f work[0]:%f 2k:%f\n",pv,workarray[0],workarray[nChr]);
      //      exit(0);
      if(pv<0)//catch underflow
	seq=log(0.0);
      else
	seq = log(pv)-aConst;//watterson
      theta_res[0].push_back(seq);
      //     ksprintf(&kb,"%s\t%d\t%f\t",header->target_name[pars->refId],pars->posi[i]+1,seq);
      theta_pos.push_back(posiToPrint[s]);
      // fprintf(stderr,"posiToPrint[s]:%d\n",posiToPrint[s]);
      if(fold==0) {
	double pairwise=0;    //Find theta_pi the pairwise
	double thL=0;    //Find thetaL sfs[i]*i;
	double thH=0;//thetaH sfs[i]*i^2
	for(size_t ii=1;ii<nChr;ii++){
	  
	  pairwise += exp(workarray[ii]+scalings[ii] );
	  double li=log(ii);
	  
	  thL += exp(workarray[ii])*ii;
	  thH += exp(2*li+workarray[ii]);
	}
	theta_res[1].push_back(log(pairwise)-aConst2);
	theta_res[2].push_back(workarray[1]);
	theta_res[3].push_back(log(thH)-aConst2);
	theta_res[4].push_back(log(thL)-aConst3);
      }else{
	double pairwise=0;    //Find theta_pi the pairwise
	for(size_t ii=1;ii<nChr+1;ii++)
	  pairwise += exp(workarray[ii]+scalings[ii] );
	theta_res[1].push_back(log(pairwise)-aConst2);
	for(int i=2;i<=4;i++)
	  theta_res[i].push_back(log(0));
      }
    }
    
#if 0 //just for printout
    for(int s=0;s<gls[0]->x;s++){
      if(arg->chooseChr==NULL)
	fprintf(stdout,"%s\t%d",curChr,posiToPrint[s]+1);
      else
	fprintf(stdout,"%s\t%d",arg->chooseChr,posiToPrint[s]+1);
      for(int i=0;i<arg->saf.size();i++)
	for(int ii=0;ii<gls[i]->y;ii++)
	  fprintf(stdout,"\t%f",log(gls[i]->mat[s][ii]));
      fprintf(stdout,"\n");
    }
#endif
    if(ret==-3&&gls[0]->x==0){//no more data in files or in chr, eith way we break;g
      //      fprintf(stderr,"breaking change of chr\n");
      break;
    }
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
    if(arg->onlyOnce)
      break;
  }
  if(theta_pos.size()>0)
    writeAllThetas(theta_dat,theta_idx,tmpChr,offs_thetas,theta_pos,theta_res,nChr);
  
  fclose(theta_idx);
  bgzf_close(theta_dat);
  if(curChr){
    //  free(curChr);//<- DRAGON leak?
    curChr=NULL;
  }
  if(tmpChr)
    free(tmpChr);
 
  delGloc(arg->saf,arg->nSites);
  destroy(gls,arg->nSites);
  delete [] prior;
  //  delete [] theta_res;
  delete [] scalings;
  delete [] posiToPrint;
  destroy_args(arg);
  return 0;
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
    //    fprintf(stderr, "\t->------------------\n\t-> ./realSFS\n\t->------------------\n");
    // fprintf(stderr,"\t-> This is the new realSFS program which works on the newer binary files from ANGSD!!\n");
    fprintf(stderr, "\t-> ---./realSFS------\n\t-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:\n\n\t-> Estimate the SFS for entire genome??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx \n");
    fprintf(stderr, "\n\t-> 1) Estimate the SFS for entire chromosome 22 ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -r chr22 \n");
    fprintf(stderr, "\n\t-> 2) Estimate the 2d-SFS for entire chromosome 22 ??\n");
    fprintf(stderr,"\t-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 \n");

    fprintf(stderr, "\n\t-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -nSites 500000000 \n");

    fprintf(stderr, "\n\t-> 4) Estimate the SFS around a gene ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 \n");
    fprintf(stderr, "\n\t-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications]\n");

    fprintf(stderr,"\n\t-> See realSFS print for possible print options\n");
    fprintf(stderr,"\t-> Use realSFS print_header for printing the header\n");
    fprintf(stderr,"\t-> Use realSFS cat for concatenating saf files\n");

    fprintf(stderr,"\n\t->------------------\n\t-> NB: Output is now counts of sites instead of log probs!!\n");
    fprintf(stderr,"\t-> NB: You can print data with ./realSFS print afile.saf.idx !!\n");
    fprintf(stderr,"\t-> NB: Higher order SFS's can be estimated by simply supplying multiple .saf.idx files!!\n");
    fprintf(stderr,"\t-> NB: Program uses accelerated EM, to use standard EM supply -m 0 \n");
    return 0;
  }
  ++argv;
  --argc;
  if(!strcasecmp(*argv,"print"))
    print<float>(--argc,++argv);
  else if(!strcasecmp(*argv,"cat"))
    saf_cat(--argc,++argv);
  else if(!strcasecmp(*argv,"fst"))
    fst(--argc,++argv);
  else if(!strcasecmp(*argv,"saf2theta"))
    saf2theta(--argc,++argv);
  else if(!strcasecmp(*argv,"print_header"))
    print_header(--argc,++argv);
  else {
    args *arg = getArgs(argc,argv);
    if(arg->saf.size()>1)
      fprintf(stderr,"\t-> Multi SFS is 'still' under development. Please report strange behaviour\n");
    if(!arg)
      return 0;

    if(isatty(fileno(stdout))){
      fprintf(stderr,"\t-> You are printing the optimized SFS to the terminal consider dumping into a file\n");
      fprintf(stderr,"\t-> E.g.: \'./realSFS");
      for(int i=0;i<argc;i++)
	fprintf(stderr," %s",argv[i]);
      fprintf(stderr," >sfs.ml.txt\'\n");   
    }
  

    main_opt<float>(arg);
    
  }
  if(bootstrap!=NULL)
    delete [] bootstrap;

  if(dumpedFiles.size()){
    fprintf(stderr,"\t-> Output filenames:\n");
    for(int i=0;i<(int)dumpedFiles.size();i++){
      fprintf(stderr,"\t\t->\"%s\"\n",dumpedFiles[i]);
      free(dumpedFiles[i]);
    }
  }
  return 0;
}
