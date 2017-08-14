/*
  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

 */

#include <cstdio>
#include <zlib.h>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "bfgs.h"
#include "realSFS.h"
#include "realSFS_saf2.cpp"
#define LENS 4096

const char*fname1 = NULL;
const char*fname2 = NULL;//only used for 2dsfs
const char*fname3 = NULL;//only used for 3dsfs
int chr1=0;
int chr2=0; //only used for 2dsfs
int chr3=0; //only used for 2dsfs

int dim =0; //will be (chr1+1)*(chr2+1); for the 2dsfs nChr1+1 for the 1dsfs

char *sfsfname=NULL;

int nThreads =4;
double tole = 1e-6;
int maxIter = 1e2;
int SIG_COND =1;
size_t nSites = 0;
int doBFGS = 0;//defaults is to use EM-optimization
int calcLike =0;//if a sfs input file has been supplied, should we just print the likelihood
int noGrad = 0;//should we use gradients for the bfgs.
int isList = 0;
pthread_t *thd=NULL;



int fexists(const char* str){
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.                             
}


size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}
size_t subfun(const char *fname,int nChr){
  size_t nbytes = fexists(fname)?fsize(fname):0;
  if(nbytes %(sizeof(double)*(nChr+1))){
    fprintf(stderr,"saffile: %s looks broken?\n",fname);
    exit(0);
  }
  return nbytes/sizeof(double)/(nChr+1);
}

size_t calcNsites(const char *fname,int nChr){
  //  fprintf(stderr,"[%s] fname: %s nChr: %d isList:%d\n",__FUNCTION__,fname,nChr,isList);
  if(isList==0)
    return subfun(fname,nChr);
  else{
    gzFile gz=getGz(fname);
    char buf[LENS];
    size_t tot=0;
    while(gzgets(gz,buf,LENS)){
      //  fprintf(stderr,"buf:%s\n",buf);
      for(char *tok = strtok(buf,"\n\t ");tok!=NULL;tok=strtok(NULL,"\n\t ")){
	if(tok[0]=='#')
	  continue;
	size_t per_file = subfun(buf,nChr);
	tot+=per_file;
      }
    }
    gzclose(gz);
    return tot;
  }
}



Matrix<double> alloc(size_t x,size_t y){
  //fprintf(stderr,"def=%f\n",def);
  Matrix<double> ret;
  ret.x=x;
  ret.y=y;
  ret.mat= new double*[x];
  for(size_t i=0;i<x;i++)
    ret.mat[i]=new double[y];
  return ret;
}
void dalloc(Matrix<double> &ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret.mat[i];
  delete [] ret.mat;
}



typedef struct emPars_t{
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
  Matrix<double> *GL1;
  Matrix<double> *GL2;
  Matrix<double> *GL3;
  int from;
  int to;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  double *grad;//used for bfgs, length is dim-1; 
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

gzFile getGz(const char *pname){
  if(fexists(pname)==0){
    fprintf(stderr,"Problem opening file: %s\n",pname);
    exit(0);
  }
  gzFile fp = Z_NULL;
  if(((fp=gzopen(pname,"r")))==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'\n",pname);
  }
  return fp;
} 



void readGL(gzFile fp,size_t nSites,int nChr,Matrix<double> &ret){
  //  fprintf(stderr,"[%s] fname:%s nSites:%d nChr:%d\n",__FUNCTION__,fname,nSites,nChr);
  ret.x=nSites;
  ret.y=nChr+1;
  size_t i;
  for(i=0;SIG_COND&&i<nSites;i++){
    int bytes_read = gzread(fp,ret.mat[i],sizeof(double)*(nChr+1));

    if(bytes_read!=0 && bytes_read<sizeof(double)*(nChr+1)){
      fprintf(stderr,"Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0)
      break;
    
    for(size_t j=0;j<nChr+1;j++)
      ret.mat[i][j] = exp(ret.mat[i][j]);
  }
  
  ret.x=i;
  if(SIG_COND==0)
    exit(0);
  
}

void readGL2(gzFile fp,size_t nSites,int nChr,Matrix<double> &ret){
  //  fprintf(stderr,"[%s] fname:%s nSites:%d nChr:%d\n",__FUNCTION__,fname,nSites,nChr);
  ret.y=nChr+1;
  char buf[LENS];
  size_t i=0;
  while(gzgets(fp,buf,LENS)){
    //strange construct for strtok, because firstcall if with buf, rest with NULL
    for(char *tok = strtok(buf,"\t\n ");tok!=NULL;tok=strtok(NULL,"\t\n ")){
      if(tok[0]=='#')
	continue;
      gzFile gz = getGz(buf);
    
      while(SIG_COND&&gzread(gz,ret.mat[i],sizeof(double)*(nChr+1))){
	for(size_t j=0;j<nChr+1;j++)
	  ret.mat[i][j] = exp(ret.mat[i][j]);
	i++;
      }
      fprintf(stderr,"Done reading file: \'%s\'\n",buf);
      gzclose(gz);
      ret.x=i;
    }
  }
  if(SIG_COND==0)
    exit(0);
  
}



std::vector<int> getPosi(const char*fname){
  char *pname = append(fname,".sfs.pos");
  if(fexists(pname)==0){
    fprintf(stderr,"Problem opening file: %s\n",pname);
    exit(0);
  }
  std::vector<int> ret;
  
  gzFile fp = gzopen(pname,"r");
  char buf[1024];
  while(gzgets(fp,buf,1024)){
    strtok(buf,"\n\t ");
    ret.push_back(atoi(strtok(NULL,"\t\n ")));
  }
  gzclose(fp);
  return ret;
}




FILE *getFILE(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}



void readSFS(const char*fname,int hint,double *ret){
  fprintf(stderr,"reading: %s\n",fname);
  FILE *fp = getFILE(fname,"r");
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
    fprintf(stderr,"\t-> Problem with size of dimension of prior %d vs %lu\n",hint,res.size());
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


void getArgs(int argc,char **argv){
  if(argc==0)
    return;
  if(argc%2){
    fprintf(stderr,"Extra args must be given as -par VALUE\n");
    exit(0);
  }
  while(*argv){
    //  fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-P"))
      nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      nSites = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-use-BFGS"))
      doBFGS = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-calcLike"))
      calcLike = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-noGrad"))
      noGrad = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-isList"))
      isList = atoi(*(++argv));
    
    
    else  if(!strcasecmp(*argv,"-start")){
      sfsfname = *(++argv);
    }else{
      fprintf(stderr,"Unknown arg: %s\n",*argv);
      exit(0);
    }
    argv++;
  }
}

double lik1(double *sfs,Matrix<double> *ret,int from,int to){
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



double lik1_master(){
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


//this sfs in this function is the transformed dim-1
double lik1_bfgs(const double *sfs_1,const void *dats){
  if(SIG_COND==0)
    exit(0);
  double p[dim];
  double ts=1;
  for(int i=0;i<dim-1;i++)
    ts +=sfs_1[i];
  emp[0].sfs[0]=1.0/ts;
  for(int i=0;i<dim-1;i++)
    emp[0].sfs[i+1]=sfs_1[i]/ts;
  return -lik1_master();
}


void lik1_grad_bfgs(double *sfs,double *grad,double fac,Matrix<double> *GL1,int from,int to){
  for(int i=0;i<dim-1;i++){
    //  fprintf(stderr,"%d:%f\n",i,sfs[i]);
    grad[i] =0.0;
  }


  for(int s=from;s<to;s++){
    double ret =GL1->mat[s][0];
    for(int i=0;i<dim-1;i++)
      ret  += sfs[i]*GL1->mat[s][i+1];
    for(int i=0;i<dim-1;i++)
      grad[i] += fac -GL1->mat[s][i+1]/ret ;
  }
  
}


void *lik1_grad_bfgs_slave(void *p){
  emPars &pars = emp[(size_t) p];
  lik1_grad_bfgs(pars.sfs,pars.grad,pars.lik,pars.GL1,pars.from,pars.to);

  return NULL;
}


//this is not the sfs but the transfomred dim-1 sfs
void lik1_grad_bfgs_master(const double *sfs,double *grad){
  if(SIG_COND==0)
    exit(0);
  double fac = 1;
  for(int i=0;i<dim-1;i++){
    
    fac += sfs[i] ;
    // fprintf(stderr,"sfs=%f\tfac:%f\n",sfs[i],fac);
  }
  fac = 1.0/fac;
  // fprintf(stderr,"fac=%f\n",fac);
  

  for(size_t i=0;i<nThreads;i++){
    memcpy(emp[i].sfs,sfs,(dim-1)*sizeof(double));
    emp[i].lik=fac;
    int rc = pthread_create(&thd[i],NULL,lik1_grad_bfgs_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);

    
  memcpy(grad,emp[0].grad,(dim-1)*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<dim-1;j++)
      grad[j] +=  emp[i].grad[j];
  }

#if 0
  double myVar = 1;
  for(int i=0;i<dim-1;i++)
    myVar += sfs[i];
  

  for(int j=0;j<dim-1;j++)
    grad[j] +=(1.0*emp[0].GL1->x)/myVar;
#endif 

#if 0
  for(int i=0;i<nThreads;i++){
    for(int j=0;j<dim-1;j++)
      fprintf(stdout,"%g ",emp[i].grad[j]);
    fprintf(stdout,"\n");
  }
  exit(0);
#endif


#if 0
  for(int j=0;j<dim-1;j++)
    fprintf(stdout,"%g ",grad[j]);
  fprintf(stdout,"\n");
  exit(0);
#endif

}








//sfs input is here the untransformed, e.g length=dim
void bfgs(double *sfs,Matrix<double> *GL1){
  double p[dim-1];
  for(int i=1;i<dim;i++)
    p[i-1] = sfs[i]/sfs[0];
  
  //initialize the stuff needed for the bfgs
  
  double para[dim-1];  /* parameters */
  double min[dim-1]; /* lower bound */
  double max[dim-1]; /* upper bound */
  int nbd[dim-1]; /* boundary constraint option; nbd[i] = 0 or 1 or 2 or 3; see bfgs.h for details */
  int noisy=-1; /* bigger value produces more detailed output */ 
  
  for (int i=0; i<dim-1; i++){
    min[i] = 0.0000001;
    max[i] = 10.0;
    nbd[i] = 2;
  }
  double fnmax ;
  if(noGrad==0)
    fnmax= findmax_bfgs( dim-1, p,NULL, lik1_bfgs, lik1_grad_bfgs_master, min, max, nbd, noisy );  
  else
    fnmax= findmax_bfgs( dim-1, p,NULL, lik1_bfgs,NULL, min, max, nbd, noisy );  
  
  //tranform back
  double ts=1;
  for(int i=0;i<dim-1;i++)
    ts  += p[i];
  sfs[0]=1.0/ts;
  for(int i=0;i<dim-1;i++)
    sfs[i+1]=p[i]/ts;


}



void emStep1(double *pre,Matrix<double> *GL1,double *post,int start,int stop){
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

  emStep1(pars.sfs,pars.GL1,pars.post,pars.from,pars.to);

  return NULL;
}


void emStep1_master(double *post){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep1_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  memcpy(post,emp[0].post,dim*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<dim;j++)
      post[j] += emp[i].post[j];
  }
  
  normalize(post,dim);

#if 0
  for(int i=0;i<nThreads;i++){
    for(int j=0;j<dim;j++)
      fprintf(stdout,"%f ",emp[i].post[j]);
    fprintf(stdout,"\n");
  }
#endif
  
}






void em1(double *sfs,Matrix<double> *GL1,double tole=0.01,int maxIter=10){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik1_master();
  else
    oldLik = lik1(sfs,GL1,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);
  fflush(stderr);

  double *tmp=new double[dim];
  
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep1_master(tmp);
    else
      emStep1(sfs,GL1,tmp,0,GL1->x);
    
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];

    if(nThreads>1)
      lik = lik1_master();
    else
      lik = lik1(sfs,GL1,0,GL1->x);

    fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  
  delete [] tmp;
}


double lik2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,size_t start,size_t stop){
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


void setThreadPars(Matrix<double> *GL1,Matrix<double> *GL2,double *sfs,int nThreads){
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
    emp[i].grad=new double[dim-1];
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


void setThreadPars(Matrix<double> *GL1,Matrix<double> *GL2,Matrix<double>*GL3,double *sfs,int nThreads){
  emp = new emPars[nThreads];
  int blockSize = GL1->x/nThreads;
  for(int i=0;i<nThreads;i++){
    emp[i].threadId = i;
    emp[i].GL1=GL1;
    emp[i].GL2=GL2;
    emp[i].GL3=GL3;
    emp[i].from =0;
    emp[i].to=blockSize;
    emp[i].sfs = sfs;
    emp[i].post=new double[dim];
    emp[i].grad=new double[dim-1];
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



double lik2_master(){
  
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




void emStep2(double *pre,Matrix<double> *GL1,Matrix<double> *GL2,double *post,int start,int stop){
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

  emStep2(pars.sfs,pars.GL1,pars.GL2,pars.post,pars.from,pars.to);

  return NULL;
}


void emStep2_master(double *post){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep2_slave,(void*) i);
    if(rc)
      fprintf(stderr,"error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  memcpy(post,emp[0].post,dim*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<dim;j++)
      post[j] += emp[i].post[j];
  }
  
  normalize(post,dim);

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





void em2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,double tole=0.01,int maxIter=10){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik2_master();
  else
    oldLik = lik2(sfs,GL1,GL2,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);

  double *tmp = new double[dim];//<- wont be cleaned up, but is only allocated once
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep2_master(tmp);
    else
      emStep2(sfs,GL1,GL2,tmp,0,GL1->x);

    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    if(nThreads>1)
      lik = lik2_master();
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

void em3(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,Matrix<double> *GL3,double tole=0.01,int maxIter=10){
  double oldLik,lik;
  if(nThreads>1)
    oldLik = lik2_master();
  else
    oldLik = lik2(sfs,GL1,GL2,0,GL1->x);
  fprintf(stderr,"startlik=%f\n",oldLik);

  double *tmp = new double[dim];//<- wont be cleaned up, but is only allocated once
  for(int it=0;SIG_COND&&it<maxIter;it++) {
    if(nThreads>1)
      emStep2_master(tmp);
    else
      emStep2(sfs,GL1,GL2,tmp,0,GL1->x);

    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    if(nThreads>1)
      lik = lik2_master();
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


Matrix<double>  merge(Matrix<double> &pop1,Matrix<double> &pop2){
  fprintf(stderr,"[%s]\n",__FUNCTION__);
  Matrix<double> ret;
  ret.x = pop1.x;
  ret.y = pop1.y*pop2.y;
  ret.mat=new double*[pop1.x];
  for(int s=0;SIG_COND&&s<pop1.x;s++){
    ret.mat[s] = new double[pop1.y*pop2.y];
    //    fprintf(stderr,"tmpdim=%zu %zu\n",pop1.y,pop2.y);
    int inc=0;
    for(int x=0;x<pop1.y;x++)
      for(int y=0;y<pop2.y;y++)
	ret.mat[s][inc++] = exp(pop1.mat[s][x]+pop2.mat[s][y]);
    delete [] pop1.mat[s];
    delete [] pop2.mat[s];
  }
  delete [] pop1.mat; delete [] pop2.mat;
  pop1.x=pop1.y=pop2.x=pop2.y=0;
  return ret;
}

void print(Matrix<double> &mat,FILE *fp){
  for(int x=0;x<mat.x;x++){
    //    fprintf(stderr,"x=%d\n",x);
    for(int y=0;y<mat.y;y++)
      fprintf(fp,"%f ",mat.mat[x][y]);
    fprintf(fp,"\n");
  }
  

}


int main_2dsfs(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./emOptim2 2dsfs pop1 pop2 nChr1 nChr2 [-start FNAME -P nThreds -tole tole -maxIter ] (only works if the two saf files covers the same region)\n");
    return 0;
  }
  argv++;
  argc--;
  fname1 = *(argv++);
  fname2 = *(argv++);
  argc -=2;
  chr1 = atoi(*(argv++));
  chr2 = atoi(*(argv++));
  argc -=2;
  getArgs(argc,argv);
  if(nSites==0){
    if(fsize(fname1)+fsize(fname2)>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n");
    //this doesnt make sense if ppl supply a filelist containing safs
     nSites=calcNsites(fname1,chr1);
  }
  fprintf(stderr,"fname1:%sfname2:%s chr1:%d chr2:%d startsfs:%s nThreads=%d tole=%f maxIter=%d nSites:%lu\n",fname1,fname2,chr1,chr2,sfsfname,nThreads,tole,maxIter,nSites);
  float bytes_req_megs = nSites*(sizeof(double)*(chr1+1) + sizeof(double)*(chr2+1)+2*sizeof(double*))/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);
  
#if 0
  //read in positions, not used, YET...
  std::vector<int> p1 = getPosi(fname1);
  std::vector<int> p2 = getPosi(fname2);
  fprintf(stderr,"nSites in pop1: %zu nSites in pop2: %zu\n",p1.size(),p2.size());
#endif

  if(nSites==0){
    if(calcNsites(fname1,chr1)!=calcNsites(fname2,chr2)){
      fprintf(stderr,"Problem with number of sites in file: %s and %s\n",fname1,fname2);
      exit(0);
    }
    nSites=calcNsites(fname1,chr1);
  }
  gzFile gz1=getGz(fname1);
  gzFile gz2=getGz(fname2);
  
  dim=(chr1+1)*(chr2+1);
  
  Matrix<double> GL1=alloc(nSites,chr1+1);
  Matrix<double> GL2=alloc(nSites,chr2+1);
  dim=GL1.y*GL2.y;
  
  double *sfs = new double[dim];
  while(1){
    if(isList ==0){
      readGL(gz1,nSites,chr1,GL1);
      readGL(gz2,nSites,chr2,GL2);
    }else{
      readGL2(gz1,nSites,chr1,GL1);
      readGL2(gz2,nSites,chr2,GL2);
    }
      
    assert(GL1.x==GL2.x);
    if(GL1.x==0)
      break;
    
    if(sfsfname!=NULL){
      readSFS(sfsfname,dim,sfs);
    }else{
      for(int i=0;i<dim;i++)
	sfs[i] = (i+1)/((double)dim);
      normalize(sfs,dim);
    }
    
    setThreadPars(&GL1,&GL2,sfs,nThreads);
    if(calcLike==0){
      if(SIG_COND) 
	em2(sfs,&GL1,&GL2,tole,maxIter);
    }
    double lik;
    if(nThreads>1)
      lik = lik1_master();
    else
      lik = lik1(sfs,&GL1,0,GL1.x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    int inc=0;
    for(int x=0;x<chr1+1;x++){
      for(int y=0;y<chr2+1;y++)
	fprintf(stdout,"%f ",log(sfs[inc++]));
      fprintf(stdout,"\n");
    }
#endif
    if(isList==1)
      break;
  }
  dalloc(GL1,nSites);
  dalloc(GL2,nSites);
  gzclose(gz1);
  gzclose(gz2);
  return 0;
}


int main_3dsfs(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./emOptim2 2dsfs pop1 pop2 nChr1 nChr2 [-start FNAME -P nThreds -tole tole -maxIter ] (only works if the two saf files covers the same region)\n");
    return 0;
  }
  argv++;
  argc--;
  fname1 = *(argv++);
  fname2 = *(argv++);
  fname3 = *(argv++);
  argc -=3;
  chr1 = atoi(*(argv++));
  chr2 = atoi(*(argv++));
  chr3 = atoi(*(argv++));
  argc -=3;
  getArgs(argc,argv);
  if(nSites==0){
    if(fsize(fname1)+fsize(fname2)+fsize(fname3)>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n");
    //this doesnt make sense if ppl supply a filelist containing safs
     nSites=calcNsites(fname1,chr1);
  }
  fprintf(stderr,"fname1:%sfname2:%sfname3:%s chr1:%d chr2:%d chr3:%d startsfs:%s nThreads=%d tole=%f maxIter=%d nSites:%lu\n",fname1,fname2,fname3,chr1,chr2,chr3,sfsfname,nThreads,tole,maxIter,nSites);
  float bytes_req_megs = nSites*(sizeof(double)*(chr1+1) + sizeof(double)*(chr2+1)+sizeof(double)*(chr3+1)+2*sizeof(double*))/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);
  
#if 0
  //read in positions, not used, YET...
  std::vector<int> p1 = getPosi(fname1);
  std::vector<int> p2 = getPosi(fname2);
  fprintf(stderr,"nSites in pop1: %zu nSites in pop2: %zu\n",p1.size(),p2.size());
#endif

  if(nSites==0){
    if((calcNsites(fname1,chr1)!=calcNsites(fname2,chr2))||((calcNsites(fname1,chr1)!=calcNsites(fname3,chr2)))){
      fprintf(stderr,"Problem with number of sites in file: %s and %s and %s\n",fname1,fname2,fname3);
      exit(0);
    }
    nSites=calcNsites(fname1,chr1);
  }
  gzFile gz1=getGz(fname1);
  gzFile gz2=getGz(fname2);
  gzFile gz3=getGz(fname3);
  
  dim=(chr1+1)*(chr2+1)*(chr3+1);
  
  Matrix<double> GL1=alloc(nSites,chr1+1);
  Matrix<double> GL2=alloc(nSites,chr2+1);
  Matrix<double> GL3=alloc(nSites,chr3+1);
  dim=GL1.y*GL2.y*GL3.y;
  
  double *sfs = new double[dim];
  while(1){
    if(isList ==0){
      readGL(gz1,nSites,chr1,GL1);
      readGL(gz2,nSites,chr2,GL2);
      readGL(gz3,nSites,chr3,GL3);
    }else{
      readGL2(gz1,nSites,chr1,GL1);
      readGL2(gz2,nSites,chr2,GL2);
      readGL2(gz3,nSites,chr3,GL3);
    }
      
    assert(GL1.x==GL2.x);
    assert(GL1.x==GL3.x);
    if(GL1.x==0)
      break;
    
    if(sfsfname!=NULL){
      readSFS(sfsfname,dim,sfs);
    }else{
      for(int i=0;i<dim;i++)
	sfs[i] = (i+1)/((double)dim);
      normalize(sfs,dim);
    }
    
    setThreadPars(&GL1,&GL2,&GL3,sfs,nThreads);
    if(calcLike==0){
      if(SIG_COND) 
	em3(sfs,&GL1,&GL2,&GL3,tole,maxIter);
    }
    double lik;
    if(nThreads>1)
      lik = lik1_master();
    else
      lik = lik1(sfs,&GL1,0,GL1.x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    int inc=0;
    for(int x=0;x<chr1+1;x++){
      for(int y=0;y<chr2+1;y++)
	fprintf(stdout,"%f ",log(sfs[inc++]));
      fprintf(stdout,"\n");
    }
#endif
    if(isList==1)
      break;
  }
  dalloc(GL1,nSites);
  dalloc(GL2,nSites);
  dalloc(GL3,nSites);
  gzclose(gz1);
  gzclose(gz2);
  gzclose(gz3);
  return 0;
}




int main_1dsfs(int argc,char **argv){
  if(argc<2){
    fprintf(stderr,"Must supply afile.saf and number of chromosomes\n");
    return 0;
  }
  fname1 = *(argv++);
  chr1 = atoi(*(argv++));
  argc-=2;
 
  getArgs(argc,argv);
  dim=chr1+1;
  //hook for new EJ banded version
  if(isNewFormat(fname1))
    return main_1dsfs_v2(fname1,chr1,nSites,nThreads,sfsfname,tole,maxIter);

  if(nSites==0){//if no -nSites is specified
    if(fsize(fname1)>getTotalSystemMemory())
      fprintf(stderr,"Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n");
    //this doesnt make sense if ppl supply a filelist containing safs
    nSites=calcNsites(fname1,chr1);
  }
  fprintf(stderr,"fname1:%s nChr:%d startsfs:%s nThreads:%d tole=%f maxIter=%d nSites=%lu\n",fname1,chr1,sfsfname,nThreads,tole,maxIter,nSites);
  float bytes_req_megs = nSites*(sizeof(double)*(chr1+1)+sizeof(double*))/1024/1024;
  float mem_avail_megs = getTotalSystemMemory()/1024/1024;//in percentile
  //  fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"The choice of -nSites will require atleast: %f megabyte memory, that is approx: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  

  Matrix<double> GL1=alloc(nSites,dim);
  gzFile gz1=getGz(fname1);  
  double *sfs=new double[dim];
  
  while(1) {
    if(isList==0)
      readGL(gz1,nSites,chr1,GL1);
    else
      readGL2(gz1,nSites,chr1,GL1);
    if(GL1.x==0)
      break;
    fprintf(stderr,"dim(GL1)=%zu,%zu\n",GL1.x,GL1.y);
   
    
  
    if(sfsfname!=NULL){
      readSFS(sfsfname,dim,sfs);
    }else{
      
      for(int i=0;i<dim;i++)
	sfs[i] = (i+1)/((double)dim);
      if(doBFGS){
	double ts=1;
	for(int i=0;i<dim-1;i++)
	  ts += 0.01/(1.0+i);
	sfs[0]=1.0/ts;
	for(int i=0;i<dim-1;i++)
	  sfs[i+1]  = (0.01/(1.0+i))/ts;
      }
      normalize(sfs,dim);
    }
    //  em2_smart(sfs2,pops,1e-6,1e3);
    setThreadPars(&GL1,NULL,sfs,nThreads);
    if(calcLike==0){
      if(doBFGS==0) 
	em1(sfs,&GL1,tole,maxIter);
      else
	bfgs(sfs,&GL1);
    }
    double lik;
    if(nThreads>1)
      lik = lik1_master();
    else
      lik = lik1(sfs,&GL1,0,GL1.x);
      
    fprintf(stderr,"likelihood: %f\n",lik);
#if 1
    for(int x=0;x<dim;x++)
      fprintf(stdout,"%f ",log(sfs[x]));
    fprintf(stdout,"\n");
    fflush(stdout);
#endif
    if(isList==1)
      break;
  }
  dalloc(GL1,nSites);
  gzclose(gz1);
  delete [] sfs;
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
    fprintf(stderr,"./emOptim2 afile.sfs nChr [-start FNAME -P nThreads -tole tole -maxIter  -nSites ]\n");
    fprintf(stderr,"nChr is the number of chromosomes. (twice the number of diploid invididuals)\n");
    return 0;
  }
  ++argv;
  argc--;

  if(isatty(fileno(stdout))){
    fprintf(stderr,"\t-> You are printing the optimized SFS to the terminal consider dumping into a file\n");
    fprintf(stderr,"\t-> E.g.: \'./emOptim2");
    for(int i=0;i<argc;i++)
      fprintf(stderr," %s",argv[i]);
    fprintf(stderr," >sfs.ml.txt\'\n");   

  }

  if(!strcasecmp(*argv,"2dsfs"))
    main_2dsfs(argc,argv);
  else if(!strcasecmp(*argv,"3dsfs"))
    main_3dsfs(argc,argv);
  else
    main_1dsfs(argc,argv);
  



  return 0;
}
