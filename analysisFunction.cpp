


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <list>
#include <assert.h>
#include <fstream>
#include "analysisFunction.h"


int angsd::getArg(const char* argName,int type,argStruct *arguments){

  int argPos = 1;
  while(argPos <arguments->argc){
    //    fprintf(stderr,"%s vs %s\n",argName,arguments->argv[argPos]);
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      //      fprintf(stderr,"HIT %s vs %s\n",argName,arguments->argv[argPos]);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atoi(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}

char* angsd::getArg(const char* argName,char* type,argStruct *arguments){
  ///fprintf(stderr,"pre %p %s \n",type,argName);
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2){
        return(strdup("-999"));
      }
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(strdup(arguments->argv[argPos+1]));  //VALGRIND says leak. don't care very small DRAGON
    }
    argPos++;
  }
  //  fprintf(stderr,"post %p\n",type);
  return(type);
  
}


char* angsd::getArg(const char* argName, const char* type,argStruct *arguments){
  //fprintf(stderr,"pre %p %s \n",type,argName);
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2){
        return(strdup("-999"));
      }
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(strdup(arguments->argv[argPos+1]));  //VALGRIND says leak. don't care very small DRAGON
    }
    argPos++;
  }
  //assert(0==1);
  //  fprintf(stderr,"post %p\n",type);
  return NULL;
  
}



float angsd::getArg(const char* argName,float type,argStruct *arguments){
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
 if(argPos==arguments->argc-1){
	fprintf(stderr,"Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atof(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}

double angsd::getArg(const char* argName,double type,argStruct *arguments){
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
 if(argPos==arguments->argc-1){
	fprintf(stderr,"Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atof(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}



double angsd::addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


double angsd::addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}


double angsd::getMax(double a,double b, double c){ 
    //get the maximum value of a, b and c
    double maxVal;// = std::max(a,std::max(b,c));
    if(a>b&&a>c)
      maxVal=a;
    else if(b>c)
      maxVal=b;
    else
      maxVal=c;
    return maxVal;
}




int angsd::fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}



angsd::Matrix<double> angsd::getMatrix(const char *name,int doBinary,int lens){
  if(!angsd::fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t";
  std::ifstream pFile(name,std::ios::in);
  
  char buffer[lens];
  std::list<double *> rows;
  int ncols =0;
  while(!pFile.eof()){
    pFile.getline(buffer,lens);
    if(strlen(buffer)==0)
      continue;
    char *tok = strtok(buffer,delims);
    std::list<double> items;
    while(tok!=NULL){
      if(doBinary)
	items.push_back(atoi(tok));
      else
	items.push_back(atof(tok));
      tok = strtok(NULL,delims);
    }
    //fprintf(stderr,"[%s] ncols:%lu\n",__FUNCTION__,items.size());
    ncols = items.size();
    double *drows = new double[items.size()];
    int i=0;
    for(std::list<double>::iterator it=items.begin();it!=items.end();it++)
      drows[i++]  = *it;
    rows.push_back(drows);
    
  }
  //  fprintf(stderr,"%s nrows:%lu\n",__FUNCTION__,rows.size());
  double **data = new double*[rows.size()];
  int i=0;
  for(std::list<double*>::iterator it=rows.begin();it!=rows.end();it++)
    data[i++]  = *it;
  
  Matrix<double> retMat;
  retMat.matrix=data;
  retMat.x = rows.size();
  retMat.y = ncols;
  return retMat;

}

void angsd::deleteMatrix(Matrix<double> mat){
  assert(mat.matrix!=NULL);
  for(int i=0;i<mat.x;i++)
    delete [] mat.matrix[i];
  delete[] mat.matrix;
  mat.matrix =NULL;
}


void angsd::printMatrix(Matrix<double> mat,FILE *file){
  fprintf(stderr,"Printing mat:%p with dim=(%d,%d)\n",mat.matrix,mat.x,mat.y);
  for(int xi=0;xi<mat.x;xi++){
    for(int yi=0;yi<mat.y;yi++)
      fprintf(file,"%f\t",mat.matrix[xi][yi]);
    fprintf(file,"\n");
  }    
}


double **angsd::get3likes(funkyPars *pars){

  double **loglike = NULL;
  loglike = new double*[pars->numSites]; 
  for(int s=0;s<pars->numSites;s++)
    loglike[s] = new double[3*pars->nInd];
  
  for(int s=0;s<pars->numSites;s++){

    if(pars->keepSites[s]==0)
      continue;
    for(int i=0;i<pars->nInd;i++){
      
      //fprintf(stderr,"mm: %d\t%d\n",pars->major[s],pars->major[s]);
      //fprintf(stderr,"%s\t%d\t%c\t%c\t",pars->sites[s].chromo,pars->sites[s].position+1,intToRef[pars->major[s]],intToRef[pars->minor[s]]);

      loglike[s][i*3+0]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]];
      loglike[s][i*3+1]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]];
      loglike[s][i*3+2]=pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]];
    }
  }
  return loglike;

}

double **angsd::get3likesRescale(funkyPars *pars){

  double **loglike = NULL;
  loglike = new double*[pars->numSites]; 
  for(int s=0;s<pars->numSites;s++)
    loglike[s] = new double[3*pars->nInd];
  
  for(int s=0;s<pars->numSites;s++){

    if(pars->keepSites[s]==0)
      continue;
    for(int i=0;i<pars->nInd;i++){
      
      //fprintf(stderr,"mm: %d\t%d\n",pars->major[s],pars->major[s]);
      //fprintf(stderr,"%s\t%d\t%c\t%c\t",pars->sites[s].chromo,pars->sites[s].position+1,intToRef[pars->major[s]],intToRef[pars->minor[s]]);

      loglike[s][i*3+0]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]];
      loglike[s][i*3+1]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]];
      loglike[s][i*3+2]=pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]];
      double mmax = loglike[s][i*3+0];
      for(int ii=1;ii<3;ii++)
	if(loglike[s][i*3+ii]>mmax)
	  mmax = loglike[s][i*3+ii];
      for(int ii=0;ii<3;ii++)
	loglike[s][i*3+ii] -=mmax;
    }
  }
  return loglike;

}


double **angsd::get3likes(funkyPars *pars,int *keepInd){
 
  int nKeep=0;
  for(int i=0;i<pars->nInd;i++){
    if(keepInd[i])
      nKeep++;
  }
 
  double **loglike = NULL;
  loglike = new double*[pars->numSites]; 
 
  for(int s=0;s<pars->numSites;s++){
     loglike[s] = new double[3*nKeep];
  }


  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)//always extract this, to avoid problems in multitrheading
      continue;
    int count=0;
    

    for(int i=0;i<pars->nInd;i++){
      if(keepInd[i]==0)
	continue;
    
      loglike[s][count*3+0]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]];
      loglike[s][count*3+1]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]];
      loglike[s][count*3+2]=pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]];
      count++;
    }
  }
 
  return loglike;

}

double **angsd::getlikes(funkyPars *pars,int *keepInd){

  int nKeep=0;
  for(int i=0;i<pars->nInd;i++){
    if(keepInd[i])
      nKeep++;
  }

  double **loglike = NULL;
  loglike = new double*[pars->numSites]; 
  for(int s=0;s<pars->numSites;s++)
    loglike[s] = new double[10*nKeep];
  
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    int count=0;
    for(int i=0;i<pars->nInd;i++){
      if(keepInd[i]==0)
	continue;
      for(int g=0;g<10;g++)
	loglike[s][count*10+g]=pars->likes[s][i*10+g];
      count++;
    }
  }
  return loglike;

}

void angsd::swapDouble (double& first, double& second)
{
        double temp = first;
        first = second;
        second = temp;
}

int angsd::matinv( double x[], int n, int m, double space[])
{
  //from rasmus nielsens code
  /* x[n*m]  ... m>=n*/
  register int i,j,k; 
  int *irow=(int*) space;
  double ee=1.0e-20, t,t1,xmax;
  double det=1.0;
  
  FOR (i,n)  {
    xmax = 0.;
    for (j=i; j<n; j++) {
      if (xmax < fabs(x[j*m+i]))  {
	xmax = fabs( x[j*m+i] );
	irow[i] = j;
      }
    }
    det *= xmax;
    if (xmax < ee)   {
      fprintf(stderr,"\nDeterminant becomes zero at %3d!\t\n", i+1);
      return(-1);
    }
    if (irow[i] != i) {
      FOR (j,m) {
	t = x[i*m+j];
	x[i*m+j] = x[irow[i] * m + j];
	x[ irow[i] * m + j] = t;
      }
    }
    t = 1./x[i*m+i];
    FOR (j,n) {
      if (j == i) continue;
      t1 = t*x[j*m+i];
      FOR(k,m)  x[j*m+k] -= t1*x[i*m+k];
      x[j*m+i] = -t1;
    }
    FOR(j,m)   x[i*m+j] *= t;
    x[i*m+i] = t;
  }                            /* i  */
  for (i=n-1; i>=0; i--) {
    if (irow[i] == i) continue;
    FOR(j,n)  {
      t = x[j*m+i];
      x[j*m+i] = x[ j*m + irow[i] ];
      x[ j*m + irow[i] ] = t;
    }
  }
  return (0);
}



void angsd::logrescale(double *ary,int len){
  int maxId = 0;
  //    fprintf(stderr,"maxid:%d maxval:%f\n",maxId,ary[maxId]);
  for(int i=1;i<len;i++){
    //if(posCounter==debug_print)
    //fprintf(stderr,"maxid:%d maxval:%f ary=%f\n",maxId,ary[maxId],ary[i]);
    if(ary[i]>ary[maxId])
      maxId=i;
  }
  double maxVal = ary[maxId];
  //if(posCounter==debug_print)
  //  fprintf(stderr,"maxval: %f\n",maxVal);
  for(int i=0;i<len;i++){
    //if(posCounter==debug_print)
    //fprintf(stderr,"%f\t%f\n",ary[i],ary[i]-maxVal);
    ary[i] = ary[i]-maxVal;
  }
  //  exit(0);
}


std::vector<char*> angsd::getFilenames(const char * name,int nInd){
 
  if(!fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(0);
  }
  const char* delims = " \t";
  std::vector<char*> ret;
  std::ifstream pFile(name,std::ios::in);

  char buffer[LENS];
  while(!pFile.eof()){
    pFile.getline(buffer,LENS);
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      if(tok[0]!='#')
	ret.push_back(strdup(buffer));
      tok = strtok(NULL,delims);
    }
  }
  if(nInd>0) {
     if(ret.size()<nInd)
      fprintf(stderr,"\t-> Number of samples is smaller than subset requested %lu vs %d\n",ret.size(),nInd);
    else{
      //   fprintf(stderr,"\t-> Will remove tail of filename list\n");
      for(int ii=nInd;ii<ret.size();ii++)
	free(ret[ii]);
      ret.erase(ret.begin()+nInd,ret.end());//we don't free this memory, it doesn't really matter
      // fprintf(stderr,"\t->  filename list now contains only: %lu\n",ret.size());
    }
#if 0
     for(size_t ii=0;ii<ret.size();ii++)
       fprintf(stderr,"%zu->%s\n",ii,ret[ii]);
     fprintf(stderr,"\n");
#endif
  }

  return ret;
}


void print_array(FILE *fp,double *ary,int len){
  for(int i=0;i<len-1;i++)
    fprintf(fp,"%f,",ary[i]);
  fprintf(fp,"%f\n",ary[len-1]);
}

void print_array(FILE *fp,int *ary,int len){
  for(int i=0;i<len-1;i++)
    fprintf(fp,"%d,",ary[i]);
  fprintf(fp,"%d\n",ary[len-1]);
}



double *angsd::readDouble(const char*fname,int hint){
  FILE *fp = NULL;
  fp = aio::getFILE(fname,"r");
  char buf[aio::fsize(fname)+1];
  if(aio::fsize(fname)!=fread(buf,sizeof(char),aio::fsize(fname),fp)){
    fprintf(stderr,"Problems reading file: %s\n will exit\n",fname);
    exit(0);
  }
  buf[aio::fsize(fname)]='\0';
  std::vector<double> res;
  res.push_back(atof(strtok(buf,"\t\n ")));
  char *tok=NULL;
  while((tok=strtok(NULL,"\t\n "))) {  
    //fprintf(stderr,"%s\n",tok);
    res.push_back(atof(tok));

  }
  //  fprintf(stderr,"size of prior=%lu\n",res.size());
  if(hint!=res.size()){
    fprintf(stderr,"problem with size of dimension of prior %d vs %lu\n",hint,res.size());
    for(size_t i=0;i<res.size();i++)
      fprintf(stderr,"%zu=%f\n",i,res[i]);
    exit(0);
  }
  double *ret = new double[res.size()];
  for(size_t i=0;i<res.size();i++)
    ret[i] = res[i];
  if(fp) fclose(fp);
  return ret;
}

int angsd::whichMax(double *d,int len){
  int r=0;
  for(int i=1;i<len;i++)
    if(d[i]>d[r])
      r=i;
  //now check if site doesnt have data.
  
  if(r==0){//only check if nothing is higher than the first
    for(int i=1;i<len;i++)
      if(d[i]!=d[0])//we see a diffrence so we have information
	return r;
    return -1;//we didnt have information 
  }else
    return r;
}


// a,c,g,t,n
// A,C,G,T,N
// 0,1,2,3,4
int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char intToRef[5] = {'A','C','G','T','N'};




int aio::fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

size_t aio::fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}

std::vector <char *> dumpedFiles;//small hack for getting a nice vector of outputfiles
FILE *aio::openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}
gzFile aio::openFileGz(const char* a,const char* b,const char *mode){
  if(0)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  gzFile fp = aio::getGz(c,mode);
  delete [] c;
  return fp;
}


BGZF *aio::openFileBG(const char* a,const char* b,const char *mode){
  if(0)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  BGZF *fp = bgzf_open(c,GZOPT);
  delete [] c;
  return fp;
}


FILE *aio::getFILE(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

gzFile aio::getGz(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  gzFile fp=Z_NULL;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}
