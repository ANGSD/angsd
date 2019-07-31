


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
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      //      fprintf(stderr,"HIT %s vs %s\n",argName,arguments->argv[argPos]);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
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
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
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
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
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
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
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
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
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

double angsd::addProtectN(double a[],int len){
  //function does: log(sum(exp(a))) while protecting for underflow
  double maxVal = a[0];

  for(int i=1;i<10;i++)
    if(maxVal<a[i])
      maxVal=a[i];

  double sumVal = 0;
  for(int i=1;i<10;i++)
    sumVal += exp(a[i]-maxVal);

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
angsd::Matrix<int> angsd::getMatrixInt(const char *name,int lens){
  if(!angsd::fexists(name)){
    fprintf(stderr,"\t-> Problems opening file: %s\n",name);
    exit(0);
  }
  const char* delims = " \t";
  std::ifstream pFile(name,std::ios::in);
  
  char buffer[lens];
  std::list<int *> rows;
  int ncols =0;
  while(!pFile.eof()){
    pFile.getline(buffer,lens);
    if(strlen(buffer)==0)
      continue;
    char *tok = strtok(buffer,delims);
    std::list<int> items;
    while(tok!=NULL){
	items.push_back(atoi(tok));
      tok = strtok(NULL,delims);
    }
    //fprintf(stderr,"[%s] ncols:%lu\n",__FUNCTION__,items.size());
    ncols = items.size();
    int *drows = new int[items.size()];
    int i=0;
    for(std::list<int>::iterator it=items.begin();it!=items.end();it++)
      drows[i++]  = *it;
    rows.push_back(drows);
    
  }
  //  fprintf(stderr,"%s nrows:%lu\n",__FUNCTION__,rows.size());
  int **data = new int*[rows.size()];
  int i=0;
  for(std::list<int*>::iterator it=rows.begin();it!=rows.end();it++)
    data[i++]  = *it;
  
  Matrix<int> retMat;
  retMat.matrix=data;
  retMat.x = rows.size();
  retMat.y = ncols;
  return retMat;

}

void angsd::deleteMatrixInt(Matrix<int> mat){
  assert(mat.matrix!=NULL);
  for(int i=0;i<mat.x;i++)
    delete [] mat.matrix[i];
  delete[] mat.matrix;
  mat.matrix =NULL;
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
      //      fprintf(stderr,"refid: %d posi:%d pars->major:%p\n",pars->refId,pars->posi[s]+1,pars->major);
      //fprintf(stderr,"mm: %d\t%d\n",pars->major[s],pars->major[s]);
      //   fprintf(stderr,"%d\t%d\t%c\t%c\t",pars->refId,pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]]);
      
      loglike[s][i*3+0]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]];
      loglike[s][i*3+1]=pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]];
      loglike[s][i*3+2]=pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]];
      double mmax = loglike[s][i*3+0];
      for(int ii=1;ii<3;ii++)
	if(loglike[s][i*3+ii]>mmax)
	  mmax = loglike[s][i*3+ii];
      for(int ii=0;(!std::isinf(mmax))&&ii<3;ii++){
	loglike[s][i*3+ii] -=mmax;
	if(std::isnan(loglike[s][i*3+ii])){
	  fprintf(stderr,"mmax: %f\n",mmax);
	  exit(0);
	}
      }
    }
  }
  return loglike;

}


double **angsd::get3likesRMlow(funkyPars *pars,int *keepInd){
 
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
   
      if(loglike[s][count*3+0] < -20 && loglike[s][count*3+1] < -20 && loglike[s][count*3+2] < -20){
	loglike[s][count*3+0] = 0;
	loglike[s][count*3+1] = 0;
	loglike[s][count*3+2] = 0;

      }
      count++;

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
//DRAGON just use std::swap
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
  for(int i=1;i<len;i++)
    if(ary[i]>ary[maxId])
      maxId=i;
  
  double maxVal = ary[maxId];
  for(int i=0;i<len;i++)
    ary[i] -= maxVal;
  
  
}


std::vector<char*> angsd::getFilenames(const char * name,int nInd){
  if(strchr(name,'\r')){
    fprintf(stderr,"\t\t-> Filelist contains carriage return. Looks like a windows file please remove hidden \'\r\' from filelist\n");
    exit(0);
  }
  
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


double angsd::sigm(double x){
  return(1/(1+exp(-x)));
}


double angsd::lbico(double n, double k){
  return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);
}

double angsd::myComb2(int k,int r, int j){
  if(j>r)
    fprintf(stderr,"%s error in k=%d r=%d j=%d\n",__FUNCTION__,k,r,j);

  double fac1= lbico(r,j)+lbico(2*k-r,2-j);
  double fac2=lbico(2*k,2);
  
  return exp(fac1-fac2);
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
    fprintf(stderr,"\t-> File: \'%s\' should contain %d values, but has %lu\n",fname,hint,res.size());
    fprintf(stderr,"\t-> If you are supplying an estimated sfs, make sure your input file is a single line (an estimate for a single region)\n");
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


//count is 5 long, A C G T N
int angsd::getRandomCount(suint *counts, int i,int depth){

  if(depth==-1){
    depth=0;
    for( int b = 0; b < 4; b++ )
      depth+=counts[b+4*i];
  }

  if(depth==0)
    return 4;

  int j = std::rand() % depth;
  int cumSum=0;
  int res=4;

  for( int b = 0; b < 4; b++ ){
    cumSum+=counts[b+4*i];
    if( cumSum > j ){
      res = b;
      break;
    }
  }
  return res;
}

// get the most frequent base, use random for tie
// depth is without N
// i is the individual
int angsd::getMaxCount(suint *counts,int i, int depth){

  if(depth==-1){
    depth=0;
    for( int b = 0; b < 4; b++ )
      depth+=counts[b+4*i];
  }

  if(depth<=0)
    return 4;

  int whichMax = 0;
  int nMax=1;  
  for(int b=1;b<4;b++){
    if (counts[b+4*i]>counts[whichMax+4*i]){
      whichMax = b;
      nMax = 1;
    }
    else if(counts[b+4*i]==counts[whichMax+4*i]){
      nMax++;
    }
  }

  if(nMax>1){ // in case of ties
    int j=0;
    int r = std::rand() % nMax;
    for(int b=1;b<4;b++){
      if(counts[b+4*i]==counts[whichMax+4*i]){
	if(r==j){
	  whichMax=b;
	  break;
	}
	j++;
      }
    }     
  }
  

  return whichMax;
}

// combine all bases to get IUPAC code
// depth is without N
// i is the individual
int angsd::getIupacCount(suint *counts,int i, double iRatio, int depth){

  if(depth==-1){
    depth=0;
    for( int b = 0; b < 4; b++ )
      depth+=counts[b+4*i];
  }

  if(depth<=0)
    return 14;

  int whichIUPAC = 0;
  double bIUPACscore = 0;
  for(int b=0;b<4;b++){
    if (double(counts[b+4*i])/double(depth)>iRatio){
      bIUPACscore = bIUPACscore + pow(b+1,2);
    }
  }
  //N
  if(bIUPACscore == 0){
    whichIUPAC = 14;
  }
  //A
  if(bIUPACscore == 1){
    whichIUPAC = 0;
  }
  //C
  if(bIUPACscore == 4){
    whichIUPAC = 1;
  }
  //G
  if(bIUPACscore == 9){
    whichIUPAC = 2;
  }
  //T
  if(bIUPACscore == 16){
    whichIUPAC = 3;
  }
  //A+G
  if(bIUPACscore == 10){
    whichIUPAC = 4;
  }
  //C+T
  if(bIUPACscore == 20){
    whichIUPAC = 5;
  }
  //G+C
  if(bIUPACscore == 13){
    whichIUPAC = 6;
  }
  //A+T
  if(bIUPACscore == 17){
    whichIUPAC = 7;
  }
  //G+T
  if(bIUPACscore == 25){
    whichIUPAC = 8;
  }
  //A+C
  if(bIUPACscore == 5){
    whichIUPAC = 9;
  }
  //C+G+T
  if(bIUPACscore == 29){
    whichIUPAC = 10;
  }
  //A+G+T
  if(bIUPACscore == 26){
    whichIUPAC = 11;
  }
  //A+C+T
  if(bIUPACscore == 21){
    whichIUPAC = 12;
  }
  //A+C+G
  if(bIUPACscore == 14){
    whichIUPAC = 13;
  }
  //A+C+G+T
  if(bIUPACscore == 30){
    whichIUPAC = 14;
  }
  return whichIUPAC;
}

//count is 4 long, A C G T
int angsd::getRandomCountTotal(suint *counts, int nInd){

  size_t totalCounts[4]={0,0,0,0};
  for(int i=0;i<4*nInd;i++)
    totalCounts[i%4] +=counts[i];   
  

  size_t depth=0;
  for( int b = 0; b < 4; b++ )
    depth+=totalCounts[b];
  

  if(depth==0)
    return 4;

  size_t j = std::rand() % depth;
  size_t cumSum=0;
  int res=4;

  for( int b = 0; b < 4; b++ ){
    cumSum+=totalCounts[b];
    if( cumSum > j ){
      res = b;
      break;
    }
  }
  return res;
}

// get the most frequent base, use random for tie
// depth is without N
// i is the individual
int angsd::getMaxCountTotal(suint *counts,int nInd){

  size_t totalCounts[4]={0,0,0,0};
  for(int i=0;i<4*nInd;i++)
    totalCounts[i%4] +=counts[i];   
  

  size_t depth=0;
  for( int b = 0; b < 4; b++ )
    depth+=totalCounts[b];
  

  if(depth==0)
    return 4;

  int whichMax = 0;
  int nMax=1;  
  for(int b=1;b<4;b++){
    if ( totalCounts[b] > totalCounts[whichMax] ){
      whichMax = b;
      nMax = 1;
    }
    else if( totalCounts[b] == totalCounts[whichMax] ){
      nMax++;
    }
  }

  if(nMax>1){ // in case of ties
    int j=0;
    int r = std::rand() % nMax;
     for(int b=1;b<4;b++){
       if( totalCounts[b] == totalCounts[whichMax] ){
	 if(r==j){
	   whichMax=b;
	   break;
	 }
	 j++;
       }

     }
      
  }


  return whichMax;
}

// combine all bases to get IUPAC code
// depth is without N
// i is the individual
int angsd::getIupacCountTotal(suint *counts,int nInd, double iRatio){

  size_t totalCounts[4]={0,0,0,0};
  for(int i=0;i<4*nInd;i++)
    totalCounts[i%4] +=counts[i];   
  

  size_t depth=0;
  for( int b = 0; b < 4; b++ )
    depth+=totalCounts[b];
  

  if(depth==0)
    return 14;

  int whichIUPAC = 0;
  double bIUPACscore = 0;
  for(int b=0;b<4;b++){
    if (double(totalCounts[b])/double(depth)>iRatio){
      bIUPACscore = bIUPACscore + pow(b+1,2);
    }
  }
  //N
  if(bIUPACscore == 0){
    whichIUPAC = 14;
  }
  //A
  if(bIUPACscore == 1){
    whichIUPAC = 0;
  }
  //C
  if(bIUPACscore == 4){
    whichIUPAC = 1;
  }
  //G
  if(bIUPACscore == 9){
    whichIUPAC = 2;
  }
  //T
  if(bIUPACscore == 16){
    whichIUPAC = 3;
  }
  //A+G
  if(bIUPACscore == 10){
    whichIUPAC = 4;
  }
  //C+T
  if(bIUPACscore == 20){
    whichIUPAC = 5;
  }
  //G+C
  if(bIUPACscore == 13){
    whichIUPAC = 6;
  }
  //A+T
  if(bIUPACscore == 17){
    whichIUPAC = 7;
  }
  //G+T
  if(bIUPACscore == 25){
    whichIUPAC = 8;
  }
  //A+C
  if(bIUPACscore == 5){
    whichIUPAC = 9;
  }
  //C+G+T
  if(bIUPACscore == 29){
    whichIUPAC = 10;
  }
  //A+G+T
  if(bIUPACscore == 26){
    whichIUPAC = 11;
  }
  //A+C+T
  if(bIUPACscore == 21){
    whichIUPAC = 12;
  }
  //A+C+G
  if(bIUPACscore == 14){
    whichIUPAC = 13;
  }
  //A+C+G+T
  if(bIUPACscore == 30){
    whichIUPAC = 14;
  }
  return whichIUPAC;
}




int ludcmp(double **a, int *indx, double &d,int n)
{
  int imax = 0;
  double big, dum, sum, temp;
  double vv[n];
  d=1;

  for (int i=0; i<n; i++){
    big=0;
    for (int j=0; j<n; j++){
      //fprintf(stderr,"%f\t",a[i][j]);
      if ((temp=fabs(a[i][j])) > big) 
	big=temp;
    }
    if(big==0){
      //fprintf(stderr,"singular matrix in ludcmp");
      return(1);
	//    assert(big!=0) ;

    }
    vv[i]=1/big;
  }

  for (int j=0; j<n; j++){
    for (int i=0; i<j; i++){
      sum = a[i][j];
      for (int k=0; k<i; k++) 
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
    }
    big=0;
    for (int i=j; i<n; i++)	{
      sum=a[i][j];
      for (int k=0; k<j; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax){
      for (int k=0; k<n; k++){
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0) 
      a[j][j] = 1.0e-20;
    if (j != n-1){
      dum = 1/(a[j][j]);
      for (int i=j+1; i<n; i++) 
	a[i][j] *= dum;
    }
  }
  return 0;
}


void lubksb(double **a, int *indx, double *b,int n)
{

  int ii=0;
  double sum;

  for (int i=0; i<n; i++){
    int ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii != 0)
      for (int j=ii-1; j<i; j++) 
	sum -= a[i][j]*b[j];
    else if (sum != 0.0) 
      ii=i+1;
    b[i]=sum;
  }
  for (int i=n-1; i>=0; i--){
    sum=b[i];
    for (int j=i+1; j<n; j++) 
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

//usefull little function to split
char *angsd::strpop(char **str,char split){
  char *tok=*str;
  while(**str){
    if(**str!=split)
      (*str)++;
    else{
      **str='\0'; (*str)++;
      break;
    }
  }
  return tok;
}




int angsd::svd_inverse(double mat[],int xLen, int yLen){
  if(xLen !=yLen){

    fprintf(stderr,"non square matrix [%s]\t[%s]\n",__FILE__,__FUNCTION__);
    exit(0);

  }
  double *col;
  double y[xLen * yLen];
  col = new double[xLen];
  double **tm;
  int *indx=new int[xLen];
  double d;
  tm = new double*[xLen];
  for (int i=0; i < xLen; i++)
    tm[i] = new double[xLen];

  for(int i=0;i<xLen;i++)
    for(int j=0;j<yLen;j++)
      tm[i][j]=mat[j*xLen+i];


  int singular=ludcmp(tm,indx,d,xLen);
  if(singular)
    return 1 ;
  
  for (int j=0; j<xLen; j++)
    {
      for (int i=0; i<xLen; i++)
	col[i]=0;
      col[j]=1;
      lubksb(tm,indx,col,xLen);
      for (int i=0; i<xLen; i++) 
	y[j*xLen+i]=col[i];
    }
  
  
  for (int j=0; j<yLen; j++)
    for (int i=0; i<xLen; i++)
      mat[j*xLen+i]=y[j*xLen+i];

  delete[] col;
  delete[] indx;
  for (int i=0; i < xLen; i++)
    delete[] tm[i];
  delete[] tm;
  return 0;
}



//function for getting density of normal distribution, has safeguards against underflow
//by emil added 24-11-2018
double angsd::dnorm(double x,double mean,double sd,int ifLog){

  double fac = 1.0/(sd*sqrt(2.0*M_PI));
  double val = exp(-(((x-mean)*(x-mean))/(2*sd*sd)));

  double lower_bound=1e-20;//emil - not for users

  if(ifLog){    
    if(val<lower_bound){
      return(log(lower_bound));
    } else{
      return (log(fac)+log(val));
    }    
  } else{
    // if val is 0 because exp(-(x-mean)*(x-mean)) is due to underflow, returns low value
    if(val<lower_bound){      
      return(lower_bound);
    } else{
      return fac*val;
    }
  }
  
}

//function for getting probability of bernoulli distribution, has safeguards against underflow
//by emil added 24-11-2018
double angsd::bernoulli(int k, double p, int ifLog){
  // if p is 0 or 1, cannot do log
  // however this because of over/underlow and p i just very close 0 or 1
  double lower_bound=1e-20;//emil - not for users
  
  if(p>1-lower_bound){
    p = 1-lower_bound;
  } else if(p<lower_bound){
    p = lower_bound;
  }
  
  if(ifLog){
    return( log(pow(p,k)*pow(1-p,1-k)) );
  } else{
    return( pow(p,k)*pow(1-p,1-k) );
  }
}



// function for getting standard derivation of a set of data
//by emil added 24-11-2018
double angsd::sd(double* phe, int size ){
  double ts = 0;
  for(int i=0;i<size;i++)
    ts += phe[i];
  double u = ts/(1.0*size);
  ts = 0;
  for(int i=0;i<size;i++)
    ts += (phe[i]-u)*(phe[i]-u);
  return ts/(1.0*(size-1.0));
}

double angsd::to_pval(Chisqdist *chisq,double f){
  return f<0?1:1-chisq->cdf(f);
}

// function for getting density of lambda function
//by emil added 12-04-2019
// from: http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
double angsd::poisson(double k,  double lambda, int ifLog) {
  if(ifLog){
    return (k * log(lambda) - lgamma(k + 1.0) - lambda);
  } else{
    return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
  }
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

char intToIupac[15] = {'A','C','G','T','R','Y','S','W','K','M','B','D','H','V','N'};


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
  FILE *fp = NULL;
  fp = fopen(c,"w");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\' check permissions\n",c);
    exit(0);
  }
  delete [] c;
  return fp;
}

BGZF *aio::openFileBG(const char* a,const char* b){

  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  dumpedFiles.push_back(strdup(c));
  BGZF *fp = bgzf_open(c,GZOPT);
  delete [] c;
  return fp;
}

htsFile *aio::openFileHts(const char* a,const char* b){

  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  dumpedFiles.push_back(strdup(c));
  htsFile *fp = hts_open(c,"w");
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

//checks that newer is newer than older
int aio::isNewer(const char *newer,const char *older){
   if (strstr(older, "ftp://") == older || strstr(older, "http://") == older)
     return 0;
  //  fprintf(stderr,"newer:%s older:%s\n",newer,older);
  // return 0;
  struct stat one;
  struct stat two;
  stat(newer, &one );
  stat(older, &two );
  
  return one.st_mtime>=two.st_mtime;
}

ssize_t aio::bgzf_write(BGZF *fp, const void *data, size_t length){
  if(length>0)
    return ::bgzf_write(fp,data,length);
  return 0;
}

void angsd::norm(double *d,size_t len){
  double ts=0;
  for(int i=0;i<len;i++)
    ts += d[i];

  for(int i=0;i<len;i++)
    d[i] /= ts;

}

int aio::tgets(gzFile gz,char**buf,int *l){
  int rlen = 0;
 neverUseGoto:
  char *tok = gzgets(gz,*buf+rlen,*l-rlen);
  if(!tok)
    return rlen;
  int tmp = tok?strlen(tok):0;
  if(tok[tmp-1]!='\n'){
    rlen += tmp;
    *l *= 2;
    *buf = (char*) realloc(*buf,*l);
    goto neverUseGoto;
  }
  rlen += tmp;
  return rlen;
}






// em freqeuncy assuming HWE
double angsd::estFreq(double *loglike,int numInds){

  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p= 0.1;
  float temp_p=p;
  double accu=0.00001;
  double accu2=0;
  float sum;
  int iter=100;

  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
     
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/numInds;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }

  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",numInds);

    for(int ii=0;1&&ii<numInds;ii++){
    
      fprintf(stderr,"1\t");
      for(int gg=0;gg<3;gg++)
	fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
     
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    //exit(0);
  }
  
  return(p);
}

int countsSample(suint *a){
  double r = drand48()*(a[0]+a[1]+a[2]+a[3]);
  if(r<a[0])
    return 0;
  else if(r>=a[0] &&r<a[1])
    return 1;
  else if(r>=a[1] &&r<a[2])
    return 2;
  else
    return 3;


}








//public domain from here http://www.johndcook.com/cpp_phi.html
double phi(double x){
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}


