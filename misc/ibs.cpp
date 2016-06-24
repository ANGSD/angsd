// g++ -o ibs ibs.cpp -lz -o3

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <limits>
#include <zlib.h>
#include <vector>
#include <pthread.h>
#include <signal.h>
#include <vector>
#include <sys/stat.h>

typedef struct
{
  int nInd;
  int nSites;
  int totalSites;
  int *keepSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
  double pi[10];
} argu;


typedef struct
{
  int nInd1;
  int nInd2;
  int nSites;
  int totalSites;
  int *keepSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
  double pi[100];
} argu2;


void info(){
 fprintf(stderr,"Arguments:\n");
 

}


double addProtectN(double a[],int len){
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

std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);

  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

int readGLF(const char* fname,double * &gls,int nInd){
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  size_t sizeMax = 1000000;
  //  double *buf = new double[sizeMax];
  gls = new double[sizeMax];
  size_t temp = gzread(fp,gls,sizeof(double)*sizeMax);
  gzclose(fp);
  temp = temp/sizeof(double)/10/nInd;
  int hej = temp;
  return hej;

}



void runEM(double *gl,argu *pars){

  float sum;
  
  int maxIter=30;
  int it;
  double W[10];
  double Wtmp[10];
  double p[10];
  for(int w=0;w<10;w++)
    p[w] = 0.1;


  for(it=0;it<maxIter;it++){
    sum=0;
    for(int w=0;w<10;w++)
      W[w]=0;
    for(int i=0;i<pars->totalSites;i++){

      if(pars->keepSites[i]==0)
	continue;
      double sum=0;
      for(int w=0;w<10;w++){
	Wtmp[w]=exp(gl[i*10+w+pars->nInd*10*pars->totalSites])*p[w];
	sum += Wtmp[w];
      }
      for(int w=0;w<10;w++)
	W[w] += Wtmp[w]/sum;
    }
    for(int w=0;w<10;w++){
      p[w] = W[w]/pars->nSites;
      //      fprintf(stdout,"%f\t",p[w]);
    }
    //        fprintf(stdout,"\n");
  }
  for(int w=0;w<10;w++)
    pars->pi[w]=p[w];

  double like=0;
  for(int i=0;i<pars->totalSites;i++){
    double sum=0;
    for(int w=0;w<10;w++){
      sum += exp(gl[i*10+w])*p[w];
    }
     like+=log(sum);

  }
    
  pars->lres = like;

}




void runEM2D(double *gl,argu2 *pars){

  float sum;
  
  int maxIter=100;
  int it;
  double W[100];
  double Wtmp[100];
  double p[100];
  for(int w=0;w<100;w++)
    p[w] = 0.01;


  for(it=0;it<maxIter;it++){
    sum=0;
    for(int w=0;w<100;w++)
      W[w]=0;
    for(int i=0;i<pars->totalSites;i++){

      if(pars->keepSites[i]==0)
	continue;
      double sum=0;
      for(int w=0;w<100;w++){
	Wtmp[w]=exp(gl[i*10+w+pars->nInd1*10*pars->totalSites]+gl[i*10+w+pars->nInd2*10*pars->totalSites])*p[w];
	sum += Wtmp[w];
      }
      for(int w=0;w<100;w++)
	W[w] += Wtmp[w]/sum;
    }
    for(int w=0;w<100;w++){
      p[w] = W[w]/pars->nSites;
      //      fprintf(stdout,"%f\t",p[w]);
    }
    //        fprintf(stdout,"\n");
  }
  for(int w=0;w<100;w++)
    pars->pi[w]=p[w];

  double like=0;
  for(int i=0;i<pars->totalSites;i++){
    double sum=0;
    for(int w=0;w<100;w++){
      sum += exp(gl[i*10+w+pars->nInd1*10*pars->totalSites]+gl[i*10+w+pars->nInd2*10*pars->totalSites])*p[w];
    }
     like+=log(sum);

  }
    
  pars->lres = like;

}


int main(int argc, char **argv){
  if(argc==1){// if no arguments, print info on program
    info();
    return 0;
  }

  const char* likeFileName = NULL;
  const char* outFileName = NULL;
  int seed =0;
  int nInd=1;
  int p1=0;
  int p2=0;
  // reading arguments
  argv++;
  while(*argv){
    if(strcmp(*argv,"-glf")==0 || strcmp(*argv,"-f")==0) likeFileName=*++argv; 
    else if(strcmp(*argv,"-outFileName")==0 || strcmp(*argv,"-o")==0) outFileName=*++argv; 
    else if(strcmp(*argv,"-nInd")==0 || strcmp(*argv,"-n")==0) nInd=atoi(*++argv);
    else if(strcmp(*argv,"-ind2")==0||strcmp(*argv,"-p2")==0) p2=atoi(*++argv);
    else if(strcmp(*argv,"-ind1")==0||strcmp(*argv,"-p1")==0) p1=atoi(*++argv);
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }



  //
  if(likeFileName==NULL){
      fprintf(stderr,"Please supply glf file: -f");
      info();
  }
  if(outFileName==NULL){
    fprintf(stderr,"Will use glf file name as prefix for output\n");
    outFileName=likeFileName;
  }

  
  FILE *flog=openFile(outFileName,".log");

  double *genoLike=NULL;

  int totalSites=readGLF(likeFileName,genoLike,nInd);
  fprintf(stdout,"read %d sites\n",totalSites);



 
  argu * myPars= new argu;
  myPars->totalSites = totalSites;
  myPars->keepSites =  new int[totalSites];

 if(p1+p2==0){
  for(int theInd=0;theInd<nInd;theInd++){
    myPars->nInd = theInd;

    int nSites=0;
    for(int i=0;i<totalSites;i++){
      myPars->keepSites[i] = 0;
      for(int w=0;w<10;w++)
	if(-genoLike[i*10+w+myPars->nInd*10*totalSites]>myPars->keepSites[i]){
	  myPars->keepSites[i] = 1;
	  nSites++;
	  break;
	}
    }  
    myPars->nSites = nSites;
    
    
    runEM(genoLike,myPars);
    

    fprintf(stdout,"%d\t%d\t%f",theInd,myPars->nSites,myPars->lres);
    for(int w=0;w<10;w++){
      fprintf(stdout,"\t%f",myPars->pi[w]);
    }
    fprintf(stdout,"\n");
  }
 }

  if(p1+p2>0){
    argu2 * myPars2D= new argu2;
    myPars2D->totalSites = totalSites;
    myPars2D->keepSites =  new int[totalSites];

    myPars2D->nInd1 = p1;
    myPars2D->nInd2 = p2;

    int nSites=0;
    for(int i=0;i<totalSites;i++){
      myPars2D->keepSites[i] = 0;
      for(int w=0;w<10;w++)
	if(-genoLike[i*10+w+myPars2D->nInd1*10*totalSites]>myPars2D->keepSites[i] && -genoLike[i*10+w+myPars2D->nInd2*10*totalSites]>myPars2D->keepSites[i]){
	  myPars2D->keepSites[i] = 1;
	  nSites++;
	  break;
	}
    }  
    myPars2D->nSites = nSites;
  
    runEM2D(genoLike,myPars2D);
    

    fprintf(stdout,"%d\t%d\t%d\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
    for(int w=0;w<100;w++){
      fprintf(stdout,"\t%f",myPars2D->pi[w]);
    }
    fprintf(stdout,"\n");

    //    runEM2D(genoLike,myPars2D);
    
  }

  ////////// clean  
  fclose(flog); 
  for(int i=0;i<dumpedFiles.size();i++){
    //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }

  delete[] genoLike;
  return 0;

}
