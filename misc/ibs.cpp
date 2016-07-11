// g++ -o ibs ibs.cpp -lz -O3

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

int model;

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
 fprintf(stderr,"Needed arguments:\n");
 fprintf(stderr,"\t-glf/-f\t\tinput GLF filename:\n");
 fprintf(stderr,"Optional arguments:\n");
 fprintf(stderr,"\t-outFileName/-o\toutput filename(prefix):\n");
 fprintf(stderr,"\t-nInd/-n\tnubmer of individuals in GLF file:\n");
 fprintf(stderr,"\t-ind1/i1\tindividuals 1:\n");
 fprintf(stderr,"\t-ind2/i2\tindividuals 2:\n");
 fprintf(stderr,"\t-allpairs/-a\tanalyse all pairs:\n");
 fprintf(stderr,"\t-maxSites/-m\tmaximum sites to analyze:\n");
 fprintf(stderr,"\t-model\t\tibs model:0 all 10 genotypes, 1 HO/HE\n");
 // fprintf(stderr,"\t-\t\t:\n");


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

int readGLF(const char* fname,double * &gls,int nInd,int maxSites){
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  int sizeMax = 100000*nInd;
  size_t nSites = 0;

  //  double *buf = new double[sizeMax];
  double *tmp_gls = new double[sizeMax];

  int temp=0;
  while( temp = gzread(fp,tmp_gls,sizeof(double)*sizeMax)){

    //    fprintf(stdout,"%d\n",temp);
    nSites += temp/sizeof(double)/10/nInd;
    if(nSites>maxSites){
      nSites = maxSites;
      break;
    }
  }

  gzclose(fp);
  delete[] tmp_gls;


  gls = new double[nSites*10*nInd];
  fp=gzopen(fname,"r");
  size_t nSites2 = gzread(fp,gls,sizeof(double)*nSites*10*nInd);
  gzclose(fp);
  nSites2 = nSites2/sizeof(double)/10/nInd;
  if(nSites!=nSites2){
    fprintf(stdout,"reading failed nSites%lu != nSites2 %lu\n",nSites,nSites2);   
  }
  return nSites2;

}



void runEM(double *gl,argu *pars){

  float sum;
  
  int maxIter=30;
  int it;
  double W[10];
  double Wtmp[10];
  double p[10];
  double pold[10];
  for(int w=0;w<10;w++)
    p[w] = 0.1;

  double tol = 0.0000001;


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

    if(model==0){ //full model. all 10 genotypes
      for(int w=0;w<10;w++)
	pold[w] = W[w]/pars->nSites;

    }
    else if(model==1){// 2 genotype model (HO/HE)
      double HO=W[0]+W[4]+W[7]+W[9];
      int w=0;
      for(int a1=0;a1<4;a1++)
	for(int a2=a1;a2<4;a2++){
	  if(a1==a2)
	    pold[w] = HO/pars->nSites/4;
	  else
	    pold[w] =(1-HO/pars->nSites)/6;
	  w++;
	}
    }
    double diff=0;
    for(int w=0;w<10;w++){
      double d = fabs(p[w]-pold[w]);
      if(d>diff)
	diff=d;
      p[w] = pold[w];
    }
    if(diff<tol){
      fprintf(stdout,"tol reached %d\t%f\n",it,diff);
      break;
    }

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
  int p1=-1;
  int p2=-1;
  int maxSites=100000;
  int all=0;
  int model=0;
  // reading arguments
  argv++;
  while(*argv){
    if(strcmp(*argv,"-glf")==0 || strcmp(*argv,"-f")==0) likeFileName=*++argv; 
    else if(strcmp(*argv,"-outFileName")==0 || strcmp(*argv,"-o")==0) outFileName=*++argv; 
    else if(strcmp(*argv,"-nInd")==0 || strcmp(*argv,"-n")==0) nInd=atoi(*++argv);
    else if(strcmp(*argv,"-ind1")==0||strcmp(*argv,"-i1")==0) p1=atoi(*++argv);
    else if(strcmp(*argv,"-ind2")==0||strcmp(*argv,"-i2")==0) p2=atoi(*++argv);
    else if(strcmp(*argv,"-allpairs")==0||strcmp(*argv,"-a")==0) all=atoi(*++argv);
    else if(strcmp(*argv,"-maxSites")==0||strcmp(*argv,"-m")==0) maxSites=atoi(*++argv);
    else if(strcmp(*argv,"-model")==0) model=atoi(*++argv);
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

  FILE *fibs;
  FILE *fibspair;
  FILE *flog=openFile(outFileName,".log");
  if(p2==-1 && all!=1)
    fibs=openFile(outFileName,".ibs");
  if((p1>=0 && p2>=0) || all)
    fibspair=openFile(outFileName,".ibspair");

  double *genoLike=NULL;

  int totalSites=readGLF(likeFileName,genoLike,nInd,maxSites);
  fprintf(stdout,"read %d sites\n",totalSites);
  fprintf(flog,"read %d sites\n",totalSites);



 
  argu * myPars= new argu;
  myPars->totalSites = totalSites;
  myPars->keepSites =  new int[totalSites];

  if(p2==-1 && all!=1){
    //header of ibs files
    fprintf(fibs,"ind\tnSites\tLlike\tpAA\tpAC\tpAG\tpAT\tpCC\tpCG\tpCT\tpGG\tpGT\tpTT\n");

    

    for(int theInd=0;theInd<nInd;theInd++){
      
      if(p1!=-1 && p1!=theInd)
	  continue;
      fprintf(stdout,"analysing individual %d\n",theInd);
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
      
      
      fprintf(fibs,"%d\t%d\t%f",theInd,myPars->nSites,myPars->lres);
      for(int w=0;w<10;w++){
	fprintf(fibs,"\t%f",myPars->pi[w]);
      }
      fprintf(fibs,"\n");
    }
  }

 if(p1>=0 && p2>=0 && all==0){


   fprintf(fibspair,"ind1\tind2\tnSites\tLlike\t");
   //paste(paste0("p",paste(rep(GENO,10),rep(GENO,each=10),sep="_")),collapse="\t")
   fprintf(fibspair,"pAA_AA\tpAC_AA\tpAG_AA\tpAT_AA\tpCC_AA\tpCG_AA\tpCT_AA\tpGG_AA\tpGT_AA\tpTT_AA\tpAA_AC\tpAC_AC\tpAG_AC\tpAT_AC\tpCC_AC\tpCG_AC\tpCT_AC\tpGG_AC\tpGT_AC\tpTT_AC\tpAA_AG\tpAC_AG\tpAG_AG\tpAT_AG\tpCC_AG\tpCG_AG\tpCT_AG\tpGG_AG\tpGT_AG\tpTT_AG\tpAA_AT\tpAC_AT\tpAG_AT\tpAT_AT\tpCC_AT\tpCG_AT\tpCT_AT\tpGG_AT\tpGT_AT\tpTT_AT\tpAA_CC\tpAC_CC\tpAG_CC\tpAT_CC\tpCC_CC\tpCG_CC\tpCT_CC\tpGG_CC\tpGT_CC\tpTT_CC\tpAA_CG\tpAC_CG\tpAG_CG\tpAT_CG\tpCC_CG\tpCG_CG\tpCT_CG\tpGG_CG\tpGT_CG\tpTT_CG\tpAA_CT\tpAC_CT\tpAG_CT\tpAT_CT\tpCC_CT\tpCG_CT\tpCT_CT\tpGG_CT\tpGT_CT\tpTT_CT\tpAA_GG\tpAC_GG\tpAG_GG\tpAT_GG\tpCC_GG\tpCG_GG\tpCT_GG\tpGG_GG\tpGT_GG\tpTT_GG\tpAA_GT\tpAC_GT\tpAG_GT\tpAT_GT\tpCC_GT\tpCG_GT\tpCT_GT\tpGG_GT\tpGT_GT\tpTT_GT\tpAA_TT\tpAC_TT\tpAG_TT\tpAT_TT\tpCC_TT\tpCG_TT\tpCT_TT\tpGG_TT\tpGT_TT\tpTT_TT\n");



   argu2 * myPars2D= new argu2;
   myPars2D->totalSites = totalSites;
   myPars2D->keepSites =  new int[totalSites];
   
   myPars2D->nInd1 = p1;
   myPars2D->nInd2 = p2;
   fprintf(stdout,"analysing pair %d %d\n",p1,p2);

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
   
   
   fprintf(fibspair,"%d\t%d\t%d\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
   for(int w=0;w<100;w++){
     fprintf(fibspair,"\t%f",myPars2D->pi[w]);
   }
   fprintf(fibspair,"\n");
   
   //    runEM2D(genoLike,myPars2D);
   
 }
 else if(all){
   fprintf(fibspair,"ind1\tind2\tnSites\tLlike\t");
   //paste(paste0("p",paste(rep(GENO,10),rep(GENO,each=10),sep="_")),collapse="\t")
   fprintf(fibspair,"pAA_AA\tpAC_AA\tpAG_AA\tpAT_AA\tpCC_AA\tpCG_AA\tpCT_AA\tpGG_AA\tpGT_AA\tpTT_AA\tpAA_AC\tpAC_AC\tpAG_AC\tpAT_AC\tpCC_AC\tpCG_AC\tpCT_AC\tpGG_AC\tpGT_AC\tpTT_AC\tpAA_AG\tpAC_AG\tpAG_AG\tpAT_AG\tpCC_AG\tpCG_AG\tpCT_AG\tpGG_AG\tpGT_AG\tpTT_AG\tpAA_AT\tpAC_AT\tpAG_AT\tpAT_AT\tpCC_AT\tpCG_AT\tpCT_AT\tpGG_AT\tpGT_AT\tpTT_AT\tpAA_CC\tpAC_CC\tpAG_CC\tpAT_CC\tpCC_CC\tpCG_CC\tpCT_CC\tpGG_CC\tpGT_CC\tpTT_CC\tpAA_CG\tpAC_CG\tpAG_CG\tpAT_CG\tpCC_CG\tpCG_CG\tpCT_CG\tpGG_CG\tpGT_CG\tpTT_CG\tpAA_CT\tpAC_CT\tpAG_CT\tpAT_CT\tpCC_CT\tpCG_CT\tpCT_CT\tpGG_CT\tpGT_CT\tpTT_CT\tpAA_GG\tpAC_GG\tpAG_GG\tpAT_GG\tpCC_GG\tpCG_GG\tpCT_GG\tpGG_GG\tpGT_GG\tpTT_GG\tpAA_GT\tpAC_GT\tpAG_GT\tpAT_GT\tpCC_GT\tpCG_GT\tpCT_GT\tpGG_GT\tpGT_GT\tpTT_GT\tpAA_TT\tpAC_TT\tpAG_TT\tpAT_TT\tpCC_TT\tpCG_TT\tpCT_TT\tpGG_TT\tpGT_TT\tpTT_TT\n");



   argu2 * myPars2D= new argu2;
   myPars2D->totalSites = totalSites;
   myPars2D->keepSites =  new int[totalSites];
   
   for(int i1=0;i1<nInd-1;i1++){
     for(int i2=i1+1;i2<nInd;i2++){

       p1=i1;
       p2=i2;
       myPars2D->nInd1 = p1;
       myPars2D->nInd2 = p2;
       fprintf(stdout,"analysing pair %d %d\n",p1,p2);
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
       
       
       fprintf(fibspair,"%d\t%d\t%d\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
       for(int w=0;w<100;w++){
	 fprintf(fibspair,"\t%f",myPars2D->pi[w]);
       }
       fprintf(fibspair,"\n");
     }
   }
   //    runEM2D(genoLike,myPars2D);
   
 }
 
 

 ////////// clean  
 fclose(flog); 
 if(p2 == -1 && all != 1)
   fclose(fibs); 
 if((p1 >= 0 && p2 >= 0) || all)
   fclose(fibspair); 
 for(int i=0;i<dumpedFiles.size();i++){
   //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
   free(dumpedFiles[i]);
 }
 
 delete[] genoLike;
 return 0;
 
}
