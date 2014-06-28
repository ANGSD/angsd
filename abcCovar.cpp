#include <cmath>

#include "shared.h"
#include "analysisFunction.h"
#include "abcCovar.h"
#include "abcFreq.h"


typedef struct{
  double *res; 
}funkyCovar;




void abcCovar::getOptions(argStruct *arguments){
  doCovar=angsd::getArg("-doCovar",doCovar,arguments);
  
  if(doCovar){
    int doCounts=0;
    int doMaf =0;
    doCounts=angsd::getArg("-doCounts",doCounts,arguments);
    doMaf=angsd::getArg("-doMaf",doMaf,arguments);
    if(doCounts==0||doMaf==0){
      fprintf(stderr,"\t-doCovar requires -doCounts 1 and -doMaf, will exit\n");
      exit(0);
    }

  }

}


void abcCovar::printArg(FILE *fp){


}

abcCovar::abcCovar(const char *outfiles,argStruct *arguments){
  
  //  fprintf(stderr,"index=%d tot_index=%d\n",index,tot_index);

  //should be changes so that it is a command options

  doCovar=0;
  outfile=NULL;
  getOptions(arguments);
  if(doCovar==0){
    shouldRun[index] = 0;
    return;
  }
  const char* postfix;
  postfix=".covar";
  outfile = aio::openFile(outfiles,postfix);
}


void abcCovar::openfile(const char *outfiles){
  //should be changes so that it is a command options
  if(doCovar==0)
    return;
  const char* suffix;
  suffix=".covar";
  outfile = aio::openFile(outfiles,suffix);

}


abcCovar::~abcCovar(){
  //  fprintf(stderr,"calling covar destructor\n");
  if(doCovar==0)
    return;

  if(outfile) fclose(outfile);
}

void abcCovar::addDefault(funkyPars *pars){
  //should be changes so that it is a command options
  if(doCovar==0)
    return;
}

void abcCovar::print(funkyPars *pars){
  funkyCovar *covStruct =(funkyCovar *) pars->extras[index];
  if(doCovar==0)
    return;
   for(int ind1=0;ind1<pars->nInd;ind1++){
    for(int ind2=0;ind2<pars->nInd;ind2++)
      fprintf(outfile,"%f\t",covStruct->res[ind1*pars->nInd+ind2]);
      //      fprintf(outfile,"%f\t",pars->covStruct->res[ind1*pars->nInd+ind2]);
    fprintf(outfile,"\n");
  }

}


void abcCovar::clean(funkyPars *pars){
  if(doCovar==0)
    return;

  //  delete [] pars->covStruct->res;//old

  //new
  funkyCovar *covStruct =(funkyCovar *) pars->extras[index];
  delete [] covStruct->res;//old
  delete covStruct;
}


void abcCovar::run(funkyPars *pars){
  if(doCovar==0)
    return;
  // fprintf(stderr,"start cover\n");
 

  double minFreq=0.20;
  double minPost=0.95;
  int minGeno=5;
  int *keepSNP;
  keepSNP=new int[pars->numSites];

 
  int *keepList;
  double *covRes;
  covRes = new double[pars->nInd*pars->nInd];
  keepList = new int[pars->nInd*pars->numSites];

 
  /////////////here only the posterior estimate for the genotype is used
  double **Egeno;
  Egeno = new double*[pars->numSites]; 

  double *avg;
  avg = new double[pars->numSites]; 

  
  for(int s=0;s<pars->numSites;s++)
    Egeno[s] = new double[pars->nInd];
  
  freqStruct *freq = (freqStruct *) pars->extras[6];
  //maybe change to posterior based on frequency
  for(int s=0;s<pars->numSites;s++){
    keepSNP[s]=0;
    avg[s]=0;
    for(int i=0;i<pars->nInd;i++){
      keepList[i*pars->nInd+s]=pars->keepSites[s];
      
      double norm = angsd::addProtect3(
				       pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]]+2*log(1-freq->freq[s]),
				       pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]]+log(2)+log(freq->freq[s])+log(1-freq->freq[s]),
				       pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]]+2*log(freq->freq[s])
				       );
      double theMax = angsd::getMax(pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]]+2*log(1-freq->freq[s]),
				    pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]]+log(2)+log(freq->freq[s])+log(1-freq->freq[s]),
				    pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]]+2*log(freq->freq[s])
				    );
      if(exp(theMax-norm)<minPost)
	keepList[i*pars->nInd+s]=0;
      else if(norm<-20)
	keepList[i*pars->nInd+s]=0;
      else if(pars->counts[s][i]+pars->counts[s][i+1]+pars->counts[s][i+2]+pars->counts[s][i+3] >8)
	keepList[i*pars->nInd+s]=0;
      else if(pars->counts[s][i]+pars->counts[s][i+1]+pars->counts[s][i+2]+pars->counts[s][i+3]<4)
	keepList[i*pars->nInd+s]=0;
      else if(pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->major[s]]]+
	 pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]]+
	 pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]]
	 >-0.0001)//anders
        keepList[i*pars->nInd+s]=0;
      if(keepList[i*pars->nInd+s])
	keepSNP[s]++;
      Egeno[s][i]=exp(pars->likes[s][i*10+angsd::majorminor[pars->major[s]][pars->minor[s]]]+log(2)+log(freq->freq[s])+log(1-freq->freq[s])-norm)+
	2*exp(pars->likes[s][i*10+angsd::majorminor[pars->minor[s]][pars->minor[s]]]+2*log(freq->freq[s])- norm);
      avg[s]+=Egeno[s][i];
    }
    avg[s]=avg[s]/pars->nInd;
  }
  ////////////


  ////////////covariance matrix
  //might try to change the order of the loops if this is too slow. Hopefully the compiler is smart enoght
  int nPolySites;

  for(int ind1=0;ind1<pars->nInd;ind1++){
    for(int ind2=0;ind2<pars->nInd;ind2++){
      covRes[ind1*pars->nInd+ind2]=0;
      nPolySites=0;
      for(int s=0;s<pars->numSites;s++){
	if(keepSNP[s]<minGeno)
	  continue;
	if(freq->freq[s]<minFreq||freq->freq[s]>1-minFreq||keepList[ind1*pars->nInd+s]==0||keepList[ind2*pars->nInd+s]==0)
	  continue;
	covRes[ind1*pars->nInd+ind2]+=
	  (Egeno[s][ind1]-avg[s])/sqrt(freq->freq[s]*(1-freq->freq[s]))*
	  (Egeno[s][ind2]-avg[s])/sqrt(freq->freq[s]*(1-freq->freq[s])); //remove sqrt at some point
	nPolySites++;
      }
      covRes[ind1*pars->nInd+ind2]=covRes[ind1*pars->nInd+ind2]/nPolySites;
    }
  }
  /*
  for(int ind1=0;ind1<pars->nInd;ind1++){
    for(int ind2=0;ind2<pars->nInd;ind2++)
      fprintf(stdout,"%f\t",covRes[ind1*pars->nInd+ind2]);
    fprintf(stdout,"\n");
  }
  */
  int sumKeep=0;
  for(int s=0;s<pars->numSites;s++)
    for(int i=0;i<pars->nInd;i++)
      sumKeep+=keepList[i*pars->nInd+s];

  int sumKeepSNP=0;
  for(int s=0;s<pars->numSites;s++)
    if(keepSNP[s]>=minGeno)
      sumKeepSNP++;
  /////////////clean up
    for(int i=0;i<pars->numSites;i++)
      delete [] Egeno[i];
    delete [] Egeno;
    delete [] avg;
    ////////////
    // fprintf(stderr,"end cover. %d genotypes %d SNPs\n",sumKeep,sumKeepSNP);


  //plot(eigen(read.table("covRes"))$vector[,1:2])
  //plot(eigen(r<-read.table("covRes"))$vector[,1:2],col=rep(1:2,c(32,73-32)))
    //pars->results->freq->pEMun=covRes;    
    
    ///pars->covStruct->res=covRes;//OLD
    
    funkyCovar *covStruct = new funkyCovar;
    covStruct->res=covRes;
    pars->extras[index] = covStruct;

}
