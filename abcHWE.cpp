/*
  thorfinn thorfinn@binf.ku.dk 19dec 2012

  Authors of this file
  filipe , anders , thorfinn

  part of angsd

 */
#include <cmath>
#include "abcFreq.h"
#include "abcHWE.h"
 

void abcHWE::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-HWE_pval\t%f\n",HWE_pval);
  fprintf(argFile,"\n");
}

void abcHWE::getOptions(argStruct *arguments){
  int tmpDoMaf =0;
  int GL =0;
  HWE_pval = angsd::getArg("-HWE_pval",HWE_pval,arguments);

  if(HWE_pval==0)
    return;
 
  doHWE=1;
  chisq = new Chisqdist(1);//<- 1degree of freedom;
  LRT_thres = chisq->invcdf(1-HWE_pval);
  fprintf(stderr,"\t-> HWE-filter using a pvalue: %.4e correspond to %.4f likelihood units\n",HWE_pval,LRT_thres);    

  
  if(doHWE==0)
    return;
  GL=angsd::getArg("-GL",GL,arguments);
  tmpDoMaf=angsd::getArg("-doMaf",tmpDoMaf,arguments);
  if(doHWE&&(tmpDoMaf==0)){
    fprintf(stderr,"You supplied -HWE_pval, you should also choose -doMaf\n");
    exit(0);
  }
  if(arguments->inputtype==INPUT_BEAGLE){
    fprintf(stderr,"Error: you cannot estimate HWE based on posterior probabilities (beagle)\n");
    exit(0);
  }
  

}

abcHWE::abcHWE(const char *outfiles,argStruct *arguments,int inputtype){
  chisq = NULL;
  outfileZ = NULL;
  doHWE=0;
  HWE_pval = 0;
  testMe=0;
  tolStop = 0.00001;


  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-HWE_pval")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doHWE==0){
    shouldRun[index] = 0;
    return;
  }
  //make output files
  const char* postfix;
  postfix=".hwe.gz";
  if(doHWE>0){
    outfileZ = aio::openFileGz(outfiles,postfix,GZOPT);
    //print header
    gzprintf(outfileZ,"Chromo\tPosition\tMajor\tMinor\tFreq\thweFreq\tF\tLRT\tp-value\n");
  }
}


abcHWE::~abcHWE(){

  if(doHWE==0)
    return;
  if(doHWE>0)
    if(outfileZ) gzclose(outfileZ);
  delete chisq;
}


void abcHWE::clean(funkyPars *pars){
  if(doHWE==0)
    return;

  funkyHWE *hweStruct =(funkyHWE *) pars->extras[index];
  delete[] hweStruct->freq;
  delete[] hweStruct->F;
  delete[] hweStruct->like0;
  delete[] hweStruct->likeF;
  delete hweStruct;
  
}

void abcHWE::print(funkyPars *pars){
  if(doHWE<=0)
    return;

  funkyHWE *hweStruct = (funkyHWE *) pars->extras[index];//new

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    freqStruct *freq = (freqStruct *) pars->extras[6];
    float lrt= 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
    float pval;
    if(lrt<0)
      pval=1;
    else
      pval=1-chisq->cdf(lrt);
    //    fprintf(stderr,"lrt:%f\n",lrt);
    gzprintf(outfileZ,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%e\t%e\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],hweStruct->freq[s],hweStruct->F[s],lrt,pval);

  }

}


void abcHWE::run(funkyPars *pars){
 
  if(doHWE==0)
    return;

  //  pars->hweStruct = new funkyHWE;//old
  funkyHWE *hweStruct = new funkyHWE;//new

  double *freq = new double[pars->numSites];
  double *F = new double[pars->numSites];
  double *like0 = new double[pars->numSites];
  double *likeF = new double[pars->numSites];
  freqStruct *freqS = (freqStruct *) pars->extras[6];

  double **loglike3;
  loglike3=angsd::get3likesRescale(pars);

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    
    //start parameter
    double x[2];
    x[0]=freqS->freq[s];
    x[1]=0.05;
    estHWE(x,loglike3[s],pars->nInd);
    freq[s]=x[0];
    F[s]=x[1];
    likeF[s] = HWE_like(x,loglike3[s],pars->nInd);
    x[1]=0.0;
    x[0]=freqS->freq[s];
    like0[s] = HWE_like(x,loglike3[s],pars->nInd);
    //    fprintf(stderr,"%f\t%f\n",x[0],x[1]);
    //fprintf(stderr,"%f\t%f\t%f\n",loglike3[s][0],loglike3[s][1],loglike3[s][2]);
    float lrt= 2*like0[s]-2*likeF[s];
    if(lrt<0)
      lrt=0;
    if(lrt<LRT_thres)
      pars->keepSites[s] =0;
  }


  hweStruct->freq=freq;
  hweStruct->F=F;
  hweStruct->like0=like0;
  hweStruct->likeF=likeF;
  pars->extras[index] = hweStruct;


  for(int s=0;s<pars->numSites;s++)
    delete[] loglike3[s];
  delete[] loglike3;

}

double abcHWE::differ(double *x,double *y){
  return (abs(x[1]-y[1]) + abs(x[0]-y[0]));
}


void abcHWE::HWE_EM(double *x,double *loglike,int nInd){
  
  double freq=x[0];
  double F=x[1];
  double p0=(pow(1-freq,2)+freq*(1-freq)*F);
  double p1=(2*freq*(1-freq)-2*freq*(1-freq)*F);
  double p2=(pow(freq,2)+freq*(1-freq)*F);

  double freq2=0;
  double F2=0;
  double norm;
  double oHO=0;
  double eHO=nInd*(pow(freq,2)+pow(1-freq,2));

  for(int i=0;i<nInd;i++){
    norm=angsd::addProtect3(log(p0)+loglike[i*3+0],log(p1)+loglike[i*3+1],log(p2)+loglike[i*3+2]);
    freq2+=exp(log(p1)+loglike[i*3+1]-norm)+exp(log(2)+log(p2)+loglike[i*3+2]-norm);
    oHO+=exp(log(p0)+loglike[i*3+0]-norm)+exp(log(p2)+loglike[i*3+2]-norm);
  //oHO+=(p0*loglike[i*3+0]+p2*loglike[i*3+2])/(2*norm);
    //fprintf(stderr,"1 %f 2 %f 3 %f\n",loglike[i*3+0],loglike[i*3+1],loglike[i*3+2]);
  }
  F2=(oHO-eHO)/(nInd-eHO);
  if(F2<0)
    F2=0;
  if(F2>1)
    F2=1;
  
  freq2=freq2/(2*nInd);
  //fprintf(stderr,"norm %f\tx[0] %f\tx[1] %f\t1 %f 2 %f 3 %f oHO %f eHO %f freq %f\n",norm,x[0],x[1],p0,p1,p2,oHO,eHO,freq2);

  if(freq2<0.0001)
    F2=0;
  //    F2=0;
  x[0]=freq2;
  x[1]=F2;
  if(freq2>1.0000001){
    fprintf(stderr,"something is wrong i HWE\t freq %f\n",freq2);
    fflush(stderr);
    exit(0);

  }
}




double abcHWE::HWE_like(double *x,double *loglike,int nInd){
  double freq=x[0];
  double F=x[1];
  double p0=(pow(1-freq,2)+freq*(1-freq)*F);
  double p1=(2*freq*(1-freq)-2*freq*(1-freq)*F);
  double p2=(pow(freq,2)+freq*(1-freq)*F);
  double totalLogLike=0;
  for(int i=0;i<nInd;i++)
    totalLogLike+=angsd::addProtect3(log(p0)+loglike[i*3+0],log(p1) + loglike[i*3+1],log(p2) + loglike[i*3+2]);
  return -totalLogLike;
}

void abcHWE::estHWE(double *x,double *loglike,int nInd){

  double l,d;
  if(testMe){
    l=HWE_like(x,loglike,nInd);
    d=l;
  }
  double *y = new double[2];
  y[0] = x[0]+2;
  y[1] = x[1]+2;
  int iter=50;
  //fprintf(stderr,"%f\n",differ(x,y));
  int printer=0;
  for(int i=0;i<iter;i++){
    HWE_EM(x,loglike,nInd);
    if(testMe){
      l=HWE_like(x,loglike,nInd);
      if(d>l+0.01){
	//   fprintf(stderr,"d %f\tl %f\n",d,l);
	printer=1;
      }
      d=l;
    }
    double dif = differ(x,y);
    // fprintf(stderr,"%f\n",dif);
    if(dif<tolStop)
      break;
    y[0] = x[0];
    y[1] = x[1];

    
  }
  if(printer && testMe ){
    x[0]=0.05;
    x[1]=0.05;
    l=HWE_like(x,loglike,nInd);
    for(int i=0;i<iter;i++){
      fprintf(stderr,"like %d\t%f\tf %f\tF %f\t%d\n",i,l,x[0],x[1],nInd);
      HWE_EM(x,loglike,nInd);
      l=HWE_like(x,loglike,nInd);
    }

    for(int ind=0;ind<nInd;ind++){
      fprintf(stdout,"likes %f, %f, %f\n",loglike[ind*3+0],loglike[ind*3+1],loglike[ind*3+2]);
    }
    fflush(stderr);
    exit(0);

  }
}
