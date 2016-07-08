/*

  test for HWE from genotype likelihhoods


  Anders albrecht@binf.ku.dk 10 April 2016

  Authors of this file:
  Anders 

  part of ANGSD: http://www.popgen.dk/angsd 

  

 */
#include <htslib/kstring.h>
#include "abcFreq.h"
#include "abcHWE.h"
 

void abcHWE::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-HWE_pval\t%f\n",HWE_pval);
  fprintf(argFile,"\n");
}

void abcHWE::getOptions(argStruct *arguments){

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

  if(arguments->inputtype==INPUT_BEAGLE||arguments->inputtype==INPUT_VCF_GP){
    fprintf(stderr,"Error: you cannot estimate HWE based on posterior probabilities \n");
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
  bufstr.s=NULL;bufstr.l=bufstr.m=0;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-HWE_pval")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  if(doHWE==0){
    shouldRun[index] = 0;
    return;
  }
  printArg(arguments->argumentFile);
  //make output files
  const char* postfix;
  postfix=".hwe.gz";
  if(doHWE>0){
    outfileZ = aio::openFileBG(outfiles,postfix);
    //print header
    const char *str = "Chromo\tPosition\tMajor\tMinor\thweFreq\tFreq\tF\tLRT\tp-value\n";
    aio::bgzf_write(outfileZ,str,strlen(str));
  }
}


abcHWE::~abcHWE(){

  if(doHWE==0)
    return;
  if(doHWE>0)
    if(outfileZ!=NULL)
      bgzf_close(outfileZ);
  delete chisq;
}


void abcHWE::clean(funkyPars *pars){
  if(doHWE==0)
    return;

  funkyHWE *hweStruct =(funkyHWE *) pars->extras[index];
  delete[] hweStruct->freq;
  delete[] hweStruct->freqHWE;
  delete[] hweStruct->F;
  delete[] hweStruct->like0;
  delete[] hweStruct->likeF;
  delete hweStruct;
  
}

void abcHWE::print(funkyPars *pars){
  if(doHWE<=0)
    return;

  funkyHWE *hweStruct = (funkyHWE *) pars->extras[index];//new
  bufstr.l=0;
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
   
    float lrt= 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
    //DRAGON lrt is sometimes nan
    float pval;
    if(std::isnan(lrt))
      pval=lrt;
    else if(lrt<0)
      pval=1;
    else			
      pval=1-chisq->cdf(lrt);
    //    fprintf(stderr,"lrt:%f\n",lrt);
    ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%e\t%e\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],hweStruct->freqHWE[s],hweStruct->freq[s],hweStruct->F[s],lrt,pval);

  }
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
}


void abcHWE::run(funkyPars *pars){
 
  if(doHWE==0)
    return;

  funkyHWE *hweStruct = new funkyHWE;

  double *freq = new double[pars->numSites];
  double *freqHWE = new double[pars->numSites];
  double *F = new double[pars->numSites];
  double *like0 = new double[pars->numSites];
  double *likeF = new double[pars->numSites];
 
  double **loglike3;
  loglike3=angsd::get3likesRescale(pars);


  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    //est under HWE
    freqHWE[s] = angsd::estFreq(loglike3[s],pars->nInd);



    //start parameter
    double x[2];
    x[1]=0.0; //F
    x[0]=freqHWE[s]; // freq

    //log like for HWE freq
    like0[s] = HWE_like(x,loglike3[s],pars->nInd);


    //start parameter for EM_F
    x[0]=freqHWE[s];
    x[1]=0.05;

    estHWE(x,loglike3[s],pars->nInd);
    freq[s]=x[0];
    F[s]=x[1];
    likeF[s] = HWE_like(x,loglike3[s],pars->nInd);


    //    fprintf(stderr,"%f\t%f\n",x[0],x[1]);
    //fprintf(stderr,"%f\t%f\t%f\n",loglike3[s][0],loglike3[s][1],loglike3[s][2]);
    float lrt= 2*like0[s]-2*likeF[s];
    if(lrt<0)
      lrt=0;
    if(lrt<LRT_thres)
      pars->keepSites[s] =0;
  }


  hweStruct->freq=freq;
  hweStruct->freqHWE=freqHWE;
  hweStruct->F=F;
  hweStruct->like0=like0;
  hweStruct->likeF=likeF;
  pars->extras[index] = hweStruct;


  for(int s=0;s<pars->numSites;s++)
    delete[] loglike3[s];
  delete[] loglike3;

}

double abcHWE::differ(double *x,double *y){
  return (fabs(x[1]-y[1]) + fabs(x[0]-y[0]));
}


void abcHWE::HWE_EM(double *x,double *loglike,int nInd){
  // freq2HWE <- function(x) {  freq <- x[2]*0.5 + x[1];   F <-	1 - x[2]/(2*freq*(1-freq)); c(freq,F);}


  

  double freq0,freq1,freq2;
  if(1){ // use F formulation instead of p0,p1,p2: compatable with the likelihood function
    double freq = x[0];
    double Fadj = freq*(1-freq)*x[1];
    freq0= pow(1-freq,2) + Fadj;
    freq1= 2*freq*(1-freq) - 2*Fadj;
    freq2= pow(freq,2) + Fadj;
  }
  else{
    freq2 = x[0];
    freq1 = x[1];
    freq0 = 1- x[0] - x[1];
  }
 
  double newFreq0=0;
  double newFreq1=0;
  double newFreq2=0;
  double norm;

  for(int i=0;i<nInd;i++){
    norm=angsd::addProtect3(log(freq0)+loglike[i*3+0],log(freq1)+loglike[i*3+1],log(freq2)+loglike[i*3+2]);
    newFreq0 += exp(log(freq0)+loglike[i*3+0] - norm);
    newFreq1 += exp(log(freq1)+loglike[i*3+1] - norm);
    newFreq2 += exp(log(freq2)+loglike[i*3+2] - norm);
  }

  newFreq0 = newFreq0/nInd;
  newFreq1 = newFreq1/nInd;
  newFreq2 = newFreq2/nInd;

  if(newFreq0 < 0.0001)
    newFreq0 = 0.0001;

  if(newFreq1 < 0.0001)
    newFreq1 = 0.0001;

  if(newFreq2 < 0.0001)
    newFreq2 = 0.0001;

  norm = newFreq0 + newFreq1 + newFreq2;
  newFreq0 = newFreq0/norm;
  newFreq1 = newFreq1/norm;
  newFreq2 = newFreq2/norm;

  
  if(newFreq0 + newFreq1 >1.00001){
    fprintf(stderr,"something is wrong i HWE\t freq %f\t F %f\n",x[0],x[1]);
    fflush(stderr);
    exit(0);
  }
 
 
  if(1){
    if(isnan(newFreq1*0.5 + newFreq2)){
      fprintf(stderr,"problems in abcHWE: \tx[0] %f\tx[1] %f\told 1 %f 2 %f 3 %f\tnew 1 %f 2 %f 3 %f \n",x[0],x[1],freq0,freq1,freq2,newFreq0,newFreq1,newFreq2);
      //      exit(0);
    }
    x[0] = newFreq1*0.5 + newFreq2;
    x[1] = 1 - newFreq1/(2*x[0]*(1-x[0])) ;
  }
  else{
    x[0] = newFreq2;
    x[1] = newFreq1;
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

  double y[2];// = new double[2];
  y[0] = x[0]+2;
  y[1] = x[1]+2;
  int iter=50;
  double dif;
  //fprintf(stderr,"%f\n",differ(x,y));
  int printer=0;
  for(int i=0;i<iter;i++){
    HWE_EM(x,loglike,nInd);
    dif = differ(x,y);
    // fprintf(stderr,"%f\n",dif);
    if(dif<tolStop)
      break;
    y[0] = x[0];
    y[1] = x[1];

  }

}
