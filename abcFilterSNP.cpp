/*
  little class that does 
1) hwe using genotype likelihoods
2) a) fisher exact
   b) GATK approach
   c) SB
   These are based on guo 2012 BMC genomics 2013 13:666
3) 2 sample wilcox rank sum for baseQ bias
   This file should be completely rewritten, much to slow and stupid. But it should give correct results

 */


#include <cmath>
#include <ctype.h>
#include "analysisFunction.h"
#include "shared.h"
#include "fet.c"
#include "chisquare.h" //<- stuff related to calculating pvalues from LRT tests
#include "abc.h"
#include "abcFilterSNP.h"
#include "abcHWE.h"
#include <htslib/kstring.h>

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




double baseQbias(tNode *tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[64];
  int minD[64];
  memset(majD,0,sizeof(int)*64);
  memset(minD,0,sizeof(int)*64);

  for(int i=0;i<nInd;i++){
    tNode &nd = tn[i];
    for(int l=0;l<nd.l;l++){
      int obB = refToInt[nd.seq[l]];

      if(obB==maj)
	majD[nd.qs[l]]++;
      else if(obB==min)
	minD[nd.qs[l]]++;
    }
  }
  
  for(int i=0;0&&i<64;i++)
    fprintf(stdout,"%d %d\n",majD[i],minD[i]);



  double U[64];
  double cumLess1 =0;
  for(int i=0;i<64;i++){
    U[i]=(double) minD[i]*(majD[i]/2.0+cumLess1);
    cumLess1 += majD[i];
  }
  double U1 =0;
  for(int i=0;i<64;i++) 
    U1+= U[i];
  //below is a check
  double cumLess2 =0;
  for(int i=0;i<64;i++){
    U[i]=(double) majD[i]*(minD[i]/2.0+cumLess2);
    cumLess2 += minD[i];
  }
  double U2 =0;
  for(int i=0;i<64;i++) 
    U2+= U[i];
  

  double mu=cumLess1*cumLess2/2.0;
  double sigu=sqrt((cumLess1*cumLess2*(cumLess1+cumLess2+1))/12.0);
  double Z=(std::min(U1,U2)-mu)/sigu;

 
  //  fprintf(stderr,"U1:%f U2:%f U1+U2:%f nObs:%f Z:%f p.value:%f\n",U1,U2,U1+U2,cumLess1*cumLess2,Z,2*phi(Z));

  return Z;
}

Chisqdist chi(1);

//guo 2012 mutat res 2012, 744(2):154-160
double sb1(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double top=b/(a+b)-d/(c+d);
  double bot=(b+d)/(a+b+c+d);
  return top/bot;
}

//the gatk way
double sb2(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double en=(b/(a+b))*(c/(c+d));
  double to=(a+c)/(a+b+c+d);
  double tre=(d/(c+d))*(a/(a+b));
  return std::max(en/to,tre/to);
}

//strandbias using fisher

double sb3(int cnts[4]){

  double left,right,twotail,prob;
  kt_fisher_exact(cnts[0], cnts[1], cnts[2], cnts[3], &left, &right, &twotail);
  return twotail;
}


void abcFilterSNP::printArg(FILE *argFile){
   fprintf(argFile,"-----BETA---------------\n%s:\n",__FILE__);
  fprintf(argFile,"doSnpStat=%d\n",doSnpStat);

}
void abcFilterSNP::run(funkyPars *pars){
  if(!doSnpStat)
    return;
  chunkyT *chk = pars->chk;
  
  if(doSnpStat==1){
    kstring_t *bufstr = new kstring_t;
    bufstr->s=NULL;bufstr->l=bufstr->m=0;
    //loop over sites;
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      //loop over samples
      int cnts[4]={0,0,0,0};
      for(int i=0;i<pars->nInd;i++){
	tNode &nd = chk->nd[s][i];
	for(int l=0;l<nd.l;l++){
	  int obB = refToInt[nd.seq[l]];
	  //	    fprintf(stderr,"%c ",nd.seq[l]);
	  int strand = (isupper(nd.seq[l])==0)<<1;
	  //  fprintf(stderr,"strand:%d\n",strand);
	  if(obB==4)
	    continue;
	  if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	    continue;
	  if(obB!=pars->major[s])
	    strand +=1;
	  //fprintf(stderr,"strand=%d\n",strand);
	  cnts[strand]++;
	}
      }
      
      ksprintf(bufstr,"%s\t%d\t%d %d %d %d\t",header->target_name[pars->refId],pars->posi[s]+1, cnts[0],cnts[1],cnts[2],cnts[3]);
      ksprintf(bufstr,"%f:%f:%f\t",sb1(cnts),sb2(cnts),sb3(cnts));
      funkyHWE *hweStruct = (funkyHWE *) pars->extras[8];//THIS IS VERY NASTY! the ordering within general.cpp is now important
      double llh = 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
      double pval = chi.cdf(-llh);
      ksprintf(bufstr,"%f:%e\t",llh,1-pval);
      
      double Z = baseQbias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(bufstr,"%f:%e\n",Z,2*phi(Z));
    }
    pars->extras[index] = bufstr;
  }
  
}

void abcFilterSNP::clean(funkyPars *fp){
  if(!doSnpStat)
    return;

  
}

void abcFilterSNP::print(funkyPars *pars){
  if(!doSnpStat)
    return;
  kstring_t *bufstr =(kstring_t*) pars->extras[index];
  gzwrite(outfileZ,bufstr->s,bufstr->l);
  free(bufstr->s);
  delete bufstr;
}


void abcFilterSNP::getOptions(argStruct *arguments){
  //default


  //from command line
  doSnpStat=angsd::getArg("-doSnpStat",doSnpStat,arguments);

  if(doSnpStat==0)
    return;
  int domajmin=0;
  //from command line
  int doHWE =0;
  domajmin=angsd::getArg("-doSNP",domajmin,arguments);
  doHWE=angsd::getArg("-doHWE",doHWE,arguments);
  if(!domajmin||!doHWE){
    fprintf(stderr,"-doSnpStat require -doSNP and -doHWE \n");
    exit(0);
  }
    
  

}


abcFilterSNP::abcFilterSNP(const char *outfiles,argStruct *arguments,int inputtype){
  doSnpStat=0;
  outfileZ = Z_NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSnpStat")||!strcasecmp(arguments->argv[1],"-doPost")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doSnpStat==0){
    shouldRun[index] =0;
    return;
  }
  

  if(doSnpStat){
    fprintf(stderr,"running doSnpStat=%d\n",doSnpStat);
    const char *postfix=".snpStat.gz";
    outfileZ = aio::openFileGz(outfiles,postfix,GZOPT);
    gzprintf(outfileZ,"Chromo\tPosition\t+Major +Minor -Major -Minor\tSB1:SB2:SB3\tHWE_LRT:HWE_pval\tbaseQ_Z:baseQ_pval\n");
  }
}

abcFilterSNP::~abcFilterSNP(){
  if(outfileZ!=Z_NULL)     gzclose(outfileZ);

}
