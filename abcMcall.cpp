#include "abc.h"
#include "analysisFunction.h"
#include "shared.h"
#include "abcMcall.h"
#include <cfloat>
void abcMcall::printArg(FILE *fp){
  fprintf(fp,"\t-> doMcall=%d\n",domcall);
  
}

double theta = 1.1e-3;
void abcMcall::run(funkyPars *pars){
  int trim=0;
  if(!domcall)
    return;
  double QS_glob[5]={0,0,0,0,0};
  double QS_ind[4][10];
  for(int i=0;i<4;i++)
    for(int j=0;j<10;j++)
      QS_ind[i][j] = 0.0;
  chunkyT *chk = pars->chk;
  for(int s=0;s<chk->nSites;s++){

    for(int i=0;i<chk->nSamples;i++) {
      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;
      
      for(int j=0;j<nd->l;j++){
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
	//filter qscore, mapQ,trimming, and always skip n/N
	if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
	  continue;
	}
	QS_ind[allele][i] += qs;
      }
      for(int i=0;0&&i<4;i++)
	fprintf(stderr,"%d) %f\n",i,QS_ind[i][0]);
    }
    for(int i=0;i<chk->nSamples;i++){
      double partsum = 0;
      for(int j=0;j<4;j++)
	partsum += QS_ind[j][i];
      for(int j=0;j<4;j++)
	QS_glob[j] = QS_ind[j][i]/partsum;
    }
    for(int i=0;0&&i<5;i++)
      fprintf(stderr,"%d) %f\n",i,QS_glob[i]);
    double liks[10*pars->nInd];//<- this will be the work array


    //switch to PL just to compare numbers
    /*
      p = 10^(-q/10)
      log10(p) = -q/10
      -10*log10(p) = q;
    */
#if 1
    for(int i=0;i<pars->nInd;i++){
      float min = FLT_MAX;
      for(int j=0;j<10;j++)
	if (min > -10*log10(exp(pars->likes[s][i*10+j])))
	  min = -10*log10(exp(pars->likes[s][i*10+j]));
      for(int j=0;j<10;j++){
	pars->likes[s][i*10+j]  = (int)(-10*log10(exp(pars->likes[s][i*10+j])) - min + .499);
	//	fprintf(stderr,"PL[%d][%d]: %f\n",i,j,pars->likes[s][i*10+j]);
      }
      //now pars->likes is in phred PL scale and integerized
      for(int j=0;j<10;j++)
	pars->likes[s][i*10+j] = log(pow(10,-pars->likes[s][10*i+j]/10));
      //now pars->likes is in logscale but we have had loss of praecision
    }
#endif
    fprintf(stderr,"pars->nInd: %d\n",pars->nInd);
    for(int i=0;i<pars->nInd;i++){
      double tsum = 0.0;
      for(int j=0;j<10;j++)
	tsum += exp(pars->likes[s][i*10+j]);
      fprintf(stderr,"tsum: %f\n",tsum);
      for(int j=0;j<10;j++){
	fprintf(stderr,"PRE[%d][%d]: %f\n",i,j,exp(pars->likes[s][i*10+j]));
	liks[i*10+j] = exp(pars->likes[s][i*10+j])/tsum;
      }
      for(int j=0;j<10;j++)
	fprintf(stderr,"lk[%d][%d]: %f\n",i,j,liks[10*i+j]);
    
    }
     exit(0);
    // Watterson factor, here aM_1 = aM_2 = 1
    double aM = 1;
    for (int i=2; i<pars->nInd*2; i++) aM += 1./i;
    theta *= aM;
    if ( theta >= 1 )
      {
	fprintf(stderr,"The prior is too big (theta*aM=%.2f), going with 0.99\n", theta);
	theta = 0.99;
      }
    theta = log(theta);
    fprintf(stderr,"theta: %f\n",theta);
    exit(0);
  }
}


void abcMcall::clean(funkyPars *fp){
  if(!domcall)
    return;

  
}

void abcMcall::print(funkyPars *fp){
  if(!domcall)
    return;
      
}


void abcMcall::getOptions(argStruct *arguments){
  fprintf(stderr,"asdfadsfadsf\n");
  //default
  domcall=0;

  //from command line
  domcall=angsd::getArg("-domcall",domcall,arguments);
  if(domcall==-999){
    domcall=0;
    printArg(stderr);
    exit(0);
  }
  if(domcall==0)
    return;

  printArg(arguments->argumentFile);

}


abcMcall::abcMcall(const char *outfiles,argStruct *arguments,int inputtype){
  getOptions(arguments);
  printArg(arguments->argumentFile);
}

abcMcall::~abcMcall(){


}
