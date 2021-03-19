#include "abc.h"
#include "analysisFunction.h"
#include "shared.h"
#include "abcMcall.h"

void abcMcall::printArg(FILE *fp){
  fprintf(fp,"\t-> doMcall=%d\n",domcall);
  
}

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

    for(int i=0;i<chk->nSamples;i++){
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
      for(int i=0;i<4;i++)
	fprintf(stderr,"%d) %f\n",i,QS_ind[i][0]);
    }
    for(int i=0;i<chk->nSamples;i++){
      double partsum = 0;
      for(int j=0;j<4;j++)
	partsum += QS_ind[j][i];
      for(int j=0;j<4;j++)
	QS_glob[j] = QS_ind[j][i]/partsum;
    }
    for(int i=0;i<5;i++)
      fprintf(stderr,"%d) %f\n",i,QS_glob[i]);
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
