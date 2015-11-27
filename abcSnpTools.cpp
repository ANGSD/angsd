#include "abc.h"
#include "abcSnpTools.h"
#include "analysisFunction.h"
#include "shared.h"


void abcSnpTools::printArg(FILE *fp){
  fprintf(fp,"doSnpTools=%d\n",doSnpTools);
  
}
void save_chr(uint16_t *p){

}
void abcSnpTools::run(funkyPars *pars){
  if(!doSnpTools)
    return;
  if(curChr!=pars->refId){
    save_chr(ebd);
  }
  double **ebds = new double*[pars->numSites];
  for(int s=0;s<pars->numSites;s++){
    double *ebd=new double[6*pars->nInd];
    for(int i=0;i<pars->nInd;i++){
      tNode *nd = pars->chk->nd[s][i];
      for (int l=0;nd&&l<nd->l;l++){


      }


    }
    ebds[s]=ebd;
  }
  
    
}


void abcSnpTools::clean(funkyPars *fp){
  if(!doSnpTools)
    return;

  
}

void abcSnpTools::print(funkyPars *fp){
  if(!doSnpTools)
    return;
      
}


void abcSnpTools::getOptions(argStruct *arguments){
  //default
  doSnpTools=0;

  //from command line
  doSnpTools=angsd::getArg("-doSnpTools",doSnpTools,arguments);
  if(doSnpTools==-999){
    doSnpTools=0;
    printArg(stderr);
    exit(0);
  }
  if(doSnpTools==0)
    return;

  printArg(arguments->argumentFile);

}


abcSnpTools::abcSnpTools(const char *outfiles,argStruct *arguments,int inputtype){
  curChr =-1;
  getOptions(arguments);
  if(doSnpTools){
    fprintf(stderr,"running doSnpTools=%d\n",doSnpTools);
    phred=new double[256];
    for(uint i=0; i<256; i++)	
      phred[i]=1-exp(-0.1*i);
  }else{
    shouldRun[index] = 0;
  }

}

abcSnpTools::~abcSnpTools(){


}
