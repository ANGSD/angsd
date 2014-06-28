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
  if(doSnpTools)
    fprintf(stderr,"running doSnpTools=%d\n",doSnpTools);
  else{
    shouldRun[index] = 0;
  }
}

abcSnpTools::~abcSnpTools(){


}
