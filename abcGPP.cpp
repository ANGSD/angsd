#include "abc.h"
#include "abcGPP.h"
#include "analysisFunction.h"
#include "shared.h"


void abcGPP::printArg(FILE *fp){
  fprintf(fp,"\t->doGPP=%d\n",doGPP);
  
}
void abcGPP::run(funkyPars *pars){
  if(!doGPP)
    return;
  
    
}


void abcGPP::clean(funkyPars *fp){
  if(!doGPP)
    return;

  
}

void abcGPP::print(funkyPars *fp){
  if(!doGPP)
    return;
      
}


void abcGPP::getOptions(argStruct *arguments){
  //default
  doGPP=0;

  //from command line
  doGPP=angsd::getArg("-doGPP",doGPP,arguments);
  if(doGPP==-999){
    doGPP=0;
    printArg(stderr);
    exit(0);
  }
  if(doGPP==0)
    return;

  printArg(arguments->argumentFile);

}


abcGPP::abcGPP(const char *outfiles,argStruct *arguments,int inputtype){
  
  getOptions(arguments);
  if(doGPP){

  }else{
    shouldRun[index] = 0;
  }

}

abcGPP::~abcGPP(){


}
