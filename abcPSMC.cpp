/*
  This is a class that dumps psmc file output

 */
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abcPSMC.h"
void abcPSMC::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-writePSMC\t%d\n",dopsmc);
  fprintf(argFile,"\t1:  (still beta, not really working)\n");
}

void abcPSMC::run(funkyPars *pars){
  if(dopsmc==0)
    return ;
   
}

void abcPSMC::clean(funkyPars *pars){
  if(dopsmc==0)
    return;
    

}

void abcPSMC::print(funkyPars *pars){
  if(dopsmc==0)
    return;

}

void abcPSMC::getOptions(argStruct *arguments){

  dopsmc=angsd::getArg("-doPSMC",dopsmc,arguments);
  if(dopsmc==0)
    return;
  int gl =0;
  gl=angsd::getArg("-gl",gl,arguments);

  if(gl==0||INPUT_BEAGLE||INPUT_VCF_GP){
    fprintf(stderr,"\nPotential problem. We need genotypes for dumping psmc\n\n");
    exit(0);
  }
  
  
}

abcPSMC::abcPSMC(const char *outfiles,argStruct *arguments,int inputtype){

  dopsmc =0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-dopsmc")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);

  if(dopsmc==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);

}


abcPSMC::~abcPSMC(){
 
}
