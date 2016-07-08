#include <ctype.h>
#include "shared.h"
#include "analysisFunction.h"
#include "abcAncestry.h"
#include <cassert>

void abcAncestry::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doAncestry\t%d \n",doAncestry);
}

void abcAncestry::getOptions(argStruct *arguments){
  doAncestry=angsd::getArg("-doAncestry",doAncestry,arguments);
  if(doAncestry==0){
    shouldRun[index]=0;
    return;
  }
}

abcAncestry::abcAncestry(const char *outfiles,argStruct *arguments,int inputtype){
  doAncestry = 0;
  outfile = NULL;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAncestry")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  

  if(doAncestry==0)
    return ;
  printArg(arguments->argumentFile);
  
  outfile = aio::openFile(outfiles,".ancestry");
}

abcAncestry::~abcAncestry(){
  if(outfile!=NULL) 
    fclose(outfile);
}

void abcAncestry::clean(funkyPars *pars){

}



void abcAncestry::print(funkyPars *pars){
  assert(pars->anc!=NULL);
  if(doAncestry==1){
  
  }
  
}


void abcAncestry::run(funkyPars *pars){


}


