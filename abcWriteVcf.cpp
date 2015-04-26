/*
  This is a class that dumps plink file output

 */
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abcCallGenotypes.h"
#include "abcWriteVcf.h"
void abcWriteVcf::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doVcf\t%d\n",doVcf);
  fprintf(argFile,"\t1: binary fam/bim/bed format (still beta, not really working)\n");
  fprintf(argFile,"\n\tNB This is a wrapper around -doGeno see more information for that option\n");
}

void abcWriteVcf::run(funkyPars *pars){
  if(doVcf==0)
    return ;
   
}

void abcWriteVcf::clean(funkyPars *pars){
  if(doVcf==0)
    return;
    

}

void abcWriteVcf::print(funkyPars *pars){
  if(doVcf==0)
    return;

}

void abcWriteVcf::getOptions(argStruct *arguments){

  doVcf=angsd::getArg("-doPlink",doVcf,arguments);
  
}

abcWriteVcf::abcWriteVcf(const char *outfiles,argStruct *arguments,int inputtype){
  fp=Z_NULL;
  doVcf =0;
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doVcf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doVcf==0){
    shouldRun[index] =0;
    return;
}
  
}


abcWriteVcf::~abcWriteVcf(){
    gzclose(fp);
}


