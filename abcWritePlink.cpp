/*
  This is a class that writes the genotype calls from -doGeno as plink
  tped/tfam files.
 */

#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "shared.h"
#include "abcCallGenotypes.h"
#include "abcWritePlink.h"
#include "aio.h"

void abcWritePlink::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doPlink\t%d\n",doPlink);
  fprintf(argFile,"\t2: tfam/tped format\n");
  fprintf(argFile,"\n\tNB This is a wrapper around -doGeno see more information for that option\n");
  fprintf(argFile,"\tNB A negative -doGeno suppresses the .geno.gz output\n");
}

void abcWritePlink::run(funkyPars *pars){
}

void abcWritePlink::clean(funkyPars *pars){
}

void abcWritePlink::print(funkyPars *pars){
  if(doPlink==0)
    return;
  genoCalls *geno =(genoCalls *) pars->extras[11];
  const char *chr = header->target_name[pars->refId];

  kstr.l=0;
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    ksprintf(&kstr,"%s %s_%d 0 %d",chr,chr,pars->posi[s]+1,pars->posi[s]+1);
    char maj = intToRef[pars->major[s]];
    char mnr = intToRef[pars->minor[s]];
    for(int i=0;i<pars->nInd;i++){
      switch(geno->dat[s][i]){
      case 0:
	ksprintf(&kstr,"\t%c %c",maj,maj);
	break;
      case 1:
	ksprintf(&kstr,"\t%c %c",maj,mnr);
	break;
      case 2:
	ksprintf(&kstr,"\t%c %c",mnr,mnr);
	break;
      default:
	kputs("\t0 0",&kstr);
      }
    }
    kputc('\n',&kstr);
  }
  if(kstr.l>0&&fwrite(kstr.s,1,kstr.l,fp)!=kstr.l){
    fprintf(stderr,"\t-> Problem writing the .tped file\n");
    exit(1);
  }
}

void abcWritePlink::getOptions(argStruct *arguments){
  doPlink=angsd::getArg("-doPlink",doPlink,arguments);
  if(doPlink==0)
    return;

  if(doPlink!=2){
    fprintf(stderr,"\t-> -doPlink %d is not supported, only -doPlink 2 (tped/tfam)\n",doPlink);
    exit(1);
  }

  int doGeno = 0;
  doGeno = angsd::getArg("-doGeno",doGeno,arguments);
  if(doGeno==0){
    fprintf(stderr,"\t-> Must supply -doGeno to write plink files (a negative value suppresses the .geno.gz output)\n");
    exit(1);
  }
}

abcWritePlink::abcWritePlink(const char *outfiles,argStruct *arguments,int inputtype){
  fp=NULL;
  doPlink =0;
  kstr.s=NULL;kstr.l=kstr.m=0;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doPlink")){
      printArg(stderr);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  if(doPlink==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);

  //family id and individual id are both the index of the individual
  FILE *fam = aio::openFile(outfiles,".tfam");
  for(int i=1;i<=arguments->nInd;i++)
    fprintf(fam,"%d %d 0 0 0 -9\n",i,i);
  fclose(fam);

  fp = aio::openFile(outfiles,".tped");
}

abcWritePlink::~abcWritePlink(){
  if(fp)
    fclose(fp);
  free(kstr.s);
}
