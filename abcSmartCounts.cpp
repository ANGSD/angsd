/*
  thorfinn thorfinn@binf.ku.dk 1april 2012

 */
#include <cmath>
#include <climits>
#include "abcSmartCounts.h"
#include "aio.h"


void abcSmartCounts::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSmartCounts\t%d\t(deprecated, write A,C,G,T counts for a single sample to .counts.bin/.counts.idx, needs -doCounts 1)\n",doSmartCounts);
  fprintf(argFile,"\n");
}

void abcSmartCounts::getOptions(argStruct *arguments){

  doSmartCounts = angsd::getArg("-doSmartCounts",doSmartCounts,arguments);

  int tmp=0;
  tmp=angsd::getArg("-doCounts",tmp,arguments);
  if(tmp==0&&doSmartCounts>0){
    fprintf(stderr,"\t-> SmartCounts require -doCounts\n");
    exit(1);
  }
  //only the counts of the first sample are stored, so refuse multi-sample input
  if(doSmartCounts>0&&arguments->nInd!=1){
    fprintf(stderr,"\t-> -doSmartCounts only supports a single sample (one bam file), but %d were given\n",arguments->nInd);
    exit(1);
  }

  if(doSmartCounts==0){
    shouldRun[index]=0;
    return;
  }  
  fprintf(stderr,"\t-> NOTE: -doSmartCounts is deprecated and will be removed in a future version of ANGSD\n");
}

//write counts for the current chromosome to the .bin file and its offset to the .idx file
void abcSmartCounts::writeChr(){
  if(curChr==-1)//no data seen yet
    return;
  int64_t retVal =bgzf_tell(fbin);
  int clen = strlen(header->target_name[curChr]);
  if(aio::bgzf_write(fbin,&clen,sizeof(int))!=(ssize_t)sizeof(int)||
     aio::bgzf_write(fbin,header->target_name[curChr],clen)!=clen||
     aio::bgzf_write(fbin,&len,sizeof(int))!=(ssize_t)sizeof(int)){
    fprintf(stderr,"\t-> Problem writing smartcounts binfile for chr: %s\n",header->target_name[curChr]);
    exit(1);
  }
  for(int i=0;i<4;i++)
    if(aio::bgzf_write(fbin,counts[i],len)!=len){
      fprintf(stderr,"\t-> Problem writing smartcounts binfile for chr: %s\n",header->target_name[curChr]);
      exit(1);
    }

  if(fwrite(&clen,sizeof(int),1,fidx)!=1||
     fwrite(header->target_name[curChr],sizeof(char),clen,fidx)!=(size_t)clen||
     fwrite(&len,sizeof(int),1,fidx)!=1||
     fwrite(&retVal,sizeof(int64_t),1,fidx)!=1){
    fprintf(stderr,"\t-> Problem writing smartcounts indexfile for chr: %s\n",header->target_name[curChr]);
    exit(1);
  }
}

void abcSmartCounts::changeChr(int newRefId){
  if(doSmartCounts==0)
    return;
  writeChr();
  curChr = newRefId;
  if(header->target_len[curChr]>INT_MAX){
    fprintf(stderr,"\t-> Chromosome '%s' is longer than %d bp, which -doSmartCounts does not support\n",header->target_name[curChr],INT_MAX);
    exit(1);
  }
  len = header->target_len[curChr];
  for(int i=0;i<4;i++){
    delete [] counts[i];
    counts[i] = new unsigned char[len];
    memset(counts[i],0,len);
  } 
}



abcSmartCounts::abcSmartCounts(const char *outfiles,argStruct *arguments,int inputtype){

  doSmartCounts =0;
  fidx=NULL;
  fbin=NULL;
  counts=NULL;
  curChr=-1;
  len =-1;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSmartCounts")){
      printArg(stderr);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  if(doSmartCounts==0)
    return;

  printArg(arguments->argumentFile);
  counts = new unsigned char*[4];
  counts[0]=counts[1]=counts[2]=counts[3]=NULL;
  
  const char * BIN= ".counts.bin";
  const char * IDX= ".counts.idx";

  fidx=aio::openFile(outfiles,IDX);
  fbin=aio::openFileBG(outfiles,BIN);

}


abcSmartCounts::~abcSmartCounts(){

  if(doSmartCounts==0)
    return;

  writeChr();
  
  for(int i=0;i<4;i++)
    delete [] counts[i];
  delete [] counts;

  if(fclose(fidx)!=0||bgzf_close(fbin)!=0){
    fprintf(stderr,"\t-> Problem closing smartcounts output files\n");
    exit(1);
  }

}


void abcSmartCounts::clean(funkyPars *pars){
  if(doSmartCounts==0)
    return;
  
}

void abcSmartCounts::print(funkyPars *pars){
  if(doSmartCounts<=0)
    return;

}


void abcSmartCounts::run(funkyPars *pars){
  if(doSmartCounts==0)
    return;
  
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;

    suint *ary=pars->counts[s];
#if 0
    int whMax = 0;
    int whMin = 0;

    for(int i=1;i<3;i++){
      if(ary[i]>whMax)
	whMax =i;
      if(ary[i]<whMin)
	whMax =i;
      
    }
#endif
    int p=pars->posi[s];
    for(int i=0;i<4;i++){
      if(ary[i]>255)
	counts[i][p] = 255;
      else
	counts[i][p] =(char) ary[i];
      //      fprintf(stderr,"%d %d\n",(unsigned char)counts[i][p],ary[i]);
    }
  }
}
