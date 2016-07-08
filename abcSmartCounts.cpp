/*
  thorfinn thorfinn@binf.ku.dk 1april 2012

 */
#include <cmath>
#include "abcSmartCounts.h"


unsigned char **counts = NULL;


void abcSmartCounts::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSmartCounts\t%d\n",doSmartCounts);
  fprintf(argFile,"\n");
}

void abcSmartCounts::getOptions(argStruct *arguments){

  doSmartCounts = angsd::getArg("-doSmartCounts",doSmartCounts,arguments);

  int tmp=0;
  tmp=angsd::getArg("-doCounts",tmp,arguments);
  if(tmp==0&&doSmartCounts>0){
    fprintf(stderr,"\t-> SmartCounts require -doCounts\n");
    exit(0);
  }

  if(doSmartCounts==0){
    shouldRun[index]=0;
    return;
  }  
}

void abcSmartCounts::changeChr(int newRefId){
  if(doSmartCounts==0)
    return;
  //  fprintf(stderr,"cur:%d new:%d\n",curChr,newRefId);
  if(curChr!=-1){
    int64_t retVal =bgzf_tell(fbin); 
    int clen = strlen(header->target_name[curChr]);
    aio::bgzf_write(fbin,&clen,sizeof(int));
    aio::bgzf_write(fbin,header->target_name[curChr],clen);
    aio::bgzf_write(fbin,&len,sizeof(int));
    for(int i=0;i<4;i++)
      aio::bgzf_write(fbin,counts[i],len);//write len of chr
    
    //write index stuff
    fprintf(stderr,"Writing index for chr: %s\n",header->target_name[curChr]);
    fwrite(&clen,sizeof(int),1,fidx);
    fwrite(header->target_name[curChr] ,sizeof(char),clen,fidx);
    fwrite(&len,sizeof(int),1,fidx);
    fwrite(&retVal,sizeof(int64_t),1,fidx);
  }
  curChr = newRefId;
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
  curChr=-1;
  len =-1;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSmartCounts")){
      printArg(stdout);
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

  int64_t retVal =bgzf_tell(fbin); 
  int clen = strlen(header->target_name[curChr]);
  aio::bgzf_write(fbin,&clen,sizeof(int));
  aio::bgzf_write(fbin,header->target_name[curChr],clen);
  aio::bgzf_write(fbin,&len,sizeof(int));
  for(int i=0;i<4;i++)
    aio::bgzf_write(fbin,counts[i],len);//write len of chr
  
  //write index stuff
  fwrite(&clen,sizeof(int),1,fidx);
  fwrite(header->target_name[curChr],sizeof(char),clen,fidx);
  fwrite(&len,sizeof(int),1,fidx);
  fwrite(&retVal,sizeof(int64_t),1,fidx);

  
  for(int i=0;i<4;i++)
    delete [] counts[i];
  delete [] counts;

  fclose(fidx);
  bgzf_close(fbin);

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
