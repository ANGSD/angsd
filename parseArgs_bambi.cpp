#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "abc.h"
#include "parseArgs_bambi.h"
#include "shared.h"
#include "analysisFunction.h"
#include "from_samtools.h"

// below is default samtools parameters
int MPLP_IGNORE_RG = 1;
int uniqueOnly = 0;
int only_proper_pairs = 1;
int remove_bads = 1;
int minMapQ =0;
int minQ = MINQ;
float downSample = 0;
int trim = 0;
int trim5 = 0;
int trim3 = 0;
int adjustMapQ =0;
int baq =0;
int checkBamHeaders = 1;
int doCheck = 1;
int MAX_SEQ_LEN = 250;
char *regfile = NULL;
char *regfiles = NULL;
char *fai_fname = NULL;
unsigned int includeflags = 2;
unsigned int discardflags = 4;
int redo_baq =0;
int cigstat =0;
void *rghash=NULL;
char *rghash_name=NULL;


int parse_region(char *extra,const bam_hdr_t *hd,int &ref,int &start,int &stop,const aMap *revMap) {
  aMap::const_iterator it;
  if(strrchr(extra,':')==NULL){//only chromosomename
    if((it = revMap->find(extra))==revMap->end()){
      fprintf(stderr,"[%s.%s():%d] Problems finding chromosome: \'%s\'\n",__FILE__,__FUNCTION__,__LINE__,extra);
      fflush(stderr);
      exit(0);
      return -1;
    }
    ref = it->second;
    start =0;
    stop = hd->target_len[ref];
    return 1;
  }
  
  char *tok = strtok(extra,":");
  if((it =revMap->find(tok))==revMap->end()){
    fprintf(stderr,"[%s.%s():%d] (-r) Problems finding chromosome: \'%s\' from header\n",__FILE__,__FUNCTION__,__LINE__,extra);
    fflush(stderr);
    return -1;
  }
  ref = it->second;
  
  start =0;
  stop = hd->target_len[ref];
  tok = extra+strlen(tok)+1;//tok now contains the rest of the string
  
  if(strlen(tok)==0)//not start and/or stop ex: chr21:
    return 1;
  
  
  if(tok[0]=='-'){//only contains stop ex: chr21:-stop
    tok =strtok(tok,"-");
    
    stop = atoi(tok);
  }else{
    //catch single point
    int isProper =0;
    for(size_t i=0;i<strlen(tok);i++)
      if(tok[i]=='-'){
	isProper=1;
	break;
      }
    //fprintf(stderr,"isProper=%d\n",isProper);
    if(isProper){
      tok =strtok(tok,"-");
      start = atoi(tok)-1;//this is important for the zero offset
      tok = strtok(NULL,"-");
      if(tok!=NULL)
	stop = atoi(tok);
    }else{
      //single point
      stop = atoi(tok);
      start =stop -1;
      
    }
    
  }
  if(stop<start){
    fprintf(stderr,"endpoint:%d is larger than startpoint:%d\n",start,stop);
    exit(0);
    
  }
  if(0){
    fprintf(stderr,"[%s] ref=%d,start=%d,stop=%d\n",__FUNCTION__,ref,start,stop);
    exit(0);
  }
  return 1;
}



void printFlagInfo(FILE *fp,unsigned int f){
  if(f&1)
    fprintf(fp,"template having multiple segments in sequencing, ");
  if(f&2)
    fprintf(fp,"each segment properly aligned according to the aligner, ");
  if(f&4)
    fprintf(fp,"segment unmapped, ");
  if(f&8)
    fprintf(fp,"next segment in the template unmapped, ");
  if(f&16)
    fprintf(fp,"SEQ being reverse complemented, ");
  if(f&32)
    fprintf(fp,"SEQ of the next segment in the templete being reverse complemented, ");
  if(f&64)
    fprintf(fp,"the first segment in the template, ");
  if(f&128)
    fprintf(fp,"the last seggment in the template, ");
  if(f&256)
    fprintf(fp,"secondary alignment, ");
  if(f&512)
    fprintf(fp,"not passing quality controls, ");
  if(f&1024)
    fprintf(fp,"PCR or optical duplicate, ");
  if(f&2048)
    fprintf(fp,"supplementary alignment, ");
  fprintf(fp,"\n");
}


void printArg(FILE *argFile,argStruct *ret){
  fprintf(argFile,"---------------\n%s: bam reader:\n",__FILE__);
  fprintf(argFile,"\t-bam/-b\t\t%s\t(list of BAM/CRAM files)\n",ret->infile);
  fprintf(argFile,"\t-i\t\t%s\t(Single BAM/CRAM file)\n",ret->infile);
  fprintf(argFile,"\t-r\t\t%s\tSupply a single region in commandline (see examples below)\n",regfile);
  fprintf(argFile,"\t-rf\t\t%s\tSupply multiple regions in a file (see examples below)\n",regfiles);
  fprintf(argFile,"\t-remove_bads\t%d\tDiscard \'bad\' reads, (flag >=256) \n",remove_bads);
  //  fprintf(argFile,"\t-type\t\t%d\n",ret->type);
  fprintf(argFile,"\t-uniqueOnly\t%d\tDiscards reads that doesn't map uniquely\n",uniqueOnly);
  fprintf(argFile,"\t-show\t\t%d\tMimic 'samtools mpileup' also supply -ref fasta for printing reference column\n",ret->show);
  fprintf(argFile,"\t-minMapQ\t%d\tDiscard reads with mapping quality below\n",minMapQ);
  fprintf(argFile,"\t-minQ\t\t%d\tDiscard bases with base quality below\n",minQ);
  fprintf(argFile,"\t-trim\t\t%d\tNumber of based to discard at both ends of the reads\n",trim);
  fprintf(argFile,"\t-trim\t\t%d\tNumber of based to discard at 5' ends of the reads\n",trim5);
  fprintf(argFile,"\t-trim\t\t%d\tNumber of based to discard at 3' ends of the reads\n",trim3);
  fprintf(argFile,"\t-only_proper_pairs %d\tOnly use reads where the mate could be mapped\n",only_proper_pairs);
  fprintf(argFile,"\t-C\t\t%d\tadjust mapQ for excessive mismatches (as SAMtools), supply -ref\n",adjustMapQ);
  fprintf(argFile,"\t-baq\t\t%d\tadjust qscores around indels (1=normal baq 2= extended(as SAMtools)), supply -ref\n",baq);
  fprintf(argFile,"\t-redo-baq\t\t%d (recompute baq, instead of using BQ tag)\n",redo_baq);
  //  fprintf(argFile,"\t-if\t\t%d\tinclude flags for each read\n",includeflags);
  // fprintf(argFile,"\t-df\t\t%d\tdiscard flags for each read\n",discardflags);
  fprintf(argFile,"\t-checkBamHeaders %d\tExit if difference in BAM headers\n",checkBamHeaders);
  fprintf(argFile,"\t-doCheck\t%d\tKeep going even if datafile is not suffixed with .bam/.cram\n",doCheck);
  fprintf(argFile,"\t-downSample\t%f\tDownsample to the fraction of original data\n",downSample);
  fprintf(argFile,"\t-nReads\t\t%d\tNumber of reads to pop from each BAM/CRAMs\n",ret->nReads);
  fprintf(argFile,"\t-minChunkSize\t%d\tMinimum size of chunk sent to analyses\n",MAX_SEQ_LEN);
  fprintf(argFile,"\t--ignore-RG\t%d\t(dev only)\n",MPLP_IGNORE_RG);
  fprintf(argFile,"\t+RG\t%s\tReadgroups to include in analysis(can be filename)\n",rghash_name);
  
  fprintf(argFile,"\n");
  fprintf(argFile,"Examples for region specification:\n");
  fprintf(argFile,"\t\tchr:\t\tUse entire chromosome: chr\n");
  fprintf(argFile,"\t\tchr:start-\tUse region from start to end of chr\n");
  fprintf(argFile,"\t\tchr:-stop\tUse region from beginning of chromosome: chr to stop\n");
  fprintf(argFile,"\t\tchr:start-stop\tUse region from start to stop from chromosome: chr\n");
  fprintf(argFile,"\t\tchr:site\tUse single site on chromosome: chr\n");
  //fprintf(argFile,"Will include read if:\n\tincludeflag:[%d] (beta)",includeflags);
  //printFlagInfo(argFile,includeflags);
  //fprintf(argFile,"Will discard read if:\n\tdiscardflag:[%d] (beta)",discardflags);
  //printFlagInfo(argFile,discardflags);
  
}

//read program parameters
void setArgsBam(argStruct *arguments){
  int seed=0;
  remove_bads = angsd::getArg("-remove_bads",remove_bads,arguments);
  uniqueOnly = angsd::getArg("-uniqueOnly",uniqueOnly,arguments);
  only_proper_pairs =angsd::getArg("-only_proper_pairs",only_proper_pairs,arguments);
   fai_fname =angsd::getArg("-f",fai_fname,arguments);
  minMapQ = angsd::getArg("-minMapQ",minMapQ,arguments);
  cigstat = angsd::getArg("-cigstat",cigstat,arguments);
  minQ = angsd::getArg("-minQ",minQ,arguments);
  downSample = angsd::getArg("-downSample",downSample,arguments);
  seed = angsd::getArg("-seed",seed,arguments);
  trim = angsd::getArg("-trim",trim,arguments);
  trim5 = angsd::getArg("-trim5",trim5,arguments);
  trim3 = angsd::getArg("-trim3",trim3,arguments);
  arguments->ref=angsd::getArg("-ref",arguments->ref,arguments);
  arguments->anc=angsd::getArg("-anc",arguments->anc,arguments);
  rghash_name= angsd::getArg("+RG",rghash_name,arguments);
  if(rghash_name&&!angsd::fexists(rghash_name))
    rghash = add_read_group_single(rghash_name);
  if(rghash_name&&angsd::fexists(rghash_name))
    rghash = add_read_groups_file(rghash_name);
  if(rghash)
    fprintf(stderr,"\t-> [READGROUP info] Number of readgroups to include: %d\n",khash_str2int_size(rghash));
  adjustMapQ = angsd::getArg("-C",adjustMapQ,arguments);
  baq = angsd::getArg("-baq",baq,arguments);
  redo_baq = angsd::getArg("-redo-baq",redo_baq,arguments);
  if(baq){
    if(baq==1)
      baq=1; //wauv
    else if(baq==2)
      baq=3;
    else{
      fprintf(stderr,"\t-> only supported options for -baq is: 1 (normal baq) and 2 (extended baq (SAMtools default)). Value supplied:%d\n",baq);
      exit(0);//ly su
    }
    if(redo_baq==1)
      baq |=4;
  }
  //  fprintf(stderr,"baq:%d redobaq:%d\n",baq,redo_baq);exit(0);
  regfile =angsd::getArg("-r",regfile,arguments);
  regfiles = angsd::getArg("-rf",regfiles,arguments);
  MAX_SEQ_LEN = angsd::getArg("-setMinChunkSize",MAX_SEQ_LEN,arguments);
  checkBamHeaders = angsd::getArg("-checkBamHeaders",checkBamHeaders,arguments);
  doCheck = angsd::getArg("-doCheck",doCheck,arguments);
  MPLP_IGNORE_RG = angsd::getArg("--ignore-RG",MPLP_IGNORE_RG,arguments);
  arguments->nReads = angsd::getArg("-nReads",arguments->nReads,arguments);
  arguments->show = angsd::getArg("-show",arguments->show,arguments);
  if(regfile && regfiles)
    fprintf(stderr,"\t-> WARNING both -r and -rf has been set \n");

  if(seed)
    srand48(seed);
  char *tmp = NULL;
  tmp = angsd::getArg("-ref",tmp,arguments);
  if(tmp==NULL && adjustMapQ!=0){
    fprintf(stderr,"\t-> Must also supply -ref for adjusting the mapping quality\n");
    exit(0);
  }
  if(tmp==NULL&&baq!=0){
    fprintf(stderr,"\t-> Must also supply -ref for adjusting base qualities (baq)\n");
    exit(0);
  }
  free(tmp);
  
  
  std::vector<char *> regionsRaw;
  if(regfiles)
    regionsRaw =  angsd::getFilenames(regfiles,0);
  
  if(regfile)
    regionsRaw.push_back(strdup(regfile));
  //  fprintf(stderr,"\t-> RegionsRaw.size():%lu hd:%p\n",regionsRaw.size(),arguments->hd);
  for(size_t i=0;i<regionsRaw.size();i++){
    regs tmpRegs;
    if(parse_region(regionsRaw[i],arguments->hd,tmpRegs.refID,tmpRegs.start,tmpRegs.stop,arguments->revMap)<0||tmpRegs.stop<tmpRegs.start){
      fprintf(stderr,"[%s] Problems with indexing: %s\n",__FUNCTION__,regionsRaw[i]);
      exit(0);
    }else
      arguments->regions.push_back(tmpRegs);
    free(regionsRaw[i]);
  }



  printArg(arguments->argumentFile,arguments);

  if(regfile)
    free(regfile);
  if(regfiles)
    free(regfiles);

}
