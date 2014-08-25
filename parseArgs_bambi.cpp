#include <vector>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "abc.h"
#include "parseArgs_bambi.h"
#include "shared.h"
#include "analysisFunction.h"
#include "indexer.h"

// below is default samtools parameters
int uniqueOnly = 0;
int only_proper_pairs = 1;
int remove_bads = 1;
int minMapQ =0;
int minQ = MINQ;
int trim = 0;
int adjustMapQ =0;
int baq =0;
int checkBamHeaders = 1;
char *regfile = NULL;
char *regfiles = NULL;

unsigned int includeflags = 2;
unsigned int discardflags = 4;


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
  fprintf(argFile,"\t-r\t\t%s\tSupply a single region in commandline (see examples below)\n",regfile);
  fprintf(argFile,"\t-rf\t\t%s\tSupply multiple regions in a file (see examples below)\n",regfiles);
  fprintf(argFile,"\t-remove_bads\t%d\tDiscard \'bad\' reads, (flag >=255) \n",remove_bads);
  //  fprintf(argFile,"\t-type\t\t%d\n",ret->type);
  fprintf(argFile,"\t-uniqueOnly\t%d\tDiscards reads that doesn't map uniquely\n",uniqueOnly);
  fprintf(argFile,"\t-show\t\t%d\tMimic 'samtools mpileup' also supply -ref fasta for printing reference column\n",ret->show);
  fprintf(argFile,"\t-minMapQ\t%d\tDiscard reads with mapping quality below\n",minMapQ);
  fprintf(argFile,"\t-minQ\t\t%d\tDiscard bases with base quality below\n",minQ);
  fprintf(argFile,"\t-trim\t\t%d\tNumber of based to discard at both ends of the reads\n",trim);
  fprintf(argFile,"\t-only_proper_pairs\t%d\tOnly use reads where the mate could be mapped\n",only_proper_pairs);
  fprintf(argFile,"\t-C\t\t%d\tadjust mapQ for excessive mismatches (as SAMtools), supply -ref\n",adjustMapQ);
  fprintf(argFile,"\t-baq\t\t%d\tadjust qscores around indels (as SAMtools), supply -ref\n",baq);
  fprintf(argFile,"\t-if\t\t%d\tinclude flags for each read\n",includeflags);
  fprintf(argFile,"\t-df\t\t%d\tdiscard flags for each read\n",discardflags);
  fprintf(argFile,"\t-checkBamHeaders\t%d\tExit if difference in BAM headers\n",checkBamHeaders);

  fprintf(argFile,"\n");
  fprintf(argFile,"Examples for region specification:\n");
  fprintf(argFile,"\t\tchr:\t\tUse entire chromosome: chr\n");
  fprintf(argFile,"\t\tchr:start-\tUse region from start to end of chr\n");
  fprintf(argFile,"\t\tchr:-stop\tUse region from beginning of chromosome: chr to stop\n");
  fprintf(argFile,"\t\tchr:start-stop\tUse region from start to stop from chromosome: chr\n");
  fprintf(argFile,"\t\tchr:site\tUse single site on chromosome: chr\n");
  fprintf(argFile,"Will include read if:\n\tincludeflag:[%d] (beta)",includeflags);
  printFlagInfo(argFile,includeflags);
  fprintf(argFile,"Will discard read if:\n\tdiscardflag:[%d] (beta)",discardflags);
  printFlagInfo(argFile,discardflags);
  
}

//read program parameters
void setArgsBam(argStruct *arguments){
  remove_bads = angsd::getArg("-remove_bads",remove_bads,arguments);
  uniqueOnly = angsd::getArg("-uniqueOnly",uniqueOnly,arguments);
  only_proper_pairs =angsd::getArg("-only_proper_pairs",only_proper_pairs,arguments);
  minMapQ = angsd::getArg("-minMapQ",minMapQ,arguments);
  minQ = angsd::getArg("-minQ",minQ,arguments);
  trim = angsd::getArg("-trim",trim,arguments);
  adjustMapQ = angsd::getArg("-C",adjustMapQ,arguments);
  baq = angsd::getArg("-baq",baq,arguments);
  regfile =angsd::getArg("-r",regfile,arguments);
  regfiles = angsd::getArg("-rf",regfiles,arguments);
  arguments->show = angsd::getArg("-show",arguments->show,arguments);

  char *tmp = NULL;
  tmp = angsd::getArg("-ref",tmp,arguments);
  if(tmp==NULL && adjustMapQ!=0){
    fprintf(stderr,"Must also supply -ref for adjusting the mapping quality\n");
    exit(0);
  }
  if(tmp==NULL&&baq!=0){
    fprintf(stderr,"Must also supply -ref for adjusting base qualities (baq)\n");
    exit(0);
  }
  free(tmp);
  
  
  std::vector<char *> regionsRaw;
  if(regfiles)
    regionsRaw =  angsd::getFilenames(regfiles,0);
  
  if(regfile)
    regionsRaw.push_back(regfile);
   
  for(size_t i=0;i<regionsRaw.size();i++){
    regs tmpRegs;
    if(parse_region(regionsRaw[i],arguments->hd,tmpRegs.refID,tmpRegs.start,tmpRegs.stop,arguments->revMap)<0||tmpRegs.stop<tmpRegs.start){
      fprintf(stderr,"[%s] problems with indexing: %s\n",__FUNCTION__,regionsRaw[i]);
      exit(0);
    }else
      arguments->regions.push_back(tmpRegs);
  }



  printArg(arguments->argumentFile,arguments);
  free(regfile);
  free(regfiles);
}
