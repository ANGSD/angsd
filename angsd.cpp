/*
  copyright (gpl)   19 sep 2012.
  Thorfinn Sand Korneliussen thorfinn@binf.ku.dk
  Anders albrechtsen albrecht@binf.ku.dk 
 
  use wiki for updated list of bugfixes and documentation
  https://github.com/ANGSD/angsd

*/

#include "version.h"
#include <iostream>//for printing time
#include <cstring> //for cstring functions
#include <cstdlib> //for exit()
#include <cstdio> //for fprintf
#include <signal.h>//for catching ctrl+c, allow threads to finish
#include <unistd.h>//for write, _exit and unlink
#include <htslib/hts.h>
#include "cigstat.h"
#include "pop1_read.h"
#include "shared.h"
#include "argStruct.h"
#include "multiReader.h"
#include "aio.h"


extern std::vector <char *> dumpedFiles;

int SIG_COND =1;//if we catch signal then quit program nicely
//set by the signal handler only, SIG_COND is also set to 0 by some errors
volatile sig_atomic_t gotSignal =0;

void printTime(FILE *fp){
  time_t rawtime;
  struct tm * timeinfo; 
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}




//this function is called from within the bamreader
void callBack_bambi(fcb *fff){

  if(fff==NULL){
    //    fprintf(stderr,"SEnding NULL this is a killswitch");
    selector(NULL);//<-send NULL which acts as a killswitch
  }else{
    funkyPars *fp = funkyPars_init();
    fp->for_callback=fff;
    fp->refId = fp->for_callback->refId;
    selector(fp);
  }
}

//the first signal stops reading new data and lets the running threads finish,
//the third exits right away. Only async-signal-safe calls are used here (write, _exit)
volatile sig_atomic_t nSignals =0;

static void sigwrite(const char *msg){
  ssize_t ret = write(STDERR_FILENO,msg,strlen(msg));
  (void) ret;
}

void handler(int s) {
  nSignals++;
  if(nSignals==1)
    sigwrite("\n\t-> Caught SIGNAL: Will try to exit nicely (no more data is read.\n\t\t\t  We will wait for the current threads to finish)\n\t-> If you really want angsd to exit uncleanly ctrl+c: 2 more times\n");
  else if(nSignals==2)
    sigwrite("\n\t-> If you really want angsd to exit uncleanly ctrl+c: 1 more time\n");
  else
    _exit(130);
  gotSignal=1;
  SIG_COND=0;
}

 //we are threading so we want make a nice signal handler for ctrl+c
void catchkill(){
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

}

//print nice info
void printProgInfo(FILE *fp){
  //the version is already printed by setArgStruct
  fprintf(fp,"\n\t-> Please use the website \"https://github.com/ANGSD/angsd\" as reference\n");
  fprintf(fp,"\t-> Use -nThreads or -P for number of threads allocated to the program\n"); 
  fprintf(fp,"Overview of methods:\n");
  fprintf(fp,"\t-GL\t\tEstimate genotype likelihoods\n");
  fprintf(fp,"\t-doCounts\tCalculate various counts statistics\n");
  fprintf(fp,"\t-doAsso\t\tPerform association study\n");
  fprintf(fp,"\t-doMaf\t\tEstimate allele frequencies\n");
  fprintf(fp,"\t-doError\tEstimate the type specific error rates\n");
  fprintf(fp,"\t-doAncError\tEstimate the errorrate based on perfect fastas\n");
  fprintf(fp,"\t-doPost\t\tCalculate genotype posteriors (shown with -doMaf)\n");
  fprintf(fp,"\t-doHWE\t\tEstimate inbreeding per site and test for HWE\n");
  fprintf(fp,"\t-HWE_pval_F\tFilter sites on a per site inbreeding test\n");
  fprintf(fp,"\t-doGeno\t\tCall genotypes\n");
  fprintf(fp,"\t-doFasta\tGenerate a fasta for a BAM file\n");
  fprintf(fp,"\t-doAbbababa\tPerform an ABBA-BABA test\n");
  fprintf(fp,"\t-doIBS\t\tSample a single base or the consensus base per individual (IBS/cov matrix)\n");
  fprintf(fp,"\t-doHaploCall\tCall haploid genotypes\n");
  fprintf(fp,"\t-sites\t\tAnalyse specific sites (can force major/minor)\n");
  fprintf(fp,"\t-doSaf\t\tEstimate the site allele frequency likelihoods (input for realSFS)\n");
  fprintf(fp,"\t-doHetPlas\tEstimate hetplasmy by calculating a pooled haploid frequency\n");
  fprintf(fp,"\n\tBelow are options that can be useful\n");
  fprintf(fp,"\t-bam\t\tOptions relating to bam reading (also -cigstat)\n\t-doMajorMinor\tInfer the major/minor using different approaches\n");  
  fprintf(fp,"\t-ref/-anc\tRead reference or ancestral genome\n");
  fprintf(fp,"\t-beagle\t\tOther input formats (-beagle, -glf, -vcf-GL, -pileup etc)\n");
  fprintf(fp,"\t-doSNPstat\tCalculate various SNPstat\n");
  fprintf(fp,"\t-doMcall\tCall variants with the bcftools multiallelic model\n");
  fprintf(fp,"\tmany others\n\n");
  fprintf(fp,"Output files:\n\t In general the specific analysis outputs specific files, but we support basic bcf output\n");
  fprintf(fp,"\t-doBcf\t\tWrapper around -dopost -domajorminor -domaf -gl -docounts\n");
  fprintf(fp,"\t-doGlf\t\tDump genotype likelihoods (shown with -GL)\n");
  fprintf(fp,"\t-doPlink\tWrite genotypes as tped/tfam\n");
  
  fprintf(fp,"For information of specific options type: \n\t./angsd METHODNAME eg \n\t\t./angsd -GL\n\t\t./angsd -doMaf\n\t\t./angsd -doAsso etc\n");
  fprintf(fp,"\t\t./angsd sites for information about indexing -sites files\n");
  fprintf(fp,"Examples:\n\tEstimate MAF for bam files in 'list'\n\t\t\'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1\'\n");
  
}




//output files are created when the analyses are initialized, remove them again if we exit due to bad arguments
void removeOutputFiles(){
  for(size_t i=0;i<dumpedFiles.size();i++)
    unlink(dumpedFiles[i]);
}

//checks that don't need the analyses, done before they make their output files (only the .arg exists)
void checkArgs(argStruct *arguments){
  for(int i=0;i<arguments->argc;i++){
    if(strcasecmp("-realSFS",arguments->argv[i])==0){
	fprintf(stderr,"-realSFS is now called -doSaf\n");
	removeOutputFiles();
	exit(1);
    }
    if(strcasecmp("-minLRT",arguments->argv[i])==0){
	fprintf(stderr,"-minLRT is now called -SNP_pval\n");
	removeOutputFiles();
	exit(1);
    }
    if(strcasecmp("-doSNP",arguments->argv[i])==0){
	fprintf(stderr,"-doSNP is now called -SNP_pval\n");
	removeOutputFiles();
	exit(1);
    }  
    if(strcasecmp("-soap",arguments->argv[i])==0){
      fprintf(stderr,"-soap is deprecated, convert soap file to bam\n");
      removeOutputFiles();
      exit(1);
    }  
    if(strcasecmp("-tglf",arguments->argv[i])==0){
      fprintf(stderr,"-tglf is deprecated, merge files and use -glf\n");
      removeOutputFiles();
      exit(1);
    }
    if(strcasecmp("-doVcf",arguments->argv[i])==0){
      fprintf(stderr,"-doVcf has been removed, use -doBcf 1 (eg. -doBcf 1 -doMajorMinor 1 -doPost 1 -gl 1 -doMaf 1 -doCounts 1 -doGeno 1)\n");
      removeOutputFiles();
      exit(1);
    }
    if(strcasecmp("-sim1",arguments->argv[i])==0){
      fprintf(stderr,"-sim1 is deprecated, use -glf filename and -isSim 1\n");
      removeOutputFiles();
      exit(1);
    }  
  }
  for(int i=0;i<arguments->argc-1;i++)
    if(angsd::isFlag(arguments->argv[i])) {//negative numbers such as -1 are values, not flags
      for(int ii=i+1;ii<arguments->argc;ii++){
	if(angsd::isFlag(arguments->argv[ii])){
	  if(0==strcasecmp(arguments->argv[i],arguments->argv[ii])){
	    fprintf(stderr,"\t-> Duplicate parameter: %s  supplied, will exit\n",arguments->argv[i]);
	    removeOutputFiles();
	    exit(1);
	  }
	}
       }
    }
}

//checks that need the analyses to have read their arguments
void parseArgStruct(argStruct *arguments){
  for(int i=1;i<arguments->argc;i++){
    if(arguments->usedArgs[i]==0){
      fprintf(stderr,"%d argument \t%s is unknown will exit\n",i,arguments->argv[i]);
      fflush(stderr);
      removeOutputFiles();
      exit(1);
    }
  }
  if(arguments->hd==NULL){
    fprintf(stderr,"\t-> Error: You must supply a -fai file such that we know the sizes of the genomes (version .4486)\n");
    removeOutputFiles();
    exit(1);
  }


}

int main(int argc, char** argv){
  if(argc>1&&!strcasecmp("sites",argv[1])){
    //from prep_sites.* used for indexing -sites files
    int main_sites(int argc,char **argv);
    return main_sites(--argc,++argv);
  }
  //users keep asking about the -fold argument in main angsd. So lets scan and exit with a proper msg
  for(int i=0;i<argc;i++)
    if(!strcasecmp(argv[i],"-fold")){
      fprintf(stderr,"\t-> The folding has been moved from main angsd to realSFS subfunction.\n\t-> Remove -fold when doing -dosaf but apply it in realSFS\n\t-> Will exit\n");
      return 1;
    }
  //same for -doThetas, which was moved to realSFS saf2theta in 0.934
  for(int i=0;i<argc;i++)
    if(!strcasecmp(argv[i],"-doThetas")){
      fprintf(stderr,"\t-> -doThetas has been removed from main angsd, thetas are now made from the saf file\n");
      fprintf(stderr,"\t-> Run -doSaf 1 without -doThetas and -pest, estimate the sfs with realSFS and then run:\n");
      fprintf(stderr,"\t->   realSFS saf2theta out.saf.idx -sfs out.sfs -outname out\n");
      fprintf(stderr,"\t->   thetaStat do_stat out.thetas.idx\n");
      fprintf(stderr,"\t-> Will exit\n");
      return 1;
    }
  
  argStruct *args=setArgStruct(argc,argv);  

  //no arguments supplied -> print info
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"--help")==0))){//if haven't been supplied with arguments, load default,print, and exit
    printProgInfo(stderr);
    return 0;
  }
  

  //print time
  clock_t t=clock();
  time_t t2=time(NULL);

   //intialize our signal handler for ctrl+c
   catchkill();
   
   if(argc==2){
     fprintf(stderr,"\t-> Analysis helpbox/synopsis information:\n");
     multiReader mr(args);
     init((void*)args);//program dies here after printing info, if a match is found,void* to avoid circular .h inclusion
     fprintf(stderr,"\nUnknown argument supplied: \'%s\'\n\n",argv[1]);
     printProgInfo(stderr);
     exit(0);//important otherwise the abc classes will try to clean up, which doesnt make sense in this context
     //     return 0;
   }

   checkArgs(args);
   multiReader *mr= new multiReader(args);

   init(args);
   parseArgStruct(args);
   //write the options to the .arg file now, otherwise it stays empty until the end of a long run
   fflush(args->argumentFile);

   //Below is main loop which will run until nomore data
   aio::doAssert(args->hd!=NULL,AT);
   aio::doAssert(args->revMap!=NULL,AT);

   extern int cigstat;
   if(cigstat)
     cigstat_init(args->outfiles);
   
   int lastRefId=-1;
   while(SIG_COND) {
     funkyPars *tmp = mr->fetch(); //YES MISTER FETCH
      if(tmp==NULL)
	break;
     if(lastRefId==-1||lastRefId!=tmp->refId){
       lastRefId=tmp->refId;
       //wait for the threads to finish the chunks of the previous chromosome before the modules change chromosome,
       //otherwise e.g. the last sites of a chromosome end up in the .saf.idx entry of the next (issue 598)
       void waiter(int refId);//located in shared.cpp
       waiter(lastRefId);
     }
    
     selector(tmp);
      
   }


   //now we have no more data lets wait and cleanup

  fprintf(stderr,"\t-> Done reading data waiting for calculations to finish\n");

  destroy_shared();

  
  //printout the filenames generated
  fprintf(stderr,"\t-> Output filenames:\n");
  for(int i=0;i<(int)dumpedFiles.size();i++){
    fprintf(stderr,"\t\t->\"%s\"\n",dumpedFiles[i]);
    fprintf(args->argumentFile,"\t\t->\"%s\"\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  // fprintf(stderr,"\n");
  fprintf(args->argumentFile,"\n");
  // fprintf(stderr,"\t");
  printTime(stderr);

  
  
  //print out nice messages
  extern size_t total_number_of_sites_unfiltred,total_number_of_sites_filtered;
  fprintf(stderr,"\t-> Arguments and parameters for all analysis are located in .arg file\n");
  //the read counts also go in the .arg file, which is what people share when something is wrong
  printReadCounts(stderr);
  printReadCounts(args->argumentFile);
  fprintf(stderr,"\t-> Total number of sites analyzed: %zu\n\t-> Number of sites retained after filtering: %zu \n",total_number_of_sites_unfiltred,total_number_of_sites_filtered);
  if(gotSignal){
    fprintf(stderr,"\t-> Interrupted by a signal, the output files are incomplete\n");
    fprintf(args->argumentFile,"\t-> Interrupted by a signal, the output files are incomplete\n");
  }else{
    fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
    fprintf(args->argumentFile, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
    fprintf(args->argumentFile, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  }
  
  //make better below
  extern int *bamSortedIds;
  delete [] bamSortedIds;
  
  delete mr;
  //check
  extern htsFormat *dingding2;
  free(dingding2);
  destroy_argStruct(args);
  //128+SIGINT, the usual exit code after ctrl+c
  return gotSignal?130:0;
}
