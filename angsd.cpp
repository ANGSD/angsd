/*
  copyright (gpl)   19 sep 2012.
  Thorfinn Sand Korneliussen thorfinn@binf.ku.dk
  Anders albrechtsen albrecht@binf.ku.dk 
 
  use wiki for updated list of bugfixes and documentation
  http://www.popgen.dk/angsd

*/

#include "version.h"
#include <cassert>
#include <iostream>//for printing time
#include <cstring> //for cstring functions
#include <cstdlib> //for exit()
#include <cstdio> //for fprintf
#include <signal.h>//for catching ctrl+c, allow threads to finish
#include <htslib/hts.h>
#include "cigstat.h"
#include "shared.h"
#include "multiReader.h"


extern std::vector <char *> dumpedFiles;

int SIG_COND =1;//if we catch signal then quit program nicely
int VERBOSE =1;

extern pthread_mutex_t mUpPile_mutex;//mutex for not deleting mUppile data untill end of data

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

int really_kill =3;

void handler(int s) {

  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want angsd to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;
  pthread_mutex_unlock(&mUpPile_mutex);
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
  fprintf(fp,"\n\t-> angsd version: %s (htslib: %s) build(%s %s)\n",ANGSD_VERSION,hts_version(),__DATE__,__TIME__); 
  fprintf(fp,"\t-> Please use the website \"http://www.popgen.dk/angsd\" as reference\n");
  fprintf(fp,"\t-> Use -nThreads or -P for number of threads allocated to the program\n"); 
  fprintf(fp,"Overview of methods:\n");
  fprintf(fp,"\t-GL\t\tEstimate genotype likelihoods\n");
  fprintf(fp,"\t-doCounts\tCalculate various counts statistics\n");
  fprintf(fp,"\t-doAsso\t\tPerform association study\n");
  fprintf(fp,"\t-doMaf\t\tEstimate allele frequencies\n");
  fprintf(fp,"\t-doError\tEstimate the type specific error rates\n");
  fprintf(fp,"\t-doAncError\tEstimate the errorrate based on perfect fastas\n");
  fprintf(fp,"\t-HWE_pval\tEst inbreedning per site or use as filter\n");
  fprintf(fp,"\t-doGeno\t\tCall genotypes\n");
  fprintf(fp,"\t-doFasta\tGenerate a fasta for a BAM file\n");
  fprintf(fp,"\t-doAbbababa\tPerform an ABBA-BABA test\n");
  fprintf(fp,"\t-sites\t\tAnalyse specific sites (can force major/minor)\n");
  fprintf(fp,"\t-doSaf\t\tEstimate the SFS and/or neutrality tests genotype calling\n");
  fprintf(fp,"\t-doHetPlas\tEstimate hetplasmy by calculating a pooled haploid frequency\n");
  fprintf(fp,"\n\tBelow are options that can be usefull\n");
  fprintf(fp,"\t-bam\t\tOptions relating to bam reading\n\t-doMajorMinor\tInfer the major/minor using different approaches\n");  
  fprintf(fp,"\t-ref/-anc\tRead reference or ancestral genome\n");
  fprintf(fp,"\t-doSNPstat\tCalculate various SNPstat\n");
  fprintf(fp,"\t-cigstat\tPrintout CIGAR stat across readlength\n");
  fprintf(fp,"\tmany others\n\n");
  fprintf(fp,"Output files:\n\t In general the specific analysis outputs specific files, but we support basic bcf output\n");
  fprintf(fp,"\t-doBcf\t\tWrapper around -dopost -domajorminor -dofreq -gl -dovcf docounts\n");
  
  fprintf(fp,"For information of specific options type: \n\t./angsd METHODNAME eg \n\t\t./angsd -GL\n\t\t./angsd -doMaf\n\t\t./angsd -doAsso etc\n");
  fprintf(fp,"\t\t./angsd sites for information about indexing -sites files\n");
  fprintf(fp,"Examples:\n\tEstimate MAF for bam files in 'list'\n\t\t\'./angsd -bam list -GL 2 -doMaf 2 -out RES -doMajorMinor 1\'\n");
  
}




void parseArgStruct(argStruct *arguments){

  //validate that there are no duplicate parameters
  for(int i=0;i<arguments->argc;i++){
    if(strcasecmp("-realSFS",arguments->argv[i])==0){
	fprintf(stderr,"-realSFS is now called -doSaf\n");
	exit(0);
    }
    if(strcasecmp("-minLRT",arguments->argv[i])==0){
	fprintf(stderr,"-minLRT is now called -SNP_pval\n");
	exit(0);
    }
    if(strcasecmp("-doSNP",arguments->argv[i])==0){
	fprintf(stderr,"-doSNP is now called -SNP_pval\n");
	exit(0);
    }  
    if(strcasecmp("-soap",arguments->argv[i])==0){
      fprintf(stderr,"-doHWE is deprecated, convert soap file to bam\n");
      exit(0);
    }  
    if(strcasecmp("-tglf",arguments->argv[i])==0){
      fprintf(stderr,"-tglf is deprecated, merge files and use -glf\n");
      exit(0);
    }
    if(strcasecmp("-sim1",arguments->argv[i])==0){
      fprintf(stderr,"-sim1 is deprecated, use -glf filename and -isSim 1\n");
      exit(0);
    }  
  }
  for(int i=0;i<arguments->argc-1;i++)
    if(arguments->argv[i][0]=='-') {
      for(int ii=i+1;ii<arguments->argc;ii++){
	if(arguments->argv[ii][0]=='-'){
	  if(0==strcasecmp(arguments->argv[i],arguments->argv[ii])){
	    fprintf(stderr,"\t-> Duplicate parameter: %s  supplied, will exit\n",arguments->argv[i]);
	    exit(0);
	  }
	}
       }
    }

  for(int i=1;i<arguments->argc;i++){
    if(arguments->usedArgs[i]==0){
      fprintf(stderr,"%d argument \t%s is unknown will exit\n",i,arguments->argv[i]);
      fflush(stderr);
      exit(0);
    }
  }
  if(arguments->hd==NULL){
    fprintf(stderr,"\t-> Error: You must supply a -fai file such that we know the sizes of the genomes (version .4486)\n");
     exit(0);
  }


}

int main(int argc, char** argv){
   fprintf(stderr,"\t-> angsd version: %s (htslib: %s) build(%s %s)\n",ANGSD_VERSION,hts_version(),__DATE__,__TIME__); 
  //no arguments supplied -> print info
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"--help")==0))){//if haven't been supplied with arguments, load default,print, and exit
    printProgInfo(stderr);
    return 0;
  }
   
   if(!strcasecmp("sites",argv[1])){
     //from prep_sites.* used for indexing -sites files
     int main_sites(int argc,char **argv);
     main_sites(--argc,++argv);
     return 0;
   }

   //print time
   clock_t t=clock();
   time_t t2=time(NULL);

   //intialize our signal handler for ctrl+c
   catchkill();

   argStruct *args=NULL;
   if(argc==2){
     fprintf(stderr,"\t-> Analysis helpbox/synopsis information:\n");
     multiReader mr(argc,argv);
     args = mr.getargs();
     
     init(args);//program dies here after printing info, if a match is found
     fprintf(stderr,"\nUnknown argument supplied: \'%s\'\n\n",argv[1]);
     printProgInfo(stderr);
     exit(0);//important otherwise the abc classes will try to clean up, which doesnt make sense in this context
     //     return 0;
   }

   multiReader *mr= new multiReader(argc,argv);
   args = mr->getargs();

   

   init(args);
   parseArgStruct(args);
   
   //Below is main loop which will run until nomore data
   assert(args->hd);
   assert(args->revMap);

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
       void changeChr(int refId);//<-ugly is located in shared.cpp, to force a change of chr, when notusing bamfiles
       changeChr(lastRefId);
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
  fprintf(stderr,"\t-> Total number of sites analyzed: %lu\n\t-> Number of sites retained after filtering: %lu \n",total_number_of_sites_unfiltred,total_number_of_sites_filtered);
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(args->argumentFile, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(args->argumentFile, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  
  //make better below
  extern int *bamSortedIds;
  delete [] bamSortedIds;
  
  delete mr;
  //check
  extern htsFormat *dingding2;
  free(dingding2);
  return 0;
}
