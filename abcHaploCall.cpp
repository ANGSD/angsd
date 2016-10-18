/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
 
    
  anders albrecht@binf.ku.dk made this.

  part of angsd
*/

#include <cmath>
#include <cstdlib>
#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "abcHaploCall.h"

void abcHaploCall::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doHaploCall\t%d\n",doHaploCall);
  fprintf(argFile,"\t(Sampling strategies)\n");
  fprintf(argFile,"\t 0:\t no haploid calling \n"); 
  fprintf(argFile,"\t 1:\t (Sample single base)\n");
  fprintf(argFile,"\t 2:\t (Concensus base)\n");
  fprintf(argFile,"\t-doCounts\t%d\tMust choose -doCount 1\n",doCount);
  fprintf(argFile,"Optional\n");
  fprintf(argFile,"\t-minMinor\t%d\tMinimum observed minor alleles\n",minMinor);
  fprintf(argFile,"\t-maxMis\t%d\tMaximum missing bases (per site)\n",maxMis);

 }




void abcHaploCall::getOptions(argStruct *arguments){

  //from command line
  doHaploCall=angsd::getArg("-doHaploCall",doHaploCall,arguments);
  if(doHaploCall==0)
    return;

  doCount=angsd::getArg("-doCounts",doCount,arguments);
  minMinor=angsd::getArg("-minMinor",minMinor,arguments);
  maxMis=angsd::getArg("-maxMis",maxMis,arguments);



  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"\t-> bam or pileup input needed for -doHaploCall \n");
    exit(0);
  }
  if( doCount==0){
    fprintf(stderr,"\t-> -doHaploCalls needs allele counts (use -doCounts 1)\n");
    exit(0);
  }

}

abcHaploCall::abcHaploCall(const char *outfiles,argStruct *arguments,int inputtype){
 
  doHaploCall=0;
  maxMis=-1;
  minMinor=0;
  doCount=0;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doHaploCall")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);


  if(doHaploCall==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);  
  //for kputs output
  bufstr.s=NULL;bufstr.l=bufstr.m=0;

  //make output files
  const char* postfix;
  postfix=".haplo.gz";

  outfileZ = aio::openFileBG(outfiles,postfix);

 

  //write header
  bufstr.l=0;
 
  ksprintf(&bufstr,"chr\tpos\tmajor");

  for(int i=0;i<arguments->nInd;i++)
      ksprintf(&bufstr,"\tind%d",i);
  
  
  ksprintf(&bufstr,"\n");
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;

}


abcHaploCall::~abcHaploCall(){

  if(doHaploCall==0)
    return; 
  fprintf(stderr,"\t-> The misc folder now contains a program to convert to plink format \'misc/haploToPlink angsdput.haplo.gz outputname\'\n");
  if(outfileZ!=NULL) bgzf_close(outfileZ);
}


void abcHaploCall::clean(funkyPars *pars){
 
  if(doHaploCall==0)
    return; 

  haploCalls *haplo =(haploCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++)
    delete[] haplo->dat[s];
  delete[] haplo->dat;
  delete[] haplo->major;

  delete haplo;
}




void abcHaploCall::printHaplo(funkyPars *pars){
  
  haploCalls *haplo =(haploCalls *) pars->extras[index];
  
  bufstr.l=0;
 
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;

    ksprintf(&bufstr,"%s\t%d\t%c\t",header->target_name[pars->refId],pars->posi[s]+1,intToRef[haplo->major[s]]);
 
    for(int i=0;i<pars->nInd;i++){
      ksprintf(&bufstr,"%c\t",intToRef[haplo->dat[s][i]]);
    }

    ksprintf(&bufstr,"\n");
  }

  if(bufstr.l>0)
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;

}

void abcHaploCall::print(funkyPars *pars){

  if(doHaploCall==0)
    return;
  
  printHaplo(pars);
}

void abcHaploCall::getHaplo(funkyPars *pars){
  
  haploCalls *haplo =(haploCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;      

    int siteCounts[5] = { 0 , 0 , 0 , 0 , 0 };

    //get random or most frequent base per individual
    for(int i=0;i<pars->nInd;i++){
      
      int dep=0;
      for( int b = 0; b < 4; b++ ){
	dep+=pars->counts[s][i*4+b];
      }

      if(dep==0){
	haplo->dat[s][i]=4;
	continue;
      }

	
      if(doHaploCall==1)//random base
	haplo->dat[s][i] = angsd::getRandomCount(pars->counts[s],i,dep);
      else if(doHaploCall==2)//most frequent base, random tie
	haplo->dat[s][i] = angsd::getMaxCount(pars->counts[s],i,dep);
      
      siteCounts[haplo->dat[s][i]]++;
    }

      
    // call major
    haplo->major[s]=4;
    int whichMax=0;
    int NnonMis=siteCounts[0];//number of A C G T for a site (non missing)
    for(int b=1;b<4;b++){
      if(siteCounts[b]>siteCounts[whichMax])
	whichMax=b;
      NnonMis += siteCounts[b];
    }
    haplo->major[s]=whichMax;
    
    //remove non polymorphic
    if( minMinor > 0 && minMinor > NnonMis - siteCounts[whichMax] )
      pars->keepSites[s] = 0;
    
    //remove sites with more than maxMis
    if(maxMis>=0 && maxMis < pars->nInd - NnonMis)
      pars->keepSites[s] = 0;
     
  }
  
}

void abcHaploCall::run(funkyPars *pars){

  if(doHaploCall==0)
    return;

  //allocate haplotype struct
  haploCalls *haplo = new haploCalls;
  haplo->dat = new int*[pars->numSites];   
  for(int s=0;s<pars->numSites;s++)
    haplo->dat[s] = new int[pars->nInd];
  haplo->major = new int[pars->numSites];

  //haploCall struct to pars
  pars->extras[index] = haplo;

  //get haplotypes
  getHaplo(pars);

}
