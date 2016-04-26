/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
 
    
  anders albrecht@binf.ku.dk made this.

  part of angsd
  ans -> anc dec 7 2013, added -ref -anc
*/

#include <cmath>
#include <cstdlib>
#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "abcIBS.h"

void abcIBS::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doIBS\t%d\n",doIBS);
  fprintf(argFile,"\t(Sampling strategies)\n");
  fprintf(argFile,"\t 0:\t no IBS \n"); 
  fprintf(argFile,"\t 1:\t (Sample single base)\n");
  fprintf(argFile,"\t 2:\t (Concensus base)\n");
  fprintf(argFile,"\t-doCounts\t%d\tMust choose -doCount 1\n",doCount);
  fprintf(argFile,"Optional\n");
  fprintf(argFile,"\t-minMinor\t%d\tMinimum observed minor alleles\n",minMinor);
  fprintf(argFile,"\t-output01\t%d\toutput 0 and 1s instead of based\n",output01);
  fprintf(argFile,"\t-maxMis\t%d\tMaximum missing bases (per site)\n",maxMis);
  fprintf(argFile,"\t-doMajorMinor\t%d\tuse input files or data to select major and minor alleles\n",majorminor);
 fprintf(argFile,"\t-makeMatrix\t%d\tprint out the ibs matrix \n",makeMatrix);

 }




void abcIBS::getOptions(argStruct *arguments){

  //from command line
  doIBS=angsd::getArg("-doIBS",doIBS,arguments);
  doCount=angsd::getArg("-doCounts",doCount,arguments);
  minMinor=angsd::getArg("-minMinor",minMinor,arguments);
  maxMis=angsd::getArg("-maxMis",maxMis,arguments);
  majorminor=angsd::getArg("-doMajorMinor",majorminor,arguments);
  output01=angsd::getArg("-output01",majorminor,arguments);
  makeMatrix=angsd::getArg("-makeMatrix",makeMatrix,arguments);
  

  if(doIBS==0)
    return;

  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"Error: bam or soap input needed for -doIBS \n");
    exit(0);
  }
  if( doCount==0){
    fprintf(stderr,"Error: -doIBSs needs allele counts (use -doCounts 1)\n");
    exit(0);
  }

}

abcIBS::abcIBS(const char *outfiles,argStruct *arguments,int inputtype){
 
  doIBS=0;
  maxMis=-1;
  minMinor=0;
  doCount=0;
  majorminor=0;
  output01=0;
  outfileMat=NULL;
  makeMatrix=1;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doIBS")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doIBS==0){
    shouldRun[index] =0;
    return;
  }
  
  //for kputs output
  bufstr.s=NULL;bufstr.l=bufstr.m=0;

  //make output files
  const char* postfix;
  postfix=".ibs.gz";

  outfileZ = aio::openFileBG(outfiles,postfix);

  const char* postfix2;
  postfix2=".ibsMat";

  if(makeMatrix)
    outfileMat =  aio::openFile(outfiles,postfix2);
 
  nInd = arguments->nInd;
  
  //write header
  bufstr.l=0;
 


  if(makeMatrix){
    ibsMatrixAll = new int[arguments->nInd*arguments->nInd];
    nonMisAll = new int[arguments->nInd*arguments->nInd];
    for(int i=0;i < arguments->nInd*arguments->nInd;i++){
      ibsMatrixAll[i] = 0;
      nonMisAll[i] = 0;
    }
  }  
 

  /// 
  if(majorminor)
    ksprintf(&bufstr,"chr\tpos\tmajor\tminor");
  else
    ksprintf(&bufstr,"chr\tpos\tmajor");

  for(int i=0;i<arguments->nInd;i++)
      ksprintf(&bufstr,"\tind%d",i);
  
  
  ksprintf(&bufstr,"\n");
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;

}


abcIBS::~abcIBS(){

  if(doIBS==0)
    return; 

  if(makeMatrix){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<nInd;j++){
	fprintf(outfileMat,"%f\t",ibsMatrixAll[i*nInd + j]*1.0/nonMisAll[i*nInd + j]);
      }
      fprintf(outfileMat,"\n");
    }
    delete [] ibsMatrixAll;
    delete [] nonMisAll;
  }


  if(outfileZ!=NULL) bgzf_close(outfileZ);
  if(outfileMat!=NULL) fclose(outfileMat);
}


void abcIBS::clean(funkyPars *pars){
 
  if(doIBS==0)
    return; 

  IBSstruct *haplo =(IBSstruct *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++)
    delete[] haplo->dat[s];
  delete[] haplo->dat;
  delete[] haplo->major;

  if(makeMatrix){
    delete[] haplo->ibsMatrix;
    delete[] haplo->nonMis;
  }

  delete haplo;

}




void abcIBS::printHaplo(funkyPars *pars){
  
  IBSstruct *haplo =(IBSstruct *) pars->extras[index];
  
  bufstr.l=0;

  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
 
    
    if(majorminor){
      for(int i=0;i<5;i++)
	intToMajorMinor[i] = -1 ;
      intToMajorMinor[pars->major[s]]=0;
      intToMajorMinor[pars->minor[s]]=1;
      ksprintf(&bufstr,"%s\t%d\t%c\t%c\t",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]]);
    }
    else
      ksprintf(&bufstr,"%s\t%d\t%c\t",header->target_name[pars->refId],pars->posi[s]+1,intToRef[haplo->major[s]]);
 
    if(output01){
     for(int i=0;i<pars->nInd;i++){
      	ksprintf(&bufstr,"%d\t",intToMajorMinor[haplo->dat[s][i]]);
	// ksprintf(&bufstr,"%c\t",haplo->dat[s][i]);
      }

    }
    else{
      for(int i=0;i<pars->nInd;i++){
	ksprintf(&bufstr,"%c\t",intToRef[haplo->dat[s][i]]);
      }
    }
    ksprintf(&bufstr,"\n");
  }

  if(bufstr.l>0)
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;


  if(makeMatrix){
    for(int i=0;i<pars->nInd;i++){
      for(int j=0;j<pars->nInd;j++){
	ibsMatrixAll[i*pars->nInd+j] += haplo->ibsMatrix[i*pars->nInd+j];
	nonMisAll[i*pars->nInd+j] += haplo->nonMis[i*pars->nInd+j];
      }
    }
  }

}

void abcIBS::print(funkyPars *pars){

  if(doIBS==0)
    return;
  
  printHaplo(pars);
}

void abcIBS::getHaplo(funkyPars *pars){
  
  IBSstruct *haplo =(IBSstruct *) pars->extras[index];

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
      
      
      if(doIBS==1)//random base
	haplo->dat[s][i] = angsd::getRandomCount(pars->counts[s],i,dep);
      else if(doIBS==2)//most frequent base, random tie
	haplo->dat[s][i] = angsd::getMaxCount(pars->counts[s],i,dep);
      
      siteCounts[haplo->dat[s][i]]++;
    }
    
    //  fprintf(stdout,"%d sfsdfsdf\n",majorminor);      
    // call major
    if(majorminor!=0)
      haplo->major[s]=pars->major[s];



    int whichMax=0;
    int NnonMis=siteCounts[0];//number of A C G T for a site (non missing)
    for(int b=1;b<4;b++){
      if(siteCounts[b]>siteCounts[whichMax])
	whichMax=b;
      NnonMis += siteCounts[b];
    }
    if(majorminor==0)
      haplo->major[s]=whichMax;
    
    //remove non polymorphic
    if( minMinor > 0 && minMinor > NnonMis - siteCounts[whichMax] )
      pars->keepSites[s] = 0;
    
    //set to missing haplotypes that are not major or minor
    if(majorminor!=0){
      if(pars->keepSites[s]==0)
	continue;      
      //get random or most frequent base per individual
      for(int i=0;i<pars->nInd;i++){
	if(haplo->dat[s][i]!=pars->major[s] && haplo->dat[s][i]!=pars->minor[s] )
	  haplo->dat[s][i] = 4;
      }
    }
    
    

    //remove sites with more than maxMis
    if(maxMis>=0 && maxMis < pars->nInd - NnonMis)
      pars->keepSites[s] = 0;
     
  }
  
}



void abcIBS::makeIBSmatrix(funkyPars *pars){
  
  IBSstruct *haplo =(IBSstruct *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;      

    
    for(int i=0;i<pars->nInd;i++){
      for(int j=0;j<pars->nInd;j++){
	if( haplo->dat[s][i]==4 || haplo->dat[s][j]==4 )
	  continue;

	haplo->nonMis[i*pars->nInd+j]++;
	if( haplo->dat[s][i] != haplo->dat[s][j] )
	  haplo->ibsMatrix[i*pars->nInd+j]++;
      }
    }

  }
 

}


void abcIBS::run(funkyPars *pars){

  if(doIBS==0)
    return;

  //allocate haplotype struct
  IBSstruct *haplo = new IBSstruct;
  haplo->dat = new int*[pars->numSites];   
  for(int s=0;s<pars->numSites;s++)
    haplo->dat[s] = new int[pars->nInd];
  haplo->major = new int[pars->numSites];

  if(makeMatrix){
    haplo->ibsMatrix = new int[pars->nInd*pars->nInd];
    haplo->nonMis = new int[pars->nInd*pars->nInd];
  }

  //haploCall struct to pars
  pars->extras[index] = haplo;

  //get haplotypes
  getHaplo(pars);


  if(makeMatrix){

    for(int i=0;i<pars->nInd*pars->nInd;i++){
      haplo->ibsMatrix[i] = 0;
      haplo->nonMis[i] = 0;
    }
    makeIBSmatrix(pars);

  }

}
