/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
 
    
  anders albrecht@binf.ku.dk made this.

  part of angsd
  ans -> anc dec 7 2013, added -ref -anc
*/

#include <cmath>
#include <cstdlib>

#include "analysisFunction.h"
#include "abcAncError.h"

void abcAncError::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAncError\t%d\n",doAncError);
  fprintf(argFile,"\t(Sampling strategies)\n");
  fprintf(argFile,"\t 0:\t no error estimation \n");
  fprintf(argFile,"\t 1:\t (Use all bases)\n");
  fprintf(argFile,"\t 2:\t (Sample single base)\n");
  fprintf(argFile,"\t 3:\t (Sample first base)\n");
  fprintf(argFile,"Required:\n\t-ref\t%s\t(fastafile containg \'perfect\' sample)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(fastafile containg outgroup)\n",ancName);
  fprintf(argFile,"\nNB: the -ref should be a fasta for a sample where you assume no errors.\nWe measure the difference between the outgroup and your -ref sample.\nThe statistic is then the excess of substitutions between your BAM file and outgroup, compared to the perfect sample. After the ANGSD run use:  Rscript R/estError.R file=angsdput.ancerror\n");
}

void abcAncError::getOptions(argStruct *arguments){

  //from command line
  doAncError=angsd::getArg("-doAncError",doAncError,arguments);

  if(doAncError==0)
    return;

  nInd=arguments->nInd;
  refName=angsd::getArg("-ref",refName,arguments);
  ancName=angsd::getArg("-anc",ancName,arguments);
 
  if(doAncError){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doAncError \n");
      exit(0);
    }
    if(refName==NULL||ancName==NULL){
      fprintf(stderr,"\t-> Must supply -ref and -anc for analysis\n");
      exit(0);
    }
  }

}

abcAncError::abcAncError(const char *outfiles,argStruct *arguments,int inputtype){
  tsk_outname=NULL;
  doAncError=0;
  currentChr=-1;
  refName=NULL;
  ancName=NULL;
  outfile = NULL;
  outfile2= NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAncError")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);


  if(doAncError==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);

  //make output files
  const char* postfix;
  postfix=".ancError";
  outfile = aio::openFile(outfiles,postfix);
  //put the name of the outputfile into tsk_outname
  tsk_outname=(char*)malloc(strlen(outfiles)+strlen(postfix)+1);
  sprintf(tsk_outname,"%s%s",outfiles,postfix);
  
  const char* postfix2;
  postfix2=".ancErrorChr";
  outfile2 = aio::openFile(outfiles,postfix2);

  //allocate allele counts
  alleleCounts = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCounts[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCounts[i][j]=0;

  alleleCountsChr = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCountsChr[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCountsChr[i][j]=0;
}


abcAncError::~abcAncError(){
  free(refName);free(ancName);

  if(doAncError==0)
    return; 
  fprintf(stderr,"\t-> To generate nice plots type in \'Rscript ANGSDDIR/R/estError.R file=\"%s\"\'\n",tsk_outname);
  fprintf(stderr,"\t-> Remember to surround the filename with \"\"\n");
  
  for(int i=0;i<nInd;i++){
    for(int j=0;j<125;j++)
      fprintf(outfile,"%lu\t",alleleCounts[i][j]);
    fprintf(outfile,"\n");
  }
  if(currentChr!=-1){
    fprintf(outfile2,"Chr: \t %s\n",header->target_name[currentChr]);
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	fprintf(outfile2,"%lu\t",alleleCountsChr[i][j]);
      fprintf(outfile2,"\n");
    }
  }
  

  for(int i=0;i<nInd;i++)
    delete[]  alleleCounts[i];
  delete [] alleleCounts; 

  for(int i=0;i<nInd;i++)
    delete[]  alleleCountsChr[i];
  delete [] alleleCountsChr; 

  if(outfile) fclose(outfile);
  if(outfile2) fclose(outfile2);
  free(tsk_outname);
}


void abcAncError::clean(funkyPars *pars){

}

void abcAncError::model1(funkyPars *pars){
    if(doAncError==1){//use all bases
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  for(int j=0;pars->chk->nd[s][i]&&j<pars->chk->nd[s][i]->l;j++){
	      alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[j]]]++;
	      alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[j]]]++;
	  }
	}
      }
    }
    if(doAncError==2){//random base
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  if(pars->chk->nd[s][i]==NULL||pars->chk->nd[s][i]->l==0)
	    continue;
	  int j = std::rand() % pars->chk->nd[s][i]->l;
	  alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[j]]]++;
	  alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[j]]]++;
	}
      }
    }
    if(doAncError==3){//first base
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	for(int i=0;i<pars->nInd;i++){
	  if(pars->chk->nd[s][i]&&pars->chk->nd[s][i]->l!=0){
	    alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[0]]]++;
	    alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i]->seq[0]]]++;
	  }
	}
      }
    }
}



void abcAncError::print(funkyPars *pars){

  if(doAncError==0)
    return;

  ///  if(doAncError==2){
  // model2(pars);
  //} else if(doAncError==1) {
    


  model1(pars);
  if(currentChr==-1)
    currentChr=pars->refId;
  if(currentChr!=pars->refId){
    fprintf(outfile2,"Chr: \t %s\n",header->target_name[currentChr]);
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	fprintf(outfile2,"%lu\t",alleleCountsChr[i][j]);
      fprintf(outfile2,"\n");
    }
    for(int i=0;i<nInd;i++)
      for(int j=0;j<256;j++)
	alleleCountsChr[i][j]=0;
    currentChr=pars->refId;
  }
    //}
   
}


void abcAncError::run(funkyPars *pars){

  if(doAncError==0)
    return;

}
