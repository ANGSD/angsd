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
#include "aio.h"

void abcAncError::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAncError\t%d\n",doAncError);
  fprintf(argFile,"\t(Sampling strategies)\n");
  fprintf(argFile,"\t 0:\t no error estimation \n");
  fprintf(argFile,"\t 1:\t (Use all bases)\n");
  fprintf(argFile,"\t 2:\t (Sample single base)\n");
  fprintf(argFile,"\t 3:\t (Sample first base)\n");
  fprintf(argFile,"Required:\n\t-ref\t%s\t(fastafile containing \'perfect\' sample)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(fastafile containing outgroup)\n",ancName);
  fprintf(argFile,"\nNB: the -ref should be a fasta for a sample where you assume no errors.\nWe measure the difference between the outgroup and your -ref sample.\nThe statistic is then the excess of substitutions between your BAM file and outgroup, compared to the perfect sample.\n");
  fprintf(argFile,"Output: .ancError (counts for the whole run) and .ancErrorChr (counts per chromosome)\n");
  fprintf(argFile,"After the ANGSD run use:  Rscript R/estError.R file=angsdput.ancError\n");
}

void abcAncError::getOptions(argStruct *arguments){

  //from command line
  doAncError=angsd::getArg("-doAncError",doAncError,arguments);

  if(doAncError==0)
    return;

  if(doAncError<0||doAncError>3){
    fprintf(stderr,"\t-> -doAncError must be 1, 2 or 3, value supplied: %d\n",doAncError);
    exit(1);
  }

  nInd=arguments->nInd;
  refName=angsd::getArg("-ref",refName,arguments);
  ancName=angsd::getArg("-anc",ancName,arguments);

  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"\t-> -doAncError needs BAM/CRAM or pileup input\n");
    exit(1);
  }
  if(refName==NULL||ancName==NULL){
    fprintf(stderr,"\t-> Must supply -ref and -anc for -doAncError\n");
    exit(1);
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
      printArg(stderr);
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
  snprintf(tsk_outname,strlen(outfiles)+strlen(postfix)+1,"%s%s",outfiles,postfix);

  const char* postfix2;
  postfix2=".ancErrorChr";
  outfile2 = aio::openFile(outfiles,postfix2);

  //allocate allele counts, () sets them to zero
  alleleCounts = new size_t *[nInd];
  alleleCountsChr = new size_t *[nInd];
  for(int i=0;i<nInd;i++){
    alleleCounts[i] = new size_t [ANCERR_NCAT]();
    alleleCountsChr[i] = new size_t [ANCERR_NCAT]();
  }
}

//write the counts of the current chromosome to .ancErrorChr and reset them
void abcAncError::printChr(){
  fprintf(outfile2,"Chr: \t %s\n",header->target_name[currentChr]);
  for(int i=0;i<nInd;i++){
    for(int j=0;j<ANCERR_NCAT;j++){
      fprintf(outfile2,"%lu\t",alleleCountsChr[i][j]);
      alleleCountsChr[i][j]=0;
    }
    fprintf(outfile2,"\n");
  }
}

abcAncError::~abcAncError(){
  free(refName);free(ancName);

  if(doAncError==0)
    return;
  fprintf(stderr,"\t-> To generate nice plots type in \'Rscript ANGSDDIR/R/estError.R file=\"%s\"\'\n",tsk_outname);
  fprintf(stderr,"\t-> Remember to surround the filename with \"\"\n");

  for(int i=0;i<nInd;i++){
    for(int j=0;j<ANCERR_NCAT;j++)
      fprintf(outfile,"%lu\t",alleleCounts[i][j]);
    fprintf(outfile,"\n");
  }
  if(currentChr!=-1)
    printChr();

  for(int i=0;i<nInd;i++){
    delete[]  alleleCounts[i];
    delete[]  alleleCountsChr[i];
  }
  delete [] alleleCounts;
  delete [] alleleCountsChr;

  if(outfile) fclose(outfile);
  if(outfile2) fclose(outfile2);
  free(tsk_outname);
}


void abcAncError::clean(funkyPars *pars){

}

//add one observed base for individual ind at site s
void abcAncError::count(int ind,funkyPars *pars,int s,char base){
  int cat = pars->anc[s]*25+pars->ref[s]*5+refToInt[(int)base];
  alleleCounts[ind][cat]++;
  alleleCountsChr[ind][cat]++;
}

void abcAncError::model1(funkyPars *pars){
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    for(int i=0;i<pars->nInd;i++){
      tNode *nd = pars->chk->nd[s][i];
      if(nd==NULL||nd->l==0)
	continue;
      if(doAncError==1){//use all bases
	for(int j=0;j<nd->l;j++)
	  count(i,pars,s,nd->seq[j]);
      }else if(doAncError==2)//random base
	count(i,pars,s,nd->seq[std::rand() % nd->l]);
      else if(doAncError==3)//first base
	count(i,pars,s,nd->seq[0]);
    }
  }
}



void abcAncError::print(funkyPars *pars){

  if(doAncError==0)
    return;

  //a new chromosome, write the counts of the previous one before counting this chunk
  if(currentChr!=-1&&currentChr!=pars->refId)
    printChr();
  currentChr=pars->refId;

  model1(pars);
}


void abcAncError::run(funkyPars *pars){

  if(doAncError==0)
    return;

}
