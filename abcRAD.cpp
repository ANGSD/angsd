/*
  DEPRECATED: -doRad is not maintained and may be removed. It was made from abcTemplate as
  a small example. Use -doCounts 1 -dumpCounts 4 for the base counts of each sample, or
  -doSnpStat 1 for the counts by strand and a strand bias test.

  -doRad 1:
  Counts the A, C, G, T and N of all samples for each strand at each site and writes them
  to the .results file. The strand is from the case of the bases (upper case is forward).
  N includes bases masked because of a low quality (-minQ).

  The counting is in ::print, so it is not threaded.

  Run command example:
  ---------------------
  ./angsd -bam bam.filelist -doRad 1 -out out
 */

#include <ctype.h> //<-used for isupper/islower
#include <cmath> //in order to construct large array/objects
#include <stdlib.h>
#include "aio.h"
#include "shared.h" //<-contains the struct defintino for funkypars
#include "analysisFunction.h" //<-contains some utility functions, usefull for parsing args
#include "abcRAD.h"//contains the analysis class defition associated with this file


//this msg is shown on screen if you type
// ./angsd -doRad
//otherwise this is printet to the file .arg
void abcRAD::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doRad\t%d (Deprecated, use -doCounts 1 -dumpCounts 4 or -doSnpStat 1)\n",doRad);
  fprintf(argFile,"\t\t1: Count and print basetypes in combination with strand (BAM/CRAM or pileup input)\n");
}

//this is the function that parses the parameters used for this analysis class
void abcRAD::getOptions(argStruct *arguments){
  //from command line
  doRad=angsd::getArg("-doRad",doRad,arguments);

  if(doRad==0){
    //nothing is run for this class (apart from the destructor)
    shouldRun[index]=0;
    return;
  }
  //-doRad 2 needed -ref but did nothing, other values made an empty file with exit code 0
  if(doRad!=1){
    fprintf(stderr,"\t-> -doRad must be 1, not %d\n",doRad);
    exit(1);
  }
  //the reads are only there for bam and pileup input, other input crashed
  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"\t-> -doRad needs BAM/CRAM or pileup input\n");
    exit(1);
  }
  fprintf(stderr,"\t-> -doRad is deprecated, use -doCounts 1 -dumpCounts 4 or -doSnpStat 1\n");
}

//constructor
abcRAD::abcRAD(const char *outfiles,argStruct *arguments,int inputtype){
  doRad = 0; //defaults= dont do analysis
  outfile = NULL;

  //first a hook for the interactive help:
  //  ./angsd -doRad
  if (arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doRad")){
      printArg(stderr);
      exit(0);
    } else
      return;
  }

  //now parse the arguments
  getOptions(arguments);
  //now print the arguments

  if(doRad==0)
    return ;

  printArg(arguments->argumentFile);
  //initalize outputfile
  outfile = aio::openFile(outfiles,".results");
  fprintf(outfile,"chromo\tposition\t+A\t+C\t+G\t+T\t+N\t-A\t-C\t-G\t-T\t-N\n");
}



//destructor
abcRAD::~abcRAD(){
  if(doRad==0)
    return;

  if(outfile!=NULL)
    fclose(outfile);
}



//this function is run, after ::run and ::print
void abcRAD::clean(funkyPars *pars){
  //we havent done any allocation so we dont need to cleanup
}



void abcRAD::print(funkyPars *pars){

  if(doRad==1){
    //count bases by strand

    //rawseqdata is in chunkyT struct (bambi_interface.h)
    chunkyT *chk = pars->chk;

    //loop over sites;
    for(int s=0;s<pars->numSites;s++){
      //sites removed by the filters were printed
      if(pars->keepSites[s]==0)
	continue;
      int bases[2][5] = {{0,0,0,0,0},{0,0,0,0,0}};

      //loop over samples
      for(int i=0;i<pars->nInd;i++){
	//all seqdata associated with single bamfile is in a tNode, NULL if the sample has no reads here (it crashed)
	tNode *nd = chk->nd[s][i];
	if(nd==NULL)
	  continue;
	//loop over the individual bases
	for(int l=0;l<nd->l;l++){
	  unsigned char c = nd->seq[l]; //this is the base
	  int strand = isupper(c)==0; //strand is defined as either small/big letters

	  //there is a lookuptable called refToInt which maps
	  //a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->4
	  bases[strand][refToInt[c]]++;
	}
      }
      //one line for each site with the counts of all samples, it was a line for each sample
      //with the counts of the samples so far
      fprintf(outfile,"%s\t%d",header->target_name[pars->refId],pars->posi[s]+1);//position is zero index internally
      for(int i=0;i<2;i++)
	for(int j=0;j<5;j++)
	  fprintf(outfile,"\t%d",bases[i][j]);
      fprintf(outfile,"\n");
    }
  }

}

void abcRAD::run(funkyPars *pars){
  if(doRad==0)
    return;
  //this functions runs for each block
}
