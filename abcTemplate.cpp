/*
  thorfinn 26june 2014 thorfinn@binf.ku.dk, part of angsd

  generally angsd will send the data for a genomic region through all analysis classes (abc*).
  the data is encapsulated in a 'funkyPars' struct (found in shared.h)
  
  For each abcClass it will do

  1) ::run
  2) ::print
  3) ::clean
  
  below are some toy example that will show
  a) how the class system works ()
  b) how to access the internal datastructure

  -doTemplate 1:
  Will simply count the number of A, C, G,T, N's for both strand

  //we will put all analysis in the ::print function. therefore it will NOT be threaded
  //the raw analysis should be put into ::run (this is threaded in angsd)
  //but we wont do this in case of simplicity
  
 */

#include <ctype.h> //<-used for isupper/islower
#include "shared.h" //<-contains the struct defintino for funkypars
#include "analysisFunction.h" //<-contains some utility functions, usefull for parsing args
#include "abcTemplate.h"//contains the analysis class defition associated with this file


//this msg is shown on screen if you type
// ./angsd -doTemplate
//otherwise this is printet to the file .arg
void abcTemplate::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doTemplate\t%d (Which analysis should we perform?)\n",doTemplate);
  fprintf(argFile,"\t\t1: Count and print basetypes in combination with strand\n");
}

//this is the function that parses the parameters used for this analysis class
void abcTemplate::getOptions(argStruct *arguments){
  //from command line
  doTemplate=angsd::getArg("-doTemplate",doTemplate,arguments);

  if(doTemplate==0){
    /*
      if this class shouldnt do any analysis,
      then setting this to zero will make sure nothing is run (apart from destructor)
    
      this could also have been accomplished by
      if(doTemplate==0) return
      in ::run ::clean ::print
    */
    
    shouldRun[index]=0;
    return;
  }
}

//constructor
abcTemplate::abcTemplate(const char *outfiles,argStruct *arguments,int inputtype){
  doTemplate = 0; //defaults= dont do analysis
  outfile = NULL;
  //first a hook for the interactive help:
  //  ./angsd -doTemplate
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doTemplate")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  //now parse the arguments
  getOptions(arguments);
  //now print the arguments
  printArg(arguments->argumentFile);
  if(doTemplate==0)
    return ;

  //initalize outputfile
  outfile = aio::openFile(outfiles,".results");
  fprintf(outfile,"Chromo\tPosition\t+A\t+C\t+G\t+T\t-A\t-C\t-G\t-T\n");
}


//destructor
abcTemplate::~abcTemplate(){
  if(outfile!=NULL) 
    fclose(outfile);
}

//this function is run, after ::run and ::print
void abcTemplate::clean(funkyPars *pars){
  //we havent done any allocation so we dont need to cleanup
}



void abcTemplate::print(funkyPars *pars){
  
  if(doTemplate==1){
    //count bases by strand

    //rawseqdata is in chunkyT struct (bambi_interface.h)
    chunkyT *chk = pars->chk;

    //loop over sites;
    for(int s=0;s<pars->numSites;s++){
      int bases[2][5] = {{0,0,0,0,0},{0,0,0,0,0}};      
      
      //loop over samples
      for(int i=0;i<pars->nInd;i++){
	//all seqdata associated with single bamfile is in a tNode
	tNode nd = chk->nd[s][i];
	//loop over the individual bases
	for(int l=0;l<nd.l;l++){
	  char c = nd.seq[l]; //this is the base
	  char q = nd.qs[l]; //this is the associated qscore, fancy shit
	  int strand = isupper(nd.seq[l])==0; //strand is defined as either small/big letters
	  
	  //there is a lookuptable called refToInt which maps
	  //a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->5
	  bases[strand][refToInt[c]]++;
	}
	//print chr and position
	fprintf(outfile,"%s\t%d",header->name[pars->refId],pars->posi[s]+1);//position is zero index internally
	
	//print the basecount
	for(int i=0;i<2;i++)
	  for(int j=0;j<5;j++)
	    fprintf(outfile,"\t%d",bases[i][j]);
	fprintf(outfile,"\n");
	

      }
      
    }
  }
  
}


void abcTemplate::run(funkyPars *pars){


}


