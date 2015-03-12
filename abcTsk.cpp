#include "shared.h"
#include <cmath>
#include "analysisFunction.h"
#include "abcTsk.h"
#include <ctype.h>

#define rlen 100
#define qslen 55

size_t mat[rlen][rlen][qslen][2][4][4];

void setZero(){
  for(int p1=0;p1<rlen;p1++)
    for(int p2=0;p2<rlen;p2++)
      for(int q=0;q<qslen;q++)
	for(int s=0;s<2;s++)
	  for(int rb=0;rb<4;rb++)
	    for(int ob=0;ob<4;ob++)
	      mat[p1][p2][q][s][rb][ob] =0;

}


void printMat(FILE *fp){
  for(int p1=0;p1<rlen;p1++)
    for(int p2=0;p2<rlen;p2++)
      for(int q=0;q<qslen;q++)
	for(int s=0;s<2;s++){
	  fprintf(fp,"%d\t%d\t%d\t%d\t",p1,p2,q,s);
	  for(int rb=0;rb<4;rb++)
	    for(int ob=0;ob<4;ob++)
	      if(s==1)
		fprintf(fp,"%lu\t",mat[p1][p2][q][s][rb][ob]);
	      else
		fprintf(fp,"%lu\t",mat[p2][p1][q][s][rb][ob]);
	  fprintf(fp,"\n");
	}

}




void abcTsk::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__); 
  fprintf(argFile,"-doThorfinn\t%d (This will simply dump a mismatch matrix)\n",doThorfinn);
  fprintf(argFile,"-ref\t%s\n",refName);

}

void abcTsk::getOptions(argStruct *arguments){
 
  doThorfinn=angsd::getArg("-doThorfinn",doThorfinn,arguments);
  refName = angsd::getArg("-ref",refName,arguments);

  if(doThorfinn==0)
    return;
  
  if(doThorfinn&&(refName==NULL)){
    fprintf(stderr,"Must supply -ref\n");
    printArg(stderr);
    exit(0);
  }



}

abcTsk::abcTsk(const char *outfiles,argStruct *arguments,int inputtype){
 doThorfinn=0;
 refName = NULL;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doThorfinn")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }



  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doThorfinn==0){
    shouldRun[index]=0;
    return;
  }
  fprintf(stderr,"Estimating calibration matrix: %s\n",outfiles);
  //make output files
  const char* postfix;
  postfix=".thorfinn";
  outfile=NULL;
  outfile = aio::openFile(outfiles,postfix);
  fprintf(outfile,"posi\tisop\tqs\tstrand");
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      fprintf(outfile,"\t%c%c",intToRef[i],intToRef[j]);
  fprintf(outfile,"\n");
  setZero();
  
}


abcTsk::~abcTsk(){
  if(doThorfinn){
    printMat(outfile);
    if(outfile) fclose(outfile);
  }
  free(refName);
}


void abcTsk::clean(funkyPars *pars){
  if(doThorfinn==0)
    return;

}

void abcTsk::print(funkyPars *pars){
  if(doThorfinn==0)
    return;
  chunkyT *chk = pars->chk;
  
  //loop over samples
  for(int i=0;i<pars->nInd;i++){
    //loop over sites;
    for(int s=0;s<pars->numSites;s++){
      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;
      for(int l=0;l<nd->l;l++){
	int refB = refToInt[pars->ref[s]];
	int obB = refToInt[nd->seq[l]];
	int strand = isupper(nd->seq[l])==0;
	if(refB==4||obB==4)
	  continue;
	mat[nd->posi[l]][nd->isop[l]][nd->qs[l]][strand][refB][obB]++;

	
      }
      

    }

  }
  
  
}


void abcTsk::run(funkyPars *pars){


}


