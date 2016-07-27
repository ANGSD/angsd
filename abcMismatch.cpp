#include <ctype.h>
#include <htslib/kstring.h>
#include "shared.h"

#include "analysisFunction.h"
#include "abcMismatch.h"


#define rlen 400
#define qslen 64

size_t ******mat=NULL;//[rlen][rlen][qslen][2][4][4];

void setZero(){
  for(int p1=0;p1<rlen;p1++)
    for(int p2=0;p2<rlen;p2++)
      for(int q=0;q<qslen;q++)
	for(int s=0;s<2;s++)
	  for(int rb=0;rb<4;rb++)
	    for(int ob=0;ob<4;ob++)
	      mat[p1][p2][q][s][rb][ob] =0;

}

void setZero2(){
  mat = new size_t*****[rlen];
  for(int p1=0;p1<rlen;p1++){
    mat[p1] = new size_t****[rlen];
    for(int p2=0;p2<rlen;p2++){
      mat[p1][p2] = new size_t***[qslen];
      for(int q=0;q<qslen;q++){
	mat[p1][p2][q] = new size_t**[2];
	for(int s=0;s<2;s++){
	  mat[p1][p2][q][s] = new size_t*[4];
	  for(int rb=0;rb<4;rb++){
	    mat[p1][p2][q][s][rb] = new size_t[4];
	    for(int ob=0;ob<4;ob++)
	      mat[p1][p2][q][s][rb][ob] =0;
	  }
	}
      }
    }
  }
}


void printMat(BGZF* fp,kstring_t &bufstr){
  bufstr.l=0;
  for(int p1=0;p1<rlen;p1++)
    for(int p2=0;p2<rlen;p2++){
      for(int q=0;q<qslen;q++){
	for(int s=0;s<2;s++){
	  for(int rb=0;rb<4;rb++){
	    size_t ts = 0;
	    for(int ob=0;ob<4;ob++){
	      if(s==1)
		ts += mat[p1][p2][q][s][rb][ob];
	      else
		ts += mat[p2][p1][q][s][rb][ob];
	    }
	    if(ts>0){
	      ksprintf(&bufstr,"%d\t%d\t%d\t%d\t%d",p1,p2,q,s,rb);
	      for(int ob=0;ob<4;ob++){
		if(s==1)
		  ksprintf(&bufstr,"\t%zu",mat[p1][p2][q][s][rb][ob]);
		else
		  ksprintf(&bufstr,"\t%zu",mat[p2][p1][q][s][rb][ob]);
		
	      }
	      ksprintf(&bufstr,"\n");
	    }
	  }
	}
      }
    }

  aio::bgzf_write(fp,bufstr.s,bufstr.l);bufstr.l=0;
}




void abcTsk::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__); 
  fprintf(argFile,"-doMisMatch\t%d (This will simply dump a mismatch matrix)\n",doMismatch);
  fprintf(argFile,"-ref\t%s\n",refName);

}

void abcTsk::getOptions(argStruct *arguments){
 
  doMismatch=angsd::getArg("-doMismatch",doMismatch,arguments);
  refName = angsd::getArg("-ref",refName,arguments);

  if(doMismatch==0)
    return;
  
  if(doMismatch&&(refName==NULL)){
    fprintf(stderr,"\t-> Must supply -ref\n");
    printArg(stderr);
    exit(0);
  }



}

abcTsk::abcTsk(const char *outfiles,argStruct *arguments,int inputtype){
 doMismatch=0;
 bufstr.s=NULL;bufstr.l=bufstr.m=0;
 refName = NULL;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doMismatch")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }



  getOptions(arguments);

  if(doMismatch==0){
    shouldRun[index]=0;
    return;
  }
  printArg(arguments->argumentFile);
  fprintf(stderr,"\t-> Estimating calibration matrix: %s\n",outfiles);
  //make output files
  const char* postfix;
  postfix=".mismatch.gz";
  outfilegz=NULL;
  outfilegz = aio::openFileBG(outfiles,postfix);
  ksprintf(&bufstr,"posi\tisop\tqs\tstrand\tRef");
  for(int j=0;j<4;j++)
    ksprintf(&bufstr,"\t%c",intToRef[j]);
  ksprintf(&bufstr,"\n");
  aio::bgzf_write(outfilegz,bufstr.s,bufstr.l);bufstr.l=0;
  setZero2();
  
}


abcTsk::~abcTsk(){
  if(doMismatch){
    printMat(outfilegz,bufstr);
    if(outfilegz!=NULL)
      bgzf_close(outfilegz);
 
    for(int p1=0;p1<rlen;p1++){
      for(int p2=0;p2<rlen;p2++){
	for(int q=0;q<qslen;q++){
	  for(int s=0;s<2;s++){
	    for(int rb=0;rb<4;rb++)
	      delete [] mat[p1][p2][q][s][rb];
	    delete [] mat[p1][p2][q][s];
	  }
	  delete [] mat[p1][p2][q];
	}
	delete [] mat[p1][p2];
      }
      delete [] mat[p1];
    }
    
    delete [] mat;
  }

  free(refName);

  
 
}


void abcTsk::clean(funkyPars *pars){
  if(doMismatch==0)
    return;

}

void abcTsk::print(funkyPars *pars){
  if(doMismatch==0)
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


