/*
  Thorfinn 31oct 2014

  refactored version of old code.
  indexing and reading of binary representation is now in prep_sites.cpp

  Code can be optimized by skipping rest of filereading for remainer of chr
osome if we have reached the last position from the keep list.

Maybe there are some memleaks. This will have to be fixed later.
*/

#include <cassert>
#include <sys/stat.h>
#include "shared.h"
#include "abc.h"
#include "analysisFunction.h"
#include "abcFilter.h"

abcFilter::~abcFilter(){
  if(fl!=NULL){
    dalloc(fl);
  }
  free(fname);
}


abcFilter::abcFilter(argStruct *arguments){
    //below if shared for all analysis classes
  shouldRun[index] = 1;
  header = arguments->hd;
  revMap = arguments->revMap;
  fl = NULL;
  //his is used by this class
  keepsChr = NULL;
  curChr = -1;
  fp = NULL;
  minInd = 0;
  fname = NULL;
  doMajorMinor =0;
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-sites")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  //get options and print them
  getOptions(arguments);
  printArg(arguments->argumentFile);

}


void abcFilter::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr regStart regStop))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos major minor))\n",fname);
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if atleast minInd of samples has data\n",minInd);
  fprintf(argFile,"\t1) You can force major/minor by -doMajorMinor 3\n\tAnd make sure file contains 4 columns (chr tab pos tab major tab minor)\n");
}

void abcFilter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-sites",fname,arguments);
  if(fname!=NULL)  
    fl = filt_read(fname);
  if(fl!=NULL)
    fprintf(stderr,"\t-> [%s] -sites is still beta, use at own risk...\n",__FILE__);


  //1=bim 2=keep
  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && fl!=NULL&& fl->hasMajMin!=1){
    fprintf(stderr,"\t-> Must supply -sites with a file containing major and minor if -doMajorMinor 3\n");
  }
  if(doMajorMinor!=3 && fl!=NULL&& fl->hasMajMin==1){
    fprintf(stderr,"\t-> Filter file contains major/minor information to use these in analysis supper \'-doMajorMinor 3\'\n");
  }
  
  minInd = angsd::getArg("-minInd",minInd,arguments);
  if(minInd >arguments->nInd){
    fprintf(stderr,"\t-> Potential problem you  filter -minInd %d but you only have %d samples?\n",minInd,arguments->nInd);
    exit(0);
  }
}


void abcFilter::run(funkyPars *p){
  //fprintf(stderr,"nsites=%d\n",p->numSites);
  p->keepSites=new int[p->numSites];
  
  for(int s=0;s<p->numSites;s++){
    p->keepSites[s]=p->nInd;
    //    p->results->freq->keepInd[s]=nInd;  
  }


  if(fl!=NULL && fl->hasMajMin==1 && doMajorMinor==3){
    p->major = new char [p->numSites];
    p->minor = new char [p->numSites];
    for(int i=0;i<p->numSites;i++){
      p->major[i] = 4;
      p->minor[i] = 4;
    }
  }
  
  if(fl!=NULL) {
    if(fl->keeps==NULL){
      for(int s=0;s<p->numSites;s++)
	p->keepSites[s] =0;
    }else{
      for(int s=0;s<p->numSites;s++){
	if(fl->keeps[p->posi[s]]==0){
	  p->keepSites[s] =0;
	}
	if(p->keepSites[s] && fl->hasMajMin==1 &&doMajorMinor==3){
	  p->major[s] = fl->major[p->posi[s]];
	  p->minor[s] = fl->minor[p->posi[s]];
	}
      }
    }
  }
  //how set the keepsites according the effective sample size persite
  //if(0!=minInd){
    if(p->chk!=NULL){
      //loop over sites;
      for(int s=0;s<p->numSites;s++){
	if(p->keepSites[s]==0)
	  continue;
	int nInfo =0;
	tNode **tn = p->chk->nd[s];
	//loop over samples;
	for(int i=0;i<p->nInd;i++){
	  if(tn[i]&&tn[i]->l!=0)
	    nInfo++;
	}
	p->keepSites[s] =nInfo;
	if(0!=minInd){
	  if(minInd>nInfo)
	    p->keepSites[s] =0;
	}
    
    }
  }
}
void abcFilter::print(funkyPars *p){
}

void abcFilter::clean(funkyPars *p){
  
}


void abcFilter::readSites(int refId) {
  //fprintf(stderr,"[%s].%s():%d -> refId:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);
  if(fl==NULL)
    return;
  filt_readSites(fl,header->target_name[refId],header->target_len[refId]);
}
