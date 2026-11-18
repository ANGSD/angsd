/*
  Thorfinn 31oct 2014

  refactored version of old code.
  indexing and reading of binary representation is now in prep_sites.cpp

  Code can be optimized by skipping rest of filereading for remainder of
  chromosome if we have reached the last position from the keep list.
*/

#include <sys/stat.h>
#include "shared.h"
#include "abc.h"
#include "analysisFunction.h"
#include "abcFilter.h"

abcFilter::~abcFilter(){
  if(fl!=NULL){
    filt_dalloc(fl);
  }
  free(fname);
}


abcFilter::abcFilter(argStruct *arguments){
  //below is shared for all analysis classes
  shouldRun[index] = 1;
  header = arguments->hd;
  revMap = arguments->revMap;
  //this is used by this class
  strict = 1;
  setMinIndDepth = 1;
  fl = NULL;
  minInd = 0;
  fname = NULL;
  doMajorMinor = 0;
  capDepth = -1;
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
  fprintf(argFile,"\t-sites\t\t%s\tFile containing sites to keep. Columns can be:\n",fname);
  fprintf(argFile,"\t\t\t\t (chr pos)\n");
  fprintf(argFile,"\t\t\t\t (chr regStart regStop)\n");
  fprintf(argFile,"\t\t\t\t (chr pos major minor)\n");
  fprintf(argFile,"\t\t\t\t (chr pos major minor af ac an)\n");
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if at least minInd of samples has data (requires BAM/CRAM or pileup input)\n",minInd);
  fprintf(argFile,"\t-setMinDepthInd\t%d\tOnly count a sample as having data if it has at least this depth\n",setMinIndDepth);
  fprintf(argFile,"\t-capDepth\t%d\tOnly use the first capDepth reads of each sample at each site, in pileup order (requires BAM/CRAM or pileup input)\n",capDepth);
  fprintf(argFile,"\t-strict\t\t%d\t(experimental) 0: -sites only removes chromosomes that are absent from the file\n",strict);
  fprintf(argFile,"\t1) You can force major/minor by -doMajorMinor 3\n\tAnd make sure file contains 4 columns (chr tab pos tab major tab minor)\n");
}

void abcFilter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-sites",fname,arguments);
  if(fname!=NULL)
    fl = filt_read(fname);

  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && fl!=NULL && fl->hasExtra==0){
    fprintf(stderr,"\t-> Must supply -sites with a file containing major and minor if -doMajorMinor 3\n");
    exit(1);
  }
  if(doMajorMinor!=3 && fl!=NULL && fl->hasExtra>0){
    fprintf(stderr,"\t-> Filter file contains major/minor information; to use these in the analysis supply \'-doMajorMinor 3\'\n");
  }
  capDepth = angsd::getArg("-capDepth",capDepth,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);
  setMinIndDepth = angsd::getArg("-setMinDepthInd",setMinIndDepth,arguments);
  strict = angsd::getArg("-strict",strict,arguments);
  if(minInd > arguments->nInd){
    fprintf(stderr,"\t-> Potential problem you filter -minInd %d but you only have %d samples?\n",minInd,arguments->nInd);
    exit(1);
  }
  if(capDepth!=-1 && capDepth<setMinIndDepth){
    fprintf(stderr,"\t-> -capDepth %d is smaller than -setMinDepthInd %d, this would discard all sites\n",capDepth,setMinIndDepth);
    exit(1);
  }
  int hasReads = arguments->inputtype==INPUT_BAM || arguments->inputtype==INPUT_PILEUP;
  if(!hasReads && (capDepth!=-1 || minInd!=0))
    fprintf(stderr,"\t-> -capDepth and -minInd require BAM/CRAM or pileup input, they are ignored for this input type\n");
}


void abcFilter::run(funkyPars *p){
  //if we are coming from BAM/CRAM then we dont need to set the below
  if(p->keepSites==NULL){
    p->keepSites=new int[p->numSites];
    for(int s=0;s<p->numSites;s++)
      p->keepSites[s]=p->nInd;
  }
  if(capDepth!=-1 && p->chk!=NULL){
    for(int s=0;s<p->numSites;s++)
      for(int i=0;i<p->nInd;i++)
	if(p->chk->nd[s][i]&&p->chk->nd[s][i]->l>capDepth)
	  p->chk->nd[s][i]->l=capDepth;
  }

  if(fl!=NULL) {
    if(fl->keeps==NULL){
      //chromosome is not in the -sites file
      for(int s=0;s<p->numSites;s++)
	p->keepSites[s] =0;
    }else if(strict){
      for(int s=0;s<p->numSites;s++)
	if(fl->keeps[p->posi[s]]==0)
	  p->keepSites[s] =0;
    }
  }

  //set keepSites to the effective sample size per site
  if(p->chk!=NULL){
    for(int s=0;s<p->numSites;s++){
      if(p->keepSites[s]==0)
	continue;
      int nInfo =0;
      tNode **tn = p->chk->nd[s];
      for(int i=0;i<p->nInd;i++)
	if(tn[i]&&tn[i]->l>=setMinIndDepth)
	  nInfo++;
      p->keepSites[s] =nInfo;
      if(minInd>nInfo)
	p->keepSites[s] =0;
    }
  }
}

void abcFilter::print(funkyPars *p){
}

void abcFilter::clean(funkyPars *p){
}


void abcFilter::readSites(int refId) {
  if(fl==NULL)
    return;
  filt_readSites(fl,header->target_name[refId],header->target_len[refId]);
}
