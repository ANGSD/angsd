/*
  thorfinn@binf.ku.dk april 15 2012
  class to call genotypes
*/

#include <assert.h>
#include <htslib/kstring.h>
#include "shared.h"

#include "analysisFunction.h"
#include "abcCallGenotypes.h"


//filename of dumped file
#define GENO ".geno.gz"

//stuff related to genotypecalling
#define GENO_MAJOR_MINOR 0x1 // write the major and the minor
#define GENO_PRINT 2 //print the called genotype 0,1,2
#define GENO_ALLELE 4 //print the called genotype AA, AC AG ...
#define GENO_ALL_POST 8 //print the all posterior
#define GENO_WRITE_POST 16 //write the post of the called genotype
#define GENO_FOR_COVAR 32 //binary dump of the posteriors.

void abcCallGenotypes::printArg(FILE *argFile){
  fprintf(argFile,"-----------------\n%s:\n\n",__FILE__);
  fprintf(argFile,"-doGeno\t%d\n",doGeno);
  fprintf(argFile,"\t1: write major and minor\n");
  fprintf(argFile,"\t2: write the called genotype encoded as -1,0,1,2, -1=not called\n");
  fprintf(argFile,"\t4: write the called genotype directly: eg AA,AC etc \n");
  fprintf(argFile,"\t8: write the posterior probability of all possible genotypes\n");
  fprintf(argFile,"\t16: write the posterior probability of called genotype\n");
  fprintf(argFile,"\t32: write the posterior probabilities of the 3 gentypes as binary\n");
  //  fprintf(argFile,"\t64: write the three posterior probability (Beagle style)\n");
  fprintf(argFile,"\t-> A combination of the above can be choosen by summing the values, EG write 0,1,2 types with majorminor as -doGeno 3\n");
  fprintf(argFile,"\t-postCutoff=%f (Only genotype to missing if below this threshold)\n",postCutoff);
  fprintf(argFile,"\t-geno_minDepth=%d\t(-1 indicates no cutof)\n",geno_minDepth);
  fprintf(argFile,"\t-geno_maxDepth=%d\t(-1 indicates no cutof)\n",geno_maxDepth);
  fprintf(argFile,"\t-geno_minMM=%f\t(minimum fraction af major-minor bases)\n",geno_minMM);
  fprintf(argFile,"\t-minInd=%d\t(only keep sites if you call genotypes from this number of individuals)\n",minInd);
  fprintf(argFile,"\n\tNB When writing the posterior the -postCutoff is not used\n");
  fprintf(argFile,"\tNB geno_minDepth requires -doCounts\n");
  fprintf(argFile,"\tNB geno_maxDepth requires -doCounts\n");
 
}


void abcCallGenotypes::getOptions(argStruct *arguments){
  //default
  //from command line
  doGeno=angsd::getArg("-doGeno",doGeno,arguments);
  if(doGeno==0)
    return;
  if(doGeno==1)
    fprintf(stderr,"\n\t-> Your choice of -doGeno will not printout called genotypes only chr pos major minor\n");
    
  int doMaf,doPost;
  doMaf=doPost= 0;

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  postCutoff=angsd::getArg("-postCutoff",postCutoff,arguments);


  geno_minMM=angsd::getArg("-geno_minMM",geno_minMM,arguments);
  geno_minDepth=angsd::getArg("-geno_minDepth",geno_minDepth,arguments);
  geno_maxDepth=angsd::getArg("-geno_maxDepth",geno_maxDepth,arguments);
  
  int doCounts=0;
  doCounts = angsd::getArg("-doCounts",doCounts,arguments);
  if((geno_minDepth!=-1 || geno_maxDepth!=-1 || geno_minMM!=-1 ) &&doCounts==0){
    fprintf(stderr,"Must supply -doCounts to use a minimum depth for GC calling\n");
    exit(0);
  }
  if(geno_maxDepth!=-1&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts to use a maximum depth for GC calling\n");
    exit(0);
  }
  

  if(arguments->inputtype!=INPUT_BEAGLE&&doPost==0){
    fprintf(stderr,"\n\t-> You need -doPost to call genotypes \n");
    exit(0);

  }

  if(doPost==1&&doMaf==0){
    fprintf(stderr,"\n\t-> You need -doMaf inorder to get posterior probabilities when using freq as prior\n");
    exit(0);
  }
  minInd=angsd::getArg("-minInd",minInd,arguments);

}

abcCallGenotypes::abcCallGenotypes(const char *outfiles,argStruct *arguments,int inputtype){
  doGeno=0;
  postCutoff= 1.0/3.0;
  outfileZ = NULL;
  geno_minMM = -1;
  geno_minDepth = -1;
  geno_maxDepth = -1;
  minInd = 0;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doGeno")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);

  if(doGeno==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);

  //make output files
  const char* postfix=GENO;
  if(doGeno>0)
    outfileZ = aio::openFileBG(outfiles,postfix);
  
}


abcCallGenotypes::~abcCallGenotypes(){
  if(outfileZ!=NULL)     
    bgzf_close(outfileZ);
  if(bufstr.s)
    free(bufstr.s);
}


void abcCallGenotypes::getGeno(funkyPars *pars){
  //pp is eiter the genotype likelihoods or the posterior probablity
  genoCalls *geno =(genoCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    pars->keepSites[s] = 0;
    for( int i =0;i<pars->nInd;i++){
      double maxPP=pars->post[s][i*3+0];
      int maxGeno=0;
      if(pars->post[s][i*3+1]>maxPP){
	maxPP=pars->post[s][i*3+1];
	maxGeno=1;
      }
      if(pars->post[s][i*3+2]>maxPP){
	maxPP=pars->post[s][i*3+2];
	maxGeno=2;
      }
      if(maxPP<=postCutoff)
	maxGeno=-1;
      if(geno_minDepth!=-1) {
	int geno_Depth = pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3];
	if(geno_Depth<geno_minDepth)
	  maxGeno=-1;
      }
      if(geno_maxDepth!=-1) {
	int geno_Depth = pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3];
	if(geno_Depth>geno_maxDepth)
	  maxGeno=-1;
      }
      if(geno_minMM!=-1){
	int geno_Depth = pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3];
	int mm_Depth = pars->counts[s][i*4+pars->major[s]] + pars->counts[s][i*4+pars->minor[s]];

	if(mm_Depth * 1.0 / geno_Depth < geno_minMM)
	  maxGeno=-1;



      }

      geno->dat[s][i]=maxGeno;
      if(geno->dat[s][i]!=-1)
	pars->keepSites[s]++;
    }
    if(pars->keepSites[s]<minInd)
      pars->keepSites[s] = 0;
  }
}



void abcCallGenotypes::printGeno(funkyPars *pars){
  genoCalls *geno =(genoCalls *) pars->extras[index];
  assert(pars->keepSites!=NULL);
  bufstr.l=0;
  int doGenoInner = abs(doGeno);
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    if(!(doGenoInner&GENO_FOR_COVAR))
      ksprintf(&bufstr,"%s\t%d\t",header->target_name[pars->refId],pars->posi[s]+1);
    if(doGenoInner&GENO_MAJOR_MINOR&&!(doGenoInner&GENO_FOR_COVAR))
      ksprintf(&bufstr,"%c\t%c\t",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
    if(doGenoInner&GENO_FOR_COVAR){
      aio::bgzf_write(outfileZ,pars->post[s],sizeof(double)*3*pars->nInd);
      continue;
    }

    for(int i=0;i<pars->nInd;i++){
      if(doGenoInner&GENO_PRINT)
	ksprintf(&bufstr,"%d\t",geno->dat[s][i]);
      if(doGenoInner&GENO_ALLELE){
	if(geno->dat[s][i]==0)
	  ksprintf(&bufstr,"%c%c\t",intToRef[pars->major[s]],intToRef[pars->major[s]]);
	else if(geno->dat[s][i]==1)
	  ksprintf(&bufstr,"%c%c\t",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==2)
	  ksprintf(&bufstr,"%c%c\t",intToRef[pars->minor[s]],intToRef[pars->minor[s]]);
	else if(geno->dat[s][i]==-1)
	  ksprintf(&bufstr,"NN\t");

      }
      if(doGenoInner&GENO_WRITE_POST){
	int which = angsd::whichMax(&pars->post[s][i*3],3);
	if(which!=-1)
	  ksprintf(&bufstr,"%f\t",pars->post[s][i*3+which]);
	else
	  ksprintf(&bufstr,"0.0\t");
	
      }


      if(doGenoInner&GENO_ALL_POST)
	for(int n=0;n<3;n++)
	  ksprintf(&bufstr,"%f\t",pars->post[s][i*3+n]);
    }
    ksprintf(&bufstr,"\n");

  }
  if(bufstr.l>0)
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
}




void abcCallGenotypes::clean(funkyPars *pars){
  if(doGeno==0)
    return;
  genoCalls *geno =(genoCalls *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++)
    delete[] geno->dat[s];
  delete[] geno->dat;

  delete geno;
}


void abcCallGenotypes::print(funkyPars *pars){
  if(doGeno<=0)
    return;
  
  genoCalls *geno =(genoCalls *) pars->extras[index];

  printGeno(pars);
  
}

void abcCallGenotypes::run(funkyPars *pars){
  
  if(doGeno==0)
    return;
  
  //allocate genoCall struct
  genoCalls *geno = new genoCalls;
  geno->dat = new int*[pars->numSites];   
  for(int s=0;s<pars->numSites;s++)
    geno->dat[s]=new int[pars->nInd];
  //genoCall struct to pars
  pars->extras[index] = geno;
  

  getGeno(pars);
}

