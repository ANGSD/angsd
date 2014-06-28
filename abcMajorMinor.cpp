/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd

  
  This class will populate the major/minor entries of the funkypars

  1) set majorminor according to the GL.
  2) set majorminor accordning to the counts of alleles
  3) wont populate the majorminor but will use the information from -sites
  4) set the major to the reference allele
  5) set the major to the ancestral allele


  models used are:
  1)
  Line Skotte, Thorfinn Sand Korneliussen, Anders Albrechtsen. Association testing for next-generation sequencing data using score statistics. Genet Epidemiol. 2012 Jul;36(5):430-7. 

  2)
  Li Y, Vinckenbosch N, Tian G, Huerta-Sanchez E, Jiang T, Jiang H, Albrechtsen A, Andersen G, Cao H, Korneliussen T, et al., 2010. Resequencing of 200 human exomes identifies an excess of low-frequency non-synonymous coding variants. Nat Genet 42:969–972. 


*/


#include <cmath> //<- for log,exp
#include <cassert>
#include <cfloat>
#include "shared.h"
#include "analysisFunction.h"
#include "abc.h"
#include "abcMajorMinor.h"

void abcMajorMinor::printArg(FILE *argFile){
  fprintf(argFile,"-------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doMajorMinor\t%d\n",doMajorMinor);
  fprintf(argFile,"\t1: Infer major and minor from GL\n");
  fprintf(argFile,"\t2: Infer major and minor from allele counts\n");
  fprintf(argFile,"\t3: use major and minor from a file (requires -sites file.txt)\n");
  fprintf(argFile,"\t4: Use reference allele as major (requires -ref)\n");
  fprintf(argFile,"\t5: Use ancestral allele as major (requires -anc)\n");
}

void abcMajorMinor::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;
  //default
  doMajorMinor=0;

  //below is used only validating if doMajorMinor has the data it needs
  int GL=0;
  int doCounts=0;

  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
   
  char *ref = NULL;
  char *anc = NULL;
  ref=angsd::getArg("-ref",ref,arguments);
  anc=angsd::getArg("-anc",anc,arguments);

  if(doMajorMinor==4&&ref==NULL){
    fprintf(stderr,"Must supply reference (-ref) when -doMajorMinor 4");
    exit(0);
  }
  if(doMajorMinor==5&&anc==NULL){
    fprintf(stderr,"Must supply ancestral (-anc) when -doMajorMinor 5");
    exit(0);
  }
  free(ref);
  free(anc);
  GL=angsd::getArg("-GL",GL,arguments);


 if(inputtype==INPUT_BEAGLE&&doMajorMinor){
   fprintf(stderr,"\t-> Potential problem: Cannot estimate the major and minor based on posterior probabilities\n");
   exit(0);
 }
 if((inputtype!=INPUT_GLF && doMajorMinor==1 && GL==0)){
   fprintf(stderr,"\t-> Potential problem: -doMajorMinor 1 is based on genotype likelihoods, you must specify a genotype likelihood model -GL \n");
   exit(0);
 }
 if(( (doMajorMinor==4 || doMajorMinor==5) && GL==0)){
   fprintf(stderr,"\t-> Potential problem: -doMajorMinor 4/5 use the genotype likelihoods to infer the minor alllee, you must specify a genotype likelihood model -GL \n");
   exit(0);
 }
 if(doMajorMinor==2&&doCounts==0){
   fprintf(stderr,"\t-> Potential problem: -doMajorMinor 2 is based on allele counts, you must specify -doCounts 1\n");
   exit(0);
   
 }
 
}

abcMajorMinor::abcMajorMinor(const char *outfiles,argStruct *arguments,int inputtype){
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doMajorMinor")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  if(doMajorMinor==0 || doMajorMinor ==3 )
    shouldRun[index] =0;

  printArg(arguments->argumentFile);
}

abcMajorMinor::~abcMajorMinor(){
}

void abcMajorMinor::clean(funkyPars *pars){
  if(doMajorMinor){
    delete [] pars->major;
    delete [] pars->minor;
    pars->major=pars->minor=NULL;
  }
}

void abcMajorMinor::print(funkyPars *pars){
}

// Function to swap major/minor reference/ancestral
void modMajorMinor(funkyPars *pars,int doMajorMinor){
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    
    if((doMajorMinor==4 &&pars->ref[s]==4)||(doMajorMinor==5 &&pars->anc[s]==4)){
      pars->keepSites[s]=0;
      continue;
    }
    int maj = pars->major[s];
    int mmin = pars->minor[s];
    

    if(doMajorMinor==4){
      if(pars->ref[s]!=maj&&pars->ref[s]!=mmin){
	//inferred major and minor is not part of the ref;
	pars->keepSites[s] =0;
	continue;
      }
      if(pars->ref[s]==mmin){
	//if reference is the inferre3d minor, then swap
	pars->major[s] = mmin;
	pars->minor[s] = maj;
      }
      //otherwise don't do anything
    }
    if(doMajorMinor==5){
      if(pars->anc[s]!=maj&&pars->anc[s]!=mmin){
	//inferred major and minor is not part of the anc;
	pars->keepSites[s] =0;
	continue;
      }
      if(pars->anc[s]==mmin){
	//if reference is the inferre3d minor, then swap
	pars->major[s] = mmin;
	pars->minor[s] = maj;
      }
      //otherwise don't do anything
    }    
  }
 
}

void majorMinorGL(funkyPars *pars,int doMajorMinor){
  
  float lmax;
  float totalLike;
  int choiceMajor;
  int choiceMinor;
  for(int s=0;s<pars->numSites;s++)  {
    if(pars->keepSites[s]==0){
      pars->major[s] = 4;
      pars->major[s] = 4;
      continue;
    }
    if(doMajorMinor==4&&refToInt[pars->ref[s]]==4){
      pars->keepSites[s]=0;
      continue;
    }
    if(doMajorMinor==5&&refToInt[pars->anc[s]]==4){
      pars->keepSites[s]=0;
      continue;
    }
    //if we have data then estimate the major/minor
    if(doMajorMinor<4) {
      lmax=-FLT_MAX;
      for(int Imajor=0;Imajor<3;Imajor++) {
	for(int Iminor=(Imajor+1);Iminor<4;Iminor++){
	  totalLike=0;
	  for(int i=0;i<pars->nInd;i++)
	    totalLike+=angsd::addProtect3(pars->likes[s][i*10+angsd::majorminor[Imajor][Imajor]]+log(0.25),
					  pars->likes[s][i*10+angsd::majorminor[Imajor][Iminor]]+log(0.5),
					  pars->likes[s][i*10+angsd::majorminor[Iminor][Iminor]]+log(0.25)
					  );
	  if(totalLike>lmax){
	    lmax=totalLike;
	    choiceMajor=Imajor;
	    choiceMinor=Iminor;
	  }
	}
//if we fix the major then skip after internal loop.
      }
      float W0;
      float W1;
      float W2;
      float sum=0;
      for(int i=0;i<pars->nInd;i++){
	W0=exp(pars->likes[s][i*10+angsd::majorminor[choiceMajor][choiceMajor]])*0.25;
	W1=exp(pars->likes[s][i*10+angsd::majorminor[choiceMajor][choiceMinor]])*0.5;
	W2=exp(pars->likes[s][i*10+angsd::majorminor[choiceMinor][choiceMinor]])*0.25;
	sum+=(W1+2*W2)/(2*(W0+W1+W2));
      }
      if(sum/pars->nInd<0.5 ){//last is because we don't want to swap if we fix majorminor
	pars->major[s]=choiceMajor;
	pars->minor[s]=choiceMinor;
      } else{
	pars->major[s]=choiceMinor;
	pars->minor[s]=choiceMajor;
      }
    }else{
      //this is the fixed major part
      lmax=-FLT_MAX;
      int Imajor=refToInt[doMajorMinor==4?pars->ref[s]:pars->anc[s]];
      for(int Iminor=0;Iminor<4;Iminor++){
	if(Imajor==Iminor)
	  continue;
	totalLike=0;
	for(int i=0;i<pars->nInd;i++)
	  totalLike+=angsd::addProtect3(pars->likes[s][i*10+angsd::majorminor[Imajor][Imajor]]+log(0.25),
					pars->likes[s][i*10+angsd::majorminor[Imajor][Iminor]]+log(0.5),
					pars->likes[s][i*10+angsd::majorminor[Iminor][Iminor]]+log(0.25)
					);
	if(totalLike>lmax){
	  lmax=totalLike;
	  choiceMajor=Imajor;
	  choiceMinor=Iminor;
	}
      }
      pars->major[s]=choiceMajor;
      pars->minor[s]=choiceMinor;

    }
    
  }
}

void majorMinorCounts(suint **counts,int nFiles,int nSites,char *major,char *minor,int *keepSites,int doMajorMinor,char *ref,char *anc) {
  assert(counts!=NULL);

  for(int s=0;s<nSites;s++){
    if(keepSites==0)
      continue;

    //part one
    //first lets get the sum of each nucleotide
    int bases[4] = {0,0,0,0};
    for(int i=0;i<nFiles;i++)
      for(size_t j=0;j<4;j++){
	bases[j] += counts[s][i*4+j];
	//	fprintf(oFile,"%d\t",bases[j]);
      }


    //now get the major/minor
    int majorA = 0;
    
    
    for(int i=1;i<4;i++)
      if (bases[i]>bases[majorA])
        majorA = i;
    //maj is now the major allele (most frequent)

    //should it be fixed to ancestral or reference?
    if(doMajorMinor==4)
      majorA =refToInt[ref[s]];
    if(doMajorMinor==4)
      majorA =refToInt[anc[s]];


    int temp=0;
    int minorA= majorA;
    for(int i=0;i<4;i++){
      if(majorA==i) //we should check the against the major allele      
        continue;
      else if (bases[i]>=temp){//always choose a minor eventhough non observed
        minorA = i;
        temp=bases[i];
      }
    }
    major[s] = majorA;
    minor[s] = minorA;
  }
}



void abcMajorMinor::run(funkyPars *pars){
  if(doMajorMinor==0 || doMajorMinor ==3 )
    return;
  if(doMajorMinor==1 && pars->likes==NULL){
    fprintf(stderr,"[%s.%s():%d] Problem here:\n",__FILE__,__FUNCTION__,__LINE__);
    exit(0);
  }

  
  //allocate and initialize
  pars->major = new char [pars->numSites];
  pars->minor = new char [pars->numSites];
  memset(pars->major,4,pars->numSites);
  memset(pars->minor,4,pars->numSites);
  
  if(doMajorMinor!=2)
    majorMinorGL(pars,doMajorMinor);
  else if(doMajorMinor==2)
    majorMinorCounts(pars->counts,pars->nInd,pars->numSites,pars->major,pars->minor,pars->keepSites,doMajorMinor,pars->ref,pars->anc);
  else
    fprintf(stderr,"[%s.%s()%d] Should never happen\n",__FILE__,__FUNCTION__,__LINE__);

  //if user has requested reference/ancestral then it is done in majorMinorGL and majorMinorCounts 0.585
  /*
  if(doMajorMinor==4||doMajorMinor==5)
    modMajorMinor(pars,doMajorMinor);
  */
}
