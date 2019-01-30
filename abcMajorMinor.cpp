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
#include "abcFilter.h"
static int isnaninf(double *d,int l){
  for(int i=0;i<l;i++)
    if(std::isnan(d[i])||std::isinf(d[i]))
      return 1;
  return 0;
}


static int isnan(double *d,int l){
  for(int i=0;i<l;i++)
    if(std::isnan(d[i]))
      return 1;
  return 0;
}

void abcMajorMinor::printArg(FILE *argFile){
  fprintf(argFile,"-------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doMajorMinor\t%d\n",doMajorMinor);
  fprintf(argFile,"\t1: Infer major and minor from GL\n");
  fprintf(argFile,"\t2: Infer major and minor from allele counts\n");
  fprintf(argFile,"\t3: use major and minor from a file (requires -sites file.txt)\n");
  fprintf(argFile,"\t4: Use reference allele as major (requires -ref)\n");
  fprintf(argFile,"\t5: Use ancestral allele as major (requires -anc)\n");
  fprintf(argFile,"\t-rmTrans: remove transitions %d\n",rmTrans);
  fprintf(argFile,"\t-skipTriallelic\t%d\n",skipTriallelic);
}

void abcMajorMinor::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;


  //below is used only validating if doMajorMinor has the data it needs
  int GL=0;
  int doCounts=0;
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);

  if(doMajorMinor==0)
    return;
  
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  rmTrans=angsd::getArg("-rmTrans",rmTrans,arguments);
   
  char *ref = NULL;
  char *anc = NULL;
  ref=angsd::getArg("-ref",ref,arguments);
  anc=angsd::getArg("-anc",anc,arguments);
  char *sites = NULL;
  sites = angsd::getArg("-sites",sites,arguments);
  if(sites==NULL&&doMajorMinor==3){
    fprintf(stderr,"\t-> You need to supply -sites for -domajorminor 3 to work. These has to have either 4 og 6 columns\n");
    exit(0);
  }
  free(sites);
  if(doMajorMinor==4&&ref==NULL){
    fprintf(stderr,"\t-> Must supply reference (-ref) when -doMajorMinor 4\n");
    exit(0);
  }
  if(doMajorMinor==5&&anc==NULL){
    fprintf(stderr,"\t-> Must supply ancestral (-anc) when -doMajorMinor 5\n");
    exit(0);
  }
  free(ref);
  free(anc);
  GL=angsd::getArg("-GL",GL,arguments);


 if(inputtype==INPUT_BEAGLE&&doMajorMinor){
   fprintf(stderr,"\t-> Potential problem: Cannot estimate the major and minor based on posterior probabilities\n");
   exit(0);
 }
 if(inputtype!=INPUT_GLF && inputtype!=INPUT_GLF3 && inputtype!=INPUT_VCF_GL && doMajorMinor==1 && GL==0 &&inputtype!=INPUT_GLF10_TEXT){
   fprintf(stderr,"\t-> Potential problem: -doMajorMinor 1 is based on genotype likelihoods, you must specify a genotype likelihood model -GL \n");
   exit(0);
 }

 if((doMajorMinor==4 || doMajorMinor==5) && GL==0){
   if(inputtype!=INPUT_GLF3 && inputtype!=INPUT_VCF_GL && inputtype!=INPUT_GLF){
     fprintf(stderr,"\t-> Potential problem: -doMajorMinor 4/5 use the genotype likelihoods to infer the minor allele, you must specify a genotype likelihood model -GL \n");
     exit(0);
   }
 }
 if(doMajorMinor==2&&doCounts==0){
   fprintf(stderr,"\t-> Potential problem: -doMajorMinor 2 is based on allele counts, you must specify -doCounts 1\n");
   exit(0);
   
 }
 doSaf = angsd::getArg("-doSaf",doSaf,arguments);
 pest = angsd::getArg("-pest",pest,arguments);
 skipTriallelic = angsd::getArg("-skipTriallelic",skipTriallelic,arguments);
 doVcf = angsd::getArg("-dobcf",doVcf,arguments);
 doGlf = angsd::getArg("-doGlf",doGlf,arguments);
}

abcMajorMinor::abcMajorMinor(const char *outfiles,argStruct *arguments,int inputtype){
  skipTriallelic=0;
  doMajorMinor = 0;
  rmTrans=0;

  doSaf = 0;
  pest = NULL;
  doVcf = 0;
  doGlf = 0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doMajorMinor")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  if(doMajorMinor==0)
    shouldRun[index] =0;
  if(doMajorMinor!=0)
  printArg(arguments->argumentFile);
}

abcMajorMinor::~abcMajorMinor(){

}

void abcMajorMinor::clean(funkyPars *pars){
  if(doMajorMinor){
    delete [] pars->major;
    delete [] pars->minor;
    pars->major=pars->minor=NULL;

    if(doVcf>0||doGlf==2){
      lh3struct *lh3 =(lh3struct*) pars->extras[index];
      for(int i=0;i<pars->numSites;i++)
	if(lh3->hasAlloced[i])
	  delete [] lh3->lh3[i];
      delete [] lh3->lh3;
      delete [] lh3->hasAlloced;
      
      delete lh3;
    }

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

void abcMajorMinor::majorMinorGL(funkyPars *pars,int doMajorMinor){
  
  float lmax;
  float totalLike;
  int choiceMajor=-1;
  int choiceMinor=-1;
  for(int s=0;s<pars->numSites;s++)  {

    if(pars->keepSites[s]==0){
      pars->major[s] = 4;
      pars->minor[s] = 4;
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
      if(choiceMajor==-1 || choiceMinor==-1){
	fprintf(stdout,"\t-> Something has gone wrong trying to infer major/minor from GLS at \'%s\' %d will discard site from analysis\n",header->target_name[pars->refId],pars->posi[s]+1);
	pars->keepSites[s]=0;
      }
      // assert(choiceMajor!=-1&&choiceMinor!=-1);
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
    //    fprintf(stderr,"AApos: %d maj:%d min:%d\n",pars->posi[s]+1,pars->major[s],pars->minor[s]);  
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
  if(doMajorMinor==0)
    return;
  if(doMajorMinor==1 && pars->likes==NULL){
    fprintf(stderr,"[%s.%s():%d] Problem here:\n",__FILE__,__FUNCTION__,__LINE__);
    exit(0);
  }
  extern abc **allMethods;
  filt *fl = ((abcFilter *) allMethods[0])->fl;
  pars->major = new char[pars->numSites];
  pars->minor = new char[pars->numSites];

  for(int i=0;i<pars->numSites;i++)
    pars->major[i]=pars->minor[i] = 4;


  //allocate and initialize
    
  //unless we want to base the majominor on counts we always use the gls
  if(doMajorMinor!=2 && doMajorMinor!=3)
    majorMinorGL(pars,doMajorMinor);
  else if(doMajorMinor==2)
    majorMinorCounts(pars->counts,pars->nInd,pars->numSites,pars->major,pars->minor,pars->keepSites,doMajorMinor,pars->ref,pars->anc);
  else if(doMajorMinor==3){
    for(int s=0;s<pars->numSites;s++){
      //      fprintf(stderr,"MMposi:%d keeps:%d\n",pars->posi[s]+1,pars->keepSites[s]);
      if(pars->keepSites[s]&&fl!=NULL&&fl->keeps[pars->posi[s]] && fl->hasExtra>0 &&doMajorMinor==3){
	pars->major[s] = fl->major[pars->posi[s]];
	pars->minor[s] = fl->minor[pars->posi[s]];
      }
      //fprintf(stderr,"MMposi2:%d keeps:%d maj:%d min:%d\n",pars->posi[s]+1,pars->keepSites[s],pars->major[s],pars->minor[s]);
    }
  }
  else
    fprintf(stderr,"[%s.%s()%d] Should never happen\n",__FILE__,__FUNCTION__,__LINE__);

  if(rmTrans){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s] == 0)
	continue;
      
      if((pars->minor[s]==0 && pars->major[s]==2) || (pars->minor[s]==2 && pars->major[s]==0 ) || (pars->minor[s]==1 && pars->major[s]==3) || (pars->minor[s]==3 && pars->major[s]==1) )
	pars->keepSites[s]=0;
      
    } 
  }
  
  if(doSaf!=0&&pest!=NULL&&skipTriallelic==1){
    //fix case of triallelic site in pest output when doing pest
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      int a=refToInt[pars->anc[s]];
      int b=refToInt[pars->major[s]];
      int c=refToInt[pars->minor[s]];
      if(a!=b&&a!=c)
	pars->keepSites[s]=0;
      
    }
    
  }
  
  //if user has requested reference/ancestral then it is done in majorMinorGL and majorMinorCounts 0.585
  /*
    if(doMajorMinor==4||doMajorMinor==5)
    modMajorMinor(pars,doMajorMinor);
  */
  if(doGlf==2||doVcf>0){
    lh3struct *lh3 = new lh3struct;
    pars->extras[index] = lh3;
    lh3->hasAlloced = new char[pars->numSites];
    memset(lh3->hasAlloced,0,pars->numSites);
    lh3->lh3 = new double *[pars->numSites];
    
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue;
      int major = pars->major[s];
      int minor = pars->minor[s];
      assert(major!=4&&minor!=4);
      
      lh3->hasAlloced[s]=1;
      lh3->lh3[s] = new double[3*pars->nInd];

      for(int i=0;i<pars->nInd;i++) {

	double val[3];
	val[0] = pars->likes[s][i*10+angsd::majorminor[major][major]];
	val[1] = pars->likes[s][i*10+angsd::majorminor[major][minor]];
	val[2] = pars->likes[s][i*10+angsd::majorminor[minor][minor]];
	if(isnan(val,3)){
	  pars->keepSites[s]=0;
	  break;
	}
	//fixed awkward case where all gls are -Inf, should only happen with -gl 6
	if(std::isinf(val[0])&&std::isinf(val[1])&&std::isinf(val[2]))
	  val[0]=val[1]=val[2]=0;
	angsd::logrescale(val,3);
	lh3->lh3[s][i*3+0]=val[0];
	lh3->lh3[s][i*3+1]=val[1];
	lh3->lh3[s][i*3+2]=val[2];
      }
      

    }
  }

}
