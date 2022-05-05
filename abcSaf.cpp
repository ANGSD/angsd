#include <cstdio>
#include <cmath>
#include <assert.h>
#include <cfloat>
#include <htslib/kstring.h>
#include "abcFreq.h"
#include "shared.h"
#include "analysisFunction.h"
#include <pthread.h>
#include "abc.h"
#include "abcSaf.h"
#include "aio.h"
#include <limits>

#define MINLIKE -1000.0 //this is for setting genotypelikelhoods to missing (EXPLAINED BELOW)

const double NEG_INF = -std::numeric_limits<double>::infinity();

int homo[4] = {0,4,7,9}; //AA,AC,AG,AT,CC,CG,CT,GG,GT,TT

namespace filipe{
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int *keepSites,realRes *r,int noTrans,int doSaf,char *major,char *minor,double *freq,double *indF,int newDim);
}

void abcSaf::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSaf\t\t%d\n",doSaf);
  fprintf(argFile,"\t   1: SAF calculation integrating over possible minor alleles\n\t   2: SAF calculation incorporating inbreeding\n\t   3: Calculate genotype probabilities using SAF (DEPRECATED; use -doPost 3)\n\t   4: SAF calculation from genotype posteriors (input is beagle text format)\n\t   5: SAF calculation conditioning on minor allele from -doMajorMinor\n");
  fprintf(argFile,"\t -underFlowProtect\t%d\n",underFlowProtect); 
  fprintf(argFile,"\t -anc\t\t%s\t(ancestral fasta)\n",anc);
  fprintf(argFile,"\t -noTrans\t%d\t(remove transitions)\n",noTrans);
  fprintf(argFile,"\t -pest\t\t%s\t(prior SFS)\n",pest);
  fprintf(argFile,"\t -isHap\t\t%d\t(samples are haploid; works with -doSaf 1 or 5)\n",isHap);
  fprintf(argFile,"\t -scoreTol\t%.1e\t(tolerance for score-limited algorithm)\n",scoreTol);
  fprintf(argFile,"\t -doPost\t%d\t(doPost 3, used for accessing SAF based variables)\n",doPost);
  fprintf(argFile,"\nNB: If -pest is supplied in addition to -doSaf then the output will be posterior probabilities of the sample allele frequency for each site\n");
  fprintf(argFile,"NB: Increasing -scoreTol will trade accuracy for reduced computation time and storage\n");
}

void abcSaf::getOptions(argStruct *arguments){
  doSaf = angsd::getArg("-doSaf",doSaf,arguments);
  doPost = angsd::getArg("-doPost",doPost,arguments);
  isHap = angsd::getArg("-isHap",isHap,arguments);
  isSim = angsd::getArg("-isSim",isSim,arguments);
  pest = angsd::getArg("-pest",pest,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);

  if(doSaf==3){ //DRAGON
    fprintf(stderr,"\t-> Please use -doPost 3 instead for -doSaf 3\n");
    exit(0);
  }

  if(doSaf>0||doPost==3){
    if(pest!=NULL){
      prior=angsd::readDouble(pest,arguments->nInd*2+1);
      int nd=arguments->nInd*2+1;

      double tts=0;
      for(int i=0;i<nd;i++)
        tts += prior[i];
      for(int i=0;i<nd;i++)
        prior[i] = log(prior[i]/tts);
    }
    lbicoTab = new double[2*arguments->nInd+1];
    tsktsktsk = 2*arguments->nInd+1;
    myComb2Tab = new double*[2*arguments->nInd+1];
    for(int i=0;i<2*arguments->nInd+1;i++){
      lbicoTab[i] = angsd::lbico(2*arguments->nInd,i);
      myComb2Tab[i] = new double[3];
      for(int j=0;j<3;j++)
        if(j<=i)
          myComb2Tab[i][j] = angsd::myComb2(arguments->nInd,i,j);
    }
    if(isHap){
      for(int i=0;i<arguments->nInd+1;i++)
        lbicoTab[i] = angsd::lbico(arguments->nInd,i);
    }
    mynchr = 2*arguments->nInd;
  }

  scoreTol = angsd::getArg("-scoreTol",scoreTol,arguments);

  noTrans = angsd::getArg("-noTrans",noTrans,arguments);

  int GL = 0;
  GL = angsd::getArg("-GL",GL,arguments);

  if(doSaf==0&&doPost!=3)
    return;

  underFlowProtect=angsd::getArg("-underFlowProtect",underFlowProtect,arguments);


  if(doSaf==0)
    return;
  if(arguments->regions.size()>1)
    fprintf(stderr,"\t-> !! You are doing -dosaf incombination with -rf, please make sure that your -rf file is sorted !!\n");
  anc = angsd::getArg("-anc",anc,arguments);
  if(doSaf && (anc==NULL&&isSim==0) ){
    if(doSaf!=3 && doSaf!=5 ){
      fprintf(stderr,"\t-> Must supply -anc for polarizing the spectrum\n");
      exit(0);
    }
  }



  if(doSaf==5){
    int doMajorMinor =0;
    doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
 
    if(doMajorMinor==0){
      
      fprintf(stderr,"\t-> for -doSaf 5 you must infer major and minor\n");
      exit(0);
    }
  }
  int ai = arguments->inputtype;
  if(GL==0 &&(ai!=INPUT_GLF && ai !=INPUT_GLF3 && ai !=INPUT_VCF_GL &&ai !=INPUT_BEAGLE&&ai!=INPUT_GLF10_TEXT)){
    fprintf(stderr,"\t-> Must supply genotype likelihoods (-GL [INT])\n");
    printArg(arguments->argumentFile);
    exit(0);
  }
  if(doSaf==2){
    fprintf(stderr,"\t-> (Using Filipe G Vieira modification of: %s)\n",__FILE__);
    int doMajorMinor =0;
    doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
    int doMaf =0;
    doMaf = angsd::getArg("-doMaf",doMaf,arguments);
    if(doMajorMinor==0||doMaf==0){
      fprintf(stderr,"\t-> Must have major/minor and MAF for using the inbreeding version of realSFS\n");
      exit(0);
    }
    char *indF_name = NULL;
    indF_name =  angsd::getArg("-indF",indF_name,arguments);
    if(indF_name==NULL){
      filipeIndF = new double[arguments->nInd];
      for(int i=0;i<arguments->nInd;i++)
	filipeIndF[i] =0;
      fprintf(stderr,"\t-> No -indF file provided will assume an inbreeding zero for all samples.\n");
      fprintf(stderr,"\t-> If no inbreeding is expected consider using -doSaf 5\n");
    }else
      filipeIndF = angsd::readDouble(indF_name,arguments->nInd);
  }
  if(doSaf==3){
    fprintf(stderr,"\t-> Will call genotypes using sample allele frequencies\n");
    int doMajorMinor =0;
    doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
    if(doMajorMinor==0){
      fprintf(stderr,"\t-> Must supply -doMajorMinor for calling genotypes\n");
      exit(0);
    }
    if(pest==NULL){
      fprintf(stderr,"\t-> You need to supply a sfs as prior (./realSFS output, -pest) to do genotypecalling with saf\n");
      exit(0);
    }
  }

}

double *abcSaf::lbicoTab = NULL;
double **abcSaf::myComb2Tab=NULL;
double *abcSaf::prior = NULL;

abcSaf::abcSaf(const char *outfiles,argStruct *arguments,int inputtype){
  tsktsktsk = 0;
  tmpChr = NULL;
  isHap = 0;
  sumBand = 0;
  minInd=0;
  //for use when dumping binary indexed saf files
  const char *SAF = ".saf.gz";
  const char *SAFPOS =".saf.pos.gz";
  const char *SAFIDX =".saf.idx";
  //for use when dumping called genotypes
  const char *GENO = ".saf.geno.gz";
  //default
  underFlowProtect = 0;
  isSim =0;
  //from command line
  anc = NULL;
  pest = NULL;
  noTrans = 0;
  prior = NULL;
  doSaf = 0;
  doPost = 0;
  outfileSAF = NULL;
  outfileSAFPOS = NULL;
  outfileSAFIDX = NULL;
  outfileGprobs = NULL;
  nnnSites = 0;
  scoreTol = 1.e-9;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSaf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);


  if(doSaf==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);  
  newDim = 2*arguments->nInd+1;
 
  if(isHap)
    newDim = arguments->nInd+1;
  if(doSaf==3){
    fprintf(stderr,"\t-> -doSaf 3 Does not work in combination with -snpstat and other snpfilters\n");
    outfileGprobs = aio::openFileBG(outfiles,GENO);
  }else if(doSaf!=0){
    outfileSAF =  aio::openFileBG(outfiles,SAF);
    outfileSAFPOS =  aio::openFileBG(outfiles,SAFPOS);
    outfileSAFIDX = aio::openFile(outfiles,SAFIDX);
    char buf[8]="safv4";
    aio::bgzf_write(outfileSAF,buf,8);
    aio::bgzf_write(outfileSAFPOS,buf,8);
    fwrite(buf,1,8,outfileSAFIDX);
    assert(bgzf_flush(outfileSAF)==0);assert(bgzf_flush(outfileSAFPOS)==0);
    offs[0] = bgzf_tell(outfileSAFPOS);
    offs[1] = bgzf_tell(outfileSAF);
    size_t tt = newDim-1;
    fwrite(&tt,sizeof(tt),1,outfileSAFIDX);
  }

}

abcSaf::~abcSaf(){
  if(doSaf&&doSaf!=3)
    writeAll();

  if(pest) free(pest);
  if(prior) delete [] prior;
  if(outfileSAF) bgzf_close(outfileSAF);;
  if(outfileSAFPOS) bgzf_close(outfileSAFPOS);
  if(outfileSAFIDX) fclose(outfileSAFIDX);

  if(outfileGprobs)  bgzf_close(outfileGprobs);
  if(lbicoTab) {
    delete [] lbicoTab;
  }
  if(myComb2Tab){
    for(int i=0;i<tsktsktsk;i++)
      delete [] myComb2Tab[i];
    delete [] myComb2Tab;
  }
  if(tmpChr) free(tmpChr);
  if(anc) free(anc);
}

void normalize_array(double *d, int len){
  double s =0;
  for(int i=0;i<len;i++)
    s+=d[i];

  for(int i=0;i<len;i++)
    d[i]=d[i]/s;
}

void normalize_array2(double *d, int len){
  double s =0;
  for(int i=0;i<len;i++)
    s+=exp(d[i]);
  s=log(s);

  for(int i=0;i<len;i++)
    d[i]=d[i]-s;
}

int isSame(double a,double b,double tolerance){
  return (fabs(a-b)<tolerance);
}

// --- banded SAF algorithm adapted from Han & Novembre 2015 Bioinformatics --- //

double logSumExp (double a, double b)
{
  // b/c angsd::addProtect2 farts out when both terms are -Inf
  double mx = a > b ? a : b;
  if (std::isinf(mx)) mx = 0.; 
  a = exp(a - mx);
  b = exp(b - mx);
  return log(a + b) + mx;
}

double saf_like_hap (double *p, double *h, int N, int j)
{
  double tmp;
  tmp = 0.;

  if(j >= 0 && j < N)
    tmp += (N-j)*p[0]*h[j];

  if(j-1 >= 0)
    tmp += j*p[1]*h[j-1];

  if(std::isnan(tmp))
  {
    fprintf(stderr, "is nan: %d\n", j);
    tmp = 0.;
  }

  return tmp;
}

void banded_saf_algo_hap (double* hj, int& lower, int& upper, double* p, const int i, const int numChr, const double tol)
{
  int mle = p[0] > p[1] ? 0 : 1; //missing: 1
  lower += mle;
  upper += mle;

  while (true)
  {
    if (lower == 0 || saf_like_hap(p, hj, numChr, lower) < tol)
      break;
    lower -= 1;
  }

  while (true)
  {
    if (upper == numChr || saf_like_hap(p, hj, numChr, upper) < tol)
      break;
    upper += 1;
  }

  for (int j=upper; j>=lower; --j)
    hj[j] = saf_like_hap (p, hj, numChr, j);
}

void saf_algo_hap (double* hj, int& lower, int& upper, double& sm, double& score_tol, double* p, const int i, const int numChr)
{
  int lower_old = lower,
      upper_old = upper;

  banded_saf_algo_hap(hj, lower, upper, p, i, numChr, score_tol);

  // normalize
  double den = 0.;
  for (int j=lower; j<=upper; ++j)
    if (hj[j] > den)
      den = hj[j];
  for (int j=lower; j<=upper; ++j)
    hj[j] /= den;

  // clean up edges
  for (int j=lower_old; j<lower; ++j) 
    hj[j] = 0.;
  for (int j=upper_old; j>upper; --j) 
    hj[j] = 0.;

  // track normalizing constant
  sm += log(den);
}

double saf_like_dip (double *p, double *h, int N, int j)
{
  double tmp;
  tmp = 0.;

  if(j >= 0 && j < N)
    tmp += (N-j)*(N-j-1)*p[0]*h[j];

  if(j-1 >= 0 && j < N)
    tmp += 2*j*(N-j)*p[1]*h[j-1];

  if(j-2 >= 0 && j <= N)
    tmp += j*(j-1)*p[2]*h[j-2];

  if(std::isnan(tmp))
  {
    fprintf(stderr, "is nan: %d\n", j);
    tmp = 0.;
  }

  return tmp;
}

void vanilla_saf_algo_dip (double* hj, int& lower, int& upper, double* p, const int i, const int numChr, const double tol)
{
  lower = 0;
  upper = 2*(i+1);
  
  for(int j=upper; j>=lower; j--)
    hj[j] = saf_like_dip(p, hj, numChr, j);
}

void banded_saf_algo_dip (double* hj, int& lower, int& upper, double* p, const int i, const int numChr, const double tol)
{
  int mle = (p[0] > p[1] && p[0] > p[2]) ? 0 : (p[1] > p[2] ? 1 : 2); //missing: 2
  lower += mle;
  upper += mle;

  while (true)
  {
    if (lower == 0 || saf_like_dip(p, hj, numChr, lower) < tol)
      break;
    lower -= 1;
  }

  while (true)
  {
    if (upper == numChr || saf_like_dip(p, hj, numChr, upper) < tol)
      break;
    upper += 1;
  }

  for (int j=upper; j>=lower; --j)
    hj[j] = saf_like_dip (p, hj, numChr, j);
}

void saf_algo_dip (double* hj, int& lower, int& upper, double& sm, double& score_tol, double* p, const int i, const int numChr)
{
  int lower_old = lower,
      upper_old = upper;

  if (p[1] < p[0] && p[1] < p[2]) // ... then update may not be unimodal, should pass over all bins
    score_tol = 0.; //this will force all future updates for this site to pass over all bins

  if (score_tol > 0.)
    banded_saf_algo_dip(hj, lower, upper, p, i, numChr, score_tol);
  else
    vanilla_saf_algo_dip(hj, lower, upper, p, i, numChr, score_tol);

  // normalize
  double den = 0.;
  for (int j=lower; j<=upper; ++j)
    if (hj[j] > den) 
      den = hj[j];
  for (int j=lower; j<=upper; ++j)
    hj[j] /= den;

  // clean up edges
  for (int j=lower_old; j<lower; ++j) 
    hj[j] = 0.;
  for (int j=upper_old; j>upper; --j) 
    hj[j] = 0.;

  // track normalizing constant
  sm += log(den);
}

int saf_sparsify_and_normalize (double* hj, int& lower, int& upper, const double tol)
{
  // rescale in logspace
  double mx = NEG_INF;
  for (int j=lower; j<=upper; ++j)
    if (hj[j] > mx)
      mx = hj[j];
  for (int j=lower; j<=upper; ++j)
    hj[j] -= mx;

  // another pass to thin band
  bool big_enough;
  int lower0 = lower,
      upper0 = upper;
  lower = upper0;
  upper = lower0;
  for (int j=lower0; j<=upper0; ++j)
  {
    big_enough = hj[j] >= log(tol);
    if (big_enough && j > upper)
      upper = j;
    if (big_enough && j < lower)
      lower = j;
  }

  if (lower == upper0 || upper == lower0)
  { // unlikely, but could happen if there's underflow
    lower = lower0;
    upper = upper0;
    // fprintf(stderr, "\t-> Banding failed\n");
    return 1;
  }
  return 0;
}

// --- different options for -doSaf --- //

// NSP 3July2020
// I did not implement banded/score-limited algo for this.
// However, I think there is a bug here -- the prior gets multiplied into the genotype likelihoods twice.
// I would not trust the output of this method until this gets investigated.
void filipe::algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int *keepSites,realRes *r,int noTrans,int doSaf,char *major,char *minor,double *freq,double *indF,int newDim) {
  //  fprintf(stderr,"liks=%p anc=%p nsites=%d nInd=%d underflowprotect=%d fold=%d keepSites=%p r=%p\n",liks,anc,nsites,numInds,underFlowProtect,fold,keepSites,r);
  assert(doSaf==2);
  int myCounter =0;

  if(anc==NULL||liks==NULL){
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }
  double sumMinors[2*numInds+1];  //the sum of the 3 different minors
  double m[numInds][3];  //HWE priors

  for(int it=0; it<nsites; it++)  {//loop over sites
    int ancB,derB;
    ancB= anc[it];
    derB=-1;

    if(ancB==4 || keepSites[it]==0){//skip if no ancestral information
      keepSites[it] =0; //
      continue;
    }
    //if the ancestral is neither major or the minor skip site
    if(ancB!=major[it]&&ancB!=minor[it]){
      keepSites[it] =0; //
      continue;
    }

    //set the resultarray to zeros
    for(int sm=0 ; sm<(2*numInds+1) ; sm++ ) sumMinors[sm] = 0;
    
    // Assign freqs to ancestral and non-ancestral alleles
    double anc_freq, non_anc_freq;
    
    anc_freq = non_anc_freq = freq[it];
    if(ancB == major[it]){
      anc_freq = 1-freq[it];
      derB = minor[it];
    }else {
      non_anc_freq = 1-freq[it];
      derB = major[it];
    }

    if(noTrans){
      if((ancB==2&&derB==0)||(ancB==0&&derB==2))
	continue;
      if((ancB==1&&derB==3)||(ancB==3&&derB==1))
	continue;
    }

    assert(ancB!=-1&&derB!=-1);
    
    double totmax = 0.0;
    
    int Aa_offset = angsd::majorminor[derB][ancB];//0-9
    int AA_offset = angsd::majorminor[derB][derB];//0-9
    int aa_offset = angsd::majorminor[ancB][ancB];//0-9
    
    double hj[2*numInds+1];
    for(int index=0;index<(2*numInds+1);index++)
      if(underFlowProtect==0)
	hj[index]=0;
      else
	hj[index]=log(0);
    
    double PAA,PAa,Paa;
    
    for(int i=0 ; i<numInds ;i++) {
      m[i][0] = pow(anc_freq,2.0) + anc_freq*non_anc_freq*indF[i];
      m[i][1] = 2.0*anc_freq*non_anc_freq - 2.0*anc_freq*non_anc_freq*indF[i];
      m[i][2] = pow(non_anc_freq,2.0) + anc_freq*non_anc_freq*indF[i];
    
      double GAA,GAa,Gaa;

      //fix issues where frequency is estimated to zero. Returns very small value
      GAA = ((m[i][2] != 0) ? (log(m[i][2])+liks[it][i*10+AA_offset]) : log(DBL_MIN));
      GAa = ((m[i][1] != 0) ? (log(m[i][1])+liks[it][i*10+Aa_offset]) : log(DBL_MIN));
      Gaa = ((m[i][0] != 0) ? (log(m[i][0])+liks[it][i*10+aa_offset]) : log(DBL_MIN));
      
      double mymax;
      if (Gaa > GAa && Gaa > GAA) mymax = Gaa;
      else if (GAa > GAA) mymax = GAa;
      else mymax = GAA;
      
      if(mymax<MINLIKE){
	Gaa = 0;
	GAa = 0;
	GAA = 0;
	totmax = totmax + mymax;
      }else{
	Gaa=Gaa-mymax;
	GAa=GAa-mymax;
	GAA=GAA-mymax;
	totmax = totmax + mymax;
      }
      
      if(underFlowProtect==0){
	PAA=exp(GAA);
	PAa=exp(GAa);
	Paa=exp(Gaa);
      }else{
	PAA =(GAA);
	PAa =(GAa);
	Paa =(Gaa);
      }
      
      
      if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa))
	fprintf(stderr,"Possible underflow  PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);

      if(i==0){
	hj[0] =Paa;
	hj[1] =PAa;
	hj[2] =PAA;
      }else{
	for(int j=2*(i+1); j>1;j--){
	  double tmp;
	  if(underFlowProtect==1)
	    tmp = angsd::addProtect3(log(m[i][2])+PAA+hj[j-2],log(m[i][1])+PAa+hj[j-1],log(m[i][0])+Paa+hj[j]);
	  else
	    tmp = m[i][2]*PAA*hj[j-2] + m[i][1]*PAa*hj[j-1] + m[i][0]*Paa*hj[j];
	  
	  if(std::isnan(tmp)){
	      fprintf(stderr,"is nan:%d\n",j );
	      hj[j] = 0;
	      break;
	  }else
	    hj[j] = tmp;
	}
	if(underFlowProtect==1){
	  hj[1] = angsd::addProtect2(log(m[i][0])+Paa+hj[1],log(m[i][1])+PAa+hj[0]);
	  hj[0] = log(m[i][0]) + Paa + hj[0];
	}
	else{
	  hj[1] = m[i][0]*Paa*hj[1] + m[i][1]*PAa*hj[0];
	  hj[0] = m[i][0]*Paa*hj[0];
	}
      }
    }

    for(int i=0;i<(2*numInds+1);i++)
      if(underFlowProtect==0){
	sumMinors[i] += hj[i];
	// sumMinors[i] += hj[i]/(1-hj[0]-hj[2*numInds]); //As in the PLoS ONE paper
      }else{
	sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]));
	//sumMinors[i] =
	//	exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-log(1-hj[0]-hj[2*numInds])));
	//As in the PLoS ONE paper
      }
     

    //sumMinors is in normal space, not log
    /*
      we do 3 things.
      1. log scaling everyting
      2. rescaling to the most likely in order to avoid underflows in the optimization
      3. we might do a fold also.
     */    

    for(int i=0;i<2*numInds+1;i++)
      sumMinors[i] = log(sumMinors[i]);
    //      angsd::logrescale(sumMinors,2*numInds+1);
    double ts=0;
    for(int i=0;i<newDim;i++)
      ts += exp(sumMinors[i]);
    ts = log(ts);
    for(int i=0;i<newDim;i++)
      sumMinors[i] = sumMinors[i]-ts;
    
    if(std::isnan(sumMinors[0]))
      r->oklist[it] = 2;
    else{
      r->oklist[it] = 1;
      r->pLikes[myCounter] =new float[2*numInds+1];
      // for now, to make this safv4
      r->pBound[myCounter] = new int[2];
      r->pBound[myCounter][0] = 0;
      r->pBound[myCounter][1] = 2*numInds+1;
      // \for now
      for(int iii=0;iii<2*numInds+1;iii++)
        r->pLikes[myCounter][iii] = sumMinors[iii];
      
      //	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
      myCounter++;
    }
  }
}

void abcSaf::algoJointPost(double **post, 
                           int nSites, 
                           int nInd, 
                           int *keepSites, 
                           realRes *r) 
{


  int counter = 0;
  int numChr = nInd*2;

  for(int s=0; s<nSites; s++)
  {
    if(keepSites[s]==0)
      continue;

    double *liks = post[s]; //we call this liks, even though it is posteriors.
    double hj[numChr+1];
    for(int i=0; i<numChr+1; i++)
      hj[i] = 0;
    double tmx = 0.;

    //initialize
    memcpy(hj, liks, 3*sizeof(double));

    int lower = 0,
        upper = 2;
    double score_tol = scoreTol;
    double p[3];
    for(int i=1; i<nInd; i++) 
    {
      p[0] = liks[i*3];
      p[1] = liks[i*3+1];
      p[2] = liks[i*3+2];
      saf_algo_dip(hj, lower, upper, tmx, score_tol, p, i, 2*(i+1));
    }

    if(saf_sparsify_and_normalize (hj, lower, upper, scoreTol))
      r->oklist[s] = 3;

    if(std::isnan(hj[lower]))
      r->oklist[s] = 2;
    else 
    {
      r->oklist[s] = 1;
      r->pLikes[counter] = new float[upper - lower + 1];
      r->pBound[counter] = new int[2];

      int k=0;
      for(int j=lower; j<=upper; j++)
        r->pLikes[counter][k++] = hj[j];

      r->pBound[counter][0] = lower;
      r->pBound[counter][1] = upper - lower + 1;

      ////debug
      //fprintf(stdout, "%u\t%u\t%u", counter, lower, upper-lower+1);
      //k=0;
      //for(int j=lower; j<=upper; ++j)
      //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
      //fprintf(stdout, "\n");
      
      counter++;
    }
  }
}


void abcSaf::algoJointHap(double **liks,
                          char *anc,int nsites,
                          int numInds,
                          int *keepSites,
                          realRes *r,
                          int noTrans) 
{
  int counter = 0;
  int numChr = numInds;

  if(anc==NULL||liks==NULL)
  {
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }

  double sumMinors[numChr+1];

  for(int it=0; it<nsites; it++) 
  {
    int major_offset = anc[it];
    if(major_offset==4 || keepSites[it]==0)
    { //skip if no ancestral information
      keepSites[it] = 0; 
      continue;
    }

    for(int j=0; j<numChr+1; j++)
      sumMinors[j] = NEG_INF;

    int lower_all = numChr,
        upper_all = 0;
    
    //loop through the 3 different minors
    for(int minor_offset=0; minor_offset<4; minor_offset++) 
    {
      if(minor_offset == major_offset)
        continue;

      if(noTrans)
      {
        if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
          continue;
        if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
          continue;
      }

      double tmx = 0.;
      int AA_offset = homo[major_offset];
      int aa_offset = homo[minor_offset];

      double p[2];
      double score_tol = scoreTol;
      double hj[numChr+1];
      for(int j=0; j<numChr+1; j++) hj[j] = 0.;
      int lower = 0,
          upper = 1;

      for(int i=0; i<numInds; i++) 
      {
        p[0] = liks[it][i*10+AA_offset];
        p[1] = liks[it][i*10+aa_offset];

        //underflow protection
        double mx = p[1] > p[0] ? p[1] : p[0];
        tmx += mx;
	  
        p[0] = mx < MINLIKE ? 0. : exp(p[0] - mx);
        p[1] = mx < MINLIKE ? 0. : exp(p[1] - mx);

        //check for underflow error, this should only occur once in a blue moon
        if(std::isnan(p[0])||std::isnan(p[1]))
          fprintf(stderr,"PAA=%f\tPaa=%f\n",p[1],p[0]);

        if(i==0)
        {
          hj[0] = p[0];
          hj[1] = p[1];
        }
        else
          saf_algo_hap(hj, lower, upper, tmx, score_tol, p, i, i+1);
      }
      
      for(int j=lower; j<=upper; j++)
        sumMinors[j] = logSumExp(sumMinors[j], log(hj[j])+tmx);

      if (lower < lower_all)
        lower_all = lower;
      if (upper > upper_all)
        upper_all = upper;
    }

    if(saf_sparsify_and_normalize (sumMinors, lower_all, upper_all, scoreTol))
      r->oklist[it] = 3;
    if(std::isnan(sumMinors[lower_all]))
      r->oklist[it] = 2;
    else
    {
      r->oklist[it] = 1;
      r->pLikes[counter] = new float[upper_all-lower_all+1];
      r->pBound[counter] = new int[2];

      int k = 0;
      for(int j=lower_all; j<=upper_all; ++j)
        r->pLikes[counter][k++] = sumMinors[j];

      r->pBound[counter][0] = lower_all;
      r->pBound[counter][1] = upper_all-lower_all+1;

      ////debug
      //fprintf(stdout, "%u\t%u\t%u", counter, lower_all, upper_all-lower_all+1);
      //k=0;
      //for(int j=lower_all; j<=upper_all; ++j)
      //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
      //fprintf(stdout, "\n");

      counter++;
    }
  }
}

void abcSaf::algoJoint(double **liks, 
                       char *anc, 
                       int nsites, 
                       int numInds,
                       int *keepSites,
                       realRes *r,
                       int noTrans) 
{
  int counter = 0;
  int numChr = 2*numInds;

  if(anc==NULL||liks==NULL)
  {
    fprintf(stderr, "\t-> problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n", __FUNCTION__, liks, anc);
    exit(0);
  }

  double sumMinors[numChr+1]; //the sum of the 3 different minors

  for(int it=0; it<nsites; it++) 
  { //loop over sites

    int major_offset = anc[it];
    if(major_offset==4 || keepSites[it]==0)
    { //skip if no ancestral information
      keepSites[it] = 0;
      continue;
    }

    { //check that fixing the ancestral produces meaningful results (anc is not major or minor)
      //AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
      double glsum[4] = {0,0,0,0};
      for(int ss=0; ss<numInds; ss++){
        glsum[0] += liks[it][10*ss];
        glsum[1] += liks[it][10*ss+4];
        glsum[2] += liks[it][10*ss+7];
        glsum[3] += liks[it][10*ss+9];
      }
      int howmanysmaller = 0;
      for(int i=0; i<4; i++){
        //fprintf(stderr,"anc: %d gl[%d]: %f\n",anc[it],i,glsum[i]);
        if(glsum[i]<=glsum[anc[it]])
          howmanysmaller++;
      }
      if(0&&howmanysmaller==2)
        continue;
      //fprintf(stderr,"How many smaller: %d\n",howmanysmaller);
    }

    for(int j=0; j<numChr+1; j++)
      sumMinors[j] = NEG_INF;

    int lower_all = numChr,
        upper_all = 0;

    
    //loop through the 3 different minors
    for(int minor_offset=0; minor_offset<4; minor_offset++) 
    {
      if(minor_offset == major_offset)
        continue;

      if(noTrans)
      {
        if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
          continue;
        if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
          continue;
      }

      double tmx = 0.; //denominator in logspace
      int Aa_offset = angsd::majorminor[minor_offset][major_offset];
      int AA_offset = angsd::majorminor[minor_offset][minor_offset];
      int aa_offset = angsd::majorminor[major_offset][major_offset];

      double p[3];
      double score_tol = scoreTol;
      double hj[numChr+1];
      for(int j=0; j<numChr+1; j++) hj[j] = 0.;
      int lower = 0,
          upper = 2;

      for(int i=0; i<numInds; i++) 
      {
        p[0] = liks[it][i*10+aa_offset];
        p[1] = liks[it][i*10+Aa_offset];
        p[2] = liks[it][i*10+AA_offset];
	
        //underflow protection
        double mx;
        if (p[2] > p[1] && p[2] > p[0]) mx = p[2];
        else if (p[1] > p[0]) mx = p[1];
        else mx = p[0];
        tmx += mx;
	  
        p[0] = mx < MINLIKE ? 1. : exp(p[0] - mx);
        p[1] = mx < MINLIKE ? 1. : exp(p[1] - mx);
        p[2] = mx < MINLIKE ? 1. : exp(p[2] - mx);

        //check for underflow error, this should only occur once in a blue moon
        if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[2]))
          fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",p[2],p[1],p[0]);

        if(i==0)
        {
          hj[0] = p[0];
          hj[1] = p[1];
          hj[2] = p[2];
        }
        else
          saf_algo_dip(hj, lower, upper, tmx, score_tol, p, i, 2*(i+1));
      }

      for(int j=lower; j<=upper; j++)
        sumMinors[j] = logSumExp(sumMinors[j], log(hj[j])+tmx);

      if (lower < lower_all)
        lower_all = lower;
      if (upper > upper_all)
        upper_all = upper;
    }

    if(saf_sparsify_and_normalize (sumMinors, lower_all, upper_all, scoreTol))
      r->oklist[it] = 3;
    if(std::isnan(sumMinors[lower_all]))
      r->oklist[it] = 2;
    else
    {
      r->oklist[it] = 1;
      r->pLikes[counter] = new float[upper_all-lower_all+1];
      r->pBound[counter] = new int[2];

      int k = 0;
      for(int j=lower_all; j<=upper_all; ++j)
        r->pLikes[counter][k++] = sumMinors[j];

      r->pBound[counter][0] = lower_all;
      r->pBound[counter][1] = upper_all-lower_all+1;

      ////debug
      //fprintf(stdout, "%u\t%u\t%u", counter, lower_all, upper_all-lower_all+1);
      //k=0;
      //for(int j=lower_all; j<=upper_all; ++j)
      //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
      //fprintf(stdout, "\n");

      counter++;
    }
  }
}

void abcSaf::algoJointMajorMinor(double **liks,
                                 int nsites,
                                 int numInds, 
                                 int *keepSites,
                                 realRes *r,
                                 char *major, 
                                 char *minor) 
{
  int counter = 0;
  int numChr = 2*numInds;

  if(liks==NULL)
  {
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p)\n", __FUNCTION__, liks);
    exit(0);
  }

  for(int it=0; it<nsites; it++) 
  {
    int major_offset = major[it];
    if(major_offset==4||keepSites[it]==0)
    { //skip if no major information
      keepSites[it] = 0;
      continue;
    }

    int minor_offset = minor[it];

    if(minor_offset == major_offset) //when would this happen?
      continue;

    if(noTrans)
    {
      if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
        continue;
      if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
        continue;
    }

    double tmx = 0.;
    int Aa_offset = angsd::majorminor[minor_offset][major_offset];
    int AA_offset = angsd::majorminor[minor_offset][minor_offset];
    int aa_offset = angsd::majorminor[major_offset][major_offset];

    double p[3];
    double score_tol = scoreTol;
    double hj[numChr+1];
    for(int j=0; j<numChr+1; j++) hj[j] = 0.;
    int lower = 0,
        upper = 2;

    for(int i=0; i<numInds; i++) 
    {
      p[0] = liks[it][i*10+aa_offset];
      p[1] = liks[it][i*10+Aa_offset];
      p[2] = liks[it][i*10+AA_offset];

      //underflow protection
      double mx;
      if (p[2] > p[1] && p[2] > p[0]) mx = p[2];
      else if (p[1] > p[0]) mx = p[1];
      else mx = p[0];
      tmx += mx;

      p[0] = mx < MINLIKE ? 0. : exp(p[0] - mx);
      p[1] = mx < MINLIKE ? 0. : exp(p[1] - mx);
      p[2] = mx < MINLIKE ? 0. : exp(p[2] - mx);

      //check for underflow error, this should only occur once in a blue moon
      if(std::isnan(p[0])||std::isnan(p[1])||std::isnan(p[2]))
        fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",p[2],p[1],p[0]);

      if(i==0)
      {
        hj[0] = p[0];
        hj[1] = p[1];
        hj[2] = p[2];
      }
      else
        saf_algo_dip(hj, lower, upper, tmx, score_tol, p, i, 2*(i+1));
    }

    for(int j=lower; j<=upper; j++)
      hj[j] = log(hj[j]) + tmx;

    if(saf_sparsify_and_normalize (hj, lower, upper, scoreTol))
      r->oklist[it] = 3;

    if(std::isnan(hj[lower]))
      r->oklist[it] = 2;
    else
    {
      r->oklist[it] = 1;
      r->pLikes[counter] = new float[upper-lower+1];
      r->pBound[counter] = new int[2];

      int k = 0;
      for(int j=lower; j<=upper; ++j)
        r->pLikes[counter][k++] = hj[j];

      r->pBound[counter][0] = lower;
      r->pBound[counter][1] = upper-lower+1;

      ////debug
      //fprintf(stdout, "%u\t%u\t%u", counter, lower, upper-lower+1);
      //k=0;
      //for(int j=lower; j<=upper; ++j)
      //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
      //fprintf(stdout, "\n");

      counter++;
    }
  } 
}

void abcSaf::algoJointMajorMinorHap(double **liks,
                                    int nsites,
                                    int numInds, 
                                    int *keepSites,
                                    realRes *r,
                                    char *major, 
                                    char *minor) 
{
  int counter = 0;
  int numChr = numInds;

  if(liks==NULL)
  {
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p)\n", __FUNCTION__, liks);
    exit(0);
  }

  for(int it=0; it<nsites; it++) 
  {
    int major_offset = major[it];
    if(major_offset==4||keepSites[it]==0)
    { //skip if no major information
      keepSites[it] = 0;
      continue;
    }

    int minor_offset = minor[it];

    if(minor_offset == major_offset) //when would this happen?
      continue;

    if(noTrans)
    {
      if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
        continue;
      if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
        continue;
    }

    double tmx = 0.;
    int AA_offset = homo[major_offset];
    int aa_offset = homo[minor_offset];

    double p[2];
    double score_tol = scoreTol;
    double hj[numChr+1];
    for(int j=0; j<numChr+1; j++) hj[j] = 0.;
    int lower = 0,
        upper = 1;

    for(int i=0; i<numInds; i++) 
    {
      p[0] = liks[it][i*10+AA_offset];
      p[1] = liks[it][i*10+aa_offset];

      //underflow protection
      double mx = p[1] > p[0] ? p[1] : p[0];
      tmx += mx;

      p[0] = mx < MINLIKE ? 0. : exp(p[0] - mx);
      p[1] = mx < MINLIKE ? 0. : exp(p[1] - mx);

      //check for underflow error, this should only occur once in a blue moon
      if(std::isnan(p[0])||std::isnan(p[1]))
        fprintf(stderr,"PAA=%f\tPaa=%f\n",p[1],p[0]);

      if(i==0)
      {
        hj[0] = p[0];
        hj[1] = p[1];
      }
      else
        saf_algo_hap(hj, lower, upper, tmx, score_tol, p, i, i+1);
    }

    for(int j=lower; j<=upper; j++)
      hj[j] = log(hj[j]) + tmx;

    if(saf_sparsify_and_normalize (hj, lower, upper, scoreTol))
      r->oklist[it] = 3;

    if(std::isnan(hj[lower]))
      r->oklist[it] = 2;
    else
    {
      r->oklist[it] = 1;
      r->pLikes[counter] = new float[upper-lower+1];
      r->pBound[counter] = new int[2];

      int k = 0;
      for(int j=lower; j<=upper; ++j)
        r->pLikes[counter][k++] = hj[j];

      r->pBound[counter][0] = lower;
      r->pBound[counter][1] = upper-lower+1;

      ////debug
      //fprintf(stdout, "%u\t%u\t%u", counter, lower, upper-lower+1);
      //k=0;
      //for(int j=lower; j<=upper; ++j)
      //  fprintf(stdout, "\t%f", r->pLikes[counter][k++]);
      //fprintf(stdout, "\n");

      counter++;
    }
  } 
}

void print_array(FILE *fp,double *ary,int len,int doLogTransform){
  //  fprintf(stderr,"Printing in logspace\n");
  if(doLogTransform){
    for (int i=0;i<len-1;i++)
      fprintf(fp,"%f\t",log(ary[i]));
    fprintf(fp,"%f\n",log(ary[len-1]));
  }else{
    for (int i=0;i<len-1;i++)
      fprintf(fp,"%f\t",(ary[i]));
    fprintf(fp,"%f\n",(ary[len-1]));
  }
}

void abcSaf::run(funkyPars  *p){
  if(p->numSites==0||(doSaf==0 ))
    return;
  //  fprintf(stderr,"-doSaf:%d -p->nind:%d\n",doSaf,p->nInd);
  if(doSaf>0&&doSaf!=3) {
    realRes *r = new realRes;
    r->oklist=new char[p->numSites];
    memset(r->oklist,0,p->numSites);
    r->pLikes=new float*[p->numSites];
    r->pBound=new int*[p->numSites];



	  if(isSim){
		for(int s=0;s<p->numSites;s++){
			if(p->keepSites[s]==0)
				continue;
			int efSize=0;
			for(int i=0;i<p->nInd;i++){
				for(int ii=1;ii<10;ii++){
					if(p->likes[s][i*10+ii]!=p->likes[s][i*10+0]){
						efSize++;
						break;
					}
				}
			}
		p->keepSites[s] = efSize;
		if(minInd!=0&&minInd>efSize){
			p->keepSites[s] = 0;
		fprintf(stderr,"\n\n%d\n\n",minInd);
		}
		}
	  }


    
    if(doSaf==1&&isHap==0)
      algoJoint(p->likes,p->anc,p->numSites,p->nInd,p->keepSites,r,noTrans);
    else if(doSaf==1&&isHap==1)
      algoJointHap(p->likes,p->anc,p->numSites,p->nInd,p->keepSites,r,noTrans);
    else if(doSaf==2){
      freqStruct *freq = (freqStruct *) p->extras[7];
      filipe::algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,p->keepSites,r,noTrans,doSaf,p->major,p->minor,freq->freq,filipeIndF,newDim);
    }else if(doSaf==4){
      algoJointPost(p->post,p->numSites,p->nInd,p->keepSites,r);
    }else if(doSaf==5&&isHap==0){
      algoJointMajorMinor(p->likes,p->numSites,p->nInd,p->keepSites,r,p->major,p->minor);
    }else if(doSaf==5&&isHap==1){
      algoJointMajorMinorHap(p->likes,p->numSites,p->nInd,p->keepSites,r,p->major,p->minor);
    }

    p->extras[index] = r;
  }
  if(doSaf==3){
    kstring_t *kbuf = new kstring_t;
    kbuf->s=NULL;kbuf->l=kbuf->m=0;
    p->extras[index] = kbuf;
    algoGeno(p->refId,p->likes,p->major,p->minor,p->numSites,p->nInd,kbuf,underFlowProtect,p->posi,p->keepSites,prior);
  }
}

void abcSaf::clean(funkyPars *p){
  if(p->numSites==0||(doSaf==0 ))
    return;
  if(doSaf==3)
    return;

  realRes *r=(realRes *) p->extras[index];
  
  int id=0;
  for(int i=0;i<p->numSites;i++)
    if(r->oklist[i]==1)
    {
      delete [] r->pLikes[id];
      delete [] r->pBound[id];
      id++;
    }

  delete [] r->pLikes;
  delete [] r->pBound;
  delete [] r->oklist;
  delete r;
}

//return value is now the number of sites used
//sumBand is now the sum of bins with data
int printFull(funkyPars *p, 
               int index, 
               BGZF *outfileSFS, 
               BGZF *outfileSFSPOS, 
               char *chr, 
               size_t &sumBand)
{
  realRes *r = (realRes *) p->extras[index];
  int counter = 0;
  for(int s=0; s<p->numSites; s++){
    if(r->oklist[s]==3)
      fprintf(stderr,"\t-> Problem with banding algorithm at site chr: \'%s\' position: %d\n",abc::header->target_name[p->refId],p->posi[s]+1);
    if(r->oklist[s]==1 && p->keepSites[s]){
      aio::bgzf_write(outfileSFS, r->pBound[counter], sizeof(int)*2);
      aio::bgzf_write(outfileSFS, r->pLikes[counter], sizeof(float)*r->pBound[counter][1]);
      sumBand += r->pBound[counter][1];
      counter++;
    }
  }

  for(int i=0; i<p->numSites; i++){
    int mypos = p->posi[i];
    if(r->oklist[i]==1&&p->keepSites[i])
      aio::bgzf_write(outfileSFSPOS, &mypos, sizeof(int));
    else if (r->oklist[i]==2){
      if(p->major!=NULL && p->minor!=NULL)
	fprintf(stderr,"PROBS at: %s\t%d anc: %d major: %d minor: %d\n",chr,p->posi[i]+1,p->anc[i],p->major[i],p->minor[i]);
      else
	fprintf(stderr,"PROBS at: %s\t%d\n",chr,p->posi[i]+1);

    }
  }
  return counter;
}

void abcSaf::print(funkyPars *p)
{
  if(p->numSites==0||(doSaf==0))
    return;
  
  if(doSaf==3)
  {
    kstring_t *buf = (kstring_t *) p->extras[index];
    aio::bgzf_write(outfileGprobs, buf->s, buf->l);
    free(buf->s);
    delete buf;
  } 
  else 
  {
    // add prior if this has been supplied
    if(prior != NULL) {

      realRes *r = (realRes *) p->extras[index];
      int counter = 0;

      for(int s=0; s<p->numSites; s++) 
      {
        float *workarray = NULL;
        int lower, nbin;

        if(r->oklist[s]==1)
        {
          workarray = r->pLikes[counter];
          lower = r->pBound[counter][0];
          nbin = r->pBound[counter][1];
          counter++;
        }
        else
          continue;
      
        double tsum = 0.;
        int k;

        k = 0;
        for(int i=lower; i<nbin; i++)
          tsum += exp(workarray[k++] + prior[i]);
        tsum = log(tsum);
      
        k = 0;
        for(int i=lower; i<nbin; i++)
          workarray[k++] += prior[i] - tsum;
    
        //pLikes now contains our posterior expectation of the different classes of frequencies
      }
    }
 
    nnnSites += printFull(p,index,outfileSAF,outfileSAFPOS,header->target_name[p->refId],sumBand);
  }   
}

// --- misc stuff --- //

// NSP 3July2020
// I left this alone when implementing score-limited algorithm. Not sure if it is even used anymore?
void abcSaf::algoGeno(int refId,double **liks,char *major,char *minor,int nsites,int numInds,kstring_t *sfsfile,int underFlowProtect,int *posi,int *keepSites,double *pest) {
  assert(pest!=NULL);
  //void algGeno(aMap &asso,int numInds,FILE *sfsfile,int underFlowProtect, double *pest) {
  //  fprintf(stderr,"UNDERFLOWPROTECT: %d\n",underFlowProtect);
  double *postp = new double[3*numInds];
#if 0
  for(int r=0;r<(2*numInds-1);r++)
    for(int j=0;j<std::min(3,r);j++){
      double res =myComb2Tab[r][j];// myComb2(numInds,r,j);
      fprintf(stderr,"myComb\t(%d,%d,%d) =%f\n",numInds,r,j,res);
    }
  exit(0);  
#endif
  
  //algorithm goes on by a site on site basis //pretty much the same as 'algo' without the prior phat
  
  //  underFlowProtect = 0;
  for(int it=0; it<nsites; it++) {//loop over sites

    int minor_offset=refToInt[minor[it]];
    int major_offset=refToInt[major[it]];
    //    fprintf(stderr,"%c %c asInt %d %d\n",major[it],minor[it],minor_offset,major_offset);
    if(keepSites[it]==0)
      continue;


    //    fprintf(sfsfile,"%s %d\t",chrnames[it->first.chromo],it->first.position);
    //1. first loop through all possible major/minor

    //    fprintf(sfsfile,"%d %d\t0 0 0 0\t",major_offset,minor_offset);
    int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
    int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
    int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
    //    fprintf(stderr,"OFFSET_INFO: min=%d\tmaj=%d aa=%d Aa=%d AA=%d\n",minor_offset,major_offset,aa_offset,Aa_offset,AA_offset);
    //part two
    double hj[2*numInds+1];
    for(int index=0;index<(2*numInds+1);index++)
      if(underFlowProtect==0)
	hj[index]=0;
      else
	hj[index]=log(0);
    double PAA,PAa,Paa;
    //    numInds =5;
    for(int i=0 ; i<numInds ;i++) {
      //	printf("AA=%f\tAa=%f\taa=%f\n",p.lk[i*3+AA_offset],p.lk[i*3+Aa_offset],p.lk[i*3+aa_offset]);
      double GAA,GAa,Gaa;
      GAA = liks[it][i*10+AA_offset];
      GAa = log(2.0)+liks[it][i*10+Aa_offset];
      Gaa = liks[it][i*10+aa_offset];

      //fprintf(stderr,"GAA=%f\tGAa=%f\tGaa=%f\n",GAA,GAa,Gaa);
      if(underFlowProtect==0){
	GAA=exp(GAA);
	GAa=exp(GAa);
	Gaa=exp(Gaa);
      }
	  
      
      PAA =(GAA);///(MAA+MAa+Maa);
      PAa =(GAa);///(MAA+MAa+Maa);
      Paa =(Gaa);///(MAA+MAa+Maa);
      
      
      //check for underflow error, this should only occur once in a blue moon
      if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	printf("PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
      }
      
      if(i==0){
	hj[0] =Paa;
	hj[1] =PAa;
	hj[2] =PAA;
      }else{
	
	for(int j=2*(i+1); j>1;j--){

	  double tmp;
	  if(underFlowProtect==1)
	    tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	  else
	    tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	      
	  if(std::isnan(tmp)){
	    fprintf(stderr,"jis nan:%d\n",j );
	    hj[j] = 0;
	    break;
	  }else
	    hj[j]  =tmp;

	}
	if(underFlowProtect==1){
	  hj[1] = angsd::addProtect2(Paa+hj[1],PAa+hj[0]);
	  hj[0] = Paa+hj[0];
	}
	else{
	  hj[1] = Paa*hj[1] + PAa*hj[0];
	  hj[0] = Paa*hj[0];
	}
	//	print_array(stdout,hj,2*numInds+1,!underFlowProtect);
      }
      
      
    }//after recursion
    //if we are underflowprotecting ht hj is in logspace
    // print_array(stdout,hj,2*numInds+1,!underFlowProtect);
    for(int i=0;i<(2*numInds+1);i++){
      //fprintf(stdout,"BICO: %f\n",log(bico(2*numInds,i)));
      if(underFlowProtect)
	hj[i] =  (hj[i]-lbicoTab[i]);
      else
	hj[i] =  exp(log(hj[i])-lbicoTab[i]);
      
    }
    //fprintf(stdout,"the full after update hj\n");
    //  print_array(stdout,hj,2*numInds+1,underFlowProtect);
    double denominator = 0;
    if(underFlowProtect)
      denominator = log(denominator);
    for(int i=0;i<=2*numInds;i++)
      if(underFlowProtect)
	denominator = angsd::addProtect2(denominator,pest[i]+angsd::addProtect2(hj[i],hj[2*numInds-i]));
      else
	denominator += exp(pest[i]) * (hj[i] +hj[2*numInds-i]);
    //    fprintf(stderr,"denom:%e\n",denominator);

    int whichGeno[numInds];
    double whichProb[numInds];

    for(int select=0;select<numInds;select++) {
      //      fprintf(stderr,"select:%d\n",select);
      double *hj = new double[2*numInds-1];
      for(int index=0;index<(2*numInds-1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,PAa,Paa;
      
      int ishadow =-1;
      for(int i=0 ; i<numInds ;i++) {
	//	printf("AA=%f\tAa=%f\taa=%f\n",p.lk[i*3+AA_offset],p.lk[i*3+Aa_offset],p.lk[i*3+aa_offset]);
	if(i!=select){
	  ishadow++;
	}else
	  continue;
	double GAA,GAa,Gaa;
	GAA = liks[it][i*10+AA_offset];
	GAa = log(2.0)+liks[it][i*10+Aa_offset];
	Gaa = liks[it][i*10+aa_offset];
	
	if(underFlowProtect==0){
	  GAA=exp(GAA);
	  GAa=exp(GAa);
	  Gaa=exp(Gaa);
	}
	  
	
	PAA =(GAA);///(MAA+MAa+Maa);
	PAa =(GAa);///(MAA+MAa+Maa);
	Paa =(Gaa);///(MAA+MAa+Maa);


	//check for underflow error, this should only occur once in a blue moon
	if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	  fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
	}
	
	if(ishadow==0){
	  hj[0] =Paa;
	  hj[1] =PAa;
	  hj[2] =PAA;
	}else{
	  
	  for(int j=2*(ishadow+1); j>1;j--){
	    //	    fprintf(stderr,"j=%d\n",j);
	    //print_array(hj,2*numInds+1);
	    double tmp;
	    if(underFlowProtect==1)
	      tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	    else
	      tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	    
	    if(std::isnan(tmp)){
	      fprintf(stderr,"jis nan:%d\n",j );

	      hj[j] = 0;
	      break;
	    }else
	      hj[j]  =tmp;
	  }
	  if(underFlowProtect==1){
	    hj[1] = angsd::addProtect2(Paa+hj[1],PAa+hj[0]);
	    hj[0] = Paa+hj[0];
	  }
	  else{
	    hj[1] = Paa*hj[1] + PAa*hj[0];
	    hj[0] = Paa*hj[0];
	  }
	}
	//ifunderflowprotect then hj is in logspace
	
      }//after recursion
      for(int i=0;i<(2*(numInds-1)+1);i++){
	//fprintf(stdout,"BICO: %f\n",log(bico(2*numInds,i)));
	if(underFlowProtect)
	  hj[i] =  (hj[i]-angsd::lbico(2*(numInds-1),i));
	else
	  hj[i] =  exp(log(hj[i])-angsd::lbico(2*(numInds-1),i));
	
      }
       
      //now do all the genocalling for individual =select
      // fprintf(stderr,"seelct=%d\tmyMaj=%d\tmyMin=%d\tAA_offset=%d Aa_offset=%d aa_offset=%d\n",select,p.major,p.minor,AA_offset,Aa_offset,aa_offset);
      double *asdf = hj;

      //print_array(stdout,asdf,2*numInds-1,!underFlowProtect);
      double *res = postp+3*select; //is always logged


      for (int j=0;j<3;j++) {//loop through 3 genotypes
	//	fprintf(stderr,"jis: %d\n",j);
	double g;
	if(j==0)
	  g=liks[it][10*select+angsd::majorminor[major[it]][major[it]]];
	else if(j==1)
	  g=liks[it][10*select+angsd::majorminor[major[it]][minor[it]]];
	else
	  g=liks[it][10*select+angsd::majorminor[minor[it]][minor[it]]];
	
	double tmp=0;
	if(underFlowProtect)
	  tmp = log(tmp);
	for(int r=0;r<2*numInds-1;r++){
	  //	  fprintf(stderr,"mycomb2:%f\tpes1: %f\tpes2: %f\n",myComb2(numInds,r+j,j),pest[r+j],pest[2*numInds-r-j]);
	  if(underFlowProtect)
	    //tmp = angsd::addProtect2(tmp, log(myComb2(numInds,r+j,j)) + (asdf[r])+pest[r+j]+pest[2*numInds-r-j]);
	    tmp = angsd::addProtect2(tmp, log(myComb2Tab[r+j][j]) + (asdf[r])+pest[r+j]+pest[2*numInds-r-j]);
	  else{
	    //fprintf(stderr,"r+j:%d\tj:%d\nold:%f\nnew:%f",r+j,j,myComb2(numInds,r+j,j), myComb2Tab[r+j][j]);
	    //exit(0);
	    //tmp += myComb2(numInds,r+j,j)*(asdf[r])*(exp(pest[r+j])+exp(pest[2*numInds-r-j]));
	    tmp += myComb2Tab[r+j][j]*(asdf[r])*(exp(pest[r+j])+exp(pest[2*numInds-r-j]));
	  }
	}
	//g is directly from glf file, tmp is in log depending on underflowprotect

	if(underFlowProtect)
	  res[j]=g+(tmp);
	else
	  res[j]=g+log(tmp);
	//	fprintf(stderr,"j:%d log(g)=%f tmp=%f g*tmp=%f res[%d]=%f logres[%d]=%f \n",j,g,log(tmp),(g*tmp),j,res[j],j,log(res[j]));

	//res is always in log
	
	

	//CODE BELOW IS A CHECK
#if 0
	double tmp1=0;
	double tmp2=0;
	if(underFlowProtect){
	  tmp1 = log(tmp1);
	  tmp2 = log(tmp2);
	}
	  
	for(int r=j;r<=2*numInds-2+j;r++) {
	  // fprintf(stderr,"tmp1: r=%d index: %d\n",r,r-j);
	    //fprintf(stderr,"logcomb: %f\n",myComb2(numInds,r,j));
	    //    fprintf(stderr,"tmp1 in preloop: %f\n",tmp1);
	    if(underFlowProtect)
	      tmp1 = angsd::addProtect2(tmp1, log(myComb2Tab[r][j])+ asdf[r-j]+pest[r]);
	    else
	      tmp1 += myComb2Tab[r][j]*(asdf[r-j])*exp(pest[r]);
	    //fprintf(stderr,"tmp1 in reloop: %f\n",tmp1);
	}

	for(int r=2-j;r<=2*numInds-j;r++){
	  //fprintf(stderr,"tmp2: r=%d index=%d\n",r,r-2+j);
	  //fprintf(stderr,"logcomb: %f\n",myComb2(numInds,r,2-j));
	  if(underFlowProtect)
	    tmp2 = angsd::addProtect2(tmp2,log(myComb2Tab[r][2-j])+ asdf[(2*(numInds-1))-(r-2+j)]+ pest[r]);
	  else
	    tmp2 += myComb2Tab[r][2-j]*(asdf[2*(numInds-1)-(r-2+j)])*exp(pest[r]);
	}
	//	fprintf(stderr,"tmp1: %f tmp2: %f\n",(tmp1),(tmp2));	
	double res2;//always log

	if(underFlowProtect)
	  res2 = g+(angsd::addProtect2(tmp1,tmp2));
	else
	  res2 = g+log(tmp1+tmp2);
	//fprintf(stdout,"%f vs %f\n",res[j],res2);
	if(1&&!isSame(res[j],res2,0.000001)){
	  fprintf(stdout,"%f vs %f\n",res[j],res2);
	}
	//check that is sums to one according the model
	 
	if(underFlowProtect)
	  fprintf(stdout,"\tASDFASDFASDF: %f\n",log((res2-denominator)));
	else
	  fprintf(stdout,"\tASDFASDFASDF: %f\n",log(exp(res2)/denominator));	

	//CODE ABOVE IS CHECK

#endif

      }//after the loop of the tree genoypes
      
      
      double shouldBeOne =0;
      if(underFlowProtect)
	shouldBeOne = exp(res[0]-(denominator))+exp(res[1]-(denominator))+exp(res[2]-(denominator));
      else
	shouldBeOne = exp(res[0]-log(denominator))+exp(res[1]-log(denominator))+exp(res[2]-log(denominator));
      //fprintf(stderr,"this should be one: %f\n",shouldBeOne);
      if((fabs(shouldBeOne-1)>0.000001)){
	fprintf(stderr,"this should be one: %f at position:(%s,%d)\n",shouldBeOne,header->target_name[refId],posi[it+1]);
	exit(0);
      }
      


      //      print_array(stdout,res,3,0);      
      
      //logrescale(res,3);
      double mySum=exp(res[0])+exp(res[1])+exp(res[2]);
      for(int i=0;i<3;i++)
	res[i] =exp(res[i])/mySum;
      int best = angsd::whichMax(res,3);
      assert(best!=-1);
      //      fprintf(stderr,"best:%d\n",best);
      whichGeno[select] = best;
      whichProb[select] = res[best];
      continue;
      if(best==0)
	whichGeno[select] = angsd::majorminor[major[it]][major[it]];
      if(best==1)
	whichGeno[select] = angsd::majorminor[major[it]][minor[it]];
      if(best==2)
	whichGeno[select] = angsd::majorminor[minor[it]][minor[it]];
      
      //print_array(sfsfile,res,3,0);      

    }//after select loop  

    //    print_array(sfsfile,whichGeno,numInds);
    // print_array(sfsfile,whichProb,numInds,0);

    ksprintf(sfsfile,"%s\t%d\t%c\t%c\t",header->target_name[refId],posi[it]+1,intToRef[major[it]],intToRef[minor[it]]);
    for(int i=0;i<numInds;i++)
      ksprintf(sfsfile,"%d ",whichGeno[i]);
    //ksprintf(sfsfile,"%d\t",whichGeno[numInds]);
    for(int i=0;i<numInds*3;i++)
      ksprintf(sfsfile,"%f ",postp[i]);
    ksprintf(sfsfile,"\n");
    //    fprintf(sfsfile,"%f\t%f\n",whichProb[numInds-1],it->second.emPhat);
  }//after all sites
  
}

void abcSaf::writeAll(){
  assert(outfileSAF!=NULL);
  assert(outfileSAFIDX!=NULL);
  assert(outfileSAFPOS!=NULL);
  //  fprintf(stderr,"nnnSites:%d\n",nnnSites);
  if(nnnSites!=0&&tmpChr!=NULL){
    size_t clen = strlen(tmpChr);
    fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
    fwrite(tmpChr,1,clen,outfileSAFIDX);
    size_t tt = nnnSites;
    fwrite(&tt,sizeof(size_t),1,outfileSAFIDX);
    fwrite(&sumBand,sizeof(size_t),1,outfileSAFIDX);
    fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
  }//else
   // fprintf(stderr,"enpty chr\n");
  //reset
  assert(bgzf_flush(outfileSAF)==0);assert(bgzf_flush(outfileSAFPOS)==0);
  offs[0] = bgzf_tell(outfileSAFPOS);
  offs[1] = bgzf_tell(outfileSAF);
  nnnSites=0;
  sumBand = 0;
}

void abcSaf::changeChr(int refId) {
  if(doSaf&&doSaf!=3)
    writeAll();

  free(tmpChr);
  tmpChr = strdup(header->target_name[refId]);
}
