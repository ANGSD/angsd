#include <cstdio>
#include <cmath>
#include <assert.h>
#include <zlib.h>
#include <cfloat>
#include <htslib/kstring.h>
#include "abcFreq.h"
#include "shared.h"
#include "analysisFunction.h"

#include "abc.h"

#include "abcSaf.h"

#define MINLIKE -1000.0 //this is for setting genotypelikelhoods to missing (EXPLAINED BELOW)

/*
  From 0.501 the realSFS method can take into account individual inbreeding coefficients
  Filipe guera
*/
namespace filipe{
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans,int doSaf,char *major,char *minor,double *freq,double *indF,int newDim);
}


void abcSaf::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSaf\t\t%d\n",doSaf);
  fprintf(argFile,"\t1: perform multisample GL estimation\n\t2: use an inbreeding version\n\t3: calculate genotype probabilities\n\t4: Assume genotype posteriors as input (still beta) \n");
  fprintf(argFile,"\t-doThetas\t\t%d (calculate thetas)\n",noTrans);
  fprintf(argFile,"\t-underFlowProtect\t%d\n",underFlowProtect); 
  fprintf(argFile,"\t-fold\t\t\t%d (deprecated)\n",fold);
  fprintf(argFile,"\t-anc\t\t\t%s (ancestral fasta)\n",anc);
  fprintf(argFile,"\t-noTrans\t\t%d (remove transitions)\n",noTrans);
  fprintf(argFile,"\t-pest\t\t\t%s (prior SFS)\n",pest);
  fprintf(argFile,"\t-isHap\t\t\t%d (is haploid beta!)\n",isHap);
  fprintf(argFile,"NB:\n\t  If -pest is supplied in addition to -doSaf then the output will then be posterior probability of the sample allelefrequency for each site\n");

}

double lbico(double n, double k){
  return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1);
}

double myComb2(int k,int r, int j){
  //fprintf(stderr,"myComb\t%d\t%d\t%d\n",k,r,j);
  //return 1.0;
  //  return 1.0/bico(2*k,2);
  if(j>r)
    fprintf(stderr,"%s error in k=%d r=%d j=%d\n",__FUNCTION__,k,r,j);
      //    return 0;
  double fac1= lbico(r,j)+lbico(2*k-r,2-j);
  double fac2=lbico(2*k,2);
  //  fprintf(stderr,"[%s]=%f\t=%f\tres=%f\n",__FUNCTION__,fac1,fac2,fac1/fac2);
  return exp(fac1-fac2);
  //return fac1;
}

void abcSaf::getOptions(argStruct *arguments){
  doSaf=angsd::getArg("-doSaf",doSaf,arguments);
  isHap=angsd::getArg("-isHap",isHap,arguments);

  int tmp=0;
  tmp=angsd::getArg("-doRealSFS",tmp,arguments);
  if(tmp){
    fprintf(stderr,"-doRealSFS is now called -doSaf (do sample allele frequency)\n");
    doSaf=tmp;
  }

  noTrans = angsd::getArg("-noTrans",noTrans,arguments);
  pest = angsd::getArg("-pest",pest,arguments);
  int GL = 0;
  GL = angsd::getArg("-GL",GL,arguments);
  doThetas= angsd::getArg("-doThetas",doThetas,arguments);


  if(doSaf==0&&doThetas==0)
    return;
  if(doThetas!=0&& doSaf==0){
    fprintf(stderr,"\t Must estimate joint llh (-doSaf) when estimating thetas\n");
    exit(0);
  }


  if(pest==NULL&& doThetas){
    fprintf(stderr,"\t Must supply a sfs used as prior for estimating thetas\n");
    exit(0);
  }

  underFlowProtect=angsd::getArg("-underFlowProtect",underFlowProtect,arguments);
  fold=angsd::getArg("-fold",fold,arguments);

  int isSim =0;
  isSim=angsd::getArg("-isSim",isSim,arguments);

  if(doSaf==0)
    return;
  anc = angsd::getArg("-anc",anc,arguments);
  if(doSaf && (anc==NULL&&isSim==0)){
    if(doSaf!=3){
      fprintf(stderr,"Must supply -anc for polarizing the spectrum\n");
      exit(0);
    }
  }
  if(pest!=NULL){
    if(fold==0)
      prior=angsd::readDouble(pest,arguments->nInd*2+1);
    else
      prior=angsd::readDouble(pest,arguments->nInd+1);
  }
  if(GL==0 &&(arguments->inputtype!=INPUT_GLF&&arguments->inputtype!=INPUT_BEAGLE&&arguments->inputtype!=INPUT_VCF_GL)){
    fprintf(stderr,"Must supply genotype likelihoods (-GL [INT])\n");
    printArg(arguments->argumentFile);
    exit(0);
  }
  if(doSaf==2){
    fprintf(stderr,"\t->(Using Filipe G Vieira modification of: %s)\n",__FILE__);
    int doMajorMinor =0;
    doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
    int doMaf =0;
    doMaf = angsd::getArg("-doMaf",doMaf,arguments);
    if(doMajorMinor==0||doMaf==0){
      fprintf(stderr,"Must have major/minor and MAF for using the inbreeding version of realSFS\n");
      exit(0);
    }
    char *indF_name = NULL;
    indF_name =  angsd::getArg("-indF",indF_name,arguments);
    if(indF_name==NULL){
      filipeIndF = new double[arguments->nInd];
      for(int i=0;i<arguments->nInd;i++)
	filipeIndF[i] =0;
      fprintf(stderr,"\t-> No -indF file provided will assume an inbreeding zero for all samples.\n");
      fprintf(stderr,"\t-> If no inbreeding is expected consider using -doSaf 1\n");
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
      fprintf(stderr,"\t-> You need to supply a sfs as prior (emOptim2 output, -pest) to do genotypecalling with saf\n");
      exit(0);
    }
  }

}

abcSaf::abcSaf(const char *outfiles,argStruct *arguments,int inputtype){
  isHap =0;
  lbicoTab = NULL;
  myComb2Tab=NULL;
  const char *SAF = ".saf";
  const char *SFSPOS =".saf.pos.gz";
  const char *THETAS =".thetas.gz";
  const char *GENO = ".saf.geno.gz";
  //default
  underFlowProtect = 0;
  fold =0;
  isSim =0;
  //from command line
  anc=NULL;
  pest=NULL;
  noTrans = 0;
  prior = NULL;
  doSaf=0;
  doThetas = 0;
  outfileSFS = NULL;
  outfileSFSPOS = Z_NULL;
  theta_fp = Z_NULL;
  outfileGprobs = Z_NULL;
  scalings = NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSaf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  printArg(arguments->argumentFile);  
  if(doSaf==0){
    shouldRun[index] =0;
    return;
  }
  lbicoTab = new double[2*arguments->nInd+1];
  myComb2Tab = new double*[2*arguments->nInd+1];
  for(int i=0;i<2*arguments->nInd+1;i++){
    lbicoTab[i] = lbico(2*arguments->nInd,i);
    myComb2Tab[i] = new double[3];
    for(int j=0;j<3;j++)
      if(j<=i)
	myComb2Tab[i][j] = myComb2(arguments->nInd,i,j);
  }
  if(isHap){
    for(int i=0;i<arguments->nInd+1;i++)
      lbicoTab[i] = lbico(arguments->nInd,i);
  }

  newDim = 2*arguments->nInd+1;
  if(fold)
    newDim = arguments->nInd+1;
  if(doSaf==3){
    outfileGprobs = aio::openFileGz(outfiles,GENO,"w6h");
  }else if(doSaf!=0&&doThetas==0){
    outfileSFS =  aio::openFile(outfiles,SAF);
    outfileSFSPOS =  aio::openFileGz(outfiles,SFSPOS,GZOPT);
  }else {
    theta_fp = aio::openFileGz(outfiles,THETAS,"w6h");;
    gzprintf(theta_fp,"#Chromo\tPos\tWatterson\tPairwise\tthetaSingleton\tthetaH\tthetaL\n");
    aConst=0;
    int nChr = 2*arguments->nInd;
    for(int i=1;i<nChr;i++)
      aConst += 1.0/i;
    aConst = log(aConst);//this is a1
    
    
    aConst2 = log((nChr*(nChr-1))/2.0);//choose(nChr,2)
    aConst3 = log((1.0*nChr-1.0));
    
    scalings = new double [nChr+1];
    for(int i=0;i<nChr+1;i++)
      scalings[i] = log(i)+log(nChr-i);
    
    
    if(fold){
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	scalings[i] = log(exp(scalings[i]) + exp(scalings[2*arguments->nInd-i]))-log(2.0);//THORFINN NEW
      
    }
  }
#if 0//just for printing the scalings
    int lim = nChr;
    if(fold)
      lim = newDim;
    for(int i=0;i<newDim;i++)
      fprintf(stderr,"SCAL[%d]\t%f\n",i,scalings[i]);

#endif

  
}


abcSaf::~abcSaf(){
  if(pest) free(pest);
  if(prior) delete [] prior;
  if(outfileSFS) fclose(outfileSFS);;
  if(outfileSFSPOS) gzclose(outfileSFSPOS);
  if(theta_fp) gzclose(theta_fp);
  if(outfileGprobs)  gzclose(outfileGprobs);
  if(lbicoTab) delete [] lbicoTab;
  if(myComb2Tab) delete [] myComb2Tab;//not cleaned properly
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

//18dec 2014 we now condition on the sites where the ancestral is either major or minor. 
void filipe::algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans,int doSaf,char *major,char *minor,double *freq,double *indF,int newDim) {
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

    if(fold) {
      //newDim is set in constructor
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	sumMinors[i] = log(sumMinors[i] + sumMinors[2*numInds-i]);
      sumMinors[newDim-1] = log(sumMinors[newDim-1])+log(2.0);
      //angsd::logrescale(sumMinors,newDim);
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
	r->pLikes[myCounter] =new double[numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(numInds+1));
	myCounter++;
      }
    }else{
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
	r->pLikes[myCounter] =new double[2*numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
	myCounter++;
      }
    }
  }
  
}

void abcSaf::algoJointPost(double **post,int nSites,int nInd,int *keepSites,realRes *r,int doFold){
  // fprintf(stderr,"[%s]\n",__FUNCTION__);
  int myCounter =0;
  for(int s=0;s<nSites;s++){
    if(keepSites[s]==0)
      continue;
    double *liks=post[s]; //we call this liks, eventhough it is posteriors.
    double *hj = new double [2*nInd+1];
    for(int index=0;index<(2*nInd+1);index++)
      hj[index]=0;
    //initalize
    memcpy(hj,liks,3*sizeof(double));
    
    for(int i=1 ; i<nInd ;i++) {
      double Paa=liks[i*3];
      double PAa=2*liks[i*3+1];
      double PAA=liks[i*3+2];
      for(int j=2*(i+1); j>1;j--)
	hj[j] = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
      hj[1] = Paa*hj[1] + PAa*hj[0];
      hj[0] = Paa*hj[0];
      
      normalize_array(hj,2*(i+1));
    }
    for(int i=0;i<(2*nInd+1);i++)
      hj[i] =  log(hj[i])-lbicoTab[i];
       
    if(doFold){
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	hj[i] = log(exp(hj[i]) + exp(hj[2*nInd-i]));
      hj[newDim-1] = hj[newDim-1]+log(2.0);
    }
    angsd::logrescale(hj,2*nInd+1);
    if(std::isnan(hj[0]))
      r->oklist[s] = 2;
    else{
      r->oklist[s] = 1;
      r->pLikes[myCounter] =hj;
      myCounter++;
    }
  }
}

//AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
int homo[4] = {0,4,7,9};


void abcSaf::algoJointHap(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans) {
  int myCounter =0;
  if(anc==NULL||liks==NULL){
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }
  double sumMinors[numInds+1];  //the sum of the 3 different minors

  for(int it=0; it<nsites; it++) {//loop over sites
    int major_offset = anc[it];
    if(major_offset==4||(keepSites[it]==0)){//skip of no ancestral information
      keepSites[it] =0; //
      //      r->oklist is zero no need to update
      continue;
    }
    //set the resultarray to zeros
    for(int sm=0 ; sm<numInds+1 ; sm++ )
      sumMinors[sm] = 0;
    
    //loop through the 3 different minors
    for(int minor_offset=0;minor_offset<4;minor_offset++) {
      if(minor_offset == major_offset)
	continue;
      if(noTrans){
	if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
	  continue;
	if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
	  continue;
      }
      double totmax = 0.0;
      //hook for only calculating one minor
      //  int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
      int AA_offset = homo[major_offset];//angsd::majorminor[minor_offset][minor_offset];//0-9
      int aa_offset = homo[minor_offset];//angsd::majorminor[major_offset][major_offset];//0-9
      //     fprintf(stderr,"%d:%d\t%d\t%d\n",major_offset,Aa_offset,AA_offset,aa_offset);
      //part two
      double hj[numInds+1];
      for(int index=0;index<(2*numInds+1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,Paa;

      for(int i=0 ; i<numInds ;i++) {
	double GAA,Gaa;

	GAA = liks[it][i*10+AA_offset];
	Gaa = liks[it][i*10+aa_offset];

	//do underlfow protection (we are in logspace here) (rasmus style)
	if(1){
	  double mymax=std::max(Gaa,GAA);
	  // fprintf(stdout,"mymax[%d]=%f\t",i,mymax);
	  
	  if(mymax<MINLIKE){
	    //	    fprintf(stderr,"\n%f %f %f\n",GAA, GAa,Gaa);
	    Gaa = 0;
	    GAA = 0;
	    totmax = totmax + mymax;
	  }else{
	    Gaa=Gaa-mymax;
	    GAA=GAA-mymax;
	    totmax = totmax + mymax;
	  }
	
	//	fprintf(stderr,"totmax=%f\n",totmax);
	//END underlfow protection (we are in logspace here) (rasmus style)
	}
	if(underFlowProtect==0){
	  PAA=exp(GAA);
	  Paa=exp(Gaa);
	}else{
	  PAA =(GAA);///(MAA+MAa+Maa);
	  Paa =(Gaa);///(MAA+MAa+Maa);
	}

	//check for underflow error, this should only occur once in a blue moon
	if(std::isnan(Paa)||std::isnan(PAA)){
	  fprintf(stderr,"PAA=%f\tPaa=%f\n",PAA,Paa);
	}
	//	fprintf(stdout,"it=%d PAA=%f\tPAa=%f\tPaa=%f\n",it,PAA,PAa,Paa);
	if(i==0){
	  hj[0] =Paa;
	  hj[1] =PAA;
	}else{
	  //fprintf(stderr,"asdf\n");
	  for(int j=i+1; j>0;j--) {
	    double tmp;
	    if(underFlowProtect==1)
	      tmp =angsd::addProtect2(PAA+hj[j-1],Paa+hj[j]);
	    else
	      tmp = PAA*hj[j-1]+Paa*hj[j];
	    
	    if(std::isnan(tmp)){
	      fprintf(stderr,"is nan:%d\n",j );
	      
	      hj[j] = 0;
	      break;
	    }else
	      hj[j]  =tmp;
	  }
	  if(underFlowProtect==1)
	    hj[0] = Paa+hj[0];
	  else
	    hj[0] = Paa*hj[0];

	}
	//ifunderflowprotect then hj is in logspace
	
      }
      
      for(int i=0;i<numInds+1;i++)
	if(underFlowProtect==0)
	  sumMinors[i] +=  exp(log(hj[i])-lbicoTab[i]+totmax);
	else
	  sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-lbicoTab[i]+totmax));
    }
    //sumMinors is in normal space, not log
    /*
      we do 3 things.
      1. log scaling everyting
      2. rescaling to the most likely in order to avoid underflows in the optimization
      3. we might do a fold also.
      
     */    
    {
      for(int i=0;i<numInds+1;i++)
	sumMinors[i] = log(sumMinors[i]);
      angsd::logrescale(sumMinors,numInds+1);
      if(std::isnan(sumMinors[0]))
	r->oklist[it] = 2;
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(numInds+1));
	myCounter++;
      }
    }
  }
  
}



void abcSaf::algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans) {
  // fprintf(stderr,"[%s]\n",__FUNCTION__);
  int myCounter =0;
  if(anc==NULL||liks==NULL){
    fprintf(stderr,"problems receiving data in [%s] will exit (likes=%p||ancestral=%p)\n",__FUNCTION__,liks,anc);
    exit(0);
  }
  double sumMinors[2*numInds+1];  //the sum of the 3 different minors

  for(int it=0; it<nsites; it++) {//loop over sites
    int major_offset = anc[it];
    if(major_offset==4||(keepSites[it]==0)){//skip of no ancestral information
      keepSites[it] =0; //
      //      r->oklist is zero no need to update
      continue;
    }
    //set the resultarray to zeros
    for(int sm=0 ; sm<(2*numInds+1) ; sm++ )
      sumMinors[sm] = 0;
    
    //loop through the 3 different minors
    for(int minor_offset=0;minor_offset<4;minor_offset++) {
      if(minor_offset == major_offset)
	continue;
      if(noTrans){
	if((major_offset==2&&minor_offset==0)||(major_offset==0&&minor_offset==2))
	  continue;
	if((major_offset==1&&minor_offset==3)||(major_offset==3&&minor_offset==1))
	  continue;
      }
      double totmax = 0.0;
      //hook for only calculating one minor
      int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
      int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
      int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
      //     fprintf(stderr,"%d:%d\t%d\t%d\n",major_offset,Aa_offset,AA_offset,aa_offset);
      //part two
      double hj[2*numInds+1];
      for(int index=0;index<(2*numInds+1);index++)
	if(underFlowProtect==0)
	  hj[index]=0;
	else
	  hj[index]=log(0);
      double PAA,PAa,Paa;

      for(int i=0 ; i<numInds ;i++) {
	//	printf("pre scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	double GAA,GAa,Gaa;
#ifdef RESCALE
	if(0){//The rescaling is now done in 'getMajorMinor()'
	double max = liks[it][i*10+0];
	for(int index=1;index<10;index++)
	  if(liks[it][i*10+index]>max)
	    max = liks[it][i*10+index];
	for(int index=0;index<10;index++)
	  liks[it][i*10+index] = liks[it][i*10+index]-max;
	}
#endif
	//printf("post scale AA=%f\tAa=%f\taa=%f\n",liks[it][i*3+AA_offset],liks[it][i*3+Aa_offset],liks[it][i*3+aa_offset]);
	GAA = liks[it][i*10+AA_offset];

	GAa = log(2.0)+liks[it][i*10+Aa_offset];
	Gaa = liks[it][i*10+aa_offset];
	//printf("[GAA] GAA=%f\tGAa=%f\tGaa=%f\n",GAA,GAa,Gaa);
	//do underlfow protection (we are in logspace here) (rasmus style)
	if(1){
	  double mymax;
	  if (Gaa > GAa && Gaa > GAA) mymax = Gaa;
	  else if (GAa > GAA) mymax = GAa;
	  else mymax = GAA;
	  // fprintf(stdout,"mymax[%d]=%f\t",i,mymax);
	  
	  if(mymax<MINLIKE){
	    //	    fprintf(stderr,"\n%f %f %f\n",GAA, GAa,Gaa);
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
	//	fprintf(stderr,"totmax=%f\n",totmax);
	//END underlfow protection (we are in logspace here) (rasmus style)
	}
	if(underFlowProtect==0){
	  PAA=exp(GAA);
	  PAa=exp(GAa);
	  Paa=exp(Gaa);
	}else{
	  PAA =(GAA);///(MAA+MAa+Maa);
	  PAa =(GAa);///(MAA+MAa+Maa);
	  Paa =(Gaa);///(MAA+MAa+Maa);
	}

	//check for underflow error, this should only occur once in a blue moon
	if(std::isnan(Paa)||std::isnan(PAa)||std::isnan(Paa)){
	  fprintf(stderr,"PAA=%f\tPAa=%f\tPaa=%f\n",PAA,PAa,Paa);
	}
	//	fprintf(stdout,"it=%d PAA=%f\tPAa=%f\tPaa=%f\n",it,PAA,PAa,Paa);
	if(i==0){
	  hj[0] =Paa;
	  hj[1] =PAa;
	  hj[2] =PAA;
	}else{
	  //fprintf(stderr,"asdf\n");
	  for(int j=2*(i+1); j>1;j--){
	    //  print_array(stdout,hj,2*numInds+1,0);
	    //print_array(hj,2*numInds+1);
	    double tmp;
	    if(underFlowProtect==1)
	      tmp =angsd::addProtect3(PAA+hj[j-2],PAa+hj[j-1],Paa+hj[j]);
	    else
	      tmp = PAA*hj[j-2]+PAa*hj[j-1]+Paa*hj[j];
	    
	    if(std::isnan(tmp)){
	      fprintf(stderr,"is nan:%d\n",j );
	      
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
	
      }
      //      fprintf(stderr,"%f %f %f\n",hj[0],hj[1],hj[2]);
      //fprintf(stdout,"\nscaledLikes=");
      for(int ii=0;0&&ii<10*numInds;ii++)
	fprintf(stdout,"%f\t",liks[it][ii]);

      //      totmax=0;
      for(int i=0;i<(2*numInds+1);i++)
	if(underFlowProtect==0)
	  sumMinors[i] +=  exp(log(hj[i])-lbicoTab[i]+totmax);
	else
	  sumMinors[i] = exp(angsd::addProtect2(log(sumMinors[i]),hj[i]-lbicoTab[i]+totmax));
    }
    //sumMinors is in normal space, not log
    /*
      we do 3 things.
      1. log scaling everyting
      2. rescaling to the most likely in order to avoid underflows in the optimization
      3. we might do a fold also.
      
     */    

    if(fold) {
      int newDim = numInds+1;
      for(int i=0;i<newDim-1;i++)// we shouldn't touch the last element
	sumMinors[i] = log(sumMinors[i] + sumMinors[2*numInds-i]);//THORFINN NEW
      sumMinors[newDim-1] = log(sumMinors[newDim-1])+log(2.0);
      angsd::logrescale(sumMinors,newDim);
      if(std::isnan(sumMinors[0]))
	r->oklist[it] = 2;
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(numInds+1));
	myCounter++;
      }
    }else{
      for(int i=0;i<2*numInds+1;i++)
	sumMinors[i] = log(sumMinors[i]);
      angsd::logrescale(sumMinors,2*numInds+1);
      if(std::isnan(sumMinors[0]))
	r->oklist[it] = 2;
      else{
	r->oklist[it] = 1;
	r->pLikes[myCounter] =new double[2*numInds+1];
	memcpy(r->pLikes[myCounter],sumMinors,sizeof(double)*(2*numInds+1));
	myCounter++;
      }
    }
  }
  
}



//Basicly the same as algoBayAll but only looping through the 3 derived given that an ancestral exists;




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
  //  fprintf(stderr,"-doSaf:%d -p->nind:%d\n",doSaf,p->nInd);exit(0);
  if(doSaf>0&&doSaf!=3){
    realRes *r = new realRes;
    r->oklist=new char[p->numSites];
    memset(r->oklist,0,p->numSites);
    r->pLikes=new double*[p->numSites];
    
    if(doSaf==1&&isHap==0)
      algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans);
    else   if(doSaf==1&&isHap==1)
      algoJointHap(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans);
    else if(doSaf==2){
      freqStruct *freq = (freqStruct *) p->extras[6];
      filipe::algoJoint(p->likes,p->anc,p->numSites,p->nInd,underFlowProtect,fold,p->keepSites,r,noTrans,doSaf,p->major,p->minor,freq->freq,filipeIndF,newDim);
    }else if(doSaf==4){
    algoJointPost(p->post,p->numSites,p->nInd,p->keepSites,r,fold);
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
  
  //  realRes *r=(realRes *) p->extras[index];
  int id=0;
  for(int i=0;i<p->numSites;i++)
    if(r->oklist[i]==1)
      delete [] r->pLikes[id++];
  delete [] r->pLikes;
  delete [] r->oklist;
  delete r;

}





void printFull(funkyPars *p,int index,FILE *outfileSFS,gzFile outfileSFSPOS,char *chr,int newDim){

  realRes *r=(realRes *) p->extras[index];
  int id=0;
  
  for(int s=0; s<p->numSites;s++){
    if(r->oklist[s]==1)
      fwrite(r->pLikes[id++],sizeof(double),newDim,outfileSFS);
  }
  
  kstring_t kbuf;kbuf.s=NULL;kbuf.l=kbuf.m=0;
  for(int i=0;i<p->numSites;i++)
    if(r->oklist[i]==1)
      ksprintf(&kbuf,"%s\t%d\n",chr,p->posi[i]+1);
    else if (r->oklist[i]==2)
      fprintf(stderr,"PROBS at: %s\t%d\n",chr,p->posi[i]+1);
  
  gzwrite(outfileSFSPOS,kbuf.s,kbuf.l);
  free(kbuf.s);
  
}


void abcSaf::calcThetas(funkyPars *pars,int index,double *prior,gzFile fpgz){
 realRes *r=(realRes *) pars->extras[index];
 int id=0;
 
 for(int i=0; i<pars->numSites;i++){
   
   if(r->oklist[i]==1){
     double *workarray = r->pLikes[id++];
    
     //First find thetaW: nSeg/a1
     double pv,seq;
     if(fold)
       pv = 1-exp(workarray[0]);
     else
       pv = 1-exp(workarray[0])-exp(workarray[2*pars->nInd]);
     
     if(pv<0)//catch underflow
       seq=log(0.0);
     else
       seq = log(pv)-aConst;//watterson
     gzprintf(fpgz,"%s\t%d\t%f\t",header->target_name[pars->refId],pars->posi[i]+1,seq);

     if(fold==0) {
       double pairwise=0;    //Find theta_pi the pairwise
       double thL=0;    //Find thetaL sfs[i]*i;
       double thH=0;//thetaH sfs[i]*i^2
       for(size_t ii=1;ii<2*pars->nInd;ii++){
	 
	 pairwise += exp(workarray[ii]+scalings[ii] );
	 double li=log(ii);
	 
	 thL += exp(workarray[ii])*ii;
	 thH += exp(2*li+workarray[ii]);
       }
       gzprintf(fpgz,"%f\t%f\t%f\t%f\n",log(pairwise)-aConst2,workarray[1],log(thH)-aConst2,log(thL)-aConst3);
     }else{
       double pairwise=0;    //Find theta_pi the pairwise
       for(size_t ii=1;ii<newDim;ii++)
	 pairwise += exp(workarray[ii]+scalings[ii] );
       
       gzprintf(fpgz,"%f\t-Inf\t-Inf\t-Inf\n",log(pairwise)-aConst2);
     }
   }else if(r->oklist[i]==2)
     fprintf(stderr,"PROBS at: %s\t%d\n",header->target_name[pars->refId],pars->posi[i]+1);
 }
}




void abcSaf::print(funkyPars *p){
  if(p->numSites==0||(doSaf==0))
    return;
  //  fprintf(stderr,"newDim:%d doSaf:%d\n",newDim,doSaf);
  if(doSaf==3){
    kstring_t *buf =(kstring_t *) p->extras[index];
    gzwrite(outfileGprobs,buf->s,buf->l);
    free(buf->s);delete buf;
  }else{
    realRes *r=(realRes *) p->extras[index];
    int id=0;
    //first addprior if this has been supplied
    if(prior!=NULL) {
      for(int s=0; s<p->numSites;s++) {
	double *workarray = NULL;
      if(r->oklist[s]==1)
	workarray = r->pLikes[id++];
      else
	continue;
      
      double tsum =exp(workarray[0] + prior[0]);
      for(int i=1;i<newDim;i++)
	tsum += exp(workarray[i]+prior[i]);
      tsum = log(tsum);
      
    for(int i=0;i<newDim;i++)
      workarray[i] = workarray[i]+prior[i]-tsum;
    
    
    //pLikes now contains our posterior expectation of the different classes of frequencies
      }
    }
 
    if(doThetas==0){
      if(isHap==0)
	printFull(p,index,outfileSFS,outfileSFSPOS,header->target_name[p->refId],newDim);
      else
	printFull(p,index,outfileSFS,outfileSFSPOS,header->target_name[p->refId],(newDim-1)/2+1);
      }
    else 
      calcThetas(p,index,prior,theta_fp);
  }   
}





//Functions below  should be added 



void abcSaf::algoGeno(int refId,double **liks,char *major,char *minor,int nsites,int numInds,kstring_t *sfsfile,int underFlowProtect,int *posi,int *keepSites,double *pest) {
  assert(pest!=NULL);
  //void algGeno(aMap &asso,int numInds,FILE *sfsfile,int underFlowProtect, double *pest) {
  //  fprintf(stderr,"UNDERFLOWPROTECT: %d\n",underFlowProtect);
  
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
	  hj[i] =  (hj[i]-lbico(2*(numInds-1),i));
	else
	  hj[i] =  exp(log(hj[i])-lbico(2*(numInds-1),i));
	
      }
       
      //now do all the genocalling for individual =select
      // fprintf(stderr,"seelct=%d\tmyMaj=%d\tmyMin=%d\tAA_offset=%d Aa_offset=%d aa_offset=%d\n",select,p.major,p.minor,AA_offset,Aa_offset,aa_offset);
      double *asdf = hj;

      //print_array(stdout,asdf,2*numInds-1,!underFlowProtect);
      double res[3]; //is always logged


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
    for(int i=0;i<numInds;i++)
      ksprintf(sfsfile,"%f ",whichProb[i]);
    ksprintf(sfsfile,"\n");
    //    fprintf(sfsfile,"%f\t%f\n",whichProb[numInds-1],it->second.emPhat);
  }//after all sites
  
}

