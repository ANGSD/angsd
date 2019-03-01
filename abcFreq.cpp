/*
  thorfinn@binf.ku.dk
  part of angsd.

  gnu license, mit, bsd whatever free license

  BUGS:
  phat doesn't put data into freq->freq. Maybe this should be fixed, in either case, nobody is using the phat estimator

*/

#include <cassert>
#include <cmath>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include "abcFreq.h"
#include "abcSaf.h"
#include "abcFilter.h"
#include "analysisFunction.h"

int abcFreq::emIter = EM_NITER;
double abcFreq::EM_start = EM_START;
double *abcFreq::indF = NULL;

//simple phat estimator from the 200danes article, one site function
double phatFun(suint *counts,int nFiles,double eps,char major,char minor) {

  double pis[nFiles];
  double wis[nFiles];

  if(major==minor)
    return 0.0;

  //if the site is variable
  for(int i=0;i<nFiles;i++){
    int ni = counts[i*4+minor];
    int nt= counts[i*4+minor] + counts[i*4+major]; 
    if(nt==0){//if we dont have any reads for individual 'i'
      pis[i] = 0;
      wis[i] = 0;
      continue;
    }
    pis[i] = (ni-eps*nt)/(nt*(1-2*eps));
    wis[i] = 2.0*nt/(nt+1.0);
  }
  
  double tmp=0;
  for(int i=0;i<nFiles;i++)
    tmp+= (pis[i]*wis[i]);
  
  double phatV = std::max(0.0,tmp/angsd::sum<double>(wis,nFiles));
  return phatV;

}

//simple phat estimator will loop over all sites in the pars
void phatLoop(funkyPars *pars,double eps,double nInd,freqStruct *freqs){
  freqs->phat = new double[pars->numSites];
  freqs->freq = new double[pars->numSites];
  for(int s=0;s<pars->numSites;s++)
    if(pars->keepSites[s]){
      freqs->phat[s] = phatFun(pars->counts[s],nInd,eps,pars->major[s],pars->minor[s]);
      freqs->freq[s] = freqs->phat[s];
    }
}

void abcFreq::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doMaf\t%d (Calculate persite frequencies \'.mafs.gz\')\n",doMaf);
  fprintf(argFile,"\t1: Frequency (fixed major and minor)\n");
  fprintf(argFile,"\t2: Frequency (fixed major unknown minor)\n");
  fprintf(argFile,"\t4: Frequency from genotype probabilities\n");
  fprintf(argFile,"\t8: AlleleCounts based method (known major minor)\n");
  fprintf(argFile,"\tNB. Filedumping is supressed if value is negative\n");
  fprintf(argFile,"-doPost\t%d\t(Calculate posterior prob 3xgprob)\n",doPost);
  fprintf(argFile,"\t1: Using frequency as prior\n");
  fprintf(argFile,"\t2: Using uniform prior\n");
  fprintf(argFile,"\t3: Using SFS as prior (still in development)\n");
  fprintf(argFile,"\t4: Using reference panel as prior (still in development), requires a site file with chr pos major minor af ac an\n");
  fprintf(argFile,"Filters:\n");
  fprintf(argFile,"\t-minMaf  \t%f\t(Remove sites with MAF below)\n",minMaf);
  fprintf(argFile,"\t-SNP_pval\t%f\t(Remove sites with a pvalue larger)\n",SNP_pval);
  fprintf(argFile,"\t-rmTriallelic\t%f\t(Remove sites with a pvalue lower)\n",rmTriallelic);
  fprintf(argFile,"Extras:\n");
  fprintf(argFile,"\t-ref\t%s\t(Filename for fasta reference)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(Filename for fasta ancestral)\n",ancName);
  fprintf(argFile,"\t-eps\t%f [Only used for -doMaf &8]\n",eps);
  
  fprintf(argFile,"\t-beagleProb\t%d (Dump beagle style postprobs)\n",beagleProb);
  fprintf(argFile,"\t-indFname\t%s (file containing individual inbreedcoeficients)\n",indFname);
  fprintf(argFile,"\t-underFlowProtect\t%d (file containing individual inbreedcoeficients)\n",underflowprotect);
  fprintf(argFile,"NB These frequency estimators requires major/minor -doMajorMinor\n");
}

//fancy little function
int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

 
void abcFreq::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;
  
  doMaf=angsd::getArg("-doMaf",doMaf,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);

  if(doMaf==0 && doPost ==0)
    return;
  
  rmTriallelic=angsd::getArg("-rmTriallelic",rmTriallelic,arguments);

  GL=angsd::getArg("-GL",GL,arguments);
  if(inputtype!=INPUT_VCF_GL && inputtype!=INPUT_GLF && inputtype!=INPUT_GLF3 && doPost && GL==0){
    fprintf(stderr,"\t-> Potential problem: You are required to choose a genotype likelihood model (-GL) for estimating genotypes posteriors.\n");
    exit(0);
  } 
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  //  fprintf(stderr,"%d %d\n",inputtype,INPUT_VCF_GL);
  if((doPost||doMaf) &&doMajorMinor==0 &&inputtype!=INPUT_VCF_GL){
    if(doMaf!=4){
      fprintf(stderr,"\t-> Potential problem: You are required to choose a major/minor estimator (-doMajorMinor)\n");
      exit(0);
    }
  } 
  if(doPost==1 &&doMaf==0){
    fprintf(stderr,"\t-> Potential problem: You are required to choose a maf estimator (-doMaf) with -doPost 1\n");
    exit(0);
  } 

  indFname = angsd::getArg("-indF",indFname,arguments);
  if(indFname !=NULL)
    indF = angsd::readDouble(indFname,arguments->nInd);





  if(doMaf==0 )//&& doPost==0?
    return;
  underflowprotect=angsd::getArg("-underFlowProtect",underflowprotect,arguments);
  chisq1 = new Chisqdist(1);
  chisq2 = new Chisqdist(2);
  chisq3 = new Chisqdist(3);



  
  minMaf=angsd::getArg("-minMaf",minMaf,arguments);
  //  assert(minMaf<=1&&minMaf>=0);

  double tmp=-1;
  tmp=angsd::getArg("-SNP_pval",tmp,arguments);
 
  if(tmp!=-1){ 
    SNP_pval = tmp;
    double pre = SNP_pval;
    //    fprintf(stderr,"ind:%f :%p: \t",SNP_pval,chisq3);
    if(SNP_pval <=1){
      if(abs(doMaf) ==2)
	SNP_pval = chisq1->invcdf(1-SNP_pval);
      if(abs(doMaf) ==1)
	SNP_pval = chisq1->invcdf(1-SNP_pval);
    }
    doSNP =1 ;
    fprintf(stderr,"\t-> SNP-filter using a pvalue: %e correspond to %f likelihood units\n",pre,SNP_pval);
  }
  if(rmTriallelic)
    SNP_pval_tri = chisq2->invcdf(1-rmTriallelic);
  refName = angsd::getArg("-ref",refName,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);

  if(abs(doMaf)&& !isPowerOfTwo((unsigned int) abs(doMaf))){
    fprintf(stderr,"\n[%s] You have selected filters for maf/lrt\n",__FILE__);
    fprintf(stderr,"If you have selected more than one MAF estimator we will choose in following order\n");
    fprintf(stderr,"\t1. knownminor EM\n");
    fprintf(stderr,"\t2. unknownminor EM\n");
    fprintf(stderr,"\t3. Posterior maf\n");
  }
  if(doSNP && abs(doMaf)&8){
    fprintf(stderr,"\t-> LRT based snp calling is not defined for count based maf estimates\n");
    exit(0);
  }
  if(doSNP&&abs(doMaf)==0){
    fprintf(stderr,"You've selected SNP-calling but no maf estimator please select -doMaf INT\n");
    exit(0);
  }

  if(doSNP==1 && abs(doMaf)==0){
    fprintf(stderr,"You've selected SNP_pval threshold please also select -doMaf \n");
    exit(0);
  }
  if(minMaf>0.0 &&(abs(doMaf)==0)){
    fprintf(stderr,"\nYou've selected minMaf but no MAF estimator, choose -doMaf\n\n");
    exit(0);
  }
  if(rmTriallelic>0.0 &&(abs(doMaf)!=1)){
    fprintf(stderr,"\n rmTriallelic only works with -doMaf 1 \n\n");
    exit(0);
  }

  int doCounts =0;
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);

  beagleProb=angsd::getArg("-beagleProb",beagleProb,arguments);
  minInd=angsd::getArg("-minInd",minInd,arguments);

  if(abs(doMaf)==0 &&doPost==0)
    return;
  if(inputtype!=INPUT_VCF_GP&&inputtype!=INPUT_BEAGLE&&inputtype!=0&&doMajorMinor==0&&inputtype!=INPUT_VCF_GL){
    fprintf(stderr,"You must specify \'-doMajorMinor\' to infer major/minor \n");
    exit(0);
  }

  if(inputtype==INPUT_BEAGLE&&abs(doMaf)!=4){
    fprintf(stderr,"Only \'-doMaf 4\' can not be performed on posterior input\n");
    exit(0);
  }
  if(inputtype!=INPUT_BEAGLE&&inputtype!=INPUT_VCF_GP&&abs(doMaf)==4){
    fprintf(stderr,"\t \'-doMaf 4\' can only be performed on genotype probabilities provided by the user (-beagle).\n");
    exit(0);
  }

  if((inputtype==INPUT_BAM || inputtype==INPUT_PILEUP )&& GL==0 &&doMaf!=8){
    fprintf(stderr,"Error: For sequence data (BAM,pileups) likehoods (-GL) must be specified for frequency estimation\n");
    exit(0);
  }

  if(doSNP&&abs(doMaf) &4){
    fprintf(stderr,"Error: -doMaf 4 cannot be used for a likelihood ratio test (doSNP) \n");
    exit(0);
  }
  if(inputtype==INPUT_BEAGLE&&doPost){
    fprintf(stderr,"Error: Cannot estimate post (doPost) based on posterior probabilites\n");
    exit(0);
  }
  if(abs(doMaf)&8&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts for MAF estimator based on counts\n");
    exit(0);
  }

  if(beagleProb && doPost==0 &&inputtype!=INPUT_VCF_GP){
    fprintf(stderr,"Must supply -doPost 1 to write beaglestyle postprobs\n");
    exit(0);
  }

  if(doPost!=0 && doMajorMinor==0 &&INPUT_VCF_GL!=inputtype){
    fprintf(stderr,"Do post requires major and minor: supply -doMajorMinor \n");
    exit(0);
  }

}

//constructor
abcFreq::abcFreq(const char *outfiles,argStruct *arguments,int inputtype){
  underflowprotect =0;
  minInd = 0;
  chisq1=chisq2=chisq3=NULL;
  inputIsBeagle =0;
  beagleProb = 0; //<-output for beagleprobs
  minMaf =-1.0;
  SNP_pval = 1;
  nInd = arguments->nInd;
  eps = 0.001;
  outfileZ2 = NULL;
  outfileZ = NULL;
  indFname = NULL;
  doMaf=0;
  rmTriallelic=0;
  GL=0;
  doSNP=0;
  doPost=0;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  //emIter=EM_NITER; //these are static see top of this file
  //EM_start = EM_START; //these are static see top of this file
  doMajorMinor=0;
  refName = NULL;
  ancName = NULL;


  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doMaf")||!strcasecmp(arguments->argv[1],"-doPost")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }



  
  getOptions(arguments);

  if(doMaf==0 && doPost==0){
    shouldRun[index]=0;
    return;
  }
  printArg(arguments->argumentFile);
  if(doMaf>0){
    //make output files
    const char* postfix;
    postfix=".mafs.gz";
    outfileZ = aio::openFileBG(outfiles,postfix);
    if(beagleProb){
      postfix=".beagle.gprobs.gz";
      outfileZ2 = aio::openFileBG(outfiles,postfix);
    }
  }else
    doMaf=abs(doMaf);
  //print header
  
  kputs("chromo\tposition\tmajor\tminor\t",&bufstr);
  if(refName!=NULL||arguments->inputtype==INPUT_PILEUP)
    kputs("ref\t",&bufstr);
  if(ancName)
    kputs("anc\t",&bufstr);
  
  if(doMaf &1)
    kputs("knownEM\t",&bufstr);
  if(doMaf &2)
    kputs("unknownEM\t",&bufstr);    
  if(doMaf &4)
    kputs("PPmaf\t",&bufstr);
  if(doMaf &8)
    kputs("phat\t",&bufstr);
  
  if(doSNP){
    if(doMaf &1)
      kputs("pK-EM\t",&bufstr);
    if(doMaf &2)
      kputs("pu-EM\t",&bufstr);
  }
  kputs("nInd\n",&bufstr);
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);
  bufstr.l=0;
  if(beagleProb){
    kputs("marker\tallele1\tallele2",&bufstr);
    for(int i=0;i<arguments->nInd;i++){
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
      kputs("\tInd",&bufstr);
      kputw(i,&bufstr);
    }
    kputc('\n',&bufstr);
    aio::bgzf_write(outfileZ2,bufstr.s,bufstr.l);
    bufstr.l=0;
  }

}


abcFreq::~abcFreq(){
  if(outfileZ!=NULL){
    bgzf_close(outfileZ);
  }
  if(outfileZ2!=NULL){
    bgzf_close(outfileZ2);
  }
  free(refName);
  free(ancName);
  delete [] indF;
  delete chisq1;
  delete chisq2;
  delete chisq3;
  free(bufstr.s);
}



void abcFreq::print(funkyPars *pars) {
  if(outfileZ==NULL&&outfileZ2==NULL)
    return;

  freqStruct *freq =(freqStruct *) pars->extras[index];

  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    //plugin chr,pos,major,minor
    kputs(header->target_name[pars->refId],&bufstr);kputc('\t',&bufstr);
    kputw(pars->posi[s]+1,&bufstr);kputc('\t',&bufstr);
    kputc(intToRef[pars->major[s]],&bufstr);kputc('\t',&bufstr);
    kputc(intToRef[pars->minor[s]],&bufstr);kputc('\t',&bufstr);


    //plugin ref, anc if exists
    if(pars->ref!=NULL)
      {kputc(intToRef[pars->ref[s]],&bufstr);kputc('\t',&bufstr);}
    if(pars->anc!=NULL)
      {kputc(intToRef[pars->anc[s]],&bufstr);kputc('\t',&bufstr);}

    
    
    if(doMaf &1)
      ksprintf(&bufstr,"%f\t",freq->freq_EM[s]);
    if(doMaf &2)
      ksprintf(&bufstr,"%f\t",freq->freq_EM_unknown[s]);
    if(doMaf &4)
      ksprintf(&bufstr,"%f\t",freq->freq[s]);
    if(doMaf &8)
      ksprintf(&bufstr,"%f\t",freq->phat[s]);
    if(doSNP){
      if(doMaf &1)
	ksprintf(&bufstr,"%e\t",angsd::to_pval(chisq1,freq->lrt_EM[s]));
      if(doMaf &2)
	ksprintf(&bufstr,"%e\t",angsd::to_pval(chisq1,freq->lrt_EM_unknown[s]));
    }

    kputw(pars->keepSites[s],&bufstr);kputc('\n',&bufstr);
  }

  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);  
  bufstr.l=0;

  if(beagleProb){
    //beagle format
    for(int s=0;s<pars->numSites;s++) {
      
      if(pars->keepSites[s]==0)
	continue;
      // fprintf(stderr,"keepsites=%d\n",pars->keepSites[s]);
      kputs(header->target_name[pars->refId],&bufstr);
      kputc('_',&bufstr);
      kputw(pars->posi[s]+1,&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->major[s],&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->minor[s],&bufstr);

      int major = pars->major[s];
      int minor = pars->minor[s];
      assert(major!=4&&minor!=4);
	
      for(int i=0;i<3*pars->nInd;i++) {
	ksprintf(&bufstr, "\t%f",pars->post[s][i]);
      }
      
      kputc('\n',&bufstr);
    
    }
    //valgrind on osx complains here check if prob on unix
    int ret=aio::bgzf_write(outfileZ2,bufstr.s,bufstr.l);
    bufstr.l=0;
    //fprintf(stderr,"ret.l:%d bufstr.l:%zu\n",ret,bufstr.l);

  }

}



void abcFreq::clean(funkyPars *pars) {
  if(doMaf==0&&doPost==0)
    return;

  freqStruct *freq =(freqStruct *) pars->extras[index];

  //cleaning
  delete [] freq->freq;
  delete [] freq->lrt;//maybe in doSNP
  delete [] freq->freq_EM;
  delete [] freq->freq_EM_unknown;
  delete [] freq->lrt_tri;
  delete [] freq->lrt_EM;
  delete [] freq->lrt_EM_unknown;
  delete [] freq->phat;
  delete freq;

  
  if(pars->post!=NULL){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->post[i];
    delete [] pars->post;
    pars->post =NULL;
  }


}

freqStruct *allocFreqStruct(){
  freqStruct *freq = new freqStruct;
  
 freq->freq_EM=NULL;
 freq->freq_EM_unknown= NULL;
 freq->lrt_tri=NULL;
 freq->lrt_EM=NULL;
 freq->lrt_EM_unknown=NULL;
 freq->freq = NULL;
 freq->lrt = NULL;
 freq->phat = NULL;
 return freq;
}

//filipes
void make_post_F(double *like,double *post,double freqEst,int nInd){
  //fprintf(stderr,"[%s]\n",__FUNCTION__);
  double *indF = abcFreq::indF;
  for(int i=0;i<nInd;i++){
    post[i*3+2]=like[i*3+2] + log(pow(freqEst,2) + (1-freqEst)*freqEst*indF[i]);
    post[i*3+1]=like[i*3+1] + log(2*(1-freqEst)*freqEst - 2*(1-freqEst)*freqEst*indF[i]);
    post[i*3+0]=like[i*3+0] + log(pow(1-freqEst,2) + (1-freqEst)*freqEst*indF[i]);
    
    double norm = angsd::addProtect3(post[i*3+0],post[i*3+1],post[i*3+2]);
    post[i*3+0]=exp(post[i*3+0]-norm);
    post[i*3+1]=exp(post[i*3+1]-norm);
    post[i*3+2]=exp(post[i*3+2]-norm);
  }
  //
}



void make_post(double *like,double *post,double freqEst,int nInd){
  if(abcFreq::indF!=NULL){
    make_post_F(like,post,freqEst,nInd);
    return;
  }
  for(int i=0;i<nInd;i++){
    //	  fprintf(stderr,"[%d]\nlik= %f %f %f\n",i,like[0][i*3+0],like[0][i*3+1],like[0][i*3+2]);
    post[i*3+2]=like[i*3+2]+2*log(freqEst);
    post[i*3+1]=like[i*3+1]+log(2)+log(1-freqEst)+log(freqEst);
    post[i*3+0]=like[i*3+0]+2*log(1-freqEst);
    //	  fprintf(stderr,"likmod= %f %f %f\n",post[0][i*3+0],post[0][i*3+1],post[0][i*3+2]);
    double norm = angsd::addProtect3(post[i*3+0],post[i*3+1],post[i*3+2]);
    post[i*3+0]=exp(post[i*3+0]-norm);
    post[i*3+1]=exp(post[i*3+1]-norm);
    post[i*3+2]=exp(post[i*3+2]-norm);
    //	  fprintf(stderr,"%f %f %f\n",post[0][i*3+0],post[0][i*3+1],post[0][i*3+2]);
  }
}


void abcFreq::run(funkyPars *pars) {
  if(doMaf==0&&doPost==0)
    return;
  freqStruct *freq = allocFreqStruct();
  pars->extras[index] = freq;
  extern abc **allMethods;
  
  //  abcFilter *abcf=(abcFilter *) allMethods[0];
  filt *fl = ((abcFilter *) allMethods[0])->fl;

  if(doMaf!=0) {

    if(doMaf&8)
      phatLoop(pars,eps,nInd,freq);
    
    if(doMaf&4)
      postFreq(pars,freq);
    if(doMaf % 4){ // doMaf =1 or doMaf 2
      likeFreq(pars,freq);
    }
    assert(freq->freq!=NULL);
    
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;

      if(freq->freq[s] < minMaf)
	pars->keepSites[s]=0;
      else if(freq->freq[s] > 1 - minMaf)
	pars->keepSites[s]=0;
     
      if(doSNP && (freq->lrt[s] < SNP_pval))
      	pars->keepSites[s]=0;
      if(rmTriallelic && (freq->lrt_tri[s] > SNP_pval_tri))
      	pars->keepSites[s]=0;
      if(minInd>0&&pars->keepSites[s]<minInd)
	pars->keepSites[s] = 0;
    }
  }
  if(doPost) {
    assert(pars->likes!=NULL);

    double **post = new double*[pars->numSites];
    double **like=angsd::get3likes(pars);
    for(int s=0;s<pars->numSites;s++){
      post[s] = new double [3*pars->nInd];

      if(pars->keepSites[s]==0){
	delete [] like[s];
	continue;
      }

      if(doPost==1)  //maf prior
	 make_post(like[s],post[s],freq->freq[s],pars->nInd);
      else if(doPost==2) {//uniform prior
	for(int i=0;i<pars->nInd;i++){
	  double norm = angsd::addProtect3(like[s][i*3+0],like[s][i*3+1],like[s][i*3+2]);
	  for(int g=0;g<3;g++)
	    post[s][i*3+g]=exp(like[s][i*3+g]-norm);
	}
      }
      else if(doPost==3){//uniform prior
	if(algoGeno(pars->likes[s],refToInt[pars->major[s]],refToInt[pars->minor[s]],pars->nInd,underflowprotect,abcSaf::prior,post[s]))
	   fprintf(stderr,"\t-> Problem calling genotypes at: (%s,%d)\n",header->target_name[pars->refId],pars->posi[s]+1);	   
      }
      else if(doPost==4){
	assert(fl!=NULL);
	//af for this pos: fl->af[pars->posi[s]]
	//an for this pos: fl->an[pars->posi[s]]
	//ac for this pos: fl->ac[pars->posi[s]]
	//	fprintf(stderr,"chr: %s pos: %d freq:%f\n",header->target_name[pars->refId],pars->posi[s]+1,fl->af[pars->posi[s]]);
	if(fl->keeps[pars->posi[s]]==0)
	  make_post(like[s],post[s],freq->freq[s],pars->nInd);
	else{
	  double exp_obs_allels_in_ngs= freq->freq[s]*pars->keepSites[s];
	  //	  fprintf(stderr,"pars->keepSite:%d\n",pars->k)
	  double aac = fl->ac[pars->posi[s]];
	  double aan = fl->an[pars->posi[s]];
	  double af = (aac+exp_obs_allels_in_ngs)/(1.0*(aan+pars->keepSites[s]));
	  fprintf(stderr,"aac:%f aan:%f nal_in_ngs:%f af:%f\n",aac,aan,exp_obs_allels_in_ngs,af);
	  make_post(like[s],post[s],af,pars->nInd);
	}

      }
      else{
	fprintf(stderr,"[%s] doPost must be 1 or 2 \n",__FUNCTION__);
	exit(0);
      
      }
      delete [] like[s];
    }
    pars->post = post; 
    delete[] like;  
  }
  
}


/*
Estimate the allele frequency from the posterior probabilities
*/
void abcFreq::postFreq(funkyPars  *pars,freqStruct *freq){
  if(freq->freq==NULL)
  freq->freq = new double[pars->numSites]; 

  for(int s=0;s<pars->numSites;s++){
    freq->freq[s]=0;
    for(int i=0;i<pars->nInd;i++){
      freq->freq[s] += pars->post[s][i*3+1]+2*pars->post[s][i*3+2];
    }
    freq->freq[s] = freq->freq[s]/(pars->nInd*2);
  }

}


void abcFreq::likeFreq(funkyPars *pars,freqStruct *freq){//method=1: bfgs_known ;method=2 em;method=4 bfgs_unknown
  
  //here only the likelihoods for the three genotypes are used. 
  double **loglike = NULL;
  if(inputIsBeagle==1)
    loglike= pars->likes;
  else
    loglike=angsd::get3likesRescale(pars);
  assert(loglike!=NULL);
  if(freq->freq == NULL)
  freq->freq = new double[pars->numSites];
  if(doSNP)
    freq->lrt = new double[pars->numSites];

  if(rmTriallelic)
    freq->lrt_tri = new double[pars->numSites];
  if(doMaf &1){ 
    freq->freq_EM = new double[pars->numSites]; 
    if(doSNP)
      freq->lrt_EM = new double[pars->numSites];
  }
  if(doMaf &2){
    freq->freq_EM_unknown =new double[pars->numSites];
    if(doSNP)
      freq->lrt_EM_unknown =new double[pars->numSites];
  }

  // number of individuals with data
  int *keepInd = pars->keepSites;
  int keepList[pars->nInd];  

  //loop though all sites and check if we have data.
  for(int s=0;s<pars->numSites;s++) {
    if(keepInd[s]==0)//if we dont have any information
      continue;
    keepInd[s]=0;//
    for(int i=0 ; i<pars->nInd ;i++) {//DRAGON CHECK THIS
      keepList[i]=1;
      //also discard if all gls are the same
      if(loglike[s][i*3+0]+loglike[s][i*3+1]+loglike[s][i*3+2]>-0.0001||((loglike[s][i*3]==loglike[s][i*3+1] )&& (loglike[s][i*3]==loglike[s][i*3+2]))){
	keepList[i]=0;
      }      
      else{
	keepInd[s]++;
      }

    }

    if(keepInd[s]==0)//if we dont have any information
      continue;

    if(doMaf &1) {
      double mstart = EM_start;
      if(freq->phat!=NULL)
	mstart = freq->phat[s];

      freq->freq_EM[s]=emFrequency(loglike[s],pars->nInd,emIter,mstart,keepList,keepInd[s]);
      if(doSNP){
	freq->lrt_EM[s]= 2*likeFixedMinor(0.0,loglike[s],pars->nInd)-2*likeFixedMinor(freq->freq_EM[s],loglike[s],pars->nInd);
	if(freq->lrt_EM[s]<0)
	  freq->lrt_EM[s]=0;
      }
    }
    
    if( doMaf &2) {
      double mstart = EM_start;
      if(freq->phat!=NULL)
	mstart = freq->phat[s];
      
      freq->freq_EM_unknown[s]=emFrequencyNoFixed(pars->likes[s],pars->nInd,emIter,mstart,keepList,keepInd[s],pars->major[s],pars->posi[s]);
      if(doSNP){
	freq->lrt_EM_unknown[s] = 2*likeNoFixedMinor(0.0,pars->likes[s],pars->nInd,pars->major[s])-2*likeNoFixedMinor(freq->freq_EM_unknown[s],pars->likes[s],pars->nInd,pars->major[s]);
	if(freq->lrt_EM_unknown[s]<0)
	  freq->lrt_EM_unknown[s]=0;
      }
    }
    if(rmTriallelic!=0){
      double lrtTri=0;
      if(freq->freq_EM[s]>0.0001)
       lrtTri=isMultiAllelic(pars->likes[s],pars->nInd,emIter,freq->freq_EM[s],keepList,keepInd[s],pars->major[s],pars->minor[s]);
      
      freq->lrt_tri[s] = lrtTri;
      
    }



  }

  if(inputIsBeagle!=1){
    for(int i=0;i<pars->numSites;i++)
      delete [] loglike[i];
    delete [] loglike;
  }
 
  for(int s=0;s<pars->numSites;s++){
    if(doMaf &2 )
      freq->freq[s]=freq->freq_EM_unknown[s];
    else if(doMaf &1 )
      freq->freq[s]=freq->freq_EM[s];
  }

  //thorfinn april 16 2012
  for(int s=0;s<pars->numSites;s++){
    if((doMaf &2) && doSNP )
      freq->lrt[s]=freq->lrt_EM_unknown[s];
    else if((doMaf &1) && doSNP)
      freq->lrt[s]=freq->lrt_EM[s];
  }

}

double abcFreq::emFrequencyNoFixed_ext(double *loglike,int numInds,int *keep,int keepInd,int major,int posi){

  return emFrequencyNoFixed(loglike,numInds, emIter,EM_start, keep,keepInd,major,posi);
}




double abcFreq::emFrequencyNoFixed_F(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int posi){

  int iMinor[3];
  int n=0;
  for(int i=0;i<4;i++){
    if(i!=major){
      iMinor[n]=i;
      n++;
    }
  }
    
  double **loglikeGeno;
  loglikeGeno =new double*[3];
  
  for(int j=0;j<3;j++){
    loglikeGeno[j] = new double[numInds*3]; 
    for(int i=0;i<numInds;i++){
      loglikeGeno[j][i*3+0]=loglike[i*10+angsd::majorminor[major][major]];
      loglikeGeno[j][i*3+1]=loglike[i*10+angsd::majorminor[major][iMinor[j]]];
      loglikeGeno[j][i*3+2]=loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]];
    }
  }
  
  float W0[3];
  float W1[3];
  float W2[3];
  start=0.4;
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;
  
  int it=0;
  float weight[3];
  float normWeight[3];
  float norm;
  int correctInd;
  
  for(it=0;it<iter;it++){
    for(int j=0;j<3;j++)
      weight[j]=-likeFixedMinor(p,loglikeGeno[j],keepInd);
    norm=angsd::addProtect3(weight[0],weight[1],weight[2]);
    for(int j=0;j<3;j++)
      normWeight[j]=exp(weight[j]-norm);
    sum=0;
    correctInd=keepInd;
    for(int i=0;i<numInds;i++){
      if(keep[i]==0)
        continue;
      for(int j=0;j<3;j++){
	W0[j]=exp(loglike[i*10+angsd::majorminor[major][major]]) * (pow(1-p,2) + (1-p)*p*indF[i]);
	W1[j]=exp(loglike[i*10+angsd::majorminor[major][iMinor[j]]]) * (2*(1-p)*p - 2*(1-p)*p*indF[i]);
	W2[j]=exp(loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]]) * (pow(p,2) + (1-p)*p*indF[i]);
	if(W0[j]+W1[j]+W2[j]<0.0000001)
	  continue;
      }
    }
    
    p=sum/correctInd;
       
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
       
    temp_p=p;
    if(std::isnan(p)){
      fprintf(stderr,"[%s] caught nan at position: %d\n",__FUNCTION__,posi);
      for(int j=0;j<3;j++)
	fprintf(stderr,"w%d %f ",j,weight[j]);
      fflush(stderr);

      for(int i=0;i<numInds;i++){
	for(int j=0;j<10;j++){
	  if(keep[i]==0)
	    continue;
	  fprintf(stderr,"%f\t",loglike[i*10+j]);
	}
	fprintf(stderr,"\n");
      }
      fprintf(stderr,"%f | %f %f %f |%f | %d  | \n",p, normWeight[0], normWeight[1], normWeight[2],angsd::addProtect3(weight[0],weight[1],weight[2]),major);
    }
  }
  
 for(int j=0;j<3;j++)
   delete[] loglikeGeno[j]; 
 delete[] loglikeGeno;

 return(p);
}



double abcFreq::emFrequencyNoFixed(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int posi){

  if(indF!=NULL)
    return emFrequencyNoFixed_F(loglike,numInds,iter,start,keep,keepInd,major,posi);
  int iMinor[3];
  int n=0;
  for(int i=0;i<4;i++){
    if(i!=major){
      iMinor[n]=i;
      n++;
    }
  }

  double **loglikeGeno;
  loglikeGeno =new double*[3];

  for(int j=0;j<3;j++){
    loglikeGeno[j] = new double[numInds*3]; 
    for(int i=0;i<numInds;i++){
      loglikeGeno[j][i*3+0]=loglike[i*10+angsd::majorminor[major][major]];
      loglikeGeno[j][i*3+1]=loglike[i*10+angsd::majorminor[major][iMinor[j]]];
      loglikeGeno[j][i*3+2]=loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]];
    }
  }

  float W0[3];
  float W1[3];
  float W2[3];
  start=0.4;
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;
  
  int it=0;
  float weight[3];
  float normWeight[3];
  float norm;
  int correctInd;

  for(it=0;it<iter;it++){
    for(int j=0;j<3;j++)
      weight[j]=-likeFixedMinor(p,loglikeGeno[j],keepInd);
    norm=angsd::addProtect3(weight[0],weight[1],weight[2]);
    for(int j=0;j<3;j++)
      normWeight[j]=exp(weight[j]-norm);
    sum=0;
    correctInd=keepInd;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      for(int j=0;j<3;j++){
	W0[j]=exp(loglike[i*10+angsd::majorminor[major][major]])*pow(1-p,2);
	W1[j]=exp(loglike[i*10+angsd::majorminor[major][iMinor[j]]])*2*p*(1-p);
	W2[j]=exp(loglike[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]])*(pow(p,2));
	if(W0[j]+W1[j]+W2[j]<0.0000001)
	  continue;
	sum+=(W1[j]+2*W2[j])/(2*(W0[j]+W1[j]+W2[j]))*normWeight[j];
      }
    }

    p=sum/correctInd;
       
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
       
    temp_p=p;
    if(std::isnan(p)){
      fprintf(stderr,"[%s] caught nan at position: %d will exit\n",__FUNCTION__,posi+1);
      for(int j=0;j<3;j++)
	fprintf(stderr,"w%d %f ",j,weight[j]);
      fflush(stderr);

      for(int i=0;i<numInds;i++){
	for(int j=0;j<10;j++){
	  if(keep[i]==0)
	    continue;
	  fprintf(stderr,"%f\t",loglike[i*10+j]);
	}
	fprintf(stderr,"\n");
      }
      fprintf(stderr,"%f | %f %f %f |%f | %d | \n",p, normWeight[0], normWeight[1], normWeight[2],angsd::addProtect3(weight[0],weight[1],weight[2]),major);

    }
  }
  
  for(int j=0;j<3;j++)
    delete[] loglikeGeno[j]; 
  delete[] loglikeGeno;
  
  return(p);
}
///////////////////////////////////////////////
double abcFreq::isMultiAllelic(double *loglike,int numInds, int iter,double freqStandard,int *keep,int keepInd,int major,int minor){

  double freq[4];
  double temp_freq[4];
  double W[10];
  float sum[4];
  for(int j=0;j<4;j++)
    freq[j] = 0.25;

  double accu=0.00001;
 
  int it=0;
  float norm;
  /* for(int i=0;i<numInds;i++){
    if(keep!=NULL && keep[i]==0)
      continue;
    for(int j=0;j<10;j++)
      fprintf(stderr,"%f ", loglike[i*10+j]);
    fprintf(stderr,"\n");
  }
  */

  for(it=0;it<iter;it++){
    for(int j=0;j<4;j++){
      temp_freq[j] = freq[j] ;
      sum[j]=0;
    }
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;

      W[0] = exp(loglike[i*10+0]) * pow(freq[0],2);     //AA 0
      W[1] = exp(loglike[i*10+1]) * 2*freq[0]*freq[1];  //AC 1
      W[2] = exp(loglike[i*10+2]) * 2*freq[0]*freq[2];  //AG 2 
      W[3] = exp(loglike[i*10+3]) * 2*freq[0]*freq[3];  //AT 3
      W[4] = exp(loglike[i*10+4]) * pow(freq[1],2);  //CC 4
      W[5] = exp(loglike[i*10+5]) * 2*freq[1]*freq[2];  //CG 5 
      W[6] = exp(loglike[i*10+6]) * 2*freq[1]*freq[3];  //CT 6
      W[7] = exp(loglike[i*10+7]) * pow(freq[2],2);       //GG 7
      W[8] = exp(loglike[i*10+8]) * 2*freq[2]*freq[3];  //GT 8
      W[9] = exp(loglike[i*10+9]) * pow(freq[3],2);        //TT 9
      norm=0;
      for(int s=0;s<10;s++){
	norm+=W[s];  
      }
      sum[0]+=(2*W[0]+W[1]+W[2]+W[3])/(2*norm);
      sum[1]+=(2*W[4]+W[1]+W[5]+W[6])/(2*norm);
      sum[2]+=(2*W[7]+W[2]+W[5]+W[8])/(2*norm);
      sum[3]+=(2*W[9]+W[3]+W[6]+W[8])/(2*norm);
       
     
    }
    //   fprintf(stderr,"%f %f %f %f\n",sum[0],sum[1],sum[2],sum[3]);

    for(int j=0;j<4;j++)
      freq[j] = sum[j] / keepInd ;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    int doBreak=1;

    for(int j=0;j<4;j++)
      if((freq[j]-temp_freq[j]>accu || temp_freq[j]-freq[j]>accu))
	doBreak=0;

    if(doBreak)
      break;
   
  }
 
  double likeAlt = likeMultiAllelic(freq,loglike,numInds,keep);
  
  double freq2[4];
  for(int j=0;j<4;j++)
    freq2[j] = 0;
  freq2[minor] = freqStandard;
  freq2[major] = 1-freqStandard;
  double likeNull = likeMultiAllelic(freq2,loglike,numInds,keep);
  //exit(0);
  // fprintf(stderr,"%f\t%f\t%f\t%f\tl=%f\t%f\tf=%f\t%d\t%d\t\n",freq[0],freq[1],freq[2],freq[3],likeAlt,likeNull,freqStandard,major,minor);
  // fprintf(stderr,"%f\t%f\t%f\t%f\tl=%f\t%f\tf=%f\t%d\t%d\t\n",freq2[0],freq2[1],freq2[2],freq2[3],likeAlt,likeNull,freqStandard,major,minor);
 
  return( 2*(likeNull - likeAlt));

}

/////////////////

double abcFreq::likeMultiAllelic(double *freq,double *loglikes,int numInds,int *keep){
  double totalLike=0;
  double W[10]; 
  for(int i=0;i<numInds;i++){
    if(keep!=NULL && keep[i]==0)
      continue;

    W[0] = loglikes[i*10+0] + log(pow(freq[0],2));     //AA 0
    W[1] = loglikes[i*10+1] + log(2*freq[0]*freq[1]);  //AC 1
    W[2] = loglikes[i*10+2] + log(2*freq[0]*freq[2]);  //AG 2 
    W[3] = loglikes[i*10+3] + log(2*freq[0]*freq[3]);  //AT 3
    W[4] = loglikes[i*10+4] + log(pow(freq[1],2));  //CC 4
    W[5] = loglikes[i*10+5] + log(2*freq[1]*freq[2]);  //CG 5 
    W[6] = loglikes[i*10+6] + log(2*freq[1]*freq[3]);  //CT 6
    W[7] = loglikes[i*10+7] + log(pow(freq[2],2));       //GG 7
    W[8] = loglikes[i*10+8] + log(2*freq[2]*freq[3]);  //GT 8
    W[9] = loglikes[i*10+9] + log(pow(freq[3],2));        //TT 9
    totalLike+=angsd::addProtectN(W,10);
  }
    return -totalLike;
}


double abcFreq::likeNoFixedMinor(double p,double *logLikes,int numInds,int major){
  //logLikes contains the 10 likelihoods for one site
  //bfgs must be used
  float partialLike[3];
  for(int j=0;j<3;j++)
    partialLike[j]=0;

  float totalLike=0;
  int iMinor[3];
  int n=0;
  for(int i=0;i<4;i++){
    if(i!=major){
      iMinor[n]=i;
      n++;
    }
  }

  for(int i=0;i<numInds;i++){
    for(int j=0;j<3;j++) 
      partialLike[j]+=
	angsd::addProtect3(
			   logLikes[i*10+angsd::majorminor[major][major]]+2*log(1-p),
			   logLikes[i*10+angsd::majorminor[major][iMinor[j]]]+log(2)+log(p)+log(1-p),
			   logLikes[i*10+angsd::majorminor[iMinor[j]][iMinor[j]]]+2*log(p));
  }
  totalLike=angsd::addProtect3(partialLike[0],partialLike[1],partialLike[2])-log(3);
  return -totalLike;

}


double abcFreq::likeFixedMinor(double p,double *logLikes,int numInds){
  double totalLike=0;
  for(int i=0;i<numInds;i++){
    totalLike+=angsd::addProtect3(logLikes[i*3+0]+log(pow(1-p,2)),
				  logLikes[i*3+1]+log(2)+log(p)+log(1-p),
				  logLikes[i*3+2]+log(pow(p,2)));
  }
    return -totalLike;
}






double abcFreq::emFrequency_ext(double *loglike,int numInds,int *keep,int keepInd){

  return (emFrequency(loglike,numInds,emIter,EM_start,keep,keepInd));
}




double abcFreq::emFrequency_F(double *loglike,int numInds, int iter,double start,int *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;


  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;

      W0=exp(loglike[i*3+0]) * (pow(1-p,2) + (1-p)*p*indF[i]);
      W1=exp(loglike[i*3+1]) * (2*(1-p)*p - 2*(1-p)*p*indF[i]);
      W2=exp(loglike[i*3+2]) * (pow(p,2) + (1-p)*p*indF[i]);
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }

  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);

    for(int ii=0;0&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	    fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;

      W0=exp(loglike[i*3+0]) * (pow(1-p,2) + (1-p)*p*indF[i]);
      W1=exp(loglike[i*3+1]) * (2*(1-p)*p - 2*(1-p)*p*indF[i]);
      W2=exp(loglike[i*3+2]) * (pow(p,2) + (1-p)*p*indF[i]);

      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    // exit(0);
  }
  
  return(p);
}


double abcFreq::emFrequency(double *loglike,int numInds, int iter,double start,int *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  if(indF!=NULL)
    return emFrequency_F(loglike,numInds,iter,start,keep,keepInd);
    
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;

  
  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      if(std::isnan(sum))
	fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f) W(%f,%f,%f) sum=%f\n",i,loglike[i*3],loglike[i*3+1],loglike[i*3+2],W0,W1,W2,sum);
    }
    
    p=sum/keepInd;
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }
  
  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);
    
    for(int ii=0;1&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1){
	//	fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
	fprintf(stderr,"\n");
      }
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"[%s.%s():%d] p=%f W %f\t%f\t%f sum=%f loglike: %f\n",__FILE__,__FUNCTION__,__LINE__,p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
      break;
    }
    p=-999;
    assert(p!=999);
    return p;
  }

  return(p);
}


int abcFreq::algoGeno(double *liks,int major_offset,int minor_offset,int numInds,int underFlowProtect,double *pest,double *postp) {

#if 0
  for(int r=0;r<(2*numInds-1);r++)
    for(int j=0;j<std::min(3,r);j++){
      double res =myComb2Tab[r][j];// myComb2(numInds,r,j);
      fprintf(stderr,"myComb\t(%d,%d,%d) =%f\n",numInds,r,j,res);
    }
  exit(0);  
#endif
  
  int Aa_offset = angsd::majorminor[minor_offset][major_offset];//0-9
  int AA_offset = angsd::majorminor[minor_offset][minor_offset];//0-9
  int aa_offset = angsd::majorminor[major_offset][major_offset];//0-9
  
  
  double hj[2*numInds+1];
  for(int index=0;index<(2*numInds+1);index++)
    if(underFlowProtect==0)
      hj[index]=0;
    else
      hj[index]=log(0);
  double PAA,PAa,Paa;
  //    numInds =5;
  for(int i=0 ; i<numInds ;i++) {
    double GAA,GAa,Gaa;
    GAA = liks[i*10+AA_offset];
    GAa = log(2.0)+liks[i*10+Aa_offset];
    Gaa = liks[i*10+aa_offset];
    
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
    }
    
      
    }//after recursion
    //if we are underflowprotecting ht hj is in logspace
    // print_array(stdout,hj,2*numInds+1,!underFlowProtect);
    for(int i=0;i<(2*numInds+1);i++){
      //fprintf(stdout,"BICO: %f\n",log(bico(2*numInds,i)));
      if(underFlowProtect)
	hj[i] =  (hj[i]-abcSaf::lbicoTab[i]);
      else
	hj[i] =  exp(log(hj[i])-abcSaf::lbicoTab[i]);
      
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
	GAA = liks[i*10+AA_offset];
	GAa = log(2.0)+liks[i*10+Aa_offset];
	Gaa = liks[i*10+aa_offset];
	
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
	double g;
	if(j==0)
	  g=liks[10*select+angsd::majorminor[major_offset][major_offset]];
	else if(j==1)
	  g=liks[10*select+angsd::majorminor[major_offset][minor_offset]];
	else
	  g=liks[10*select+angsd::majorminor[minor_offset][minor_offset]];
	
	double tmp=0;
	if(underFlowProtect)
	  tmp = log(tmp);
	for(int r=0;r<2*numInds-1;r++){
	  //	  fprintf(stderr,"mycomb2:%f\tpes1: %f\tpes2: %f\n",myComb2(numInds,r+j,j),pest[r+j],pest[2*numInds-r-j]);
	  if(underFlowProtect)
	    //tmp = angsd::addProtect2(tmp, log(myComb2(numInds,r+j,j)) + (asdf[r])+pest[r+j]+pest[2*numInds-r-j]);
	    tmp = angsd::addProtect2(tmp, log(abcSaf::myComb2Tab[r+j][j]) + (asdf[r])+pest[r+j]+pest[2*numInds-r-j]);
	  else{
	    //fprintf(stderr,"r+j:%d\tj:%d\nold:%f\nnew:%f",r+j,j,myComb2(numInds,r+j,j), myComb2Tab[r+j][j]);
	    //exit(0);
	    //tmp += myComb2(numInds,r+j,j)*(asdf[r])*(exp(pest[r+j])+exp(pest[2*numInds-r-j]));
	    tmp += abcSaf::myComb2Tab[r+j][j]*(asdf[r])*(exp(pest[r+j])+exp(pest[2*numInds-r-j]));
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
	return 1;
      }
      


      //      print_array(stdout,res,3,0);      
      
      //logrescale(res,3);
      double mySum=exp(res[0])+exp(res[1])+exp(res[2]);
      for(int i=0;i<3;i++)
	res[i] =exp(res[i])/mySum;
      
      //print_array(sfsfile,res,3,0);      
      delete [] hj;
    }//after select loop
    return 0;
}
