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
#include "abcFreq.h"

int abcFreq::emIter = EM_NITER;
double abcFreq::EM_start = EM_START;
double *abcFreq::indF = NULL;

double abcFreq::to_pval(Chisqdist *chisq,double f){
  return f<0?1:1-chisq->cdf(f);
}

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
  fprintf(argFile,"NB These frequency estimators requires major/minor -doMajorMinor\n");
}

//fancy little function
int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

 
void abcFreq::getOptions(argStruct *arguments){
  int inputtype = arguments->inputtype;
  

  rmTriallelic=angsd::getArg("-rmTriallelic",rmTriallelic,arguments);

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
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
      if(rmTriallelic)
	SNP_pval_tri = chisq2->invcdf(1-rmTriallelic);
    }
    doSNP =1 ;
    fprintf(stderr,"\t-> SNP-filter using a pvalue: %e correspond to %f likelihood units\n",pre,SNP_pval);
  }

  refName = angsd::getArg("-ref",refName,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);

  if(abs(doMaf)&& !isPowerOfTwo((unsigned int) abs(doMaf))){
    fprintf(stderr,"\n[%s] You have selected filters for maf/lrt\n",__FILE__);
    fprintf(stderr,"If you have selected more than one MAF estimator we will choose in following order\n");
    fprintf(stderr,"\t1. knownminor EM\n");
    fprintf(stderr,"\t2. unknownminor EM\n");
    fprintf(stderr,"\t3. Posterior maf\n");
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
  chisq1=chisq2=chisq3=NULL;
  inputIsBeagle =0;
  beagleProb = 0; //<-output for beagleprobs
  minMaf =-1.0;
  SNP_pval = 1;
  nInd = arguments->nInd;
  eps = 0.001;
  outfileZ2 = Z_NULL;
  outfileZ = Z_NULL;
  indFname = NULL;
  doMaf=0;
  rmTriallelic=0;
  GL=0;
  doSNP=0;
  doPost=0;

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
  printArg(arguments->argumentFile);
  if(doMaf==0 && doPost==0){
    shouldRun[index]=0;
    return;
  }
  if(doMaf>0){
    //make output files
    const char* postfix;
    postfix=".mafs.gz";
    outfileZ = aio::openFileGz(outfiles,postfix,GZOPT);
    if(beagleProb){
      postfix=".beagle.gprobs.gz";
      outfileZ2 = aio::openFileGz(outfiles,postfix,GZOPT);
    }
  }else
    doMaf=abs(doMaf);
  //print header
  kstring_t bufstr;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  kputs("chromo\tposition\tmajor\tminor\t",&bufstr);
  if(refName!=NULL)
    kputs("ref\t",&bufstr);
  if(ancName)// inputtyp=3 is tglf, if tglf then -posi has been supplied
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
  gzwrite(outfileZ,bufstr.s,bufstr.l);
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
    gzwrite(outfileZ2,bufstr.s,bufstr.l);
  }

  free(bufstr.s);
}


abcFreq::~abcFreq(){
  if(outfileZ!=Z_NULL){
    fprintf(stderr,"Closing file1 Z_NULL\n");
    gzclose(outfileZ);
    outfileZ=Z_NULL;
  }
  if(outfileZ2!=Z_NULL){
    fprintf(stderr,"Closing file2 Z_NULL\n");
    gzclose(outfileZ2);
    outfileZ2=Z_NULL;
  }
  free(refName);
  free(ancName);
  delete [] indF;
  delete chisq1;
  delete chisq2;
  delete chisq3;
}



void abcFreq::print(funkyPars *pars) {
  if(outfileZ==Z_NULL&&outfileZ2==Z_NULL)
    return;
  kstring_t bufstr;
  bufstr.s=NULL; bufstr.l=bufstr.m=0;

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
	ksprintf(&bufstr,"%e\t",to_pval(chisq1,freq->lrt_EM[s]));
      if(doMaf &2)
	ksprintf(&bufstr,"%e\t",to_pval(chisq1,freq->lrt_EM_unknown[s]));
    }

    kputw(pars->keepSites[s],&bufstr);kputc('\n',&bufstr);
  }

  gzwrite(outfileZ,bufstr.s,bufstr.l);  
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
    int ret=gzwrite(outfileZ2,bufstr.s,bufstr.l);
    //fprintf(stderr,"ret.l:%d bufstr.l:%zu\n",ret,bufstr.l);
    bufstr.l=0;
  }
  free(bufstr.s);
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
      }if(doPost==1)  //maf prior
	 make_post(like[s],post[s],freq->freq[s],pars->nInd);
      else if(doPost==2){//uniform prior
	for(int i=0;i<pars->nInd;i++){
	  double norm = angsd::addProtect3(like[s][i*3+0],like[s][i*3+1],like[s][i*3+2]);
	  for(int g=0;g<3;g++)
	    post[s][i*3+g]=exp(like[s][i*3+g]-norm);
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
  //fprintf(stderr,"keepSites[0] %d\n",pars->keepSites[0]);
  for(int s=0;s<pars->numSites;s++) {
    if(keepInd[s]==0)//if we dont have any information
      continue;
    keepInd[s]=0;//
    for(int i=0 ; i<pars->nInd ;i++) {//DRAGON CHECK THIS
      //fprintf(stderr,"size %d\nind %d\t loglike:%f\t%f\t%f\n",s,i,loglike[s][i*3+0],loglike[s][i*3+1],loglike[s][i*3+2]);
      keepList[i]=1;
      if(loglike[s][i*3+0]+loglike[s][i*3+1]+loglike[s][i*3+2]>-0.0001){
	//	fprintf(stderr,"size %d\nind %d\t loglike:%f\t%f\t%f\n",s,i,loglike[s][i*3+0],loglike[s][i*3+1],loglike[s][i*3+2]);
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

    for(int ii=0;1&&ii<numInds;ii++){
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
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    exit(0);
  }
  
  return(p);
}

