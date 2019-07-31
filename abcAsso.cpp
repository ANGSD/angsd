/*
  emil emil.jorsboe@bio.ku.dk jan9 2019
  has implemented new methods -doAsso 4, 5 and 6

  thorfinn thorfinn@binf.ku.dk dec17 2012 
  has modified so no sites[].chromo etc are used
  
  anders albrecht@binf.ku.dk made this.

  logistic regression changed 612new/613

  part of angsd

*/

#include <cmath>
#include <htslib/kstring.h>
#include "shared.h"
#include "analysisFunction.h"

#include "abcFreq.h"
#include "abcAsso.h"

void abcAsso::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAsso\t%d\n",doAsso);
  fprintf(argFile,"\t1: Frequency Test (Known Major and Minor)\n");
  fprintf(argFile,"\t2: Score Test\n");
  fprintf(argFile,"\t4: Latent genotype model\n");
  fprintf(argFile,"\t5: Score Test with latent genotype model - hybrid test\n");
  fprintf(argFile,"\t6: Dosage regression\n");
  fprintf(argFile,"\t7: Latent genotype model (wald test)\n");
  fprintf(argFile,"  Frequency Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\t\n\n",yfile1);
  fprintf(argFile,"  Score, Latent, Hybrid and Dosage Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\n",yfile1);

  fprintf(argFile,"\t-yCount\t\t%s\t(File containing count phenotypes)\n",yfile2);
  fprintf(argFile,"\t-yQuant\t\t%s\t(File containing phenotypes)\n",yfile3);
  fprintf(argFile,"\t-cov\t\t%s\t(File containing additional covariates)\n",covfile);
  fprintf(argFile,"\t-model\t%d\n",model);
  fprintf(argFile,"\t1: Additive/Log-Additive (Default)\n");
  fprintf(argFile,"\t2: Dominant\n");
  fprintf(argFile,"\t3: Recessive\n\n");
  fprintf(argFile,"\t-minHigh\t%d\t(Require atleast minHigh number of high credible genotypes)\n",minHigh);
  fprintf(argFile,"\t-minCount\t%d\t(Require this number of minor alleles, estimated from MAF)\n",minCount);

  fprintf(argFile,"\t-assoThres\t%f\tThreshold for logistic regression\n",assoThres);
  fprintf(argFile,"\t-assoIter\t%d\tNumber of iterations for logistic regression\n",assoIter);
  fprintf(argFile,"\t-emThres\t%f\tThreshold for convergence of EM algorithm in doAsso 4 and 5\n",emThres);
  fprintf(argFile,"\t-emIter\t%d\tNumber of max iterations for EM algorithm in doAsso 4 and 5\n\n",emIter);
  fprintf(argFile,"\t-doPriming\t%d\tPrime EM algorithm with dosage derived coefficients (0: no, 1: yes - default) \n\n",doPriming);
  
  fprintf(argFile,"  Hybrid Test Options:\n");
  fprintf(argFile,"\t-hybridThres\t\t%f\t(p-value value threshold for when to perform latent genotype model)\n",hybridThres);
  fprintf(argFile,"Examples:\n\tPerform Frequency Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 1 -GL 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\tPerform Score Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 2 -GL 1 -doPost 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\n");
}


void abcAsso::getOptions(argStruct *arguments){


  doAsso=angsd::getArg("-doAsso",doAsso,arguments);
  if(doAsso==3){
    fprintf(stderr,"\t-> -doAsso 3 is deprecated from version 0.615 \n");
    exit(0);
  }
  if(doAsso==0)
    return;
  
  doMaf=angsd::getArg("-doMaf",doMaf,arguments);

  adjust=angsd::getArg("-adjust",adjust,arguments);
  model=angsd::getArg("-model",model,arguments);
  minCov=angsd::getArg("-minCov",minCov,arguments);
  dynCov=angsd::getArg("-dynCov",dynCov,arguments);
  minHigh=angsd::getArg("-minHigh",minHigh,arguments);
  doPrint=angsd::getArg("-doPrint",doPrint,arguments);
  minCount=angsd::getArg("-minCount",minCount,arguments);
  sitePerm=angsd::getArg("-sitePerm",sitePerm,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  covfile=angsd::getArg("-cov",covfile,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  hybridThres=angsd::getArg("-hybridThres",hybridThres,arguments);
  assoThres=angsd::getArg("-assoThres",assoThres,arguments);
  assoIter=angsd::getArg("-assoIter",assoIter,arguments);
  emThres=angsd::getArg("-emThres",emThres,arguments);
  emIter=angsd::getArg("-emIter",emIter,arguments);
  doPriming=angsd::getArg("-doPriming",doPriming,arguments);

  yfile1=angsd::getArg("-yBin",yfile1,arguments);  
  if(yfile1!=NULL)
    isBinary=1;

  yfile2=angsd::getArg("-yCount",yfile2,arguments);  
  if(yfile2!=NULL)
    isCount=1;    
    
  yfile3=angsd::getArg("-yQuant",yfile3,arguments);
  if(yfile3!=NULL)
    isQuant=1;    

  
  if(doPrint)
    fprintf(stderr,"finished [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doAsso && doMaf==0){
    fprintf(stderr,"Error: you must estimate the maf (-doMaf) in order to perform association \n");
    exit(0);

  }

  if(doAsso && yfile1==NULL && yfile2==NULL && yfile3==NULL){
    fprintf(stderr,"Error: you must provide a phenotype file (-yBin or -yQuant) to perform association \n");
    exit(0);
   }
 if(doAsso && arguments->inputtype==INPUT_BEAGLE&&doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on posterior input\n");
    exit(0);
  }
  if(doAsso && arguments->inputtype!=INPUT_BEAGLE&&(doAsso==2)&&doPost==0){
    fprintf(stderr,"Error: For doAsso=2 you must estimate the posterior probabilites for the genotypes (-doPost 1) \n");
    exit(0);
  }
  if((isBinary+isCount+isQuant)>1){
    fprintf(stderr,"Error: Only one of -yBin, yQuant and -yCount \n");
    exit(0);
  }

}
 
abcAsso::abcAsso(const char *outfiles,argStruct *arguments,int inputtype){
  multiOutfile=NULL;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  model=1;
  doPrint=0;
  doAsso=0;
  GL=0;
  doPost=0;
  isBinary=0;
  isCount=0;
  isQuant=0;
  sitePerm=0;//not for users
  covfile=NULL;
  yfile1=NULL;
  yfile2=NULL;
  yfile3=NULL;
  minHigh=10;
  minCount=10;
  dynCov=0;//not for users
  minCov=5;//not for users
  adjust=1;//not for users  
  doMaf=0;  
  assoIter=100;
  assoThres=1e-06;
  emIter=40;
  emThres=1e-04;  
  hybridThres=0.05;
  doPriming=1;
  //from command line
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAsso")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);

  if(doAsso==0){
    shouldRun[index] =0;
    return;
  }

  printArg(arguments->argumentFile);

  //read phenotype - binary or count (poisson)
  if(isBinary)
    ymat = angsd::getMatrix(yfile1,1,100000);  
  else if(isCount)
    ymat = angsd::getMatrix(yfile2,1,100000);  
  else
    ymat = angsd::getMatrix(yfile3,0,100000);

  //read covariates 
  if(covfile!=NULL)
    covmat = angsd::getMatrix(covfile,0,100000);
  else{
    covmat.x=0;
    covmat.y=0;
    covmat.matrix=NULL;
  }
  
  if(covfile!=NULL&&(covmat.x!=ymat.x)){
    fprintf(stderr,"The number of covariates (%d) does not match the number of phenotypes (%d)\n",covmat.x,ymat.x);
    exit(0);
  }
 
  if(!isBinary && doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on quantitative traits\n");
    exit(0);
  }

  // check cov and ymat
  check_pars(covmat,ymat,isBinary);

  //make output files  
  multiOutfile = new BGZF*[ymat.y];
  const char* postfix;
  postfix=".lrt";

  if(covfile!=NULL&&minCov>0){
    int keepList[ymat.x];
    for(int i=0 ; i < ymat.x;i++) {
      keepList[i]=1;
      for(int yi=0;yi<ymat.y;yi++) {
	if(ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
      }
      for(int ci=0;ci<covmat.y;ci++) {
	if(covmat.matrix[i][ci]==-999)
	  keepList[i]=0;
      }
    }
    int nCov=0;
    int count[covmat.y];
    for(int ci=0;ci<covmat.y;ci++) {
      count[ci]=0;
      for(int i=0 ; i < ymat.x;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][ci]!=0){
	  count[ci]++;
	}
      }
  
      if(count[ci]<minCov){
	fprintf(stderr,"Error: Cov #%d only has %d non zero entries\n",ci,count[ci]);
      }
      else
	nCov++;
      
    }
    if(!dynCov&&covmat.y!=nCov){
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      exit(0);

    }
    else if(covmat.y!=nCov){
      //      angsd::printMatrix(covmat,stderr);
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      angsd::Matrix<double> newmat;
      newmat.x=covmat.x;
      newmat.y=nCov;
      newmat.matrix=new double*[covmat.x];
      for(int xi=0;xi<covmat.x;xi++){
	newmat.matrix[xi] = new double[nCov];
	int tempCount=0;
	for(int ci=0;ci<covmat.y;ci++){
	  if(count[ci]>minCov){
	    newmat.matrix[xi][tempCount]=covmat.matrix[xi][ci];
	    tempCount++;
	  }
	}
      }
      angsd::deleteMatrix(covmat);
      covmat=newmat;
      
    }
  }

  //open outfiles
  for(int i=0;i<ymat.y;i++){
    char ary[5000];
    snprintf(ary,5000,"%s%d.gz",postfix,i);
    multiOutfile[i] = NULL;
    multiOutfile[i] = aio::openFileBG(outfiles,ary);
  }

  //print header
  bufstr.l=0;


  // hybrid model score model first, then if significant use EM asso  
  if(doAsso==5){
    ksprintf(&bufstr,"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRTscore\thigh_WT/HE/HO\tLRTem\tbeta\tSE\temIter\n");
    //EM asso or dosage model
  } else if(doAsso==4 or doAsso==7){
    ksprintf(&bufstr,"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\tbeta\tSE\thigh_WT/HE/HO\temIter\n");
  } else if(doAsso==6){
    ksprintf(&bufstr,"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\tbeta\tSE\thigh_WT/HE/HO\n");    
  } else if(doAsso==2)
    ksprintf(&bufstr,"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\thigh_WT/HE/HO\n");
  else
    ksprintf(&bufstr,"Chromosome\tPosition\tMajor\tMinor\tFrequency\tLRT\n");
  for(int yi=0;yi<ymat.y;yi++)
    aio::bgzf_write(multiOutfile[yi],bufstr.s,bufstr.l);
  bufstr.l=0;
}

// desctructor
abcAsso::~abcAsso(){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doAsso==0)
    return;
  for(int i=0;i<ymat.y;i++)
    if(multiOutfile[i]!=NULL)
      bgzf_close(multiOutfile[i]);
  delete [] multiOutfile;

  if(covfile!=NULL)
    angsd::deleteMatrix(covmat);
  angsd::deleteMatrix(ymat);

}


void abcAsso::clean(funkyPars *pars){
  if(doAsso==0)
    return;

  assoStruct *assoc =(assoStruct*) pars->extras[index];

  if(assoc->stat!=NULL)
    for(int yi=0;yi<ymat.y;yi++)
      delete[] assoc->stat[yi];

  delete[] assoc->stat;

  if(doAsso==5){
    delete[] assoc->highWt;
    delete[] assoc->highHe;
    delete[] assoc->highHo;
    delete[] assoc->betas;
    delete[] assoc->SEs;
    for(int yi=0;yi<ymat.y;yi++)
      delete[] assoc->statOther[yi];
    delete[] assoc->statOther;
  } else if(doAsso==2){
    delete[] assoc->highWt;
    delete[] assoc->highHe;
    delete[] assoc->highHo;    
  } else if(doAsso==4 or doAsso==6 or doAsso==7){
    delete[] assoc->highWt;
    delete[] assoc->highHe;
    delete[] assoc->highHo;    
    delete[] assoc->betas;
    delete[] assoc->SEs;
  }

  if(doAsso==4 or doAsso==5 or doAsso==7){
    delete[] assoc->emIter;
  }
  
  if(assoc->keepInd!=NULL)
    for( int yi =0;yi<ymat.y;yi++)
      delete[] assoc->keepInd[yi];
  delete[] assoc->keepInd;
  
  delete assoc;
}



void abcAsso::print(funkyPars *pars){
  if(doAsso==0)
    return;

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  printDoAsso(pars);
}

assoStruct *allocAssoStruct(){
  assoStruct *assoc = new assoStruct;

  assoc->SEs=NULL;
  assoc->betas=NULL;
  assoc->stat=NULL;
  assoc->statOther=NULL; 
  assoc->keepInd=NULL;
  assoc->highWt = NULL;
  assoc->highHe = NULL;
  assoc->highHo = NULL;
  assoc->emIter = NULL;
  
  return assoc;
}


void abcAsso::run(funkyPars *pars){

  if(doAsso==0)
    return;
  
  assoStruct *assoc = allocAssoStruct();
  pars->extras[index] = assoc;

  if(doAsso==2 or doAsso==4 or doAsso==5 or doAsso==6 or doAsso==7){
    assoc->highWt=new int[pars->numSites];
    assoc->highHe=new int[pars->numSites];
    assoc->highHo=new int[pars->numSites];
  }
  
  if(doAsso==1){
    frequencyAsso(pars,assoc);
  } else if(doAsso==2){
    scoreAsso(pars,assoc);
  } else if(doAsso==4){
    assoc->betas=new double[pars->numSites];
    assoc->SEs=new double[pars->numSites];
    assoc->emIter=new int[pars->numSites];
    emAsso(pars,assoc);    
  } else if(doAsso==6){
    assoc->betas=new double[pars->numSites];
    assoc->SEs=new double[pars->numSites];
    dosageAsso(pars,assoc);    
  } else if(doAsso==5){
    assoc->betas=new double[pars->numSites];
    assoc->SEs=new double[pars->numSites];
    assoc->emIter=new int[pars->numSites];
    hybridAsso(pars,assoc);
  } else if(doAsso==7){
    assoc->betas=new double[pars->numSites];
    assoc->SEs=new double[pars->numSites];
    assoc->emIter=new int[pars->numSites];
    emAssoWald(pars,assoc);    
  }


}

void abcAsso::check_pars(angsd::Matrix<double> &cov, angsd::Matrix<double> &phe, int isBinary){

  // check if all values of first column in cov file is 1
  int interceptChecker = 0;
  // check if phenotype 0 or 1 - for quantitative
  int isBinaryQuan = 1;
  for(int j=0;j<phe.y;j++){
    for(int i=0;i<phe.x;i++){      
      //if logistic regression check if phenotypes are 0 or 1
      if(isBinary){
	isBinaryQuan = 0;

	if(phe.matrix[i][j]!=0 and phe.matrix[i][j]!=1 and phe.matrix[i][j]!=-999){
	  fprintf(stderr,"Phenotypes are not binary (0 or 1) for logistic model!\n");
	  exit(1);
	}

      } else if(isCount){
	isBinaryQuan = 0;
	
	if(ceil(phe.matrix[i][j])!=floor(phe.matrix[i][j])){
	  fprintf(stderr,"Phenotypes are not count data (must be integer) for poisson model!\n");
	  exit(1);
	}

      } else{
	if(phe.matrix[i][j]!=0 and phe.matrix[i][j]!=1){
	  isBinaryQuan = 0;
	}
      }
    }
    
    if(isBinaryQuan){
      fprintf(stderr,"\n");
      fprintf(stderr,"#######################################################################################\n");
      fprintf(stderr,"#######################################################################################\n");
      fprintf(stderr,"WARNING: Phenotype nr. %i for linear model appears to be binary, run a logistic model instead!\n",j);
      fprintf(stderr,"#######################################################################################\n");
      fprintf(stderr,"#######################################################################################\n");
      fprintf(stderr,"\n");
    }
    isBinaryQuan=1;
    
  }

  for(int j=0;j<cov.y;j++){
    for(int i=0;i<cov.x;i++){

      if(cov.matrix[i][j]==1){
	interceptChecker++;
      }
    }

    // p->len is equal to number of indis in plink file
    // if these are same length, means user has column of 1s in .cov file
    if(interceptChecker==cov.x){
      fprintf(stderr,"Column of 1s for intercept should not be included in covariate,\n");
      fprintf(stderr,"as this will be added automatically.\n");
      exit(1);
    }

    interceptChecker=0;
  }
 
}


			 
void abcAsso::frequencyAsso(funkyPars  *pars,assoStruct *assoc){

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);

    fflush(stderr);
    exit(0);
  }

  if(ymat.y!=1){
    fprintf(stderr,"Only one phenotype allowed for -doAsso 1\n");
    fflush(stderr);
    exit(0);
  }

  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++)
    stat[yi] = new double[pars->numSites];
    
  int y[pars->nInd];
  for(int i=0;i<pars->nInd;i++)
    y[i]=ymat.matrix[i][0];

  double **like0;//genotype likelihood for controls
  double **like1;//genotype likelihood for cases
  double **likeAll;//genotype likelihood for cases and controls
  int Ncases=0; //number of cases
  int Ncontrols=0; //number of cases
  int Nall=0; //number of cases and controls
  int cases[pars->nInd];
  int controls[pars->nInd];
  int all[pars->nInd];

  for(int i=0;i<pars->nInd;i++){
    cases[i]=0;
    controls[i]=0;
    all[i]=0;
    if((int)y[i]==1){
      Ncases++;
      cases[i]=1;
      Nall++;
      all[i]=1;
    }
    if((int)y[i]==0){
      Ncontrols++;
      controls[i]=1;
      Nall++;
      all[i]=1;
    }
    if(doPrint)
      fprintf(stderr,"all, case, control: %d %d %d\n",all[i],cases[i],controls[i]);
  }
  if(doPrint)
    fprintf(stderr,"count complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doAsso==1){
    like0=angsd::get3likesRMlow(pars,controls);
    like1=angsd::get3likesRMlow(pars,cases);
    likeAll=angsd::get3likesRMlow(pars,all);
  }

 if(doPrint)
    fprintf(stderr,"like complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  for(int s=0;s<pars->numSites;s++){//loop overs sites
    stat[0][s]=-999;
    if(pars->keepSites[s]==0)
      continue;
    if(doAsso==1){

      if(doPrint)
	fprintf(stderr,"do freq [%s]\t[%s]\n",__FILE__,__FUNCTION__);

      double score0=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(like0[s],Ncontrols,NULL,Ncontrols),like0[s],Ncontrols);
      //likelihood for the cases
      double score1=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(like1[s],Ncases,NULL,Ncases),like1[s],Ncases);
      //likelihood for all individuals
      double scoreNull=abcFreq::likeFixedMinor(abcFreq::emFrequency_ext(likeAll[s],Nall,NULL,Nall),likeAll[s],Nall);
      //likelhood ratio statistics \sim chi^2
      double LRT=-2*(score0+score1-scoreNull);
      stat[0][s]=LRT;
    }
  }
  assoc->stat=stat;
  for(int s=0;s<pars->numSites;s++){
    delete[] like0[s];
    delete[] like1[s];
    delete[] likeAll[s];
  }
   delete[] like0;
   delete[] like1;
   delete[] likeAll;

 if(doPrint)
    fprintf(stderr,"finish [%s]\t[%s]\n",__FILE__,__FUNCTION__);

}

void abcAsso::scoreAsso(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }

  int **keepInd  = new int*[ymat.y];
  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
        
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;

    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }

	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }  
      double *y = new double[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 
      freqStruct *freq = (freqStruct *) pars->extras[6];
      stat[yi][s]=doAssociation(pars,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc);
      
      //cleanup
      delete [] y;
      delete [] keepList;
      
    } //phenotypes end
    
    delete [] keepListAll;
  } // sites end

  assoc->stat=stat;
  assoc->keepInd=keepInd;
}


///////////////////////////////////////////////////////////////////////////////
// INSERTED BY EMIL 
///////////////////////////////////////////////////////////////////////////////

double abcAsso::standardError(double* start, angsd::Matrix<double> *design, angsd::Matrix<double> *postAll, double *y, double *post, int isBinary, int isCount,int *keepList, int nInd){

  //problem is that not all y values might be included since some might have missing information
  std::vector<double> ynew(design->x);         
  
  //only need individuals with valid phenotypic information
  int count=0;
  for(int i=0;i<nInd;i++){
    if(keepList[i]){      
      ynew[count]=y[i];
      count++;
    }
  }
  
  if(postAll!=NULL){

    double ret;
    for(int i=0;i<design->x;i++){
      double m = 0;
      for(int j=0;j<design->y;j++){
	m += design->matrix[i][j]*start[j];
      }

      // because y is indis long, each entry has to be used 3 times - design has 3*indis rows
      size_t indexy = (size_t)floor(i/3);
      
      if(isBinary){
		
	double prob = exp(m)/(exp(m)+1.0);
	// now handling if p is 0 or 1 (makes it very small or large)
	m = angsd::bernoulli(ynew[indexy],prob,1);	

      } else if(isCount){
	double lambda = exp(m);
	m = angsd::poisson(ynew[indexy],lambda,1);
	
      } else{
	// density function of normal distribution
	m = angsd::dnorm(ynew[indexy],m,start[design->y],1);
      }
      
      double tmp = m + log(post[i]);
      postAll->matrix[(size_t)floor(i/3)][i % 3] = tmp;
    }
    
    // x is number of indis, y is 3
    postAll->x = (size_t) design->x/3;
    postAll->y = 3;

    // to get rowSums
    std::vector<double> postTmp(postAll->x);  
  
    for(int i=0;i<postAll->x;i++){
      double tmp = 0;
      double maxval = postAll->matrix[i][0];
      for(int j=0;j<postAll->y;j++){
	// find max - part of trick for doing log(p1+p2+p3)
	maxval = std::max(maxval,postAll->matrix[i][j]);     
      }
      // trick to avoid over/underflow - log(exp(log(val1)-log(max)) + ...) + log(max) = (exp(log(val1))/exp(log(max)))*(max) + ...
      // same (exp(log(val1))/exp(log(max)))*(max)
      postTmp[i] = log(exp(postAll->matrix[i][0]-maxval)+exp(postAll->matrix[i][1]-maxval)+exp(postAll->matrix[i][2]-maxval)) + maxval;  
    }
  
    // divide each entry of a row with sum of that row
    int df = postAll->x - design->y;
    for(int i=0;i<postAll->x;i++){
      double tmp = 0;
      for(int j=0;j<postAll->y;j++){          
	postAll->matrix[i][j] -= postTmp[i];
      }    
    }
    
    //need to flatten the weights, which is p(s|y,G,phi,Q,f)
    // so first four values is for first indi, and so on...
    double weights[postAll->x*postAll->y];
    int a = 0;
    for(int i=0;i<postAll->x;i++)
      for(int j=0;j<postAll->y;j++){
	weights[a++] = exp(postAll->matrix[i][j]);
	//check if issue with weights
	if(exp(postAll->matrix[i][j])!=exp(postAll->matrix[i][j]) or std::isinf(exp(postAll->matrix[i][j]))){
	  fprintf(stderr,"Issue with weights being nan or inf\n");
	  return(-9);
	}
      }
    
    std::vector<double> yTilde(design->x);         
    
    for(int i=0;i<design->x;i++)
      for(int x=0;x<design->y;x++){
	yTilde[i] += design->matrix[i][x]*start[x];
      }

    // has to be expected value - mu, therefore using link function exp/(1+exp)
    if(isBinary){
      for(int i=0;i<design->x;i++)
	yTilde[i] = exp(yTilde[i]) / (1+exp(yTilde[i]));
    }

    // link function exp
    if(isCount){
      for(int i=0;i<design->x;i++)
	yTilde[i] = exp(yTilde[i]);
    }
    
    int count=0;

    // this has colSums for each col for each indi (3 rows)
    std::vector<double> summedDesign((design->x/3)*design->y);   
    std::vector<double> tmp(design->x);
        
    for(int i=0;i<design->x;i++){
      // because y is indis long, each entry has to be used 3 times - design has 3*indis rows
      size_t indexy = (size_t)floor(i/3);      
      // getting the residuals
        
      if(isBinary || isCount){
	tmp[i] = (ynew[indexy]-yTilde[i]) * weights[i];
      } else{
	tmp[i] = ((ynew[indexy]-yTilde[i]) / (start[design->y]*start[design->y]))*weights[i];
      }    
      
      for(int x=0;x<design->y;x++){                  
	if(i % 3 == 2){
	  // for each col take 3 rows for indis and take sum
	  summedDesign[count*design->y+x]=design->matrix[i-2][x]*tmp[i-2]+design->matrix[i-1][x]*tmp[i-1]+design->matrix[i][x]*tmp[i];	 
	}	
      }
      if(i % 3 == 2){
	count++;
      }      
    }
    
    double invSummedDesign[design->y*design->y];
    for(int i=0;i<design->y*design->y;i++)
      invSummedDesign[i]=0;
    
    // this is doing the matrix product of (X)^T*W*X 
    for(int x=0;x<design->y;x++)
      for(int y=0;y<design->y;y++)
	for(int i=0;i<(int)floor(design->x/3);i++){
	  invSummedDesign[x*design->y+y]+=summedDesign[i*design->y+x]*summedDesign[i*design->y+y];
	}

    // doing inv((X)^T*W*X)
    int singular=angsd::svd_inverse(invSummedDesign, design->y, design->y);
    
    return(sqrt(invSummedDesign[0]));
    
  } else{
    
    std::vector<double> yTilde(design->x);         
    
    for(int i=0;i<design->x;i++)
      for(int x=0;x<design->y;x++){
	yTilde[i] += design->matrix[i][x]*start[x];
      }

    // has to be expected value - mu, therefore using link function exp/(1+exp)
    if(isBinary){
      for(int i=0;i<design->x;i++)
	yTilde[i] = exp(yTilde[i]) / (1+exp(yTilde[i]));
    }

    // link function exp
    if(isCount){
      for(int i=0;i<design->x;i++)
	yTilde[i] = exp(yTilde[i]);
    }        

    std::vector<double> summedDesign(design->x*design->y); 
    std::vector<double> tmp(design->x);

    //design has 3*indis rows
    for(int i=0;i<design->x;i++){

      //ynew is indis long, so this is to be able to use that each value 3 times
      size_t indexy = (size_t)floor(i/3);
      
      if(isBinary || isCount){
	tmp[i] = (ynew[indexy]-yTilde[i]);
      } else{
	tmp[i] = ((ynew[indexy]-yTilde[i]) / (start[design->y]*start[design->y]));
      }

      for(int x=0;x<design->y;x++){                  	
	summedDesign[i*design->y+x]=design->matrix[i][x]*tmp[i];	 		
      }
    }             
    
    double invSummedDesign[design->y*design->y];
    for(int i=0;i<design->y*design->y;i++)
      invSummedDesign[i]=0;
    
    // this is doing the matrix product of (X)^T*W*X 
    for(int x=0;x<design->y;x++)
      for(int y=0;y<design->y;y++)
	for(int i=0;i<design->x;i++){
	  invSummedDesign[x*design->y+y]+=summedDesign[i*design->y+x]*summedDesign[i*design->y+y];
	}

    // doing inv((X)^T*W*X)
    int singular=angsd::svd_inverse(invSummedDesign, design->y, design->y);
    
    return(sqrt(invSummedDesign[0]));
  }
}


double abcAsso::dosageAssoc(funkyPars *p,angsd::Matrix<double> *design,angsd::Matrix<double> *designNull,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s,assoStruct *assoc,int model, int isBinary, int isCount, double* start, int fullModel){

  int maf0=1;

  // if too rare do not run analysis
  if(freq>0.995||freq<0.005){
    maf0=0;
  }

  int highWT=0;
  int highHE=0;
  int highHO=0;
 
  double covMatrix[(covmat.y+2)*keepInd];
  double covMatrixNull[(covmat.y+1)*keepInd];
  double yfit[keepInd];
    
  double yNull[keepInd];
  int count = 0;
  double post[keepInd*3];
  
  designNull->x=keepInd;
  //covars + intercept
  designNull->y=covmat.y+1;
     
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      yNull[count]=yOrg[i];
      post[count*3+0]=postOrg[i*3+0];
      post[count*3+1]=postOrg[i*3+1];
      post[count*3+2]=postOrg[i*3+2];

      covMatrixNull[count] = 1;
      
      designNull->matrix[count][0] = 1;
      for(int c=0;c<covmat.y;c++){
	designNull->matrix[count][c+1] = covmat.matrix[i][c];
	covMatrixNull[(c+1)*keepInd+count] = covmat.matrix[i][c];
      }
      count++;    
    }
  }

  // need intercept 
  double startNull[covmat.y+2];

  for(int i=0;i<(covmat.y+2);i++){
    startNull[i] = drand48()*2-1;
  }

  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }
  
  if(isBinary){
    getFitBin(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  } else if(isCount){
    getFitPois(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  } else{      
    getFit(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  }

  double llhNull = logLike(startNull,yNull,designNull,post,isBinary,isCount,0);
    
  start[0] = drand48()*2-1;

  double y[keepInd];
  count=0;
      
  design->x=keepInd;
  design->y=covmat.y+2;
  
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      y[count]=yOrg[i];      
      // 1: add, 2: dom, 3: rec
      if(model==3){
	design->matrix[count][0] = post[count*3+2];
	covMatrix[count] = post[count*3+2];
      } else if(model==2){
	design->matrix[count][0] = post[count*3+2] + post[count*3+1];
	covMatrix[count] = post[count*3+2] + post[count*3+1];
      } else{
	// getting dosage of the minor allele or second allele of beagle file
	design->matrix[count][0] = 2*post[count*3+2] + post[count*3+1];
	covMatrix[count] = 2*post[count*3+2] + post[count*3+1];
      }

       //if rec model, WT+HE->WT, HO->HE
      if(model==3){
if((post[count*3+0]+post[count*3+1])>0.90)
	  highWT++;
	if(post[count*3+2]>0.90)
	  highHE++;
      } else if(model==2){
	 //if dom model, WT->WT, HE+HO->HE, HO->0
	if(post[count*3+0]>0.90)
	  highWT++;
	if((post[count*3+1]+post[count*3+2])>0.90)
	  highHE++;
      } else{
	if(post[count*3+0]>0.90)
	  highWT++;
	if(post[count*3+1]>0.90)
	  highHE++;
	if(post[count*3+2]>0.90)
	  highHO++;
      }
           
      covMatrix[keepInd+count] = 1;
      design->matrix[count][1] = 1;
      for(int c=0;c<covmat.y;c++){
	// design matrix has 1 rows per indi with dosage
	design->matrix[count][2+c] = covmat.matrix[i][c];
	covMatrix[(c+2)*keepInd+count] = covmat.matrix[i][c];
      }
      count++;
    }

    
  }
  
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
  
  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;
  
  if(nGeno<2)
    return(-999);//set_snan(lrt);

  //basically the same, just renamed for normScoreEnv and binomScoreEnv functions
  int numInds = keepInd;
  //freq*numInds*2 is the expected number of minor alleles
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
    return(-999);//set_snan(lrt);      set_snan(lrt);
  
  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }
   
  if(isBinary){
    getFitBin(yfit,y,covMatrix,keepInd,design->y,start);
  } else if(isCount){
    getFitPois(yfit,y,covMatrix,keepInd,design->y,start);
  } else{
    getFit(yfit,y,covMatrix,keepInd,design->y,start);  
  }
 
  double llh = logLike(start,y,design,post,isBinary,isCount,0);

  // likelihood ratio - chi square distributed according to Wilk's theorem
  double LRT = -2*(llh-llhNull);

  return(LRT);  
}
 

void abcAsso::dosageAsso(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }

  // covariates + geno + intercept + sd(y)
  double *start = new double[covmat.y+3];
  int **keepInd  = new int*[ymat.y];
  // we can also get coef of genotype
  double **stat = new double*[ymat.y*2];
 
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  angsd::Matrix<double> designNull;
  designNull.x=pars->nInd;
  //covars + intercept
  designNull.y=covmat.y+1;
  designNull.matrix=new double*[pars->nInd];
    
  angsd::Matrix<double> design; 
  design.x=3*pars->nInd;
  design.y=covmat.y+2;
  design.matrix=new double*[pars->nInd];
   
  for(int xi=0;xi<pars->nInd;xi++){    
    designNull.matrix[xi] = new double[covmat.y+1];
    design.matrix[xi] = new double[covmat.y+2];            
  }

  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
        
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;
    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }
	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }  
      double *y = new double[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 
      freqStruct *freq = (freqStruct *) pars->extras[6];

      stat[yi][s]=dosageAssoc(pars,&design,&designNull,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc,model,isBinary,isCount,start,1);
      
      //if not enough, WT, HE or HO or ind to run test
      if(stat[yi][s]< -900){
	assoc->betas[s]=NAN;
	assoc->SEs[s]=NAN;
      } else{
	assoc->betas[s]=start[0];
	assoc->SEs[s]=standardError(start,&design,NULL,y,pars->post[s],isBinary,isCount, keepList, pars->nInd);
      }
      
      //cleanup
      delete [] y;
      delete [] keepList;
      
    } //phenotypes end
    
    delete [] keepListAll;
  } // sites end

  delete [] start;

  angsd::deleteMatrix(design);
  angsd::deleteMatrix(designNull);

  assoc->stat=stat;
  assoc->keepInd=keepInd;
}



// fits a linear model using weighted least squares, 3 observations per individual -
// up to 3 possible genotypes, when only seq data is known
int abcAsso::getFitWLS(double* start, double* y, double** covMatrix, double* weights, int nInd3, int nEnv, int df){

   /*
     linear regression. Fits a linear model with weights
     y is the responce (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the posterior p(G|X) - function as weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)
   */


   std::vector<double> yw(nInd3);
   std::vector<double> xw(nInd3*nEnv);
   
   //double yw[nInd3]; //<-stripped phenos scaled by stripped weights
   //double xw[nInd3*nEnv]; //<-stripped designs

   int cnt=0;
   for(int i=0;i<nInd3;i++){     
     yw[i] = y[i]*sqrt(weights[i]);
     for(int j=0;j<nEnv;j++){       
       xw[i*nEnv+j] = covMatrix[i][j] * sqrt(weights[i]);       
     }        
   }
    
   double XtX[nEnv*nEnv];
   for(int i=0;i<nEnv*nEnv;i++)
     XtX[i]=0;

   // this is doing the matrix product of (X)^T*W*X 
   for(int x=0;x<nEnv;x++)
     for(int y=0;y<nEnv;y++)
       for(int i=0;i<nInd3;i++)
	 XtX[x*nEnv+y]+=xw[i*nEnv+x]*xw[i*nEnv+y];

#if 0
   //print before inversion
   fprintf(stderr,"BEFORE:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
#endif

   // doing inv((X)^T*W*X)
   int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
   if(singular)
     return 1;
   
#if 0
   //print after inversion
   fprintf(stderr,"AFTER:\n");
   for(int i=0;i<nEnv;i++){
     for(int j=0;j<nEnv;j++)
       fprintf(stderr,"%f ",XtX[i*nEnv+j]);
     fprintf(stderr,"\n");
   } 
#endif

   std::vector<double> Xt_y(nEnv);
   std::vector<double> invXtX_Xt_y(nEnv);
   //double Xt_y[nEnv];
   //double invXtX_Xt_y[nEnv]; 
   for(int i=0;i<nEnv;i++)
     Xt_y[i]=0;
   for(int i=0;i<nEnv;i++)
     invXtX_Xt_y[i]=0;

   // doing (X)^T*W*Y
   for(int x=0;x<nEnv;x++)
     for(int i=0;i<nInd3;i++)
       Xt_y[x]+=xw[i*nEnv+x]*yw[i];

   // calculating the coefs: inv((X)^T*W*X)*((X)^T*W*Y)
   for(int x=0;x<nEnv;x++)
     for(int y=0;y<nEnv;y++)
       invXtX_Xt_y[x] += XtX[y*nEnv+x]*Xt_y[y];

   for(int x=0;x<nEnv;x++){
     start[x]=invXtX_Xt_y[x];
   }

   std::vector<double> yTilde(nInd3);
   //double yTilde[nInd3];
   for(int i=0;i<nInd3;i++)
     yTilde[i] = 0;

   for(int i=0;i<nInd3;i++)
     for(int x=0;x<nEnv;x++)
       yTilde[i] += covMatrix[i][x]*start[x];
       
   double ts=0;
   for(int i=0;i<nInd3;i++){
     // getting the residuals
     double tmp = y[i]-yTilde[i];    
     ts += tmp*tmp*weights[i];
   }
        
   if(df==-1){
     start[nEnv] = sqrt(ts/(1.0*(nInd3-nEnv)));
   } else{
     start[nEnv] = sqrt(ts/(1.0*df));
   }
   return 0;
        
 }


int abcAsso::getFitWLSBin(double* start, double* y, double** covMatrix, double* weights, int nInd3, int nEnv, int df){

  //emil - made this settable
  double tol = assoThres;
  
   /*
     Fits logistic model using iterativly weighted least squares (IWLS)
     linear regression. Fits a logistic model with weights
     y is the binary response (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)

     ASSUMES ALL INDIVIDUALS HAVE VALID WEIGHTS
   */
  
   std::vector<double> yw(nInd3);
   std::vector<double> xw(nInd3*nEnv);
   
   //double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   //double ww[nIndW]; //<-stripped weights
   //double xw[nIndW*nEnv]; //<-stripped designs

   for(int i=0;i<nInd3;i++){
     yw[i] = y[i];
     for(int j=0;j<nEnv;j++){
       xw[i*nEnv+j] = covMatrix[i][j];       
     }     
   }

   std::vector<double> mustart(nInd3);
   //double mustart[nInd3];
   // if weights do not exists, have all weights be 1     
   for(int i=0;i<nInd3;i++){
     mustart[i] = (weights[i]*yw[i] + 0.5)/(weights[i] + 1);
   }

   std::vector<double> eta(nInd3);
   //double eta[nInd3];
   for(int i=0;i<nInd3;i++){
     //link
     eta[i] = log(mustart[i]/(1-mustart[i]));
   }

   std::vector<double> mu(nInd3);
   //double mu[nInd3];
   for(int i=0;i<nInd3;i++){
     //linkinv
     mu[i] = 1 / (1 + exp(-eta[i]));
   }

   std::vector<double> mu0(nInd3);
   std::vector<double> muetaval(nInd3);
   std::vector<double> z(nInd3);
   std::vector<double> w(nInd3);
   std::vector<double> Xt_y(nEnv);
   std::vector<double> invXtX_Xt_y(nEnv);   
   double XtX[nEnv*nEnv];

   for(int i=0;i<nEnv*nEnv;i++)
     XtX[i]=0;
   
   //   double mu0[nInd3];
   //double muetaval[nInd3];
   //double z[nInd3];
   //double w[nInd3];
   //   double XtX[nEnv*nEnv] = {0};
   //double Xt_y[nEnv] = {0};
   //double invXtX_Xt_y[nEnv] = {0}; 
   
   // we run this 20 times...
   for(int t=0;t<assoIter;t++){
     
     for(int i=0;i<nInd3;i++){
       // because mu is same as linkinv(eta)
       muetaval[i] = mu[i]*(1-mu[i]); 
     }

     for(int i=0;i<nInd3;i++){
       // can be calculated together as do not depent on eachother
       z[i] = eta[i]+(yw[i]-mu[i])/muetaval[i];       
       w[i] = sqrt( (weights[i]*(muetaval[i]*muetaval[i])) / (mu[i]*(1-mu[i])) );       
     }
      
     // this is doing the matrix product of (X)^T*W*X
     // takes all columns of second matrix and puts on first column of first matrix - stores at first 0 to nEnv-1 values of XtX
     // takes all columns of second matrix and puts on second column of first matrix - stores at nEnv to 2*nEnv-1 values of XtX
     // thereby the same as taking rows of transposed first matrix (cols of org) and putting it on all columns of second matrix
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 for(int i=0;i<nInd3;i++){
	   // t(xw) %*% xw
	   XtX[x*nEnv+y]+=xw[i*nEnv+x]*w[i]*xw[i*nEnv+y]*w[i];
	 }
       }
     }     

#if 0
     //print before inversion
     fprintf(stderr,"BEFORE:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif

     //fprintf(stderr,"XtX %f %f %f %f\n",XtX[0],XtX[1],XtX[2],XtX[3]);
     
     int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
     if(singular)
       return 1;

     //fprintf(stderr,"XtX %f %f %f %f\n",XtX[0],XtX[1],XtX[2],XtX[3]);
        
#if 0
     //print after inversion
     fprintf(stderr,"AFTER:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif 
    
     // this is doing the matrix product of (X)^T*W*Y
     // takes first column of first matrix and puts on first and only column of second matrix
     // takes second column of first matrix and puts on first and only column of second matrix
     for(int x=0;x<nEnv;x++){
       for(int i=0;i<nInd3;i++){	 
	 Xt_y[x] += xw[i*nEnv+x]*w[i]*z[i]*w[i];
       }       
     }
     
     // calculating coefficients so inv((X)^T*W*X)*((X)^T*W*Y)
     // the coefficients are stored in start
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 invXtX_Xt_y[x] += XtX[y*nEnv+x]*Xt_y[y];
       }
       start[x]=invXtX_Xt_y[x];
     }
          
     // eta
     for(int i=0;i<nInd3;i++){
       // clear values of eta
       eta[i]=0;
       for(int x=0;x<nEnv;x++)
	 eta[i] += xw[i*nEnv+x]*start[x];
     }
     double diff = 0;
     
     // mu
     for(int i=0;i<nInd3;i++){
       //linkinv
       mu[i] = 1 / (1 + exp(-eta[i]));
       // we cannot do this for first iteration diff between previous (mu0) and current iteration (mu)
       diff += fabs(mu[i]-mu0[i]);
       // mu0 has values of previous iteration
       mu0[i] = 1 / (1 + exp(-eta[i]));    
     }

     if(diff<tol & t>0){
       break;
     }
     
     // clear those that have +=
     // much faster than setting values 0 than in for loop
     memset(XtX, 0, sizeof(XtX));
     Xt_y.assign(nEnv,0);
     invXtX_Xt_y.assign(nEnv,0);
     //memset(Xt_y, 0, sizeof(Xt_y));
     //memset(invXtX_Xt_y, 0, sizeof(invXtX_Xt_y));
     
   }
   
   return 0;
}



int abcAsso::getFitWLSPois(double* start, double* y, double** covMatrix, double* weights, int nInd3, int nEnv, int df){

  //emil - made this settable
  double tol = assoThres;
  
   /*
     Fits poisson model using iterativly weighted least squares (IWLS)
     linear regression. Fits a poisson model with weights
     y is the count data response (nInd)
     covMatrix is the design matrix [nInd][nEnv]
     weights is the weights (nInd*nEnv)
     int df is the number of degrees of freedom
     nInd is the number of individuals
     nEnv is the number of predictors (including the intercept)

     ASSUMES ALL INDIVIDUALS HAVE VALID WEIGHTS
   */
  
   std::vector<double> yw(nInd3);
   std::vector<double> xw(nInd3*nEnv);
   
   //double yw[nIndW]; //<-stripped phenos scaled by stripped weights
   //double ww[nIndW]; //<-stripped weights
   //double xw[nIndW*nEnv]; //<-stripped designs

   for(int i=0;i<nInd3;i++){
     yw[i] = y[i];
     for(int j=0;j<nEnv;j++){
       xw[i*nEnv+j] = covMatrix[i][j];       
     }     
   }

   std::vector<double> mustart(nInd3);
   //double mustart[nInd3];
   // if weights do not exists, have all weights be 1     
   for(int i=0;i<nInd3;i++){
     mustart[i] = (weights[i]*yw[i] + 0.5)/(weights[i] + 1);
   }

   std::vector<double> eta(nInd3);
   //double eta[nInd3];
   for(int i=0;i<nInd3;i++){
     //link
     eta[i] = log(mustart[i]);
   }

   std::vector<double> mu(nInd3);
   //double mu[nInd3];
   for(int i=0;i<nInd3;i++){
     //linkinv
     mu[i] = exp(eta[i]);
   }

   std::vector<double> mu0(nInd3);
   std::vector<double> muetaval(nInd3);
   std::vector<double> z(nInd3);
   std::vector<double> w(nInd3);
   std::vector<double> Xt_y(nEnv);
   std::vector<double> invXtX_Xt_y(nEnv);   
   double XtX[nEnv*nEnv];

   for(int i=0;i<nEnv*nEnv;i++)
     XtX[i]=0;
   
   //   double mu0[nInd3];
   //double muetaval[nInd3];
   //double z[nInd3];
   //double w[nInd3];
   //   double XtX[nEnv*nEnv] = {0};
   //double Xt_y[nEnv] = {0};
   //double invXtX_Xt_y[nEnv] = {0}; 
   
   // we run this 20 times...
   for(int t=0;t<assoIter;t++){
     
     for(int i=0;i<nInd3;i++){
       // because mu is same as linkinv(eta)
       muetaval[i] = mu[i]; 
     }

     for(int i=0;i<nInd3;i++){
       // can be calculated together as do not depent on eachother
       z[i] = eta[i]+(yw[i]-mu[i])/muetaval[i];       
       w[i] = sqrt( (weights[i]*(muetaval[i]*muetaval[i])) / mu[i] );       
     }
      
     // this is doing the matrix product of (X)^T*W*X
     // takes all columns of second matrix and puts on first column of first matrix - stores at first 0 to nEnv-1 values of XtX
     // takes all columns of second matrix and puts on second column of first matrix - stores at nEnv to 2*nEnv-1 values of XtX
     // thereby the same as taking rows of transposed first matrix (cols of org) and putting it on all columns of second matrix
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 for(int i=0;i<nInd3;i++){
	   // t(xw) %*% xw
	   XtX[x*nEnv+y]+=xw[i*nEnv+x]*w[i]*xw[i*nEnv+y]*w[i];
	 }
       }
     }     

#if 0
     //print before inversion
     fprintf(stderr,"BEFORE:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif

     //fprintf(stderr,"XtX %f %f %f %f\n",XtX[0],XtX[1],XtX[2],XtX[3]);
     
     int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
     if(singular)
       return 1;

     //fprintf(stderr,"XtX %f %f %f %f\n",XtX[0],XtX[1],XtX[2],XtX[3]);
        
#if 0
     //print after inversion
     fprintf(stderr,"AFTER:\n");
     for(int i=0;i<nEnv;i++){
       for(int j=0;j<nEnv;j++)
	 fprintf(stderr,"%f ",XtX[i*nEnv+j]);
       fprintf(stderr,"\n");
     } 
#endif 
    
     // this is doing the matrix product of (X)^T*W*Y
     // takes first column of first matrix and puts on first and only column of second matrix
     // takes second column of first matrix and puts on first and only column of second matrix
     for(int x=0;x<nEnv;x++){
       for(int i=0;i<nInd3;i++){	 
	 Xt_y[x] += xw[i*nEnv+x]*w[i]*z[i]*w[i];
       }       
     }
     
     // calculating coefficients so inv((X)^T*W*X)*((X)^T*W*Y)
     // the coefficients are stored in start
     for(int x=0;x<nEnv;x++){
       for(int y=0;y<nEnv;y++){
	 invXtX_Xt_y[x] += XtX[y*nEnv+x]*Xt_y[y];
       }
       start[x]=invXtX_Xt_y[x];
     }
          
     // eta
     for(int i=0;i<nInd3;i++){
       // clear values of eta
       eta[i]=0;
       for(int x=0;x<nEnv;x++)
	 eta[i] += xw[i*nEnv+x]*start[x];
     }
     double diff = 0;
     
     // mu
     for(int i=0;i<nInd3;i++){
       //linkinv
       mu[i] = exp(eta[i]);
       // we cannot do this for first iteration diff between previous (mu0) and current iteration (mu)
       diff += fabs(mu[i]-mu0[i]);
       // mu0 has values of previous iteration
       mu0[i] = exp(eta[i]);    
     }

     if(diff<tol & t>0){
       break;
     }
     
     // clear those that have +=
     // much faster than setting values 0 than in for loop
     memset(XtX, 0, sizeof(XtX));
     Xt_y.assign(nEnv,0);
     invXtX_Xt_y.assign(nEnv,0);
     //memset(Xt_y, 0, sizeof(Xt_y));
     //memset(invXtX_Xt_y, 0, sizeof(invXtX_Xt_y));
     
   }
   
   return 0;
}



//pat[add/rec][g]
int pat[3][3] = 
  // add with geno
  {{0,1,2},
    // dom 
   {0,1,1},
   // rec 
   {0,0,1}};

  
double abcAsso::logLike(double *start,double* y,angsd::Matrix<double> *design,double *post,int isBinary,int isCount,int fullModel){
  double ret = 0;
 
  double tmp=0;
  for(int i=0;i<design->x;i++){
    double m=0;
    for(int j=0;j<design->y;j++){
      m += design->matrix[i][j]*start[j];
    }

    if(isBinary){
      m = angsd::bernoulli(y[i],(exp(m)/(exp(m)+1)),0);
      
    } else if(isCount){
      double lambda = exp(m);
      // now handling if p is 0 or 1 (makes it very small or large)
      m = angsd::poisson(y[i],lambda,0);	
      
    } else{
      // design only has n param
      m = angsd::dnorm(y[i],m,start[design->y],0);

    }
    
    if(fullModel){
      tmp += m*post[i];
      if((i % 3 )==2){
	// same as taking log of all values and then summing and then flipping sign (plus to minus)
	ret -= log(tmp);    
	tmp = 0;
      }	  	  
    } else{
      //likelihood for null model
      ret -= log(m);
    }    
  }  
  return ret;
}


//double updateEM(funkyPars *p){
double abcAsso::logupdateEM(double* start,angsd::Matrix<double> *design,angsd::Matrix<double> *postAll,double* y,int keepInd,double* post,int isBinary,int isCount,int fullModel, int iter){
  
  double meanPheno = 0;
  if(iter){
    for(int i=0;i<design->x;i++){
      meanPheno += y[i];
    }
    meanPheno = meanPheno/design->x;    
  }
    
  for(int i=0;i<design->x;i++){
    double m = 0;   
    for(int j=0;j<design->y;j++){
      m += design->matrix[i][j]*start[j];
    }

    if(isBinary){
      double prob;      
      if(iter==0 and not doPriming){
	prob = exp(meanPheno)/(exp(meanPheno)+1.0);	
      } else{
	prob = exp(m)/(exp(m)+1.0);	
      }
      // now handling if p is 0 or 1 (makes it very small or large)
      m = angsd::bernoulli(y[i],prob,1);            
   } else if(isCount){
      double lambda;
      if(iter==0 and not doPriming){
	lambda = exp(meanPheno);	      
      } else{
	lambda = exp(m);
      }            
      // now handling if p is 0 or 1 (makes it very small or large)
      m = angsd::poisson(y[i],lambda,1);            

    } else{      
      // density function of normal distribution
      if(iter==0 and not doPriming){
	m = angsd::dnorm(y[i],meanPheno,start[design->y],1);
      } else{
	m = angsd::dnorm(y[i],m,start[design->y],1);
      }
    }
    
    double tmp = m + log(post[i]);
    postAll->matrix[(size_t)floor(i/3)][i % 3] = tmp;
  }
  
  // x is number of indis, y is 3
  postAll->x = (size_t) design->x/3;
  postAll->y = 3;

  // to get rowSums
  std::vector<double> postTmp(postAll->x);  
  //  double postTmp[postAll->x];
  
  for(int i=0;i<postAll->x;i++){
    double tmp = 0;
    double maxval = postAll->matrix[i][0];
    for(int j=0;j<postAll->y;j++){
      // find max - part of trick for doing log(p1+p2+p3)
      maxval = std::max(maxval,postAll->matrix[i][j]);     
    }
    // trick to avoid over/underflow - log(exp(log(val1)-log(max)) + ...) + log(max) = (exp(log(val1))/exp(log(max)))*(max) + ...
    // same (exp(log(val1))/exp(log(max)))*(max)
    postTmp[i] = log(exp(postAll->matrix[i][0]-maxval)+exp(postAll->matrix[i][1]-maxval)+exp(postAll->matrix[i][2]-maxval)) + maxval;  
  }
   
  // divide each entry of a row with sum of that row
  int df = postAll->x - design->y;
  for(int i=0;i<postAll->x;i++){
    double tmp = 0;
    for(int j=0;j<postAll->y;j++){          
      postAll->matrix[i][j] -= postTmp[i];
    }    
  }
  
  //need to flatten the weights, which is p(s|y,G,phi,Q,f)
  // so first four values is for first indi, and so on...
  double weigths[postAll->x*postAll->y];
  int a = 0;
  for(int i=0;i<postAll->x;i++)
    for(int j=0;j<postAll->y;j++){
      weigths[a++] = exp(postAll->matrix[i][j]);
      //check if issue with weights
      if(exp(postAll->matrix[i][j])!=exp(postAll->matrix[i][j]) or std::isinf(exp(postAll->matrix[i][j]))){
	  fprintf(stderr,"Issue with weights being nan or inf\n");
	  return(-9);
      }
    }   
  
  if(isBinary){    
    getFitWLSBin(start,y,design->matrix,weigths,keepInd*3,design->y,df);    
  } else if(isCount){
    getFitWLSPois(start,y,design->matrix,weigths,keepInd*3,design->y,df);
  } else{
    //double resi[nInd];
    getFitWLS(start,y,design->matrix,weigths,keepInd*3,design->y,df);
  }
     
  return logLike(start,y,design,post,isBinary,isCount,fullModel);
}




double abcAsso::doEMasso(funkyPars *p,angsd::Matrix<double> *design,angsd::Matrix<double> *designNull,angsd::Matrix<double> *postAll,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s,assoStruct *assoc,int model, int isBinary, int isCount, double* start, int fullModel){

  int maf0=1;

  // if too rare do not run analysis
  if(freq>0.995||freq<0.005){
    maf0=0;
  }

  int highWT=0;
  int highHE=0;
  int highHO=0;
  
  double yNull[keepInd];
  int count = 0;
  double post[keepInd*3];

  double covMatrixNull[(covmat.y+1)*keepInd];
  double yfit[keepInd];

  //for doing dosage regression to prime EM (give good starting coefs)
  double covMatrixDose[(covmat.y+2)*keepInd];
  
  designNull->x=keepInd;
  //covars + intercept
  designNull->y=covmat.y+1;
     
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      yNull[count]=yOrg[i];
      designNull->matrix[count][0] = 1;
      covMatrixNull[count] = 1;

      if(model==3){
	covMatrixDose[count] = postOrg[i*3+2];
      } else if(model==2){
	covMatrixDose[count] = postOrg[i*3+2] + postOrg[i*3+1];
      } else{
	// getting dosage of the minor allele or second allele of beagle file
	covMatrixDose[count] = 2*postOrg[i*3+2] + postOrg[i*3+1];
      }      
      covMatrixDose[keepInd+count] = 1;
            
      for(int c=0;c<covmat.y;c++){
	designNull->matrix[count][c+1] = covmat.matrix[i][c];
	covMatrixNull[(c+1)*keepInd+count] = covmat.matrix[i][c];
	covMatrixDose[(c+2)*keepInd+count] = covmat.matrix[i][c];
      }
      count++;    
    }
  }

  // both need intercept and sd(y)
  double startNull[covmat.y+2];

  for(int i=0;i<(covmat.y+1);i++){
    startNull[i] = drand48()*2-1;
  }

  startNull[covmat.y+1]=angsd::sd(yNull,keepInd); 
  
  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }

  if(isBinary){
    getFitBin(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  } else if(isCount){
    getFitPois(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  } else{      
    getFit(yfit,yNull,covMatrixNull,keepInd,designNull->y,startNull);
  }
  
  double llhNull = logLike(startNull,yNull,designNull,post,isBinary,isCount,0);

  if(doPriming){
    
    int dimDose=designNull->y+1;    
    //do dosage regression first to prime EM
    if(isBinary){
      getFitBin(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    } else if(isCount){
      getFitPois(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    } else{      
      getFit(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    }
    
  } else{
    //first element of start is coef of genotype    
    start[0] = drand48()*2-1;
    
    for(int i=0;i<(covmat.y+2);i++){
      start[i+1] = startNull[i];      
    }
  }
  
  double y[3*keepInd];
  count=0;
  
  design->x=3*keepInd;
  design->y=covmat.y+2;
  
  postAll->x=3*keepInd;
  postAll->y=3;
  
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      
      for(int j=0;j<3;j++){
	y[count*3+j]=yOrg[i];	
	post[count*3+j]=postOrg[i*3+j];
		
	// minus 1 (0-indexed), model parameter is 1: add, 2: dom, 3: rec
	design->matrix[count*3+j][0] = pat[model-1][j];
	design->matrix[count*3+j][1] = 1;
	for(int c=0;c<covmat.y;c++){
	  // design matrix has 3 rows per indi - possible genotypes 
	  design->matrix[count*3+j][2+c] = covmat.matrix[i][c];
	}
      }
      
      //if rec model, WT+HE->WT, HO->HE
      if(model==3){
	if((post[count*3+0]+post[count*3+1])>0.90)
	  highWT++;
	if(post[count*3+2]>0.90)
	  highHE++;
      } else if(model==2){
	 //if dom model, WT->WT, HE+HO->HE, HO->0
	if(post[count*3+0]>0.90)
	  highWT++;
	if((post[count*3+1]+post[count*3+2])>0.90)
	  highHE++;
      } else{
	if(post[count*3+0]>0.90)
	  highWT++;
	if(post[count*3+1]>0.90)
	  highHE++;
	if(post[count*3+2]>0.90)
	  highHO++;
      }
      
      count++;
    }
  }
   
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;

  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;

  if(nGeno<2){
    return(-999);//set_snan(lrt);
  }

  //basically the same, just renamed for normScoreEnv and binomScoreEnv functions
  int numInds = keepInd;
  //freq*numInds*2 is the expected number of minor alleles
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount){
    return(-999);//set_snan(lrt);      set_snan(lrt);
  }
  
  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }  
  
  double pars0[design->y+1]; 
  memcpy(pars0,start,sizeof(double)*(design->y+1));
  int nIter=0;
  
  double llh0 = logLike(start,y,design,post,isBinary,isCount,fullModel);  
  for(int i=0;i<emIter;i++){
    nIter++;
    double llh1 = logupdateEM(start,design,postAll,y,keepInd,post,isBinary,isCount,fullModel,i);
            
    if(fabs(llh1-llh0)<emThres and llh1 > 0){	 
      // Converged
      assoc->emIter[s] = nIter;
      break;
    } else if(llh0<llh1 & (not doPriming | not (i==0))){
      // Fit caused increase in likelihood, will roll back to previous step
      memcpy(start,pars0,sizeof(double)*(design->y+1));
      assoc->emIter[s] = nIter;                  
      break;
      // if we problems with the weights then break
    } else if(llh1<-4){
      for(int i=0;i<design->y+1;i++){
	start[i]=NAN;
      }
      llh0=NAN;
      llhNull=NAN;
      assoc->emIter[s] = nIter;    
      break;
    }
    
    llh0=llh1;
    memcpy(pars0,start,sizeof(double)*(design->y+1));
        
  }

  
        
  // likelihood ratio - chi square distributed according to Wilk's theorem
  double LRT = -2*(llh0-llhNull);
  return(LRT);  
}

void abcAsso::emAsso(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }
  
  double *start = new double[covmat.y+3];
  int **keepInd  = new int*[ymat.y];
  // we can also get coef of genotype
  double **stat = new double*[ymat.y*2];

  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  
  angsd::Matrix<double> designNull;
  designNull.x=pars->nInd;
  //covars + intercept
  designNull.y=covmat.y+1;
  designNull.matrix=new double*[pars->nInd];
    
  angsd::Matrix<double> design; 
  design.x=3*pars->nInd;
  design.y=covmat.y+2;
  design.matrix=new double*[3*pars->nInd];
  
  angsd::Matrix<double> postAll;
  postAll.x=3*pars->nInd;
  postAll.y=3;
  postAll.matrix=new double*[3*pars->nInd];
  
  for(int xi=0;xi<pars->nInd;xi++){    
    designNull.matrix[xi] = new double[covmat.y+1];
    for(int i=0;i<3;i++){
      design.matrix[3*xi+i] = new double[covmat.y+2];    
      postAll.matrix[3*xi+i] = new double[postAll.y];
    }
  }
      
  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
        
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;
    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }

	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }
      
      double *y = new double[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 
      freqStruct *freq = (freqStruct *) pars->extras[6];
  
      stat[yi][s]=doEMasso(pars,&design,&designNull,&postAll,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc,model,isBinary,isCount,start,1);
            
      //if not enough, WT, HE or HO or ind to run test
      if(stat[yi][s] < -900){
	assoc->betas[s]=NAN;
	assoc->SEs[s]=NAN;
      } else{
	assoc->betas[s]=start[0]; // giving coefs
	assoc->SEs[s]=standardError(start,&design,&postAll,y,pars->post[s],isBinary,isCount,keepList,pars->nInd);
      }
	
      //cleanup
      delete [] y;
      delete [] keepList;
      
    } //phenotypes end
    
    delete [] keepListAll;
  } // sites end

  delete [] start;

  designNull.x=pars->nInd;
  design.x=3*pars->nInd;
  postAll.x=3*pars->nInd;
  
  angsd::deleteMatrix(postAll);
  angsd::deleteMatrix(design);
  angsd::deleteMatrix(designNull);
  
  assoc->stat=stat;
  assoc->keepInd=keepInd;

}



double abcAsso::doEMassoWald(funkyPars *p,angsd::Matrix<double> *design,angsd::Matrix<double> *postAll,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s,assoStruct *assoc,int model, int isBinary, int isCount, double* start, int fullModel){

  int maf0=1;

  // if too rare do not run analysis
  if(freq>0.995||freq<0.005){
    maf0=0;
  }

  int highWT=0;
  int highHE=0;
  int highHO=0;
  
  int count = 0;
  double post[keepInd*3];
  
  double yfit[keepInd];
  double yNull[keepInd];

  //for doing dosage regression to prime EM (give good starting coefs)
  double covMatrixDose[(covmat.y+2)*keepInd];
       
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      yNull[count]=yOrg[i];
      if(model==3){
	covMatrixDose[count] = postOrg[i*3+2];
      } else if(model==2){
	covMatrixDose[count] = postOrg[i*3+2] + postOrg[i*3+1];
      } else{
	// getting dosage of the minor allele or second allele of beagle file
	covMatrixDose[count] = 2*postOrg[i*3+2] + postOrg[i*3+1];
      }      
      covMatrixDose[keepInd+count] = 1;
            
      for(int c=0;c<covmat.y;c++){
	covMatrixDose[(c+2)*keepInd+count] = covmat.matrix[i][c];
      }
      count++;    
    }
  }
  
  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }

  //how many parameters or columns in design matrix - covars + intercept + geno

  
  if(doPriming){
    
    int dimDose = covmat.y+2;
    
    //do dosage regression first to prime EM
    if(isBinary){
      getFitBin(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    } else if(isCount){
      getFitPois(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    } else{      
      getFit(yfit,yNull,covMatrixDose,keepInd,dimDose,start);
    }

  }  else{
    //first element of start is coef of genotype    
    start[0] = drand48()*2-1;
    
    for(int i=0;i<(covmat.y+2);i++){
      start[i+1] = drand48()*2-1;
    }
  }
    
  double y[3*keepInd];
  count=0;
      
  design->x=3*keepInd;
  design->y=covmat.y+2;

  postAll->x=3*keepInd;
  postAll->y=3;
  
  // WLS (weighted) model with genotypes
  for(int i=0;i<p->nInd;i++){
    if(keepList[i]){
      
      for(int j=0;j<3;j++){
	y[count*3+j]=yOrg[i];	
	post[count*3+j]=postOrg[i*3+j];
		
	// minus 1 (0-indexed), model parameter is 1: add, 2: dom, 3: rec
	design->matrix[count*3+j][0] = pat[model-1][j];
	design->matrix[count*3+j][1] = 1;
	for(int c=0;c<covmat.y;c++){
	  // design matrix has 3 rows per indi - possible genotypes 
	  design->matrix[count*3+j][2+c] = covmat.matrix[i][c];
	}
      }
      
      //if rec model, WT+HE->WT, HO->HE
      if(model==3){
	if((post[count*3+0]+post[count*3+1])>0.90)
	  highWT++;
	if(post[count*3+2]>0.90)
	  highHE++;
      } else if(model==2){
	 //if dom model, WT->WT, HE+HO->HE, HO->0
	if(post[count*3+0]>0.90)
	  highWT++;
	if((post[count*3+1]+post[count*3+2])>0.90)
	  highHE++;
      } else{
	if(post[count*3+0]>0.90)
	  highWT++;
	if(post[count*3+1]>0.90)
	  highHE++;
	if(post[count*3+2]>0.90)
	  highHO++;
      }
      
      count++;
    }
  }
  
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;

  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;

  if(nGeno<2){
    return(-999);//set_snan(lrt);
  }

  //basically the same, just renamed for normScoreEnv and binomScoreEnv functions
  int numInds = keepInd;
  //freq*numInds*2 is the expected number of minor alleles
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount){
    return(-999);//set_snan(lrt);      set_snan(lrt);
  }
  
  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }  

  double pars0[design->y+1]; 
  memcpy(pars0,start,sizeof(double)*(design->y+1));
  int nIter=0;
  
  double llh0 = logLike(start,y,design,post,isBinary,isCount,fullModel);
  
  for(int i=0;i<emIter;i++){
    nIter++;
    double llh1 = logupdateEM(start,design,postAll,y,keepInd,post,isBinary,isCount,fullModel,i);
    
    if(fabs(llh1-llh0)<emThres and llh1 > 0){	 
      // Converged
      assoc->emIter[s] = nIter;
      break;
    } else if(llh0<llh1 & (not doPriming | not (i==0))){
      // Fit caused increase in likelihood, will roll back to previous step
      memcpy(start,pars0,sizeof(double)*(design->y+1));
      assoc->emIter[s] = nIter;                  
      break;
      // if we problems with the weights then break
    } else if(llh1<-4){
      for(int i=0;i<design->y+1;i++){
	start[i]=NAN;
      }
      llh0=NAN;
      assoc->emIter[s] = nIter;    
      break;
    }
    
    llh0=llh1;
    memcpy(pars0,start,sizeof(double)*(design->y+1));
        
  }

        
  // just returns 0 for success
  return(0);  
}

// em asso with wald test
void abcAsso::emAssoWald(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }
  
  double *start = new double[covmat.y+3];
  int **keepInd  = new int*[ymat.y];
  // we can also get coef of genotype
  double **stat = new double*[ymat.y*2];
    
  angsd::Matrix<double> design; 
  design.x=3*pars->nInd;
  design.y=covmat.y+2;
  design.matrix=new double*[3*pars->nInd];
  
  angsd::Matrix<double> postAll;
  postAll.x=3*pars->nInd;
  postAll.y=3;
  postAll.matrix=new double*[3*pars->nInd];
  
  for(int xi=0;xi<pars->nInd;xi++){    
    for(int i=0;i<3;i++){
      design.matrix[3*xi+i] = new double[covmat.y+2];    
      postAll.matrix[3*xi+i] = new double[postAll.y];
    }
  }
      
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
        
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;
    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }

	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }  
      double *y = new double[pars->nInd];

      double tmp = 0;
      
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 
      freqStruct *freq = (freqStruct *) pars->extras[6];
  
      tmp=doEMassoWald(pars,&design,&postAll,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc,model,isBinary,isCount,start,1);
            
      //if not enough, WT, HE or HO or ind to run test
      if(tmp < -900){
	assoc->betas[s]=NAN;
	assoc->SEs[s]=NAN;
      } else{
	assoc->betas[s]=start[0]; // giving coefs
	assoc->SEs[s]=standardError(start,&design,&postAll,y,pars->post[s],isBinary,isCount,keepList,pars->nInd);
      }

      stat[yi][s]=(start[0]*start[0])/((assoc->SEs[s])*(assoc->SEs[s]));
	
      //cleanup
      delete [] y;
      delete [] keepList;
      
    } //phenotypes end
    
    delete [] keepListAll;
  } // sites end

  delete [] start;

  design.x=3*pars->nInd;
  postAll.x=3*pars->nInd;
  
  angsd::deleteMatrix(postAll);
  angsd::deleteMatrix(design);
  
  assoc->stat=stat;
  assoc->keepInd=keepInd;

}


void abcAsso::hybridAsso(funkyPars  *pars,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }

  double *start = new double[covmat.y+3];
  int **keepInd  = new int*[ymat.y];
  // we can also get coef of genotype
  double **stat = new double*[ymat.y*2];
  double **statOther = new double*[ymat.y*2];

  //chisq distribution with df=1
  Chisqdist *chisq1 = new Chisqdist(1);
  
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    statOther[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  angsd::Matrix<double> designNull;
  designNull.x=pars->nInd;
  //covars + intercept
  designNull.y=covmat.y+1;
  designNull.matrix=new double*[pars->nInd];
    
  angsd::Matrix<double> design; 
  design.x=3*pars->nInd;
  design.y=covmat.y+2;
  design.matrix=new double*[3*pars->nInd];
  
  angsd::Matrix<double> postAll;
  postAll.x=3*pars->nInd;
  postAll.y=3;
  postAll.matrix=new double*[3*pars->nInd];
  
  for(int xi=0;xi<pars->nInd;xi++){    
    designNull.matrix[xi] = new double[covmat.y+1];
    for(int i=0;i<3;i++){
      design.matrix[3*xi+i] = new double[covmat.y+2];    
      postAll.matrix[3*xi+i] = new double[postAll.y];
    }
  }
  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
        
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;

    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }
	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }
      
      double *y = new double[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 
      freqStruct *freq = (freqStruct *) pars->extras[6];
      stat[yi][s]=doAssociation(pars,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc);

      // cutoff has been changed to p-value
      if(angsd::to_pval(chisq1,stat[yi][s])<hybridThres){
	// add extra params in stat for storing LRT and beta
	statOther[yi][s]=doEMasso(pars,&design,&designNull,&postAll,pars->post[s],y,keepInd[yi][s],keepList,freq->freq[s],s,assoc,model,isBinary,isCount,start,1);
	assoc->SEs[s]=standardError(start,&design,&postAll,y,pars->post[s],isBinary,isCount,keepList,pars->nInd);
	assoc->betas[s]=start[0];
      } else{
	assoc->emIter[s]=0;
	assoc->betas[s]=NAN;
	assoc->SEs[s]=NAN;
	statOther[yi][s]=NAN;
      }
      
      //cleanup
      delete [] y;
      delete [] keepList;
      
    } //phenotypes end
    
    delete [] keepListAll;
  } // sites end

  delete [] start;
  delete chisq1;
  
  designNull.x=pars->nInd;
  design.x=3*pars->nInd;
  postAll.x=3*pars->nInd;
  
  angsd::deleteMatrix(postAll);
  angsd::deleteMatrix(design);
  angsd::deleteMatrix(designNull);

  assoc->statOther=statOther;
  assoc->stat=stat;
  assoc->keepInd=keepInd;
}

int abcAsso::getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv,double *start){

  //emil added start as an argument passed to the function - which gets estimated coefs
  
  /*
    linear regression. Fits a linear model. 
    res is the predicted values  (eta) = X%*%coef
    Y is the responce
    covMatrix is the design matrix (nInd x nEnv)
    nInd is the number of individuals
    nEnv is the number of predictors (including the intersept)
  */

  
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  for(int i=0;i<nEnv;i++)
    Xt_y[i]=0;
  for(int i=0;i<nEnv;i++)
    invXtX_Xt_y[i]=0;
 
  //get t(X)%*%y
  for(int x=0;x<nEnv;x++)//col X
    for(int i=0;i<nInd;i++)
      Xt_y[x]+=covMatrix[x*nInd+i]*Y[i];

  //get inv(t(X)%*%X)
  double XtX[nEnv*nEnv];
  for(int i=0;i<nEnv*nEnv;i++)
    XtX[i]=0;

  for(int x=0;x<nEnv;x++)//col X
    for(int y=0;y<nEnv;y++)//row Xt
      for(int i=0;i<nInd;i++){
	XtX[x*nEnv+y]+=covMatrix[y*nInd+i]*covMatrix[x*nInd+i];	
      }

  //  double workspace[2*nEnv];
  //  angsd::matinv(XtX, nEnv, nEnv, workspace);
  int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
  if(singular)
    return 1 ;

  //get (inv(t(X)%*%X))%*%(t(X)%*%y) //this is the coef!
  for(int x=0;x<nEnv;x++)//col X
    for(int y=0;y<nEnv;y++)//row Xt
      invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];
  
  //get X%*%(inv(t(X)%*%X))%*%(t(X)%*%y)
  for(int j=0;j<nInd;j++){//row Xt
    res[j]=0;
    //storing coefs

  }

  for(int x=0;x<nEnv;x++){//row Xt
      start[x]=invXtX_Xt_y[x];
  }

  double ts = 0;
  for(int j=0;j<nInd;j++){//row Xt
    for(int x=0;x<nEnv;x++){
      res[j]+=covMatrix[x*nInd+j]*invXtX_Xt_y[x];
    }
    ts += (Y[j]-res[j])*(Y[j]-res[j]);
  }
       
  start[nEnv] = sqrt(ts/(1.0*(nInd-nEnv)));
  
  return 0;
}

int abcAsso::getFitBin(double *res,double *Y,double *covMatrix,int nInd,int nEnv, double *start){

  //emil - tested that tolerance makes logistic regression run much faster
  //double tol = 1e-6;
  double tol = assoThres;
  
  /*
    logistic regression. Fits a logistic regression model. 
    res is the estimated coefficients
    Y is the responce
    covMatrix is the design matrix (nInd x nEnv)
    nInd is the number of individuals
    nEnv is the number of predictors (including the intersept)
    // R code
    getFitBin<-function(y,X){
        b<-rep(0,ncol(X))#estimates
	for(i in 1:20){ 
	    eta  <- as.vector(1/(1+exp(-(X%*%b))))
	    change<-solve(t(X*eta*(1-eta)) %*% X ) %*% t(X) %*% ( y - eta )
	    b<-b+ change
	    if(sum(abs(change))<1e-6)
	        break
	    }

	b
    }
    ///////
    coef<-getFitBin(y,X)
    yTilde <- sigm(X%*%coef)
  */
  
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);
  
  double coef[nEnv];
  double eta[nInd];
  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  double XtX[nEnv*nEnv];

  for(int x=0;x<nEnv;x++)//col X
    coef[x]=0;

  //emil - tested that tolerance makes logistic regression run much faster
  //for(int iter=0;iter<100;iter++){
  for(int iter=0;iter<assoIter;iter++){

    //set to zero
    for(int x=0;x<nEnv;x++){
      Xt_y[x]=0;
      invXtX_Xt_y[x]=0;
      for(int y=0;y<nEnv;y++)
	XtX[x*nEnv + y]=0;
    }
    for(int i=0;i<nInd;i++){
      eta[i]=0;
    }

    //eta <- 1/(1+exp(-(X%*%b)))
    for(int i=0;i<nInd;i++){
      for(int x=0;x<nEnv;x++)//col X
      	eta[i]+=coef[x]*covMatrix[x*nInd+i];
      eta[i] = 1.0/(1+exp(-eta[i])); 
    }
     
    //get t(X) %*% ( y - eta )
    for(int x=0;x<nEnv;x++)//col X
      for(int i=0;i<nInd;i++)
	Xt_y[x]+=covMatrix[x*nInd+i]*(Y[i]-eta[i]);

    //get solve(t(X*eta*(1-eta)) %*% X )
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	for(int i=0;i<nInd;i++)
	  XtX[x*nEnv+y]+=covMatrix[y*nInd+i] * eta[i] * covMatrix[x*nInd+i];

    // double workspace[2*nEnv];
    //angsd::matinv(XtX, nEnv, nEnv, workspace);
    int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
    if(singular)
      return 1 ;
  
    //S = svd_inverse(S,flag);     
    //get (inv(t(X)%*%X))%*%(t(X)%*%y)
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];

    double diff = 0;
    for(int x=0;x<nEnv;x++)
      diff += fabs(invXtX_Xt_y[x]);
  
     for(int x=0;x<nEnv;x++)
       coef[x]+=invXtX_Xt_y[x];

     if(diff<tol)
       break;
  }
  
  //yTilde <- X%in%coef
  for(int j=0;j<nInd;j++){//row Xt
    res[j]=0;   
  }

  for(int x=0;x<nEnv;x++){//row Xt
    start[x]=coef[x];
  }

  for(int j=0;j<nInd;j++)//row Xt
    for(int x=0;x<nEnv;x++)
      res[j]+=covMatrix[x*nInd+j]*coef[x];
  for(int j=0;j<nInd;j++)
    res[j]= angsd::sigm(res[j]);

  //   for(int x=0;x<nEnv;x++)
  //   fprintf(stdout,"%f\t",coef[x]);
  //fprintf(stdout,"\n");

  return 0;
}


int abcAsso::getFitPois(double *res,double *Y,double *covMatrix,int nInd,int nEnv, double *start){

  //emil - tested that tolerance makes poisson regression run much faster
  //double tol = 1e-6;
  double tol = assoThres;
  
  /*
    poisson regression. Fits a poisson regression model. 
    res is the estimated coefficients
    Y is the responce
    covMatrix is the design matrix (nInd x nEnv)
    nInd is the number of individuals
    nEnv is the number of predictors (including the intersept)
    // R code
    getFitBin<-function(y,X){
        b<-rep(0,ncol(X))#estimates
	for(i in 1:20){ 
	    eta  <- as.vector(1/(1+exp(-(X%*%b))))
	    change<-solve(t(X*eta*(1-eta)) %*% X ) %*% t(X) %*% ( y - eta )
	    b<-b+ change
	    if(sum(abs(change))<1e-6)
	        break
	    }

	b
    }
    ///////
    coef<-getFitBin(y,X)
    yTilde <- sigm(X%*%coef)
  */
  
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);
  
  double coef[nEnv];
  double eta[nInd];
  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  double XtX[nEnv*nEnv];

  for(int x=0;x<nEnv;x++)//col X
    coef[x]=0;

  //emil - tested that tolerance makes logistic regression run much faster
  //for(int iter=0;iter<100;iter++){
  for(int iter=0;iter<assoIter;iter++){

    //set to zero
    for(int x=0;x<nEnv;x++){
      Xt_y[x]=0;
      invXtX_Xt_y[x]=0;
      for(int y=0;y<nEnv;y++)
	XtX[x*nEnv + y]=0;
    }
    for(int i=0;i<nInd;i++){
      eta[i]=0;
    }

    //eta <- 1/(1+exp(-(X%*%b)))
    for(int i=0;i<nInd;i++){
      for(int x=0;x<nEnv;x++)//col X
      	eta[i]+=coef[x]*covMatrix[x*nInd+i];
      eta[i] = exp(eta[i]); 
    }
     
    //get t(X) %*% ( y - eta )
    for(int x=0;x<nEnv;x++)//col X
      for(int i=0;i<nInd;i++)
	Xt_y[x]+=covMatrix[x*nInd+i]*(Y[i]-eta[i]);

    //get solve(t(X*eta*(1-eta)) %*% X )
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	for(int i=0;i<nInd;i++)
	  XtX[x*nEnv+y]+=covMatrix[y*nInd+i] * eta[i] * covMatrix[x*nInd+i];

    // double workspace[2*nEnv];
    //angsd::matinv(XtX, nEnv, nEnv, workspace);
    int singular=angsd::svd_inverse(XtX,nEnv,nEnv);
    if(singular)
      return 1 ;
  
    //S = svd_inverse(S,flag);     
    //get (inv(t(X)%*%X))%*%(t(X)%*%y)
    for(int x=0;x<nEnv;x++)//col X
      for(int y=0;y<nEnv;y++)//row Xt
	invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];

    double diff = 0;
    for(int x=0;x<nEnv;x++)
      diff += fabs(invXtX_Xt_y[x]);
  
     for(int x=0;x<nEnv;x++)
       coef[x]+=invXtX_Xt_y[x];

     if(diff<tol)
       break;
  }
  
  //yTilde <- X%in%coef
  for(int j=0;j<nInd;j++){//row Xt
    res[j]=0;   
  }

  for(int x=0;x<nEnv;x++){//row Xt
    start[x]=coef[x];
  }

  for(int j=0;j<nInd;j++)//row Xt
    for(int x=0;x<nEnv;x++)
      res[j]+=covMatrix[x*nInd+j]*coef[x];
  for(int j=0;j<nInd;j++)
    res[j]=exp(res[j]);

  //   for(int x=0;x<nEnv;x++)
  //   fprintf(stdout,"%f\t",coef[x]);
  //fprintf(stdout,"\n");

  return 0;
}


///////////////////////////////////////////////////////////////////////////////
// EMIL END
///////////////////////////////////////////////////////////////////////////////

double abcAsso::doAssociation(funkyPars *pars,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s,assoStruct *assoc){
  if(doPrint)
    fprintf(stderr,"Staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);
  
  double covMatrix[(covmat.y+1)*keepInd];
  double y[keepInd];
  double post[keepInd*3];
  double start[keepInd];
  
  int count=0;
  for(int i=0;i<pars->nInd;i++){
    if(keepList[i]){
      y[count]=yOrg[i];
      for(int g=0;g<3;g++)
	post[count*3+g]=postOrg[i*3+g];
      count++;
    }
  }

  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }

  int nEnv;
  if(adjust==1)
    nEnv=(covmat.y+1);
  else
    nEnv=1;

  int num;
  for(int j=0;j<keepInd;j++)
    covMatrix[j]=1;
  if(covmat.matrix!=NULL){
    num=0;
    for(int j=0;j<covmat.x;j++){
      if(keepList[j]==0)
	continue;
      for(int i=1;i<nEnv;i++)
	covMatrix[i*keepInd+num]=covmat.matrix[j][i-1];   
      num++;
    }
  }

  //fprintf(stderr,"number of inds %d\n",keepInd);

  // permutation
  if(sitePerm){
    if((covmat.y+1)==1){
      for(int i=0 ; i<keepInd ;i++){	
	int j = rand() % (keepInd);
	angsd::swapDouble(y[j],y[i]); 
      }
    }
    else{
      int col0=0; 
      for(int i=0 ; i<covmat.x ;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][0]<0.5)
	  col0++;
      }
      for(int i=0 ; i<col0 ;i++) {
	int j = rand() % (col0);
	angsd::swapDouble(y[j],y[i]);
      }
      for(int i=0 ; i<keepInd-col0 ;i++) {
	int j = rand() % (keepInd-col0);
	angsd::swapDouble(y[j+col0],y[i+col0]);
      }
      if(col0<500||keepInd-col0<500){
	fprintf(stderr,"colTrouble %d %d\n",col0,keepInd-col0);
      }
    }

  }

  //
 
  double *yfit = new double[keepInd];
  if(nEnv==1){
    double mean=0;
    for(int i=0;i<keepInd;i++)
      mean+=y[i];
    mean=mean/keepInd;
    for(int i=0;i<keepInd;i++)
      yfit[i]=mean;
  }
  else{
    if(isBinary){
      if(getFitBin(yfit,y,covMatrix,keepInd,nEnv,start))
	return -999;
    } else if(isCount){
      if(getFitPois(yfit,y,covMatrix,keepInd,nEnv,start))
	return -999;      
    } else
      if(getFit(yfit,y,covMatrix,keepInd,nEnv,start))
	return -999;
  }
  
  //for(int i=0;i<keepInd;i++)
  //   fprintf(stdout,"%f\t",y[i]);
  // exit(0);

  if(model==2){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+1]+=post[i*3+2];
      post[i*3+2]=0;
    }
  }
  if(model==3){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+0]+=post[i*3+1];
      post[i*3+1]=post[i*3+2];
      post[i*3+2]=0;
    }
  }

  double stat;
  if(isBinary)
    stat = binomScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);
  else if(isCount)
    stat = poisScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);
  else
    stat = normScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,assoc,s);

  delete[] yfit;
  return stat;

}

double abcAsso::normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  int rankProb=0;
  double sum=0;
  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  double U=0;
  int highWT=0;
  int highHE=0;
  int highHO=0;
  
  double sumEx=0;
  
  for(int i=0;i<numInds;i++)
    sum+=pow(y[i]-ytilde[i],2);
  double var=sum/(numInds-nEnv);
  
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;
 
  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i])/var;
    
    //Vaa<-Vaa+1/var*Xe[tal,]%*%t(Xe[tal,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= (1/var)*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }
    
    for(int x=0;x<nEnv;x++)
      Vab[x]+= (1/var)*Ex[i]*cov[x*numInds+i];
    
    //Vab<-Vab+1/var*Ex[tal]*cbind(Xe[tal,])
    
    if(post[i*3+0]>0.90)
      highWT++;
    if(post[i*3+1]>0.90)
      highHE++;
    if(post[i*3+2]>0.90)
      highHO++;
  }//recursion done
  
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
  
  for(int i =0; i<numInds;i++)
    sumEx+=Ex[i];
    
  //  double Vab=sumEx/var;
  double Vbb=0;
  for(int i =0; i<numInds;i++)
    Vbb+=(1/var-pow(y[i]-ytilde[i],2)/pow(var,2))*Ex2[i]+pow(y[i]-ytilde[i],2)/pow(var,2)*pow(Ex[i],2);
  
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  double workspace[2*nEnv];
  rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
  
  double I=0;
  
  //inv(Vaa)%*%Vab
  double invVaa_Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    invVaa_Vab[x]=0;
  
  //NB! Vaa is now the inverse 
  for(int Nx=0;Nx<nEnv;Nx++)
    for(int Ny=0;Ny<nEnv;Ny++)
      invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];
  
  //I<-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  for(int x=0;x<nEnv;x++)
    I+=Vab[x]*invVaa_Vab[x];
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  //s  fprintf(stderr,"tVab_invVaa_Vab: %f\n",I);  
  I=Vbb-I;


  //the observed varians of the dispersion 
  double Vbs=0;
  for(int i =0; i<numInds;i++)
    Vbs+=Ex[i]*(y[i]-ytilde[i])/pow(var,2);
  
  double Vss=0;
  for(int i =0; i<numInds;i++)
    Vss+=pow(y[i],2)+2*y[i]*ytilde[i];
  Vss=Vss*pow(var,-3)-numInds/(4*M_PI*pow(var,2));

  //fprintf(stderr,"Vbs %f Vss %f\n",Vbs,Vss);
  I=I-pow(Vbs,2)/Vss;

  
  double lrt =pow(U,2)/I;

  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;

  if(nGeno<2)
    lrt=-999;//set_snan(lrt);
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
    lrt=-999;//set_snan(lrt);      set_snan(lrt);
  if(rankProb!=0)
    lrt=-99;
  //      fprintf(lrtfile,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde,var,U,Vaa,Vab,Vbb,I,lrt);
  //fprintf(lrtfile,"\t%d\t%f",numInds,lrt);
  /*
  if(verbose>0)
    fprintf(fp,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde[0],var,U,Vaa[0],Vab[0],Vbb,I,lrt);
    else
    fprintf(fp,"\t%d\t%f",numInds,lrt);
  */

  if((0||lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
    for(int i=0 ; i<numInds ;i++) {
      fprintf(stderr,"y: %f\t  yfit: %f \t post %f %f %f\tEx %f %f\tU %f cov: ",y[i],ytilde[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      for(int j=0;j<nEnv;j++)
	fprintf(stderr,"%f\t",cov[j*numInds+i]);
      fprintf(stderr,"\n");
 
    }
    for(int j=0;j<pow(nEnv,2);j++)
      fprintf(stderr,"Vaa: %f\t",Vaa[j]); 
    fprintf(stderr,"rank %d\tlrt: %f\t",rankProb,lrt); 
    fprintf(stderr,"\n");                                                                                                            
 
  
    for(int j=0;j<nEnv;j++)
      fprintf(stderr,"Vab: %f\t",Vab[j]);
    fprintf(stderr,"\n");

    fflush(stderr);
    fflush(stdout);
    fprintf(stderr,"The test is unreliable. You should increase -minHigh\n");


      //    exit(0);
  }
  
  delete []  Ex;
  delete []  Ex2;
  return lrt;

}
  


double abcAsso::binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  int rankProb=0;

  double U=0;
  int highWT=0;
  int highHE=0;
  int highHO=0;
  double sumEx=0;
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;


  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i]);   

  // Vaa<-Vaa+yTilde[i]*(1-yTilde[i])*A[i,]%*%t(A[i,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= ytilde[i]*(1-ytilde[i])*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }

    for(int x=0;x<nEnv;x++)
      Vab[x]+= ytilde[i]*(1-ytilde[i])*Ex[i]*cov[x*numInds+i];

    //Vba<-Vba+yTilde[i]*(1-yTilde[i])*A[i,]*Ex[i]

    if(post[i*3+0]>0.9)
      highWT++;
    if(post[i*3+1]>0.9)
      highHE++;
    if(post[i*3+2]>0.9)
      highHO++;
  }//recursion done
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
 

    for(int i =0; i<numInds;i++)
      sumEx+=Ex[i];
    // double Vaa=ytilde[0]*(1-ytilde[0])*numInds;
    //    double Vab=ytilde[0]*(1-ytilde[0])*sumEx;
    double Vbb=0;
    for(int i =0; i<numInds;i++)
      Vbb+=(ytilde[i]*(1-ytilde[i])-pow(y[i]-ytilde[i],2))*Ex2[i]+pow(y[i]-ytilde[i],2)*pow(Ex[i],2);

    double workspace[2*nEnv];
    rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
 

    double I =0;

    double invVaa_Vab[nEnv];
    for(int x=0;x<nEnv;x++)
      invVaa_Vab[x]=0;

    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++)
	invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];

    for(int x=0;x<nEnv;x++)
      I+=Vab[x]*invVaa_Vab[x];
    I=Vbb-I;

    double lrt =pow(U,2)/I;

    int nGeno=0;
    if(highWT >= minHigh)
      nGeno++;
    if(highHE >= minHigh)
      nGeno++;
    if(highHO >= minHigh)
      nGeno++;
    if(nGeno<2)
      lrt=-999;
    //freq*numInds*2 is the expected number of minor alleles
    if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
      lrt=-999;
    if(rankProb!=0)
      lrt=-99;
    //dispersion matrix has zero corners

    if((lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
      for(int i=0 ; i<numInds ;i++) {
	fprintf(stderr,"y: %f\t  post %f %f %f\tEx %f %f\tU %f\n",y[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      }
      //      exit(0);
      fprintf(stderr,"The test is unreliable. You should increase -minHigh\n");
      
      
    }

    delete []  Ex;
    delete []  Ex2;
   
    return lrt;
}


// emil my try so far
double abcAsso::poisScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  int rankProb=0;

  double U=0;
  int highWT=0;
  int highHE=0;
  int highHO=0;
  double sumEx=0;
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;


  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i]);   

    // emil - variance of poisson
    double lambda = ytilde[i];
    
    // Vaa<-Vaa+yTilde[i]*(1-yTilde[i])*A[i,]%*%t(A[i,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	// emil - lambda is varinace of poisson
	Vaa[Nx*nEnv+Ny]+= lambda*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }

    for(int x=0;x<nEnv;x++){
      // emil - lambda is varinace of poisson
      Vab[x]+= lambda*Ex[i]*cov[x*numInds+i];
    }
    //Vba<-Vba+yTilde[i]*(1-yTilde[i])*A[i,]*Ex[i]

    if(post[i*3+0]>0.9)
      highWT++;
    if(post[i*3+1]>0.9)
      highHE++;
    if(post[i*3+2]>0.9)
      highHO++;
  }//recursion done
  assoc->highWt[s] = highWT;
  assoc->highHe[s] = highHE;
  assoc->highHo[s] = highHO;
  
  
  for(int i =0; i<numInds;i++)
    sumEx+=Ex[i];
  // double Vaa=ytilde[0]*(1-ytilde[0])*numInds;
  //    double Vab=ytilde[0]*(1-ytilde[0])*sumEx;
  double Vbb=0;
  for(int i =0; i<numInds;i++){
    // emil - variance of poisson
    double lambda = ytilde[i];
    
    Vbb+=(lambda-pow(y[i]-ytilde[i],2))*Ex2[i]+pow(y[i]-ytilde[i],2)*pow(Ex[i],2);
  }
  double workspace[2*nEnv];
  rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
  
  double I=0;
  
  double invVaa_Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    invVaa_Vab[x]=0;
  
  for(int Nx=0;Nx<nEnv;Nx++)
    for(int Ny=0;Ny<nEnv;Ny++)
      invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];
  
  for(int x=0;x<nEnv;x++)
    I+=Vab[x]*invVaa_Vab[x];
  I=Vbb-I;
  
  double lrt = pow(U,2)/I;
  
  int nGeno=0;
  if(highWT >= minHigh)
    nGeno++;
  if(highHE >= minHigh)
    nGeno++;
  if(highHO >= minHigh)
    nGeno++;
  if(nGeno<2)
    lrt=-999;
  //freq*numInds*2 is the expected number of minor alleles
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
    lrt=-999;
  if(rankProb!=0)
    lrt=-99;
  //dispersion matrix has zero corners
  
  if((lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
    for(int i=0 ; i<numInds ;i++) {
      fprintf(stderr,"y: %f\t  post %f %f %f\tEx %f %f\tU %f\n",y[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
    }
    //      exit(0);
    fprintf(stderr,"The test is unreliable. You should increase -minHigh\n");
    
    
  }
  
  delete []  Ex;
  delete []  Ex2;
  
  return lrt;
}


void abcAsso::printDoAsso(funkyPars *pars){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  
  freqStruct *freq = (freqStruct *) pars->extras[6];
  assoStruct *assoc= (assoStruct *) pars->extras[index];
  for(int yi=0;yi<ymat.y;yi++){
    bufstr.l=0;
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0){//will skip sites that have been removed      
	continue;
     }
 
      if(doAsso==5){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d/%d/%d\t%f\t%f\t%f\t%i\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s],assoc->statOther[yi][s],assoc->betas[s],assoc->SEs[s],assoc->emIter[s]);
      } else if(doAsso==4 or doAsso==7){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%f\t%f\t%d/%d/%d\t%i\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->betas[s],assoc->SEs[s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s],assoc->emIter[s]);	 
       } else if(doAsso==6){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%f\t%f\t%d/%d/%d\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->betas[s],assoc->SEs[s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s]);
       } else if(doAsso==2){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d/%d/%d\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->keepInd[yi][s],assoc->stat[yi][s],assoc->highWt[s],assoc->highHe[s],assoc->highHo[s]);
      } else{
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],freq->freq[s],assoc->stat[yi][s]);

      }
    }
    aio::bgzf_write(multiOutfile[yi],bufstr.s,bufstr.l);bufstr.l=0;
  }
}

