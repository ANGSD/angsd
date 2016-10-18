#include <list>
#include <cmath>
#include <cassert>
#include "bfgs.h"
#include "analysisFunction.h"
#include "abcError.h"
#include "abcFreq.h"

//from analysisMaf.cpp
double phatFun(suint *counts,int nFiles,double eps,char major,char minor) ;

//the errorestimation allocates errormatrix defined by 3 genotypes and the and 4 different nucleotides
//it also uses the nSites x 3*nInd with the likelihoods.
//it aloso used logfactorial

// To avoid the spurious alloc/dealloc these are only allocated once and removed once per optimization.
//the following struct is everything that the optimzation need
typedef struct{
  float *logfact;
  double **errors;//<- these are the ones that are optimized
  int numSites;
  int nInd;
  int *major;  int *minor;
  suint **counts;
  
  double ****errorProbs;//huge matrix with errorprobs
  double **loglike3;//matrix with likelihoods

  float EM_start;
  int emIter;
} error_struct;

void abcError::printArg(FILE *argFile){
  fprintf(argFile,"---------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doError\t%d\n",doError);
  fprintf(argFile,"\t1: SYK method, joint typespecific errors (Multisample)\n");
  fprintf(argFile,"\t-minSites\t%d\n",minSites);
  fprintf(argFile,"\t-errors\t\t%s\t(Filename for starterrors)\n",errorFname);
  fprintf(argFile,"\t-emIter\t\t%d\n",emIter);
  fprintf(argFile,"\t-minPhat\t%f\t(Minimum phat)\n",minPhat);
  fprintf(argFile,"\t-eps\t\t%f\t(Estimate of errorrate)\n",eps);
  fprintf(argFile,"\tNB this method requires -doMajorMinor 2\n");
  fprintf(argFile,"\n");

}


void abcError::consolidate(funkyPars *p) {
  for(int s=0;s<p->numSites;s++)
    if(p->keepSites[s]){
      point pp;
      pp.major = p->major[s];
      pp.minor = p->minor[s];
      pp.aLine = new suint[p->nInd*4];
      memcpy(pp.aLine,p->counts[s],4*sizeof(suint)*p->nInd);
      aList.push_back(pp);
    }
  
  //  fprintf(stderr,"sizeof alist=%lu\n",aList.size());
  if(aList.size()>minSites||(p->killSig==1)){
    fprintf(stderr,"\t Triggering error estimation with nsites=%lu\n",aList.size());
    int numSites = aList.size();
    int *major = new int[aList.size()];
    int *minor = new int[aList.size()];
    suint **counts = new suint*[aList.size()];
    int pos = 0;
    for(myList::iterator it=aList.begin();it!=aList.end();it++){
      major[pos] = it->major;
      minor[pos] = it->minor;
      counts[pos] = it->aLine;
      pos++;
    }
    aList.clear();
    //trigger optimiztion
    likeFixedMinorError_bfgs_tsk2(errors,numSites,p->nInd,major,minor,counts,emIter,EM_START); //this should be the same but faster

  }
}



//tsk 30may no allocs
void getLikesFullError3Genotypes_tsk(int numSites,int nInd, int* major,int *minor,suint **counts,double ****errorProbs,float *logfactorial,double **loglikes){

  double *logError;
  for(int s=0;s<numSites;s++){
    if(major[s]==minor[s]){
      if(major[s]>0)
	minor[s]=0;
      else
	minor[s]=1;
    }

    for(int g=0;g<3;g++){
      logError=errorProbs[g][major[s]][minor[s]];
      for(int i=0;i<nInd;i++){
	loglikes[s][i*3+g]=logfactorial[counts[s][i*4+0]+counts[s][i*4+1]+counts[s][i*4+2]+counts[s][i*4+3]];
        for(int j=0;j<4;j++)
          loglikes[s][i*3+g]+=-logfactorial[counts[s][i*4+j]]+counts[s][i*4+j]*logError[j];
      }
    }

    ///printing stuff
    if(0&&s==2){//print infomation about counts and computed likelihood
      for(int g=0;g<3;g++){
	for(int i=0;i<4;i++){
	  logError=errorProbs[g][major[s]][minor[s]];
	fprintf(stderr,"%f\t",logError[i]);
	}
	fprintf(stderr,"\n");
	}
      for(int i=0;i<nInd;i++)
      fprintf(stderr,"%d\t%d\t%f\t%f\t%f\t\n",counts[s][i*4+major[s]],counts[s][i*4+minor[s]],loglikes[s][i*3],loglikes[s][i*3+1],loglikes[s][i*3+2]);
      
    }
  }
  //  delete []logfactorial;
  //  return loglikes;
}


//tsk 30 may noallocs
void generateErrorPointers_tsk(double **error,int nG,int nA,double ****errorProbs) {
  //estimates the probabilities for the four bases given the major and minor
  //The function is only called once and is of no importance speedwise
  //HOWEVER for errorestimation thi is run in every mainfunction call
  //The "errorProbs" is global 

  for(int iGeno=0;iGeno<nG;iGeno++) {
    //    errorProbs[iGeno]= new double**[nA];
    for(int iMajor=0;iMajor<nA;iMajor++){
      //      errorProbs[iGeno][iMajor]= new double*[nA];
      for(int iMinor=0;iMinor<nA;iMinor++){
	//        errorProbs[iGeno][iMajor][iMinor]= new double[nA];	
	if(iMajor==iMinor)
	  continue;
        int other1=-1;
        int other2=-1;
        for(int i=0;i<3;i++){
          other1=i;
          if(i!=iMajor&&i!=iMinor)
            break;
        }
        for(int i=0;i<4;i++){
          other2=i;
          if(i!=iMajor&&i!=iMinor&&i!=other1)
            break;
        }

	errorProbs[iGeno][iMajor][iMinor][iMajor]=angsd::addProtect2(log(1-iGeno*0.5)+log(1-error[iMajor][iMinor]-error[iMajor][other1]-error[iMajor][other2]),log(iGeno)-log(2)+log(error[iMinor][iMajor]));
//prob Major
	errorProbs[iGeno][iMajor][iMinor][iMinor]=angsd::addProtect2(log(iGeno)-log(2)+log(1-error[iMinor][iMajor]-error[iMinor][other1]-error[iMinor][other2]),log(1-iGeno*0.5)+log(error[iMajor][iMinor]));
//prob Minor
	errorProbs[iGeno][iMajor][iMinor][other1]=angsd::addProtect2(log(iGeno)-log(2)+log(error[iMinor][other1]),log(1-iGeno*0.5)+log(error[iMajor][other1]));
//prob other1
	errorProbs[iGeno][iMajor][iMinor][other2]=angsd::addProtect2(log(iGeno)-log(2)+log(error[iMinor][other2]),log(1-iGeno*0.5)+log(error[iMajor][other2]));
//prob other2
	double sum=0;
	for(int ia=0;ia<4;ia++)
	  sum+=exp(errorProbs[iGeno][iMajor][iMinor][ia]);
	if((sum>1.00005||sum<0.999999)&&iMajor!=iMinor){//never happens, just a check
	  for(int j=0;j<4;j++)
	    fprintf(stderr,"%f\t ",errorProbs[iGeno][iMajor][iMinor][j]);
	  fprintf(stderr,"G %d M %d m %d o1 %d o2 %d sum %f %f\n",iGeno,iMajor,iMinor,other1,other2,sum,angsd::addProtect2(log((1-iGeno*0.5)*(1-error[iMajor][iMinor]-error[iMajor][other1]-error[iMajor][other2])),log(iGeno*0.5*error[iMajor][iMinor])));
	  
	}
	
	
      }
    }
  }
  
  //  return errorProbs;
}


//30 may no allocs
double likeFixedMinorError_wrapper_tsk(const double *para,const void *dats){ //est error
  //plug in para into a our matrix structure
  error_struct *loc_pars = (error_struct *) dats;
  int num=0;
  for(int i=0;i<4;i++) for(int j=0;j<4;j++)
      if(i!=j)
	loc_pars->errors[i][j]=para[num++];
  
  double **loglike3 = loc_pars->loglike3;//this points to allocated but empty data



  double ****errorProbs = loc_pars->errorProbs;

  //compute the probability of the four alleles given major, minor and errors
  
  //error numGeno,numAllele
  generateErrorPointers_tsk(loc_pars->errors,3,4,loc_pars->errorProbs);

  //loglike now contains all genotype likelihoods. dim=nsites times (3*nind)
  //  double **loglike3 = getLikesFullError3Genotypes(loc_pars,errorProbs);
  
  //  fprintf(stderr,"loc_adsf=%f\n",loc_pars->logfact[10]);
  getLikesFullError3Genotypes_tsk(loc_pars->numSites,loc_pars->nInd,loc_pars->major,loc_pars->minor,loc_pars->counts,errorProbs,loc_pars->logfact, loglike3);

  // we now just call em directly instead of going through estFreq
  //get the likelihood for all sites
  double totalLogLike=0;
  for(int s=0;s<loc_pars->numSites;s++){
    int keepInd=0;
    int keeps[loc_pars->nInd];
    for(int i=0;i<loc_pars->nInd;i++)
      if((loglike3[s][i*3+0]+loglike3[s][i*3+1]+loglike3[s][i*3+2])!=0){
	float max=angsd::getMax(loglike3[s][i*3+0],loglike3[s][i*3+1],loglike3[s][i*3+2]);
	if(max<-20)
	  keeps[i] =0;
	else{
	  keepInd++;
	  keeps[i]=1;
	}
      }else{
	keeps[i]=0;
      }
    if(keepInd>0)
      totalLogLike += abcFreq::likeFixedMinor(abcFreq::emFrequency(loglike3[s],loc_pars->nInd,loc_pars->emIter,loc_pars->EM_start,keeps,keepInd),loglike3[s],loc_pars->nInd);
  }

  return totalLogLike;
}


double ****allocErrorProbs(int nG, int nA){
  
  double ****errorProbs = new double***[3];
  for(int iGeno=0;iGeno<nG;iGeno++){
    errorProbs[iGeno]= new double**[nA];
    for(int iMajor=0;iMajor<nA;iMajor++){
      errorProbs[iGeno][iMajor]= new double*[nA];
      for(int iMinor=0;iMinor<nA;iMinor++){
	errorProbs[iGeno][iMajor][iMinor]= new double[nA];
	
      }
    }
  }
  return errorProbs;
}


float *abcError::logfact(int len ) 
/* log(n!) */ 
{
  if(len<2){
    fprintf(stderr,"\t-> You should give a larger value\n");
    exit(0);
  }
  float *res = new float[len];
  res[0]=res[1]=0.0;
  for (int i=2; i<len; i++){
    res[i]=res[i-1]+log((double)i);
    //    fprintf(fpp2,"%f\n",res[i]);
  }
  return res;
}


void abcError::killGlobalErrorProbs(double ****errorProbs){
  if(errorProbs==NULL)
    return;
  int nG=3; //genotypes
  int nA=4; //the four alleles

  for(int iGeno=0;iGeno<nG;iGeno++){
    for(int iMajor=0;iMajor<nA;iMajor++){
      for(int iMinor=0;iMinor<nA;iMinor++)
	delete [] errorProbs[iGeno][iMajor][iMinor];
      delete [] errorProbs[iGeno][iMajor];
    }
    delete [] errorProbs[iGeno];
  }
  
  delete [] errorProbs;
}


void abcError::getOptions(argStruct *arguments){

  doError=angsd::getArg("-doError",doError,arguments);
  if(doError==0)
    return;

  int doMajorMinor = 0;
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  emIter=angsd::getArg("-emIter",emIter,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  minSites=angsd::getArg("-minSites",minSites,arguments);
  minPhat=angsd::getArg("-minPhat",minPhat,arguments);
  eps=angsd::getArg("-eps",eps,arguments);
  
  errorFname = angsd::getArg("-errors",errorFname,arguments);

  int it = arguments->inputtype;
  if(doError && doCounts==0){
    fprintf(stderr,"\t-> Must supply -doCounts to use syk error estimator\n");
    exit(0);
  }

  if((doError &&(it!=INPUT_BAM) )&& (doError && it!=INPUT_PILEUP) ){
    fprintf(stderr,"Error rates can be estimated based on either \n 1) BAM files or pileup files\n");
    exit(0);
  }
  
  if(doError && doMajorMinor!=2){
    fprintf(stderr,"You must infer the major and minor genotypes from the sequencing data (-doMajorMinor 2)\n");
    exit(0);

  }


}

abcError::abcError(const char *outfiles,argStruct *arguments,int inputtype){
  //  outfile1 = outfile2 = NULL;
  outfile2 = NULL;
  minPhat = 0.005;
  eps = 0.001;
  doError=0;
  emIter=100;
  //not used in this class for things other than checking
  doCounts =0;
  minSites = 10000;

  //from command line

  errorFname = NULL;

  EM_start = EM_START;
  emIter = EM_NITER;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doError")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  
  if(doError==0){
    shouldRun[index] = 0;
    return;
  }
  printArg(arguments->argumentFile);  
  void readError(double **errors,const char *fname);
  
  double errorsDefault[4][4]={{0       ,0.00031 , 0.00373 , 0.000664},
			      {0.000737,   0    , 0.000576, 0.001702},
			      {0.001825,0.000386,    0    , 0.000653},
			      {0.00066 ,0.003648, 0.000321,    0    },
  };
  
  errors = new double *[4];
  for(int i=0;i<4;i++){
    errors[i] = new double[4];
    for(int j=0;j<4;j++)
      errors[i][j] = errorsDefault[i][j];
  }

  if(errorFname!=NULL)
    readError(errors,errorFname);


  //make output files
  //  const char* postfix1=".error.chunkordered";
  const char* postfix2=".error.chunkunordered";
  //  outfile1 = aio::openFile(outfiles,postfix1);
  outfile2 = aio::openFile(outfiles,postfix2);

}


abcError::~abcError(){
  if(doError==0)
    return;
  //  if(outfile1) fclose(outfile1);
  if(outfile2) fclose(outfile2);
}



void abcError::print(funkyPars *pars){
  if(doError==0)
    return;

  /*
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      fprintf(outfile1,"%f\t",pars->errors[i][j]);
  fprintf(outfile1,"\n");
  fflush(outfile1);
  */

}


void abcError::clean(funkyPars *pars){
  if(doError==0)
    return;

  /*
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      fprintf(outfile1,"%f\t",pars->errors[i][j]);
  fprintf(outfile1,"\n");
  fflush(outfile1);
  */

}


void abcError::run(funkyPars *pars){
  if(doError==0)
    return;
  assert(pars->counts!=NULL);
    

  //  pars->opt->emIter=emIter;
  if(pars->major==NULL||pars->minor==NULL){
    pars->major = new char[pars->numSites];
    pars->minor = new char[pars->numSites];
    setMajorMinor(pars->counts,pars->major,pars->minor,pars->nInd,pars->numSites);
  }
  for(int i=0;i<pars->numSites;i++){
    if(pars->keepSites[i])
      if(phatFun(pars->counts[i],pars->nInd,eps,pars->major[i],pars->minor[i])<minPhat)
	pars->keepSites[i] =0;
    //fprintf(stderr,"%d %d\n",pars->major[i],pars->minor[i]);
  }
  

 
  consolidate(pars); //this function will merge data, and maybe trigger optimization
 

}

//the errors parameter is not used for now but can be used for supplying the start point
void abcError::likeFixedMinorError_bfgs_tsk2(double **errors,int numSites,int nInd,int *major,int *minor,suint **counts,int emIter,double EM_start){ //est error
  double start[12];
  double lbound[12];
  double ubound[12];
  int lims[12];
  for(int i=0;i<12;i++){
    start[i]=0.01;
    lbound[i]=0.0000001;
    ubound[i]=0.1;
    lims[i]=2;
  }
  //prepare struct that are sent to optimziation
  error_struct *es = new error_struct;
  es->logfact = logfact(LOGMAX);
  es->errors =new double *[4]; //errors;//<-this will be overridden by 'start' above
  for(int i=0;i<4;i++)
    es->errors[i] = new double[4];
  es->numSites = numSites;
  es->nInd = nInd;
  es->major=major;
  es->minor=minor;
  es->counts=counts;
  es->errorProbs = allocErrorProbs(3,4);
  es->loglike3 = new double*[es->numSites];
  for(int i=0;i<es->numSites;i++)
    es->loglike3[i] = new double[3*es->nInd];

  
  es->EM_start = EM_START;
  es->emIter = EM_NITER;

  double res=findmax_bfgs(12,start,(void *)es,likeFixedMinorError_wrapper_tsk,NULL,lbound,ubound,lims,-1);
  //  fprintf(stderr,"[%s] \t like %f isoptim is: %f\n",__FUNCTION__,res,start[0]);	

  int num=0;
  //plug the optimized start value back into our matrix form
  for(int i=0;i<4;i++){//printout
    for(int j=0;j<4;j++){
	if(i==j)
	  es->errors[i][j] = 0.0;
	else
	  es->errors[i][j] = start[num++];

    }
  }
  
  
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      fprintf(outfile2,"%f\t",es->errors[i][j]);
  fprintf(outfile2,"\n");
  fflush(outfile2);
  //cleanup
  killGlobalErrorProbs(es->errorProbs);
  for(int i=0;i<3;i++)
    delete [] es->loglike3[i];
  delete [] es->loglike3;
  delete [] es->logfact;
  for(int i=0;i<4;i++)
    delete [] es->errors[i];
  delete [] es->errors;
  delete es;
}




void abcError::setMajorMinor(suint **counts,char *major,char *minor,int nInd,int nSites){
  for(int s=0;s<nSites;s++){//loop over sites
    //get basecounts
    int bases[4] = {0,0,0,0};
    for(int i=0;i<nInd;i++)
      for(int j=0;j<4;j++)
	bases[j] += counts[s][i*4+j];
    
    //find major
    int maj = 0;
    for(int i=1;i<4;i++)
      if (bases[i]>bases[maj])
	maj = i;
    //maj is now the major allele (most frequent)
    
    //find minor
    int temp=0;
    int mino= maj;
    for(int i=0;i<4;i++){
      if(maj==i) //we should check the against the major allele      
        continue;
      else if (bases[i]>temp){
        mino = i;
        temp=bases[i];
      }
    }
    major[s] = maj;
    minor[s] = mino;
    
  }


}

double **** abcError::generateErrorPointers(double **error,int nG,int nA) {
  //  if(errorProbs==NULL)//
  double ****errorProbs= new double***[3];
  //estimates the probabilities for the four bases given the major and minor
  //The function is only called once and is of no importance speedwise
  //The "errorProbs" is global 

  for(int iGeno=0;iGeno<nG;iGeno++) {
    errorProbs[iGeno]= new double**[nA];
    for(int iMajor=0;iMajor<nA;iMajor++){
      errorProbs[iGeno][iMajor]= new double*[nA];
      for(int iMinor=0;iMinor<nA;iMinor++){
        errorProbs[iGeno][iMajor][iMinor]= new double[nA];	
	if(iMajor==iMinor)
	  continue;
        int other1=-1;
        int other2=-1;
        for(int i=0;i<3;i++){
          other1=i;
          if(i!=iMajor&&i!=iMinor)
            break;
        }
        for(int i=0;i<4;i++){
          other2=i;
          if(i!=iMajor&&i!=iMinor&&i!=other1)
            break;
        }

	errorProbs[iGeno][iMajor][iMinor][iMajor]=angsd::addProtect2(log(1-iGeno*0.5)+log(1-error[iMajor][iMinor]-error[iMajor][other1]-error[iMajor][other2]),log(iGeno)-log(2)+log(error[iMinor][iMajor]));
//prob Major
	errorProbs[iGeno][iMajor][iMinor][iMinor]=angsd::addProtect2(log(iGeno)-log(2)+log(1-error[iMinor][iMajor]-error[iMinor][other1]-error[iMinor][other2]),log(1-iGeno*0.5)+log(error[iMajor][iMinor]));
//prob Minor
	errorProbs[iGeno][iMajor][iMinor][other1]=angsd::addProtect2(log(iGeno)-log(2)+log(error[iMinor][other1]),log(1-iGeno*0.5)+log(error[iMajor][other1]));
//prob other1
	errorProbs[iGeno][iMajor][iMinor][other2]=angsd::addProtect2(log(iGeno)-log(2)+log(error[iMinor][other2]),log(1-iGeno*0.5)+log(error[iMajor][other2]));
//prob other2
	    double sum=0;
	    for(int ia=0;ia<4;ia++)
	      sum+=exp(errorProbs[iGeno][iMajor][iMinor][ia]);
	    if((sum>1.00005||sum<0.999999)&&iMajor!=iMinor){//never happens, just a check
	      for(int j=0;j<4;j++)
		fprintf(stderr,"%f\t ",errorProbs[iGeno][iMajor][iMinor][j]);
	      fprintf(stderr,"G %d M %d m %d o1 %d o2 %d sum %f %f\n",iGeno,iMajor,iMinor,other1,other2,sum,angsd::addProtect2(log((1-iGeno*0.5)*(1-error[iMajor][iMinor]-error[iMajor][other1]-error[iMajor][other2])),log(iGeno*0.5*error[iMajor][iMinor])));
	      
	    }


      }
    }
  }

  return errorProbs;
}
