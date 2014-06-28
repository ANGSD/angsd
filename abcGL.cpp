/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd

  
  This class will calculate the GL in 4 differnt ways

  1) SAMtools 0.1.16+ version 
  2) Simple GATK model
  3) SOAPsnp
  4) SYK

  4 different output formats are supplied
  
  1) binary 10xdouble persample
  2) beagle output (requires estimation of major/minor)\
  3) binary beagle
  4) text output of the 10 llhs persample

*/



#include <cmath>
#include <zlib.h>
#include <assert.h>
#include "kstring.h"//<-used for buffered output when dumping beagle 0.204
#include "bfgs.h"
#include "analysisFunction.h"
#include "abc.h"

#include "abcGL.h"
#include "abcError.h"
extern int refToInt[256];

static float *logfactorial=NULL;

void readError(double **errors,const char *fname){
  fprintf(stderr,"will try to read errorestimates from file:%s\n",fname);
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  
  char buf[LENS];
  double res[16];
  for(int i=0;i<16;i++) res[i] = 0;

  int nLines =0;
  while(fgets(buf,LENS,fp)){
    res[0] += atof(strtok(buf," \t\n"));
    for(int i=1;i<16;i++)
      res[i] += atof(strtok(NULL," \t\n"));
    nLines ++;
  }
  for(int j=0;j<16;j++)
    fprintf(stderr,"%f\t",res[j]);
  fprintf(stderr,"\nEstimating errors using nChunks:%d\n",nLines);
  int pos =0;
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      errors[i][j] = res[pos++]/(1.0*nLines);
  if(fp) fclose(fp);
}


void abcGL::printArg(FILE *argFile){
  fprintf(argFile,"---------------------\n%s:\n",__FILE__);

  fprintf(argFile,"\t-GL=%d: \n",GL);
  fprintf(argFile,"\t1: SAMtools\n");
  fprintf(argFile,"\t2: GATK\n");
  fprintf(argFile,"\t3: SOAPsnp\n");
  fprintf(argFile,"\t4: SYK\n");

  fprintf(argFile,"\t-trim\t\t%d\t\t(zero means no trimming)\n",trim);
  fprintf(argFile,"\t-tmpdir\t\t%s/\t(used by SOAPsnp)\n",angsd_tmpdir);
  fprintf(argFile,"\t-errors\t\t%s\t\t(used by SYK)\n",errorFname);
  fprintf(argFile,"\t-minInd\t\t%d\t\t(0 indicates no filtering)\n",minInd);
  fprintf(argFile,"\n");
  fprintf(argFile,"Filedumping:\n");
  fprintf(argFile,"\t-doGlf\t%d\n",doGlf);
  fprintf(argFile,"\t1: binary glf (10 log likes)\t%s\n",postfix);
  fprintf(argFile,"\t2: beagle likelihood file\t%s\n",beaglepostfix);
  fprintf(argFile,"\t3: binary 3 times likelihood\t%s\n",postfix);
  fprintf(argFile,"\t4: text version (10 log likes)\t%s\n",postfix);
  fprintf(argFile,"\n");

}

void abcGL::getOptions(argStruct *arguments){

  //parse all parameters that this class could use
  GL=angsd::getArg("-GL",GL,arguments);
  trim = angsd::getArg("-trim",trim,arguments);
  angsd_tmpdir = angsd::getArg("-tmpdir",angsd_tmpdir,arguments);
  doGlf=angsd::getArg("-doGlf",doGlf,arguments);
  errorFname = angsd::getArg("-errors",errorFname,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);


  int doCounts=0;
  int doMajorMinor =0;
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(arguments->inputtype==INPUT_GLF&&GL!=0){
    fprintf(stderr,"Can't calculate genotype likelihoods from glf files\n");
    exit(0);
  }
  if(arguments->inputtype==INPUT_GLF)
    return;
  if(doGlf&&GL==0){
    fprintf(stderr,"\t-> You need to choose a genotype likelihood model -GL for dumping genotype likelihoods\n");
    exit(0);
  }
  if(GL==0&&doGlf==0){
    shouldRun[index] =0;
    return;
  }
  if(GL==0)
    return;

  if(GL<0||GL>4){
    fprintf(stderr,"\t-> You've choosen a GL model=%d, only 1,2,3,4 are implemented\n",GL);
    exit(0);
  }
  if(GL==4&&(doCounts==0)){
    fprintf(stderr,"\t-> Must supply -doCounts for SYK model\n");
    exit(0);
  }
  if((doGlf==2||doGlf==3) && doMajorMinor==0){
    fprintf(stderr,"\t-> For dumping beaglestyle output you need to estimate major/minor: -doMajorMinor\n");
    exit(0);
  }
  if(arguments->inputtype==INPUT_BEAGLE&&doGlf){
    fprintf(stderr,"\t-> cannot output likelihoods (doGlf) when input is beagle\n");
    exit(0);
  }
 
  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"Error: Likelihoods can only be estimated based on SOAP input and uppile input\n");
    exit(0);
  }


  printArg(arguments->argumentFile);
}

abcGL::abcGL(const char *outfiles,argStruct *arguments,int inputtype){
  errors = NULL;
  postfix = ".glf.gz";
  beaglepostfix = ".beagle.gz";
  

  trim =0;
  GL=0;
  doGlf=0;
  errorFname = NULL;
  errorProbs = NULL;
  GL=0;
  minInd=0;
  angsd_tmpdir = strdup("angsd_tmpdir");
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-GL")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  printArg(arguments->argumentFile);

  //  if(GL==0)
  //  return;
  if(GL==1)
    bam_likes_init();
  else if(GL==2)
    gatk_init();
  else if(GL==3){
    soap.init(arguments->nInd,angsd_tmpdir);
    if(soap.doRecal)
      fprintf(stderr,"[%s] Will calculate recalibration matrices, please don't do any other analysis\n",__FILE__);
    else
      fprintf(stderr,"[%s] Will use precalculated calibration matrices\n",__FILE__);

  }else if(GL==4) {
    //default errormatrix
    double errorsDefault[4][4]={{0       ,0.00031 , 0.00373 , 0.000664},
				{0.000737,   0    , 0.000576, 0.001702},
				{0.001825,0.000386,    0    , 0.000653},
				{0.00066 ,0.003648, 0.000321,    0    },
    };
    //allocate and plug in default values
    errors = new double *[4];
    for(int i=0;i<4;i++){
      errors[i] = new double[4];
      for(int j=0;j<4;j++)
	errors[i][j] = errorsDefault[i][j];
    }
    if(errorFname!=NULL)
      readError(errors,errorFname);
    errorProbs = abcError::generateErrorPointers(errors,3,4);
  }
  
  gzoutfile = Z_NULL;
  bufstr.s=NULL; bufstr.l=bufstr.m=0;// <- used for buffered output 
  bufstr.l=0;
  if(doGlf){
 
    if(doGlf!=2)
      gzoutfile = aio::openFileGz(outfiles,postfix,GZOPT);
    else{
      gzoutfile = aio::openFileGz(outfiles,beaglepostfix,GZOPT);
      
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
      gzwrite(gzoutfile,bufstr.s,bufstr.l);
    }
 
  }

}


abcGL::~abcGL(){

  free(angsd_tmpdir);
  
  if(GL==0&&doGlf==0)
    return;
  else if(GL==1)
    bam_likes_destroy();
  else if(GL==2)
    gatk_destroy();
  else if(GL==4)
    abcError::killGlobalErrorProbs(errorProbs);
  if(doGlf)    gzclose(gzoutfile);
    
  
  if(bufstr.s!=NULL)
    free(bufstr.s);

  if(errors){
    for(int i=0;i<4;i++)
      delete [] errors[i];
    delete [] errors;
  }
  delete [] logfactorial;
}

void abcGL::clean(funkyPars *pars){

  if(pars->likes!=NULL){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->likes[i];
    delete [] pars->likes;
  }
  
}


void abcGL::print(funkyPars *pars){
  if(doGlf)
    printLike(pars);
}


void abcGL::run(funkyPars *pars){
  assert(pars!=NULL);
  
  if(GL==0)
    return;
  //assert(pars->chk!=NULL);
  double **likes = NULL;
  if(soap.doRecal!=1)
    likes = new double*[pars->chk->nSites];
  
  if(GL==1)
    call_bam(pars->chk,likes,trim);
  else if(GL==2)
    call_gatk(pars->chk,likes,trim);
  else if(GL==3){
    soap.run(pars->chk,likes,pars->ref,trim);
    //we dont estimate GL but make a calibration matrix
    if(soap.doRecal==1)
      return;

  }else if(GL==4)
    getLikesFullError10Genotypes(pars->numSites,pars->nInd,pars->counts,errorProbs,pars->keepSites,likes);

  pars->likes = likes;
  


  /*
    if trimming has been requested, then some site might not contain data,
    we therefore set keepsites to zero for these sites
    while we are at it, lets also count the effective sample size persite
  */
  if(1){
    for(int s=0;s<pars->numSites;s++){
      //     fprintf(stderr,"keepSites[%d]=%d\n",s,pars->keepSites[s]);
      if(pars->keepSites[s]==0)
	continue;
      int efSize=0;
      for(int i=0;i<pars->nInd;i++){
	for(int ii=1;ii<10;ii++){
	  if(pars->likes[s][i*10+ii]!=pars->likes[s][i*10+0]){
	    efSize++;
	    //	    break;
	  }
	}
      }
      pars->keepSites[s] = efSize;
 
      if(minInd!=0&&minInd>efSize)
	pars->keepSites[s] = 0;
    }
  }

  //rescale the genotype likelihoods to loglike ratios.
  if(1){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
 
      for(int i=0;i<pars->nInd;i++)
	angsd::logrescale(pars->likes[s] +i*10,10);

    }
  }


}

void abcGL::getLikesFullError10Genotypes(int numSites,int nInd,suint **counts,double ****errorProbs,int *keepSites,double **loglikes) {

  //only calculate this once

  if(logfactorial==NULL)//dont bother populating if exists.
    logfactorial=abcError::logfact(LOGMAX); //calculate log factorials

  double *logError;

  for(int s=0;s<numSites;s++){
    loglikes[s] = new double [10*nInd];
    if(keepSites[s]==0)
      continue;
    for(int allele1=0;allele1<4;allele1++) {
      for(int allele2=allele1;allele2<4;allele2++){
	int Gindex=angsd::majorminor[allele1][allele2];
	int geno=0;
	int m2=allele2;
	if(allele1!=allele2)
	  geno++;
	else{//total grimt must redo
	m2++;
	 if(m2>3)
	  m2=0;
	}
	logError=errorProbs[geno][allele1][m2];
	for(int i=0;i<nInd;i++){
	  loglikes[s][i*10+Gindex]=logfactorial[counts[s][i*4+0]+counts[s][i*4+1]+counts[s][i*4+2]+counts[s][i*4+3]]; //should be computed before these loops for faster implimentation
	  for(int j=0;j<4;j++)
	    loglikes[s][i*10+Gindex]+=-logfactorial[counts[s][i*4+j]]+counts[s][i*4+j]*logError[j];
	}
      }
    }
  }

}



void abcGL::printLike(funkyPars *pars) {
  assert(pars->likes!=NULL);

  
  if(doGlf==1){
    //glffinn format
    for(int i=0;i<pars->numSites;i++){
      if(pars->keepSites[i]==0)
	continue;
      gzwrite(gzoutfile,pars->likes[i],sizeof(double)*10*pars->nInd);
    }
  }
  else if(doGlf==2){
    //beagle format
    bufstr.l = 0; //set tmpbuf beginning to zero
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue;
      
      kputs(header->name[pars->refId],&bufstr);
      kputc('_',&bufstr);
      kputw(pars->posi[s]+1,&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->major[s],&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->minor[s],&bufstr);

      int major = pars->major[s];
      int minor = pars->minor[s];
      assert(major!=4&&minor!=4);
	
      for(int i=0;i<pars->nInd;i++) {
	
	double norm=exp(pars->likes[s][i*10+angsd::majorminor[major][major]])+exp(pars->likes[s][i*10+angsd::majorminor[major][minor]])+exp(pars->likes[s][i*10+angsd::majorminor[minor][minor]]);
	double val1 = exp(pars->likes[s][i*10+angsd::majorminor[major][major]])/norm;
	double val2 = exp(pars->likes[s][i*10+angsd::majorminor[major][minor]])/norm;
	double val3 = exp(pars->likes[s][i*10+angsd::majorminor[minor][minor]])/norm;
	ksprintf(&bufstr, "\t%f",val1);
	ksprintf(&bufstr, "\t%f",val2);
	ksprintf(&bufstr, "\t%f",val3);
      }
      
      kputc('\n',&bufstr);

    }
    gzwrite(gzoutfile,bufstr.s,bufstr.l);
  }
  else if(doGlf==3) { //FGV v0.208 Aug,28
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0) //TSK 0.441 sep 25
	continue;
      int major = pars->major[s];
      int minor = pars->minor[s] ;
      assert(major!=4&&minor!=4);

      for(int i=0;i<pars->nInd;i++) {
	double dump[3];
	dump[0] = pars->likes[s][i*10+angsd::majorminor[major][major]] ;
	dump[1] = pars->likes[s][i*10+angsd::majorminor[major][minor]] ;
	dump[2] = pars->likes[s][i*10+angsd::majorminor[minor][minor]] ;
	gzwrite(gzoutfile,dump,3*sizeof(double));
      }
    }
  } else if(doGlf==4){
    bufstr.l=0;
    //otherwise print textoutput
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      kputs(header->name[pars->refId],&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->posi[s]+1,&bufstr);
      for(int i=0;i<10*pars->nInd;i++)      
	ksprintf(&bufstr, "\t%f",pars->likes[s][i]);

      kputc('\n',&bufstr);
    }
    gzwrite(gzoutfile,bufstr.s,bufstr.l);
  }


}

