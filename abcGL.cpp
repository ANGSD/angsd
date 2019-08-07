/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd

  
  This class will calculate the GL in 4 differnt ways

  1) SAMtools 0.1.16+ version 
  2) Simple GATK model
  3) SOAPsnp
  4) SYK
  5) new fancy method
  6) simple sample gls
  4 different output formats are supplied
  
  1) binary 10xdouble persample
  2) beagle output (requires estimation of major/minor)\
  3) binary beagle
  4) text output of the 10 llhs persample

*/



#include <cmath>
#include <assert.h>
#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "abc.h"

#include "abcGL.h"
#include "abcError.h"
#include "phys_likes.h"
#include "abcMajorMinor.h"

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
  fprintf(argFile,"\t5: phys\n");
  fprintf(argFile,"\t6: Super simple sample an allele type GL. (1.0,0.5,0.0)\n");
  fprintf(argFile,"\t7: outgroup gls\n");
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
  fprintf(argFile,"\t5: binary saf files (usefull for realSFS)\t%s\n",postfix);
  fprintf(argFile,"\n");

}

void abcGL::getOptions(argStruct *arguments){

  //parse all parameters that this class could use
  GL=angsd::getArg("-GL",GL,arguments);

  if(0&&GL==0)//DRAGON
    return;
  
  trim = angsd::getArg("-trim",trim,arguments);
  angsd_tmpdir = angsd::getArg("-tmpdir",angsd_tmpdir,arguments);
  doGlf=angsd::getArg("-doGlf",doGlf,arguments);
  errorFname = angsd::getArg("-errors",errorFname,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);

  // should parse a list of quals e.g. 0,20,30. Meaning three bins: 0-19; 20-29; 30+

  int doCounts=0;
  int doMajorMinor =0;
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(GL!=0 && (arguments->inputtype==INPUT_GLF || arguments->inputtype==INPUT_GLF3 || arguments->inputtype==INPUT_VCF_GL)){
    fprintf(stderr,"Can't calculate genotype likelihoods from -glf/-glf3/VCF files\n");
    exit(0);
  }
  if(arguments->inputtype==INPUT_GLF||arguments->inputtype==INPUT_GLF3||arguments->inputtype==INPUT_VCF_GL||arguments->inputtype==INPUT_GLF10_TEXT)
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
  if(GL==7){
    fprintf(stderr,"\t-> GL model=%d (outgroup gls) is BETA\n", GL);
    // return;
  }
  
  if(( GL<0||GL>7 )) {
    fprintf(stderr,"\t-> You've choosen a GL model=%d, only 1,2,3,4,5,6,7 are implemented\n",GL);
    exit(0);
  }
  if(GL==4&&(doCounts==0)){
    fprintf(stderr,"\t-> Must supply -doCounts 1 for SYK model\n");
    exit(0);
  }
  if(GL==6&&(doCounts==0)){
    fprintf(stderr,"\t-> Must supply -doCounts 1 for -gl 6\n");
    exit(0);
  }
  /*
  if(doGlf==2){
    fprintf(stderr,"\t-> BEAGLE format 3.0 is deprecated\n\t-> Consider using -doVCF 1\n");
    for(int j=3;j>0;j--){
      fprintf(stderr,"\t-> Program will continue in %d seconds    \n",j);fflush(stderr);
      sleep(1);
    }

  }
  */
  if((doGlf==2||doGlf==3) && doMajorMinor==0){
    fprintf(stderr,"\t-> For dumping beaglestyle output you need to estimate major/minor: -doMajorMinor\n");
    exit(0);
  }
  if(arguments->inputtype==INPUT_BEAGLE&&doGlf){
    fprintf(stderr,"\t-> cannot output likelihoods (doGlf) when input is beagle\n");
    exit(0);
  }
 
  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"Error: Likelihoods can only be estimated based on BAM input and uppile input\n");
    exit(0);
  }


  printArg(arguments->argumentFile);
}

abcGL::abcGL(const char *outfiles,argStruct *arguments,int inputtype){
  nnnSites =0;
  outfileSAF = NULL;
  outfileSAFPOS = NULL;
  outfileSAFIDX = NULL;
  errors = NULL;
  postfix = ".glf.gz";
  beaglepostfix = ".beagle.gz";
  tmpChr = NULL;
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

  //  if(GL==0) //why?
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
  }else if(GL==5){
      phys_init(arguments->nams);
  }else if(GL==6){
    simple_init();
  }else if(GL==7){
    ancestral_lik.init(arguments->nInd, angsd_tmpdir);
    if(ancestral_lik.doRecal){
      fprintf(stderr, "\t-> [%s] Will generate the counts file, please do not run other analyses\n", __FILE__);
    } else {
      fprintf(stderr, "\t-> [%s] Will use error matrix already estimated with python script on the counts matrix generated in the previous run\n", __FILE__);
    } // ancestral_init();
  }
  
  gzoutfile = gzoutfile2 = NULL;
  bufstr.s=NULL; bufstr.l=bufstr.m=0;// <- used for buffered output 
  const char *SAF = ".saf.gz";
  const char *SAFPOS =".saf.pos.gz";
  const char *SAFIDX =".saf.idx";

  if(doGlf==5){
    outfileSAF =  aio::openFileBG(outfiles,SAF);
    outfileSAFPOS =  aio::openFileBG(outfiles,SAFPOS);
    outfileSAFIDX = aio::openFile(outfiles,SAFIDX);
    char buf[8]="safv3";
    aio::bgzf_write(outfileSAF,buf,8);
    aio::bgzf_write(outfileSAFPOS,buf,8);
    fwrite(buf,1,8,outfileSAFIDX);
    offs[0] = bgzf_tell(outfileSAFPOS);
    offs[1] = bgzf_tell(outfileSAF);
    size_t tt = 9;
    fwrite(&tt,sizeof(tt),1,outfileSAFIDX);
  }
  
  
  if(doGlf){
 
    if(doGlf!=2){
      gzoutfile = aio::openFileBG(outfiles,postfix);
      if(doGlf==3)
	gzoutfile2 = aio::openFileBG(outfiles,".glf.pos.gz");
    }else{
      gzoutfile = aio::openFileBG(outfiles,beaglepostfix);
      
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
      aio::bgzf_write(gzoutfile,bufstr.s,bufstr.l);bufstr.l=0;
    }
 
  }

}


abcGL::~abcGL(){
  if(GL>0&&doGlf==5){
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
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
    }//else
    // fprintf(stderr,"enpty chr\n");
    //reset
    offs[0] = bgzf_tell(outfileSAFPOS);
    offs[1] = bgzf_tell(outfileSAF);
    nnnSites=0;
  }
  
  free(angsd_tmpdir);
  
  if(GL==0&&doGlf==0)
    return;
  else if(GL==1)
    bam_likes_destroy();
  else if(GL==2)
    gatk_destroy();
  else if(GL==4)
    abcError::killGlobalErrorProbs(errorProbs);
  else if(GL==5)
    phys_destroy();
  else if(GL==6)
    simple_destroy();

  if(doGlf)    bgzf_close(gzoutfile);
    
  if(gzoutfile!=NULL)
    bgzf_close(gzoutfile2);

  if(bufstr.s!=NULL)
    free(bufstr.s);

  if(errors){
    for(int i=0;i<4;i++)
      delete [] errors[i];
    delete [] errors;
  }
  delete [] logfactorial;
  if(doGlf==5){
    if(outfileSAF) bgzf_close(outfileSAF);;
    if(outfileSAFPOS) bgzf_close(outfileSAFPOS);
    if(outfileSAFIDX) fclose(outfileSAFIDX);
  }
}

void abcGL::clean(funkyPars *pars){

  if(pars->likes!=NULL){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->likes[i];
    delete [] pars->likes;
    pars->likes=NULL;
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
  if(soap.doRecal!=1 || ancestral_lik.doRecal!=1)
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
  else if(GL==5)
    call_phys(pars->chk,likes,trim);
  else if(GL==6){
    call_simple(pars->counts,pars->keepSites,likes,pars->numSites,pars->nInd);
  }else if(GL==7){
    ancestral_lik.run(pars->chk, likes, pars->ref, pars->anc, pars->keepSites, trim);
    if(ancestral_lik.doRecal==1){
      return;
    }
  }
  pars->likes = likes;
  


  /*
    if trimming has been requested, then some site might not contain data,
    we therefore set keepsites to zero for these sites
    while we are at it, lets also count the effective sample size persite
  */
  if(1){
    for(int s=0;s<pars->numSites;s++){
      
      if(pars->keepSites[s]==0)
	continue;
      int efSize=0;
      for(int i=0;i<pars->nInd;i++){

	for(int ii=1;ii<10;ii++){
	  if(pars->likes[s][i*10+ii]!=pars->likes[s][i*10+0]){
	    efSize++;
	    break;
	  }
	}
      }
      pars->keepSites[s] = efSize;
 
      if(minInd!=0&&minInd>efSize)
	pars->keepSites[s] = 0;
      //      fprintf(stderr,"posi:%d keepSites2[%d]=%d\n",pars->posi[s]+1,s,pars->keepSites[s]);
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
      aio::bgzf_write(gzoutfile,pars->likes[i],sizeof(double)*10*pars->nInd);
    }
  }
  else if(doGlf==2){
    //beagle format
    bufstr.l = 0; //set tmpbuf beginning to zero
    for(int s=0;s<pars->numSites;s++) {
      lh3struct *lh3 = (lh3struct*) pars->extras[index+1];
      if(pars->keepSites[s]==0||lh3->hasAlloced[s]==0)
	continue;
      
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

     
      for(int i=0;i<pars->nInd;i++) {
	double val[3];
	val[0]= exp(lh3->lh3[s][i*3+0]);
	val[1]= exp(lh3->lh3[s][i*3+1]);
	val[2]= exp(lh3->lh3[s][i*3+2]);
	angsd::norm(val,3);
	ksprintf(&bufstr, "\t%f",val[0]);
	ksprintf(&bufstr, "\t%f",val[1]);
	ksprintf(&bufstr, "\t%f",val[2]);
	assert(!std::isnan(val[0]));
	assert(!std::isnan(val[1]));
	assert(!std::isnan(val[2]));
      }

      if(bufstr.l!=0)
	kputc('\n',&bufstr);

    }
    aio::bgzf_write(gzoutfile,bufstr.s,bufstr.l);bufstr.l=0;
  }
  else if(doGlf==3) { //FGV v0.208 Aug,28
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0) //TSK 0.441 sep 25
	continue;
      char major = pars->major[s];
      char minor = pars->minor[s] ;
      assert(major!=4&&minor!=4);

      for(int i=0;i<pars->nInd;i++) {
	double dump[3];
	dump[0] = pars->likes[s][i*10+angsd::majorminor[major][major]] ;
	dump[1] = pars->likes[s][i*10+angsd::majorminor[major][minor]] ;
	dump[2] = pars->likes[s][i*10+angsd::majorminor[minor][minor]] ;
	aio::bgzf_write(gzoutfile,dump,3*sizeof(double));
      }
      bufstr.l=0;
      ksprintf(&bufstr,"%s\t%d\t",header->target_name[pars->refId],pars->posi[s]+1);
      ksprintf(&bufstr,"%c\t%c\n",intToRef[major],intToRef[minor]);
      aio::bgzf_write(gzoutfile2,bufstr.s,bufstr.l);bufstr.l=0;
    }
  } else if(doGlf==4){
    bufstr.l=0;
    //otherwise print textoutput
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      kputs(header->target_name[pars->refId],&bufstr);
      kputc('\t',&bufstr);
      kputw(pars->posi[s]+1,&bufstr);
      for(int i=0;i<10*pars->nInd;i++)      
	ksprintf(&bufstr, "\t%f",pars->likes[s][i]);

      kputc('\n',&bufstr);
    }
    aio::bgzf_write(gzoutfile,bufstr.s,bufstr.l);bufstr.l=0;
  }else if(doGlf==5){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      nnnSites++;
      float tmp[10*pars->nInd];
      for(int i=0;i<10;i++)
	tmp[i] = pars->likes[s][i];
      aio::bgzf_write(outfileSAF,tmp,sizeof(float)*10*pars->nInd);
      int mypos = pars->posi[s];
      aio::bgzf_write(outfileSAFPOS,&mypos,sizeof(int));
    }


  }


}

void abcGL::changeChr(int refId) {
  fprintf(stderr,"Charnge chr:%d\n",refId);
  if(GL>0&&doGlf==5){
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
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
    }//else
    // fprintf(stderr,"enpty chr\n");
    //reset
    offs[0] = bgzf_tell(outfileSAFPOS);
    offs[1] = bgzf_tell(outfileSAF);
    nnnSites=0;
    
    free(tmpChr);
    tmpChr = strdup(header->target_name[refId]);


  }
}
