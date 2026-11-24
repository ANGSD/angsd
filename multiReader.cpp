#include "abc.h"
#include "shared.h"
#include "multiReader.h"
#include "parseArgs_bambi.h"
#include "version.h"
#include "sample.h"
#include "aio.h"

//returns NULL if the fai file is malformed
bam_hdr_t *getHeadFromFai(const char *fname){
  std::vector<char *> chrs;
  std::vector<uint32_t> lengths;
  FILE *fp = aio::getFILE(fname,"r");
  char *buf = NULL;
  size_t bufLen = 0;
  int lineNo = 0;
  int bad = 0;
  while(getline(&buf,&bufLen,fp)!=-1){
    lineNo++;
    char *tok = strtok(buf,"\t \r\n");
    if(!tok)//skip empty lines
      continue;
    char *name = tok;
    tok = strtok(NULL,"\t \r\n");
    char *end = NULL;
    unsigned long len = tok?strtoul(tok,&end,10):0;
    if(!tok||end==tok||*end!='\0'||len>UINT32_MAX){
      if(len>UINT32_MAX)
	fprintf(stderr,"\t-> fai file: \'%s\' line %d: reference \'%s\' is longer than %u bp, which is not supported\n",fname,lineNo,name,UINT32_MAX);
      else
	fprintf(stderr,"\t-> fai file: \'%s\' looks malformed at line %d (reference: \'%s\'), expected <name><TAB><length>\n",fname,lineNo,name);
      bad = 1;
      break;
    }
    chrs.push_back(strdup(name));
    lengths.push_back((uint32_t) len);
  }
  free(buf);
  fclose(fp);
  if(!bad&&chrs.size()==0){
    fprintf(stderr,"\t-> fai file: \'%s\' contains no references\n",fname);
    bad = 1;
  }
  if(bad){
    for(size_t i=0;i<chrs.size();i++)
      free(chrs[i]);
    return NULL;
  }
  bam_hdr_t *ret = bam_hdr_init();
  ret->l_text = strlen(fname);
  ret->text =(char*)  malloc(strlen(fname)+1);
  ret->text = strcpy(ret->text,fname);
  ret->n_targets = chrs.size();
  ret->target_len = (uint32_t*) malloc(sizeof(uint32_t)*chrs.size());
  ret->target_name = (char**) malloc(sizeof(char*)*chrs.size());
  for(size_t i=0;i<chrs.size();i++){
    ret->target_len[i] = lengths[i];
    ret->target_name[i] = chrs[i];//header takes ownership
  }
  return ret;
}

aMap *buildRevTable(const bam_hdr_t *hd){
  aio::doAssert(hd!=NULL,AT);
  aMap *ret = new aMap;
  for(int i=0;i<hd->n_targets;i++){
    char *name = strdup(hd->target_name[i]);
    if(!ret->insert(std::pair<char *,int>(name,i)).second){
      fprintf(stderr,"\t-> Reference name: \'%s\' occurs more than once in the header/fai file. Will exit\n",name);
      exit(1);
    }
  }
  return ret;
}


void multiReader::printArg(FILE *argFile,argStruct *args){
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
  fprintf(argFile,"\t-nLines\t%d\t(Number of lines to read)\n",nLines);
  fprintf(argFile,"\t-beagle\t%s\t(Beagle Filename (can be .gz), read as genotype probabilities)\n",fname);
  fprintf(argFile,"\t-beagle-GL\t%s\t(Beagle Filename (can be .gz), read as genotype likelihoods, e.g. from -doGlf 2)\n",fname);
  fprintf(argFile,"\t-bgen\t%s\t(Bgen Filename)\n",fname);
  fprintf(argFile,"\t-vcf-GL\t%s\t(vcf Filename (can be bcf compressed or uncompressed))\n",fname);
  fprintf(argFile,"\t-vcf-PL\t%s\t(vcf Filename (can be bcf compressed or uncompressed))\n",fname);
  fprintf(argFile,"\t-vcf-GP\t%s\t(vcf Filename (can be bcf compressed or uncompressed))\n",fname);
  fprintf(argFile,"\t-glf\t%s\t(glf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-glf3\t%s\t(glf3 Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-glf10_text\t%s\t(text glf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-pileup\t%s\t(pileup Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-intName %d\t(Assume First column is chr_position)\n",intName);
  fprintf(argFile,"\t-isSim\t%d\t(Simulated data assumes ancestral is A)\n",isSim);
  fprintf(argFile,"\t-nInd\t%d\t\t(Number of individuals)\n",nInd);
  fprintf(argFile,"\t-minQ\t%d\t(minimum base quality; only used in pileupreader)\n",minQ);
  fprintf(argFile,"\t-fai\t%s\t(fai file)\n",args->fai);
}


void multiReader::getOptions(argStruct *arguments){
  nLines=angsd::getArg("-nLines",nLines,arguments);
  arguments->nLines = nLines;
  fname=angsd::getArg("-beagle",fname,arguments);
  fname=angsd::getArg("-beagle-GL",fname,arguments);
  fname=angsd::getArg("-bgen",fname,arguments);
  fname=angsd::getArg("-pileup",fname,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  fname=angsd::getArg("-glf",fname,arguments);
  fname=angsd::getArg("-glf3",fname,arguments);
  fname=angsd::getArg("-vcf-GL",fname,arguments);
  fname=angsd::getArg("-vcf-GP",fname,arguments);
  fname=angsd::getArg("-vcf-pl",fname,arguments);
  fname=angsd::getArg("-glf10_text",fname,arguments);
  intName=angsd::getArg("-intName",intName,arguments);
  isSim=angsd::getArg("-isSim",isSim,arguments);
  nInd=angsd::getArg("-nInd",nInd,arguments);
  arguments->nInd = nInd;
  //read fai if suppplied (used by other than native bam readers)

  char *tmptmp=NULL;
  tmptmp=angsd::getArg("-fai",arguments->fai,arguments);
  if(tmptmp){
    free(arguments->fai);
    arguments->fai=tmptmp;
  }
  
  printArg(arguments->argumentFile,arguments);
}


multiReader::multiReader(argStruct *args_arg){
  gz=Z_NULL;
  myglf=NULL;myglf_text=NULL;myvcf=NULL;mpil=NULL;bglObj=NULL;bgenObj=NULL;
  revMap=NULL;
  pl_or_gl=-1;
  
  nLines=50;
  
  fname=NULL;
  intName=1;
  minQ = MINQ;

  nInd =0;
  isSim =0;

  args=NULL;
  args = args_arg;
  void printTime(FILE *fp);
  printTime(args->argumentFile); 

  //type = args->inputtype;

  if(args->argc==2) {
    if((!strcasecmp(args->argv[1],"-beagle")) ||
       (!strcasecmp(args->argv[1],"-beagle-GL")) ||
       (!strcasecmp(args->argv[1],"-bgen")) ||
       (!strcasecmp(args->argv[1],"-glf")) ||
       (!strcasecmp(args->argv[1],"-glf3")) ||
       (!strcasecmp(args->argv[1],"-pileup")) ||
       (!strcasecmp(args->argv[1],"-vcf-GL")) ||
       (!strcasecmp(args->argv[1],"-vcf-pl")) ||
       (!strcasecmp(args->argv[1],"-glf10_text")) ||
       (!strcasecmp(args->argv[1],"-vcf-GP"))) {
      printArg(stderr,args);
      exit(0);
    }else if ((!strcasecmp(args->argv[1],"-bam")) ||
	      (!strcasecmp(args->argv[1],"-b"))){
      setArgsBam(args);
      exit(0);
    }else
      return;
  }
  
  getOptions(args);
  if(args->fai==NULL){
    int printAndExit =0;
    switch(args->inputtype)
      {
      case INPUT_GLF:
	printAndExit=1;
	break;
      case INPUT_GLF10_TEXT:
	printAndExit=1;
	break;
      case INPUT_GLF3:
	printAndExit=1;
	break;
      case INPUT_BEAGLE:
	printAndExit=1;
	break;
      case INPUT_BEAGLE_GL:
	printAndExit=1;
	break;
      case INPUT_BGEN:
	printAndExit=1;
	break;
      case INPUT_PILEUP:
	printAndExit=1;
	break;
      }
    if(printAndExit){
      fprintf(stderr,"\t-> Must supply -fai file\n");
      exit(1);
    }
  }
  if(args->fai){
    if(!(args->hd=getHeadFromFai(args->fai)))
      exit(1);
  }else{
    if(args->nams.size()==0){
      fprintf(stderr,"\t-> Must choose inputfile -bam/-glf/-glf3/-pileup/-i/-vcf-gl/-vcf-gp/-vcf-pl/-glf10_text filename\n");
      exit(1);
    }
    if(args->inputtype==INPUT_BAM){
      htsFile *in=sam_open(args->nams[0],"r");
      if(in==NULL){
	fprintf(stderr,"\t-> Problems opening file: \'%s\' will exit\n",args->nams[0]);
	exit(1);
      }
      args->hd= sam_hdr_read(in);
      hts_close(in);
      if(args->hd==NULL){
	fprintf(stderr,"\t-> Problems reading header from file: \'%s\' will exit\n",args->nams[0]);
	exit(1);
      }
    }
  }

  if((args->inputtype==INPUT_PILEUP||args->inputtype==INPUT_GLF||args->inputtype==INPUT_GLF3||args->inputtype==INPUT_GLF10_TEXT)){
    if(nInd==0){
      fprintf(stderr,"\t-> Must supply -nInd when using -glf/-glf3/-pileup/-glf10_text files\n");
      exit(1);
    }
    for(int i=0;i<nInd;i++){
      char bufbuf[256];
      snprintf(bufbuf,256,"sim%d",i);
      bam_smpl_add(args->sm, bufbuf, NULL);
    }
  }else
    args->nInd = args->nams.size();

  if(args->inputtype==INPUT_VCF_GP||args->inputtype==INPUT_VCF_GL){
    char *tmptmp=NULL;
    tmptmp=angsd::getArg("-vcf-pl",tmptmp,args);
    if(tmptmp!=NULL)
      pl_or_gl =0;
    free(tmptmp);tmptmp=NULL;
    
    tmptmp=angsd::getArg("-vcf-gl",tmptmp,args);
    if(tmptmp!=NULL)
      pl_or_gl =1;
    free(tmptmp);tmptmp=NULL;
    
    tmptmp=angsd::getArg("-vcf-gp",tmptmp,args);
    if(tmptmp!=NULL)
      pl_or_gl =2;
    free(tmptmp);tmptmp=NULL;



    
    myvcf = new vcfReader(args->infile,NULL,pl_or_gl,&args->regions);
    //the contigs are taken from the vcf header, also when -fai is given
    bam_hdr_destroy(args->hd);
    args->hd=myvcf->bamhdr;
    args->nInd = myvcf->hs->nsamples;
  }
  //make args->hd
  revMap = buildRevTable(args->hd);
  args->revMap = revMap;
  setArgsBam(args);
  if(args->inputtype==INPUT_BAM){
    //read bamfiles
    extern int checkBamHeaders;
    extern int doCheck;
    extern char *fai_fname;
    bufReader *initializeBufReaders2(const std::vector<char *> &vec,int exitOnError,int doCheck,char *fai_fname,bam_sample_t *sm);
    args->rd = initializeBufReaders2(args->nams,checkBamHeaders,doCheck,fai_fname,args->sm);
    fprintf(stderr, "[%s] %d samples in %lu input files\n", __func__, args->sm->n, args->nams.size());
  }

  
  //readers below that read through a gzFile
  int usesGz = args->inputtype==INPUT_PILEUP||args->inputtype==INPUT_GLF||args->inputtype==INPUT_GLF3||args->inputtype==INPUT_GLF10_TEXT||args->inputtype==INPUT_BEAGLE||args->inputtype==INPUT_BEAGLE_GL;
  if(fname!=NULL&&usesGz){
    gz=gzopen(fname,"r");
    if(gz==Z_NULL){
      fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
      exit(1);
    }
  }
  switch(args->inputtype){
  case INPUT_PILEUP:{
    mpil = new mpileup(args->nInd,gz,args->revMap,minQ);
    break;
  }
  case INPUT_GLF:{
    myglf = new glfReader(args->nInd,gz,10);
    break;
  }
  case INPUT_GLF3:{

    isSim = 1; //Added by FGV on 22/02/2015: GLF3 is always simulated data until a better alternative can be found
    myglf = new glfReader(args->nInd,gz,3);
    break;
  }
  case INPUT_GLF10_TEXT:{
    myglf_text = new glfReader_text(args->nInd,gz,args->revMap);
    break;
  }

  case INPUT_BEAGLE:{
    bglObj = new beagle_reader(gz,args->revMap,intName,args->nInd,0);
    break;
  }
  case INPUT_BEAGLE_GL:{
    bglObj = new beagle_reader(gz,args->revMap,intName,args->nInd,1);
    break;
  }
    
  case INPUT_BGEN:{
    bgenObj = new bgenReader(args->infile,args->revMap,intName,args->nInd);
    break;
  }
    
  default:{
    break;
  }

  }
  if(args->inputtype==INPUT_VCF_GL||args->inputtype==INPUT_VCF_GP){
    fprintf(stderr,"\t-> VCF still beta. Remember that\n");
    fprintf(stderr,"\t   1. indels are are discarded\n");
    fprintf(stderr,"\t   2. will use chrom, pos PL,GL,GP columns and tags\n");
    fprintf(stderr,"\t   3. GL/PL tags are interpreted as log10 and are scaled to ln\n");
    fprintf(stderr,"\t   4. GP tags are interpreted directly as unscaled post probs)\n");
    fprintf(stderr,"\t   5. FILTER column is currently NOT used (not sure what concensus is)\n");
    fprintf(stderr,"\t   6. -sites does NOT work with vcf input but -r does\n");
    fprintf(stderr,"\t   7. vcffilereading is still BETA, please report strange behaviour\n");
    fprintf(stderr,"\t   8. Consider adding \'-doMajorMinor 1\'\n");
  }
}

multiReader::~multiReader(){
  if(revMap){
    for(aMap::iterator it=revMap->begin();it!=revMap->end();++it)
      free((char*)it->first);
  }
 
  delete revMap;

  if(mpil)
    delete mpil;
  if(myvcf)
    delete myvcf;
  if(myglf)
    delete myglf;
  if(myglf_text)
    delete myglf_text;
  if(bglObj)
    delete bglObj;
  if(bgenObj)
    delete bgenObj;

  if(gz!=Z_NULL)
    gzclose(gz);
  free(fname);

 
  
}

funkyPars *multiReader::fetch(){
  funkyPars *fp = NULL;
  switch(args->inputtype){
  case INPUT_PILEUP:{
    fp = mpil->fetch(nLines);
    break;
  }
  case INPUT_GLF:
  case INPUT_GLF3:{
    fp = myglf->fetch(nLines); 
    break;
  }
  case INPUT_GLF10_TEXT:{
    fp = myglf_text->fetch(nLines); 
    break;
  }
    
  case INPUT_VCF_GL:{
    fp = myvcf->fetch(nLines); 
    break;
  }
  case INPUT_VCF_GP:{
    fp = myvcf->fetch(nLines); 
    break;
  }
  case INPUT_BEAGLE:
  case INPUT_BEAGLE_GL:{
    fp = bglObj->fetch(nLines); 
    break;
  }
  case INPUT_BGEN:{
    fp = bgenObj->fetch(nLines); 
    break;
  }
  case INPUT_BAM:{
    bammer_main(args);
    break;
    
  }
    
  default:{
    fprintf(stderr,"Unknown inputformat: %d\n",args->inputtype);
    exit(1);
  }
  }
  if(isSim==1&&fp!=NULL){
    fp->anc = new char[fp->numSites];
    memset(fp->anc,0,fp->numSites);
  }
  
  if(fp!=NULL && fp->refId==-1){
    fprintf(stderr,"\n\t-> Unknown reference found (not in header/fai file). Will exit\n");
    exit(1);
  }
  return fp;

}
