#include <cassert>
#include "abc.h"
#include "shared.h"
#include "multiReader.h"
#include "parseArgs_bambi.h"
#include "version.h"
#include "sample.h"
#include "aio.h"

bam_hdr_t *getHeadFromFai(const char *fname){
  std::vector<char *> chrs;
  std::vector<int> lengths;
  FILE *fp = aio::getFILE(fname,"r");
  char buf[1024];
  while(fgets(buf,1024,fp)){ 
    char *tok = strtok(buf,"\t \n");
    //fprintf(stderr,"tok: %s\n",tok);
    chrs.push_back(strdup(tok));//<-strdup so don't clean here
    tok = strtok(NULL,"\t \n");
    if(!tok){
      fprintf(stderr,"\t-> fai file looks malformed? last reference: %s\n",chrs.back());
      fclose(fp);
      return NULL;
    }
    // fprintf(stderr,"tok: %s\n",tok);
    lengths.push_back(atoi(tok));
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
    ret->target_name[i] =strdup(chrs[i]);
  }
 
  for(uint i=0;i<chrs.size();i++)
    free(chrs[i]);
  fclose(fp);
  return ret;
}

aMap *buildRevTable(const bam_hdr_t *hd){
  assert(hd);
  aMap *ret = new aMap;
  for(int i=0;i<hd->n_targets;i++){
    ret->insert(std::pair<char *,int>(strdup(hd->target_name[i]),i));
  }
  for(aMap::iterator it= ret->begin();0&&it!=ret->end();++it)
    fprintf(stderr,"%s %d\n",it->first,it->second);
  return ret;
}



void multiReader::printArg(FILE *argFile,argStruct *args){
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
  fprintf(argFile,"\t-nLines\t%d\t(Number of lines to read)\n",nLines);
  fprintf(argFile,"\t-beagle\t%s\t(Beagle Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-vcf-GL\t%s\t(vcf Filename (can be bcf compressed or uncompressed))\n",fname);
  fprintf(argFile,"\t-vcf-PL\t%s\t(vcf Filename (can be bcf compressed or uncompressed))\n",fname);
  fprintf(argFile,"\t-vcf-GP\t%s\t(vcf Filename (can be bcf compressed or uncompressed))(*not used)\n",fname);
  fprintf(argFile,"\t-glf\t%s\t(glf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-pileup\t%s\t(pileup Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-intName %d\t(Assume First column is chr_position)\n",intName);
  fprintf(argFile,"\t-isSim\t%d\t(Simulated data assumes ancestral is A)\n",isSim);
  fprintf(argFile,"\t-nInd\t%d\t\t(Number of individuals)\n",nInd);
  fprintf(argFile,"\t-minQ\t%d\t(minimum base quality; only used in pileupreader)\n",minQ);
  fprintf(argFile,"\t-fai\t%s\t(fai file)\n",args->fai);
  fprintf(argFile,"\t-minQ\t%d\t(minimum base quality; only used in pileupreader)\n",minQ);

  
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
}


void multiReader::getOptions(argStruct *arguments){

  nLines=angsd::getArg("-nLines",nLines,arguments);
  arguments->nLines = nLines;
  fname=angsd::getArg("-beagle",fname,arguments);
  fname=angsd::getArg("-pileup",fname,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  fname=angsd::getArg("-glf",fname,arguments);
  fname=angsd::getArg("-glf3",fname,arguments);
  fname=angsd::getArg("-vcf-GL",fname,arguments);
  fname=angsd::getArg("-vcf-GP",fname,arguments);
  fname=angsd::getArg("-vcf-pl",fname,arguments);
  fname=angsd::getArg("-glf10_text",fname,arguments);
  intName=angsd::getArg("-intName",intName,arguments);
  isSim=angsd::getArg("-isSim",intName,arguments);
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
  myglf=NULL;myvcf=NULL;mpil=NULL;bglObj=NULL;
  
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
       (!strcasecmp(args->argv[1],"-glf")) ||
       (!strcasecmp(args->argv[1],"-glf3")) ||
       (!strcasecmp(args->argv[1],"-pileup")) ||
       (!strcasecmp(args->argv[1],"-vcf-GL")) ||
       (!strcasecmp(args->argv[1],"-vcf-pl")) ||
       (!strcasecmp(args->argv[1],"-glf10_text")) ||
       (!strcasecmp(args->argv[1],"-vcf-GP"))) {
      printArg(stdout,args);
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
      case INPUT_PILEUP:
	printAndExit=1;
	break;
      }
    if(printAndExit){
      fprintf(stderr,"\t-> Must supply -fai file\n");
      exit(0);
    }
  }
  
  if(args->fai){
    if(!(args->hd=getHeadFromFai(args->fai)))
      exit(0);
  }else{
    if(args->nams.size()==0){
      fprintf(stderr,"\t-> Must choose inputfile -bam/-glf/-glf3/-pileup/-i/-vcf-gl/-vcf-gp/-vcf-pl/-glf10_text filename\n");
      exit(0);
    }
    if(args->inputtype==INPUT_BAM){
      htsFile *in=sam_open(args->nams[0],"r");
      assert(in);
      args->hd= sam_hdr_read(in);
      hts_close(in);
    }
  }
 
  if(!(INPUT_VCF_GL||INPUT_VCF_GP)){
    if(args->hd==NULL){
      fprintf(stderr,"For non-bams you should include -fai arguments\n");
      exit(0);
    }
  }

  if((args->inputtype==INPUT_PILEUP||args->inputtype==INPUT_GLF||args->inputtype==INPUT_GLF3||args->inputtype==INPUT_GLF10_TEXT)){
    if(nInd==0){
      fprintf(stderr,"\t-> Must supply -nInd when using -glf/-glf3/-pileup/-glf10_text files\n");
      exit(0);
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
    myvcf = new vcfReader(args->infile,NULL,pl_or_gl,&args->regions);
    args->hd=myvcf->bamhdr;
    args->nInd = myvcf->hs->nsamples;
  }
  //make args->hd
  
  revMap = buildRevTable(args->hd);
  args->revMap = revMap;
  setArgsBam(args);
  if(args->inputtype==INPUT_VCF_GL){
    if(args->regions.size()==1){
      char tmp[1024];
      int start=args->regions[0].start;
      int stop=args->regions[0].stop;
      int ref=args->regions[0].refID;
      snprintf(tmp,1024,"%s:%d-%d",args->hd->target_name[ref],start+1,stop);
      //    fprintf(stderr,"tmp:%s\n",tmp);
      //    exit(0);
      myvcf->seek(tmp);
    }
    
  }
  
  if(args->inputtype==INPUT_BAM){
    //read bamfiles
    extern int checkBamHeaders;
    extern int doCheck;
    extern char *fai_fname;
    bufReader *initializeBufReaders2(const std::vector<char *> &vec,int exitOnError,int doCheck,char *fai_fname,bam_sample_t *sm);
    args->rd = initializeBufReaders2(args->nams,checkBamHeaders,doCheck,fai_fname,args->sm);
    fprintf(stderr, "[%s] %d samples in %lu input files\n", __func__, args->sm->n, args->nams.size());
  }

  
  if(fname==NULL)
    return;

  gz=Z_NULL;
  gz=gzopen(fname,"r");   
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
    exit(0);
  }

  switch(args->inputtype){
  case INPUT_PILEUP:{
    mpil = new mpileup(args->nInd,gz,args->revMap,minQ);
    break;
  }
  case INPUT_GLF:{
    myglf = new glfReader(args->nInd,gz,10,isSim);
    break;
  }
  case INPUT_GLF3:{
    isSim = 1; //Added by FGV on 22/02/2015: GLF3 is always simulated data until a better alternative can be found
    myglf = new glfReader(args->nInd,gz,3,isSim);
    break;
  }
  case INPUT_GLF10_TEXT:{
    myglf_text = new glfReader_text(args->nInd,gz,args->revMap);
    break;
  }

  case INPUT_BEAGLE:{
    bglObj = new beagle_reader(gz,args->revMap,intName,args->nInd);
    break;
  }
  default:{
    break;
  }
  }
  if(args->inputtype==INPUT_VCF_GL||args->inputtype==INPUT_VCF_GL){
    fprintf(stderr,"\t-> VCF still beta. Remember that\n");
    fprintf(stderr,"\t   1. indels are are discarded\n");
    fprintf(stderr,"\t   2. will use chrom, pos PL columns\n");
    fprintf(stderr,"\t   3. GL tags are interpreted as log10 and are scaled to ln (NOT USED)\n");
    fprintf(stderr,"\t   4. GP tags are interpreted directly as unscaled post probs (spec says phredscaled...) (NOT USED)\n");
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
  if(bglObj)
    delete bglObj;

  if(gz!=Z_NULL)
    gzclose(gz);
  free(fname);

 
  
}

funkyPars *multiReader::fetch(){
  //  fprintf(stderr,"fetching\n");`
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
  case INPUT_BEAGLE:{
    fp = bglObj->fetch(nLines); 
    break;
  }
  case INPUT_BAM:{
    bammer_main(args);
    break;
    
  }
    
  default:{
    fprintf(stderr,"Unknown inputformat: %d\n",args->inputtype);
#if 1
  fprintf(stderr,"bam:%d\n",INPUT_BAM);
  fprintf(stderr,"glf:%d\n",INPUT_GLF);
  fprintf(stderr,"glf3:%d\n",INPUT_GLF3);
  fprintf(stderr,"blg:%d\n",INPUT_BEAGLE);
  fprintf(stderr,"plp:%d\n",INPUT_PILEUP);
  fprintf(stderr,"vcf_GP:%d\n",INPUT_VCF_GP);
  fprintf(stderr,"vcf_GL:%d\n",INPUT_VCF_GL);
#endif

    exit(0);
    break;
  }
  }
  if(fp&&0)
    fprintf(stderr,"numSites:%d\n",fp->numSites);
  if(fp!=NULL && fp->refId==-1){
    fprintf(stderr,"\n\t-> Unkown refid found will close program refid:%d\n",fp->refId);
    return NULL;
  }
  return fp;

}
