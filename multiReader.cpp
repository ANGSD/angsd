#include <libgen.h>//for checking if output dir exists 'dirname'
#include <cassert>
#include "abc.h"
#include "shared.h"
#include "multiReader.h"
#include "parseArgs_bambi.h"
#include "version.h"

#define ARGS ".arg"
void checkIfDir(char *fname){
  char *dirNam2 = strdup(fname);
  char *dirNam=dirname(dirNam2);
  
  if(strlen(dirNam)>1 &&!aio::fexists(dirNam)){
    fprintf(stderr,"\t Folder: \'%s\' doesn't exist, please create\n",dirNam);
    exit(0);
  }
  //free(dirNam); VALGRIND on osx wants this line
  free(dirNam2);
  
}


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
  aMap *ret = new aMap;
  for(int i=0;i<hd->n_targets;i++){
    ret->insert(std::pair<char *,int>(strdup(hd->target_name[i]),i));
  }
  for(aMap::iterator it= ret->begin();0&&it!=ret->end();++it)
    fprintf(stderr,"%s %d\n",it->first,it->second);
  return ret;
}


void setInputType(argStruct *args){
#if 0
  fprintf(stderr,"bam:%d\n",INPUT_BAM);
  fprintf(stderr,"glf:%d\n",INPUT_GLF);
  fprintf(stderr,"glf3:%d\n",INPUT_GLF3);
  fprintf(stderr,"blg:%d\n",INPUT_BEAGLE);
  fprintf(stderr,"plp:%d\n",INPUT_PILEUP);
  fprintf(stderr,"vcf_gp:%d\n",INPUT_VCF_GP);
  fprintf(stderr,"vcf_gl:%d\n",INPUT_VCF_GL);
#endif


  char *tmp =NULL;
  tmp = angsd::getArg("-glf",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp = NULL;
  tmp = angsd::getArg("-glf3",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF3;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-vcf-GP",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GP;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    char *tmp_fai = NULL;
    tmp_fai = angsd::getArg("-fai",tmp_fai,args);
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-vcf-GL",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GL;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-vcf-PL",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GL;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-pileup",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_PILEUP;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-beagle",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_BEAGLE;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  tmp=NULL;
  tmp = angsd::getArg("-i",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_BAM;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  int nInd = 0;
  nInd = angsd::getArg("-nInd",nInd,args);
  tmp=NULL;
  tmp = angsd::getArg("-bam",tmp,args);
  tmp = angsd::getArg("-b",tmp,args);
  if(tmp!=NULL && args->argc>2){
    args->inputtype=INPUT_BAM;
    args->infile = tmp;
    args->nams = angsd::getFilenames(args->infile,nInd);
    
    return;
  }

}



argStruct *setArgStruct(int argc,char **argv) { 

  argStruct *arguments = new argStruct;
  arguments->argumentFile = stderr;
  arguments->outfiles =NULL;
  arguments->hd=NULL;
  arguments->revMap= NULL;
  arguments->argc=argc;
  arguments->argv=argv;
  arguments->usedArgs= new int[argc+1];//well here we allocate one more than needed, this is only used in the ./angsd -beagle version
  for(int i=0;i<argc;i++)
    arguments->usedArgs[i]=0;
  
  arguments->inputtype=-1;
  arguments->infile = NULL;
  arguments->show =0;

  //check output filename
  arguments->outfiles = angsd::getArg("-out",arguments->outfiles,arguments);
  if(arguments->outfiles==NULL){
    if(argc!=2)
      fprintf(stderr,"\t-> No \'-out\' argument given, output files will be called \'angsdput\'\n");
    arguments->outfiles = strdup("angsdput");
  }
  checkIfDir(arguments->outfiles);
  

  //print arguments into logfile
  if(argc>2)
    arguments->argumentFile=aio::openFile(arguments->outfiles,ARGS);
   
  setInputType(arguments);
  //  fprintf(stderr,"setintputtype:%d\n",arguments->inputtype);exit(0);
  return arguments;
}


void multiReader::printArg(FILE *argFile){
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
  fprintf(argFile,"\t-nLines\t%d\t(Number of lines to read)\n",nLines);
  fprintf(argFile,"\t-bpl\t%d (bytesPerLine)\n",bytesPerLine);
  fprintf(argFile,"\t-beagle\t%s\t(Beagle Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-vcf-GL\t%s\t(vcf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-vcf-GP\t%s\t(vcf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-glf\t%s\t(glf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-pileup\t%s\t(pileup Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-intName %d\t(Assume First column is chr_position)\n",intName);
  fprintf(argFile,"\t-isSim\t%d\t(Simulated data assumes ancestral is A)\n",isSim);
  fprintf(argFile,"\t-nInd\t%d\t\t(Number of individuals)\n",nInd);
  fprintf(argFile,"\t-minQ\t%d\t(minimum base quality; only used in pileupreader)\n",minQ);
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
}


void multiReader::getOptions(argStruct *arguments){

  nLines=angsd::getArg("-nLines",nLines,arguments);
  arguments->nLines = nLines;
  bytesPerLine=angsd::getArg("-bpl",bytesPerLine,arguments);
  fname=angsd::getArg("-beagle",fname,arguments);
  fname=angsd::getArg("-pileup",fname,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  fname=angsd::getArg("-glf",fname,arguments);
  fname=angsd::getArg("-glf3",fname,arguments);
  fname=angsd::getArg("-vcf-GL",fname,arguments);
  fname=angsd::getArg("-vcf-GP",fname,arguments);
  fname=angsd::getArg("-vcf-pl",fname,arguments);
  intName=angsd::getArg("-intName",intName,arguments);
  isSim=angsd::getArg("-isSim",intName,arguments);
  nInd=angsd::getArg("-nInd",nInd,arguments);
  arguments->nInd = nInd;

  //read fai if suppplied (used by other than native bam readers)
  
  fai=angsd::getArg("-fai",fai,arguments);
  
  printArg(arguments->argumentFile);
}


multiReader::multiReader(int argc,char**argv){
  gz=Z_NULL;
  myglf=NULL;myvcf=NULL;mpil=NULL;bglObj=NULL;
  fai = NULL;
  bytesPerLine = 33554432 ;//2^15 about 33 megs perline/persites should be enough
  nLines=50;
  
  fname=NULL;
  intName=1;
  minQ = MINQ;

  nInd =0;
  isSim =0;
  args=NULL;
  args = setArgStruct(argc,argv);
  fprintf(args->argumentFile,"\t-> Command: \n");
  for(int i=0;i<argc;i++)
    fprintf(args->argumentFile,"%s ",argv[i]);
  //  fprintf(args->argumentFile,"\n\n");
  if(args->argumentFile!=stderr)
    fprintf(args->argumentFile,"\n\t-> angsd version: %s (htslib: %s) build(%s %s)\n",ANGSD_VERSION,hts_version(),__DATE__,__TIME__); 
  void printTime(FILE *fp);
  printTime(args->argumentFile); 


  type = args->inputtype;

  if(args->argc==2) {
    if((!strcasecmp(args->argv[1],"-beagle")) ||
       (!strcasecmp(args->argv[1],"-glf")) ||
       (!strcasecmp(args->argv[1],"-glf3")) ||
       (!strcasecmp(args->argv[1],"-pileup")) ||
       (!strcasecmp(args->argv[1],"-vcf-GL")) ||
       (!strcasecmp(args->argv[1],"-vcf-pl")) ||
       (!strcasecmp(args->argv[1],"-vcf-GP"))) {
      printArg(stdout);
      exit(0);
    }else if ((!strcasecmp(args->argv[1],"-bam")) ||
	      (!strcasecmp(args->argv[1],"-b"))){
      setArgsBam(args);
      exit(0);
    }else
      return;
  }
  
  getOptions(args);
  if(fai==NULL){
    int printAndExit =0;
    switch(args->inputtype)
      {
      case INPUT_GLF:
	printAndExit=1;
	break;
      case INPUT_GLF3:
	printAndExit=1;
	break;
      case INPUT_VCF_GP:
	printAndExit=1;
	break;
      case INPUT_VCF_GL:
	printAndExit=1;
	break;
      case INPUT_BEAGLE:
	printAndExit=1;
	break;
      }
    if(printAndExit){
      fprintf(stderr,"\t-> Must supply -fai file\n");
      exit(0);

    }
  }

  
  if(fai){
    if(!(hd=getHeadFromFai(fai)))
      exit(0);
  }else{
    if(args->nams.size()==0){
      fprintf(stderr,"\t-> Must choose inputfile -bam/-glf/-glf3/-pileup/-i/-vcf-gl/-vcf-gp/-vcf-pl filename\n");
      exit(0);
    }
    htsFile *in=sam_open(args->nams[0],"r");
    assert(in);
    hd= sam_hdr_read(in);
    hts_close(in);
  }
  args->hd = hd;
  
  if(args->hd==NULL){
    fprintf(stderr,"For non-bams you should include -fai arguments\n");
    exit(0);
  }  



  if((type==INPUT_PILEUP||type==INPUT_GLF||type==INPUT_GLF3||type==INPUT_VCF_GP||type==INPUT_VCF_GL)){
    if(nInd==0){
      fprintf(stderr,"\t-> Must supply -nInd when using -glf/-glf3/-pileup/-vcf-GL/-vcf-GP files\n");
      exit(0);
    }
  }else
    args->nInd = args->nams.size();

  revMap = buildRevTable(args->hd);

  args->revMap = revMap;

  setArgsBam(args);
  if(fname==NULL)
    return;

  gz=Z_NULL;
  gz=gzopen(fname,"r");   
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
    exit(0);
  }
  switch(type){
  case INPUT_PILEUP:{
    mpil = new mpileup(args->nInd,gz,bytesPerLine,args->revMap,minQ);
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
  case INPUT_VCF_GP:{
    myvcf = new vcfReader(args->nInd,gz,bytesPerLine,args->revMap);
    break;
  }
  case INPUT_VCF_GL:{
    fprintf(stderr,"input_vcf\n");
    myvcf = new vcfReader(args->nInd,gz,bytesPerLine,args->revMap);
    break;
  }
  case INPUT_BEAGLE:{
    bglObj = new beagle_reader(bytesPerLine,gz,args->revMap,intName,args->nInd);
    break;
  }
    
  default:{
    //    fprintf(stderr,"\t[%s] assuming what ?\n",__FUNCTION__);
    //pl = new pileups (faifile,lStart,lStop,filenames,type);
    break;
  }
  }
  if(args->inputtype==INPUT_VCF_GL||args->inputtype==INPUT_VCF_GL){
    fprintf(stderr,"\t-> VCF still beta. Remember that\n");
    fprintf(stderr,"\t   1. indels and sites that are non diallelic are discarded\n");
    fprintf(stderr,"\t   2. will use chrom, pos, ref, alt columns\n");
    fprintf(stderr,"\t   3. GL tags are interpreted as log10 and are scaled to ln\n");
    fprintf(stderr,"\t   4. GP tags are interpreted directly as unscaled post probs (spec says phredscaled...)\n");
    fprintf(stderr,"\t   5. FILTER column is currently NOT used (not sure what concensus is)\n");

  }
}

multiReader::~multiReader(){
  if(revMap){
    for(aMap::iterator it=revMap->begin();it!=revMap->end();++it)
      free((char*)it->first);
  }
 
  delete revMap;

  

  free(fai);

  switch(type){
  case INPUT_PILEUP:{
    delete mpil;
    break;
  }
  case INPUT_VCF_GP:{
    delete myvcf;
    break;
  }
  case INPUT_VCF_GL:{
    delete myvcf;
    break;
  }
  case INPUT_GLF:
  case INPUT_GLF3:{
    delete myglf;
    break;
  }    
  case INPUT_BEAGLE:{
    delete bglObj;
    break;
  }    
  default:{
    break;
  }
  }
  if(gz!=Z_NULL)
    gzclose(gz);
  free(fname);
  
  for(unsigned i=0;i<args->nams.size();i++)
    free(args->nams[i]);
  
  delete []   args->usedArgs;
  free(args->outfiles);
  free(args->infile);
  bam_hdr_destroy(args->hd);
  if(args->argumentFile!=stderr) fclose(args->argumentFile);
  delete args;
  
}
funkyPars *multiReader::fetch(){
  //  fprintf(stderr,"fetching\n");`
  funkyPars *fp = NULL;
  switch(type){
  case INPUT_PILEUP:{
    fp = mpil->fetch(nLines);
    break;
  }
  case INPUT_GLF:
  case INPUT_GLF3:{
    fp = myglf->fetch(nLines); 
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
    fprintf(stderr,"Unkown inputformat: %d\n",type);
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
    fprintf(stderr,"\t-> Unkown refid found will close program\n");
    return NULL;
  }
  return fp;

}
