#include "abc.h"
#include "shared.h"
#include "multiReader.h"
#include "parseArgs_bambi.h"
#include <libgen.h>//for checking if output dir exists 'dirname'
#define ARGS ".arg"
void checkIfDir(char *fname){
  char *dirNam2 = strdup(fname);
  char *dirNam=dirname(dirNam2);
  
  if(strlen(dirNam)>1 &&!aio::fexists(dirNam)){
    fprintf(stderr,"\t Folder: \'%s\' doesn't exist, please create\n",dirNam);
    exit(0);
  }
  free(dirNam2);
  
}


aHead *getHeadFromFai(const char *fname){
  std::vector<char *> chrs;
  std::vector<int> lengths;
  FILE *fp = aio::getFILE(fname,"r");
  char buf[1024];
  while(fgets(buf,1024,fp)){ 
    chrs.push_back(strdup(strtok(buf,"\t \n")));//<-strdup so don't clean here
    lengths.push_back(atoi(strtok(NULL,"\t \n")));
  }
  
  aHead *ret = new aHead;
  ret->l_text = strlen(fname);
  ret->text = new char[strlen(fname)+1];
  ret->text = strcpy(ret->text,fname);
  ret->n_ref = chrs.size();
  ret->l_name = new int[chrs.size()];
  ret->l_ref = new int [chrs.size()];
  ret->name = new char*[chrs.size()];
  for(size_t i=0;i<chrs.size();i++){
    ret->l_name[i] = strlen(chrs[i]);
    ret->l_ref[i] = lengths[i];
    //    ret->name[i] = chrs[i];
    ret->name[i] =(char*) malloc(strlen(chrs[i])+1);
    strcpy(ret->name[i],chrs[i]);
  }

  for(uint i=0;i<chrs.size();i++)
    free(chrs[i]);
  fclose(fp);
  return ret;
}

aMap *buildRevTable(const aHead *hd){
  aMap *ret = new aMap;
  for(int i=0;i<hd->n_ref;i++){
    ret->insert(std::pair<char *,int>(strdup(hd->name[i]),i));
  }
  for(aMap::iterator it= ret->begin();0&&it!=ret->end();++it)
    fprintf(stderr,"%s %d\n",it->first,it->second);
  return ret;
}


void setInputType(argStruct *args){
#if 0
  fprintf(stderr,"bam:%d\n",INPUT_BAM);
  fprintf(stderr,"glf:%d\n",INPUT_GLF);
  fprintf(stderr,"blg:%d\n",INPUT_BEAGLE);
  fprintf(stderr,"plp:%d\n",INPUT_PILEUP);
#endif


  char *tmp =NULL;
  tmp = angsd::getArg("-glf",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF;
    args->infile = tmp;
    args->nams.push_back(args->infile);
    return;
  }

  tmp = angsd::getArg("-pileup",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_PILEUP;
    args->infile = tmp;
    args->nams.push_back(args->infile);
    return;
  }
  tmp = angsd::getArg("-beagle",tmp,args);
  if(tmp!=NULL){
    char *tmp_fai = NULL;
    tmp_fai = angsd::getArg("-fai",tmp_fai,args);
    if(tmp_fai==NULL){
      fprintf(stderr,"\t-> You must supply a fai file (-fai) when using -beagle input\n");
      exit(0);

    }
    args->inputtype=INPUT_BEAGLE;
    args->infile = tmp;
    args->nams.push_back(args->infile);
    return;
  }
  tmp = angsd::getArg("-i",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_BAM;
    args->infile = tmp;
    args->nams.push_back(args->infile);
    return;
  }
  int nInd = 0;
  nInd = angsd::getArg("-nInd",nInd,args);
  tmp = angsd::getArg("-bam",tmp,args);
  tmp = angsd::getArg("-b",tmp,args);
  if(tmp!=NULL&&args->argc>2){
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
  return arguments;
}


void multiReader::printArg(FILE *argFile){
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
  fprintf(argFile,"-nLines=%d\n",nLines);
  fprintf(argFile,"-bytesPerLine=%d\n",bytesPerLine);
  fprintf(argFile,"\t-beagle\t%s\t(Beagle Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-glf\t%s\t(glf Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-pileup\t%s\t(pileup Filename (can be .gz))\n",fname);
  fprintf(argFile,"\t-intName=%d\t(Assume First column is chr_position)\n",intName);
  fprintf(argFile,"\t-isSim=%d\t(Assume First column is chr_position)\n",isSim);
  fprintf(argFile,"\t-from=%d\t(Assume First column is chr_position)\n",from);
  fprintf(argFile,"\t-to=%d\t(Assume First column is chr_position)\n",to);
  fprintf(argFile,"\t-nInd=%d\t(Assume First column is chr_position)\n",nInd);
  fprintf(argFile,"\t-minQ=%d\t(minium q only used in pileupreader)\n",minQ);
  fprintf(argFile,"----------------\n%s:\n",__FILE__); 
}
void multiReader::getOptions(argStruct *arguments){

  nLines=angsd::getArg("-nLines",nLines,arguments);
  arguments->nLines = nLines;
  bytesPerLine=angsd::getArg("-bytesPerLine",bytesPerLine,arguments);
  fname=angsd::getArg("-beagle",fname,arguments);
  fname=angsd::getArg("-pileup",fname,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  fname=angsd::getArg("-glf",fname,arguments);
  intName=angsd::getArg("-intName",intName,arguments);
  isSim=angsd::getArg("-isSim",intName,arguments);
  nInd=angsd::getArg("-nInd",nInd,arguments);
  from=angsd::getArg("-from",from,arguments);
  to=angsd::getArg("-to",to,arguments);
  arguments->nInd = nInd;
  //read fai if suppplied (used by other than native bam readers)
  
  fai=angsd::getArg("-fai",fai,arguments);
  
  printArg(arguments->argumentFile);
}


multiReader::multiReader(int argc,char**argv){

  fai = NULL;
  bytesPerLine = 33554432 ;//2^15 about 33 megs perline/persites should be enough
  nLines=50;
  
  fname=NULL;
  intName=1;
  minQ = MINQ;

  nInd =0;
  nInd2=-1;
  from =to=-1;
  isSim =0;
  args=NULL;
  args = setArgStruct(argc,argv);
  fprintf(args->argumentFile,"\t-> Command: \n");
  for(int i=0;i<argc;i++)
    fprintf(args->argumentFile,"%s ",argv[i]);
  //  fprintf(args->argumentFile,"\n\n");
  extern double VERS;
  fprintf(args->argumentFile,"\n\t-> angsd version: %.3f\t build(%s %s)\n",VERS,__DATE__,__TIME__); 
  void printTime(FILE *fp);
  printTime(args->argumentFile);  


  type = args->inputtype;

  if(args->argc==2){
    if((!strcasecmp(args->argv[1],"-beagle"))||!strcasecmp(args->argv[1],"-glf")||(!strcasecmp(args->argv[1],"-pileup"))){
      printArg(stdout);
      exit(0);
    }else if ((!strcasecmp(args->argv[1],"-bam"))|| (!strcasecmp(args->argv[1],"-b"))){
  
      setArgsBam(args);
      exit(0);
    }else
      return;
  }
  getOptions(args);

  if(fai){
    hd=getHeadFromFai(fai);
  }else{
    if(args->nams.size()==0){
      fprintf(stderr,"\t-> Must choose inputfile -bam/-glf/-pileup/-i filename\n");
      exit(0);
    }
    hd= getHd_andClose(args->nams[0]);
  }
  args->hd = hd;
  
  if(args->hd==NULL){
    fprintf(stderr,"For non-bams you should include -fai arguments\n");
    exit(0);
  }  



  if((type==INPUT_PILEUP||type==INPUT_GLF)){
    if(nInd==0){
      fprintf(stderr,"\t-> Must supply -nInd when using raw GLF files\n");
      exit(0);
    }
  }else
    args->nInd = args->nams.size();

  revMap = buildRevTable(args->hd);

  args->revMap = revMap;
  setArgsBam(args);
  gz=gzopen(fname,"r");   

  switch(type){
  case INPUT_PILEUP:{
    mpil = new mpileup(args->nInd,gz,bytesPerLine,args->revMap,minQ);
    break;
  }
  case INPUT_GLF:{
    myglf = new glfReader(args->nInd,nInd2,from,to,gz,isSim);
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

}

multiReader::~multiReader(){
  if(revMap){
    for(aMap::iterator it=revMap->begin();it!=revMap->end();++it)
      free((char*)it->first);
  }
  dalloc(hd);
  delete revMap;

  

  free(fai);

  switch(type){
  case INPUT_PILEUP:{
    delete mpil;
    break;
  }
  case INPUT_GLF:{
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
  case INPUT_GLF:{
    fp = myglf->fetch(nLines); 
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
    fprintf(stderr,"whats up\n");
    exit(0);
    break;
  }
  }
  if(fp&&0)
    fprintf(stderr,"numSites:%d\n",fp->numSites);
    
  return fp;

}
