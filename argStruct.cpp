#include <libgen.h>//for checking if output dir exists 'dirname'
#include <fstream>
#include <sstream>
#include <string>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include "version.h"
#include "argStruct.h"
#include "aio.h"

#define ARGS ".arg"

void whatIsTheInput(int type){
  switch(type){
  case INPUT_BAM:
    fprintf(stderr,"\t-> Inputtype is BAM/CRAM\n");
    break;
  case INPUT_GLF:
    fprintf(stderr,"\t-> Inputtype is GLF\n");
    break;
  case INPUT_GLF3:
    fprintf(stderr,"\t-> Inputtype is GLF3\n");
    break;
  case INPUT_BEAGLE:
    fprintf(stderr,"\t-> Inputtype is beagle\n");
    break;
  case INPUT_BGEN:
    fprintf(stderr,"\t-> Inputtype is bgen\n");
    break;
  case INPUT_PILEUP:
    fprintf(stderr,"\t-> Inputtype is pileup\n");
    break;
  case INPUT_VCF_GP:
    fprintf(stderr,"\t-> Inputtype is vcf/bcf gp\n");
    break;
  case INPUT_VCF_GL:
    fprintf(stderr,"\t-> Inputtype is vcf/bcf gl\n");
    break;
  case INPUT_GLF10_TEXT:
    fprintf(stderr,"\t-> Inputtype is GLF10 text\n");
    break;
  default:{
    fprintf(stderr,"\t-> Unknown input type, must supply one of -b/-bam/-i/-glf/-glf3/-glf10_text/-vcf-GL/-vcf-GP/-vcf-PL/-pileup/-beagle/-bgen\n");
    exit(1);
   }
  }
}

//returns 1 if the string can be fully parsed as a number (e.g. "-1", "-0.5", "1e-6")
static int isNumber(const char *s){
  if(s==NULL||*s=='\0')
    return 0;
  char *end=NULL;
  strtod(s,&end);
  return *end=='\0';
}

//a flag is something like "-out", but not "-" (stdin/stdout) and not a negative number like "-1"
int angsd::isFlag(const char *s){
  return s!=NULL && s[0]=='-' && s[1]!='\0' && !isNumber(s);
}

/*
  Finds argName in argv and returns the position of its value.
  Returns -1 if the flag is not present.
  Returns  0 if argc==2 (the "angsd -doXXX" help mode), callers then return their -999 sentinel.
  Every flag in angsd takes exactly one value, so the value following a flag is skipped
  and can never be mistaken for a flag itself.
*/
static int findArg(const char *argName,argStruct *arguments){
  int argPos = 1;
  while(argPos<arguments->argc){
    const char *cur = arguments->argv[argPos];
    if(strcasecmp(cur,argName)==0){
      if(arguments->argc==2)
        return 0;
      arguments->usedArgs[argPos]=1;
      if(argPos==arguments->argc-1){
        fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
        exit(1);
      }
      const char *val = arguments->argv[argPos+1];
      if(angsd::isFlag(val)){
        fprintf(stderr,"\t-> Must supply a parameter for: %s (found the flag \'%s\' instead)\n",argName,val);
        exit(1);
      }
      arguments->usedArgs[argPos+1]=1;
      return argPos+1;
    }
    argPos += angsd::isFlag(cur) ? 2 : 1;
  }
  return -1;
}

static long parseLong(const char *argName,const char *val){
  char *end=NULL;
  errno=0;
  long v=strtol(val,&end,10);
  if(end==val||*end!='\0'||errno==ERANGE||v<INT_MIN||v>INT_MAX){
    fprintf(stderr,"\t-> Parameter for %s must be an integer, got: \'%s\'\n",argName,val);
    exit(1);
  }
  return v;
}

static double parseDouble(const char *argName,const char *val){
  char *end=NULL;
  errno=0;
  double v=strtod(val,&end);
  if(end==val||*end!='\0'||errno==ERANGE){
    fprintf(stderr,"\t-> Parameter for %s must be a number, got: \'%s\'\n",argName,val);
    exit(1);
  }
  return v;
}

//the possible input types where a single filename is supplied
static const struct{
  const char *flag;
  int type;
}singleInputs[]={
  {"-glf",INPUT_GLF},
  {"-glf3",INPUT_GLF3},
  {"-glf10_text",INPUT_GLF10_TEXT},
  {"-vcf-GP",INPUT_VCF_GP},
  {"-vcf-GL",INPUT_VCF_GL},
  {"-vcf-PL",INPUT_VCF_GL},//PL is read by the same vcf reader as GL, multiReader converts PL to GL
  {"-pileup",INPUT_PILEUP},
  {"-beagle",INPUT_BEAGLE},
  {"-bgen",INPUT_BGEN},
  {"-i",INPUT_BAM},
};

void setInputType(argStruct *args){
  if(args->fai){
    free(args->fai);
    args->fai=NULL;
  }
  args->fai = angsd::getArg("-fai",args->fai,args);

  const char *usedFlag=NULL;
  for(size_t i=0;i<sizeof(singleInputs)/sizeof(singleInputs[0]);i++){
    char *tmp = angsd::getArg(singleInputs[i].flag,(char*)NULL,args);
    if(tmp==NULL)
      continue;
    if(usedFlag){
      fprintf(stderr,"\t-> Only one input can be supplied, both %s and %s were given\n",usedFlag,singleInputs[i].flag);
      exit(1);
    }
    usedFlag = singleInputs[i].flag;
    args->inputtype = singleInputs[i].type;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
  }

  //a filelist of bam/cram files
  char *bam = angsd::getArg("-bam",(char*)NULL,args);
  char *b = angsd::getArg("-b",(char*)NULL,args);
  if(bam&&b){
    fprintf(stderr,"\t-> Both -bam and -b were given, these are the same option, please only use one\n");
    exit(1);
  }
  char *list = bam?bam:b;
  if(list==NULL)
    return;
  if(usedFlag){
    fprintf(stderr,"\t-> Only one input can be supplied, both %s and %s were given\n",usedFlag,bam?"-bam":"-b");
    exit(1);
  }
  int nInd = 0;
  nInd = angsd::getArg("-nInd",nInd,args);
  args->inputtype=INPUT_BAM;
  args->infile = list;
  args->nams = angsd::getFilenames(args->infile,nInd);
}

void checkIfDir(char *fname){
  char *dirNam2 = strdup(fname);
  char *dirNam=dirname(dirNam2);
  if(strlen(dirNam)>1 &&!aio::fexists(dirNam)){
    fprintf(stderr,"\t Folder: \'%s\' doesn't exist, please create\n",dirNam);
    exit(1);
  }
  free(dirNam2);
}

argStruct *setArgStruct(int argc,char **argv) {
  //value-initialization, all pointers and counters not set below start out as NULL/0
  argStruct *arguments = new argStruct();
  arguments->argumentFile = stderr;
  arguments->argc=argc;
  arguments->argv=argv;
  arguments->nReads = 50;
  arguments->sm=bam_smpl_init();
  ASSERT(arguments->sm);

  //one more than needed, findArg marks the value after a flag (argPos+1) as used
  arguments->usedArgs= new int[argc+1]();

  arguments->inputtype=-1;

  kstring_t kstr;
  kstr.s=NULL;kstr.l=kstr.m=0;
  ksprintf(&kstr,"angsd version: %s (htslib: %s) build(%s %s)\n",ANGSD_VERSION,hts_version(),__DATE__,__TIME__);
  arguments->version = strdup(kstr.s);
  kstr.l=0;
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  ksprintf (&kstr, "%s", asctime (timeinfo) );
  //asctime ends with a newline, which we don't want in eg. the vcf/bcf header
  if(kstr.l>0 && kstr.s[kstr.l-1]=='\n')
    kstr.s[--kstr.l]='\0';
  arguments->datetime = strdup(kstr.s);
  kstr.l = 0;
  for(int i=0;i<argc;i++)
    ksprintf(&kstr,"%s ",argv[i]);
  arguments->cmdline = kstr.s;
  fprintf(stderr,"\t-> %s",arguments->version);
  fprintf(stderr,"\t-> %s\n",arguments->cmdline);
  //check output filename
  arguments->outfiles = angsd::getArg("-out",arguments->outfiles,arguments);
  if(arguments->outfiles==NULL){
    if(argc>2)
      fprintf(stderr,"\t-> No \'-out\' argument given, output files will be called \'angsdput\'\n");
    arguments->outfiles = strdup("angsdput");
  }
  checkIfDir(arguments->outfiles);

  if(argc<=2)
    return arguments;
  //print arguments into logfile
  arguments->argumentFile=aio::openFile(arguments->outfiles,ARGS);
  fprintf(arguments->argumentFile,"\t-> %s",arguments->version);
  fprintf(arguments->argumentFile,"\t-> %s\n",arguments->cmdline);
  setInputType(arguments);
  whatIsTheInput(arguments->inputtype);
  return arguments;
}



int angsd::getArg(const char* argName,int type,argStruct *arguments){
  int pos = findArg(argName,arguments);
  if(pos==0)
    return -999;
  if(pos<0)
    return type;
  return (int) parseLong(argName,arguments->argv[pos]);
}

float angsd::getArg(const char* argName,float type,argStruct *arguments){
  int pos = findArg(argName,arguments);
  if(pos==0)
    return -999;
  if(pos<0)
    return type;
  return (float) parseDouble(argName,arguments->argv[pos]);
}

double angsd::getArg(const char* argName,double type,argStruct *arguments){
  int pos = findArg(argName,arguments);
  if(pos==0)
    return -999;
  if(pos<0)
    return type;
  return parseDouble(argName,arguments->argv[pos]);
}

//the returned string is always newly allocated (the caller must free it), except when the flag isn't found where 'type' is returned as is
char* angsd::getArg(const char* argName,char* type,argStruct *arguments){
  int pos = findArg(argName,arguments);
  if(pos==0)
    return strdup("-999");
  if(pos<0)
    return type;
  return strdup(arguments->argv[pos]);
}

//as above, but since the default is const the default is returned as a newly allocated copy (or NULL)
char* angsd::getArg(const char* argName, const char* type,argStruct *arguments){
  int pos = findArg(argName,arguments);
  if(pos==0)
    return strdup("-999");
  if(pos<0)
    return type?strdup(type):NULL;
  return strdup(arguments->argv[pos]);
}



/*
  Reads a filelist, one or more filenames per line separated by spaces/tabs.
  Everything after a '#' is treated as a comment.
  If nInd>0 only the first nInd files are used.
*/
std::vector<char*> angsd::getFilenames(const char * name,int nInd){
  if(strchr(name,'\r')){
    fprintf(stderr,"\t\t-> Filelist contains carriage return. Looks like a windows file please remove hidden \'\\r\' from filelist\n");
    exit(1);
  }

  if(!aio::fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(1);
  }
  std::ifstream pFile(name,std::ios::in);
  if(!pFile.is_open()){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(1);
  }

  std::vector<char*> ret;
  std::string line;
  size_t lineNo=0;
  while(std::getline(pFile,line)){
    lineNo++;
    if(line.find('\r')!=std::string::npos){
      fprintf(stderr,"\t\t-> Filelist: \'%s\' contains carriage return at line %zu. Looks like a windows file please remove hidden \'\\r\' from: \'%s\'\n",name,lineNo,name);
      exit(1);
    }
    size_t hash = line.find('#');
    if(hash!=std::string::npos)
      line.erase(hash);
    std::istringstream iss(line);
    std::string tok;
    while(iss>>tok)
      ret.push_back(strdup(tok.c_str()));
  }

  if(ret.empty()){
    fprintf(stderr,"\t-> Filelist: \'%s\' doesn't contain any filenames\n",name);
    exit(1);
  }

  if(nInd>0) {
    if(ret.size()<(size_t)nInd){
      fprintf(stderr,"\t-> Number of samples in \'%s\' is smaller than the requested -nInd: %zu vs %d\n",name,ret.size(),nInd);
      exit(1);
    }
    for(size_t ii=nInd;ii<ret.size();ii++)
      free(ret[ii]);
    ret.erase(ret.begin()+nInd,ret.end());
  }

  return ret;
}

void destroy_argStruct(argStruct *args){
  free(args->cmdline);
  free(args->version);
  free(args->datetime);
  if(args->sm)
    bam_smpl_destroy(args->sm);
  for(unsigned i=0;i<args->nams.size();i++)
    free(args->nams[i]);
  free(args->fai);
  delete [] args->usedArgs;
  free(args->outfiles);
  free(args->infile);
  free(args->ref);
  free(args->anc);
  //args->rd is freed in bammer_main, args->revMap is owned by multiReader
  bam_hdr_destroy(args->hd);
  if(args->argumentFile!=stderr)
    fclose(args->argumentFile);
  delete args;
}
