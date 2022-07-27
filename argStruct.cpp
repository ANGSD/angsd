#include <cassert>
#include <libgen.h>//for checking if output dir exists 'dirname'
#include <fstream>
#include "version.h"
#include "argStruct.h"
#include "aio.h"

#define ARGS ".arg"
#define LENS 10000

void whatIsTheInput(int type){
  switch(type){
  case INPUT_BAM:
    fprintf(stderr,"\t-> Inputtype is BAM/CRAM\n");
    break;
  case INPUT_GLF:
    fprintf(stderr,"\t-> Inputtype is GLF\n");
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
    fprintf(stderr,"\t-> Unknown input type\n");
    exit(0);
   }
  }
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
  fprintf(stderr,"glf10_text:%d\n",INPUT_GLF10_TEXT);
  fprintf(stderr,"bgen:%d\n",INPUT_BGEN);
#endif
  if(args->fai){
    free(args->fai);
    args->fai=NULL;
  }
  args->fai = angsd::getArg("-fai",args->fai,args);
  char *tmp =NULL;
  tmp = angsd::getArg("-glf",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp = NULL;
  tmp = angsd::getArg("-glf3",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF3;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp = NULL;
  tmp = angsd::getArg("-glf10_text",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_GLF10_TEXT;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-vcf-GP",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GP;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-vcf-GL",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GL;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-vcf-PL",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_VCF_GL;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-pileup",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_PILEUP;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-beagle",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_BEAGLE;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-bgen",tmp,args);
  if(tmp!=NULL){
    args->inputtype=INPUT_BGEN;
    args->infile = tmp;
    args->nams.push_back(strdup(args->infile));
    return;
  }
  free(tmp); 
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
  free(tmp);
  tmp=NULL;
  tmp = angsd::getArg("-bam",tmp,args);
  tmp = angsd::getArg("-b",tmp,args);
  if(tmp!=NULL && args->argc>2){
    args->inputtype=INPUT_BAM;
    args->infile = tmp;
    args->nams = angsd::getFilenames(args->infile,nInd);
    
    return;
  }
  if(args->inputtype==INPUT_BAM){
    

  }
}

void checkIfDir(char *fname){
  char *dirNam2 = strdup(fname);
  char *dirNam=dirname(dirNam2);
  if(strlen(dirNam)>1 &&!aio::fexists(dirNam)){
    fprintf(stderr,"\t Folder: \'%s\' doesn't exist, please create\n",dirNam);
    exit(0);
  }
  free(dirNam2);
}

argStruct *setArgStruct(int argc,char **argv) { 

  argStruct *arguments = new argStruct;
  arguments->cmdline = NULL;
  arguments->version = NULL;
  arguments->argumentFile = stderr;
  arguments->outfiles =NULL;
  arguments->fai=arguments->anc=arguments->ref=NULL;
  arguments->hd=NULL;
  arguments->revMap= NULL;
  arguments->argc=argc;
  arguments->argv=argv;
  arguments->nReads = 50;
  arguments->sm=NULL;
  arguments->sm=bam_smpl_init();
  assert(arguments->sm);
 
  arguments->usedArgs= new int[argc+1];//well here we allocate one more than needed, this is only used in the ./angsd -beagle version
  for(int i=0;i<argc;i++)
    arguments->usedArgs[i]=0;
  
  arguments->inputtype=-1;
  arguments->infile = NULL;
  arguments->show =0;
  arguments->fai = NULL;

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
    if(argc!=2)
      fprintf(stderr,"\t-> No \'-out\' argument given, output files will be called \'angsdput\'\n");
    arguments->outfiles = strdup("angsdput");
  }
  checkIfDir(arguments->outfiles);
  
  if(argc<=2)
    return arguments;
  //print arguments into logfile
  if(argc>2)
    arguments->argumentFile=aio::openFile(arguments->outfiles,ARGS);
  fprintf(arguments->argumentFile,"\t-> %s",arguments->version);
  fprintf(arguments->argumentFile,"\t-> %s\n",arguments->cmdline);
  setInputType(arguments);
  whatIsTheInput(arguments->inputtype);
  return arguments;
}



int angsd::getArg(const char* argName,int type,argStruct *arguments){

  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      //      fprintf(stderr,"HIT %s vs %s\n",argName,arguments->argv[argPos]);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atoi(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}

char* angsd::getArg(const char* argName,char* type,argStruct *arguments){
  ///fprintf(stderr,"pre %p %s \n",type,argName);
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2){
        return(strdup("-999"));
      }
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(strdup(arguments->argv[argPos+1]));  //VALGRIND says leak. don't care very small DRAGON
    }
    argPos++;
  }
  //  fprintf(stderr,"post %p\n",type);
  return(type);
  
}


char* angsd::getArg(const char* argName, const char* type,argStruct *arguments){
  //fprintf(stderr,"pre %p %s \n",type,argName);
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2){
        return(strdup("-999"));
      }
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
      if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(strdup(arguments->argv[argPos+1]));  //VALGRIND says leak. don't care very small DRAGON
    }
    argPos++;
  }
  //assert(0==1);
  //  fprintf(stderr,"post %p\n",type);
  return NULL;
  
}



float angsd::getArg(const char* argName,float type,argStruct *arguments){
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
 if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atof(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}

double angsd::getArg(const char* argName,double type,argStruct *arguments){
  int argPos = 1;
  while(argPos <arguments->argc){
    if (strcasecmp(arguments->argv[argPos],argName)==0){
      if(arguments->argc==2)
        return(-999);
      arguments->usedArgs[argPos]=1;
      arguments->usedArgs[argPos+1]=1;
 if(argPos==arguments->argc-1){
	fprintf(stderr,"\t-> Must supply a parameter for: %s\n",argName);
	exit(0);
      }
      return(atof(arguments->argv[argPos+1]));  
    }
    argPos++;
  }
  return(type);
}



std::vector<char*> angsd::getFilenames(const char * name,int nInd){
  if(strchr(name,'\r')){
    fprintf(stderr,"\t\t-> Filelist contains carriage return. Looks like a windows file please remove hidden \'\r\' from filelist\n");
    exit(0);
  }
  
  if(!aio::fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(0);
  }
  const char* delims = " \t";
  std::vector<char*> ret;
  std::ifstream pFile(name,std::ios::in);

  char buffer[LENS];
  while(!pFile.eof()){
    pFile.getline(buffer,LENS);
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      if(tok[0]!='#')
	ret.push_back(strdup(buffer));
      tok = strtok(NULL,delims);
    }
  }
  if(nInd>0) {
     if(ret.size()<nInd)
      fprintf(stderr,"\t-> Number of samples is smaller than subset requested %lu vs %d\n",ret.size(),nInd);
    else{
      //   fprintf(stderr,"\t-> Will remove tail of filename list\n");
      for(int ii=nInd;ii<ret.size();ii++)
	free(ret[ii]);
      ret.erase(ret.begin()+nInd,ret.end());//we don't free this memory, it doesn't really matter
      // fprintf(stderr,"\t->  filename list now contains only: %lu\n",ret.size());
    }
     for(size_t ii=0;ii<ret.size();ii++)
       if(strchr(ret[ii],'\r')){
	 fprintf(stderr,"\t\t-> Filelist: \'%s\' contains carriage return. Looks like a windows file please remove hidden \'\r\' from: \'%s\'\n",name,ret[ii]);
	 exit(0);
       }
#if 0
     for(size_t ii=0;ii<ret.size();ii++)
       fprintf(stderr,"%zu->%s\n",ii,ret[ii]);
     fprintf(stderr,"\n");
#endif
  }

  return ret;
}

 void destroy_argStruct(argStruct *args){
   if(args->cmdline)
     free(args->cmdline);
   if(args->version)
     free(args->version);
   if(args->datetime)
     free(args->datetime);
   if(args->sm)
     bam_smpl_destroy(args->sm);
   
   for(unsigned i=0;i<args->nams.size();i++)
     free(args->nams[i]);
   if(args->fai)
     free(args->fai);
   delete []   args->usedArgs;
   free(args->outfiles);
   free(args->infile);
   if(args->anc)
     free(args->anc);
   bam_hdr_destroy(args->hd);
   if(args->argumentFile!=stderr) fclose(args->argumentFile);
   delete args;
 }
