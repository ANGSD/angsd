#include <ctype.h>
#include "shared.h"
#include "analysisFunction.h"
#include "abcScounts.h"
#include <cassert>



typedef unsigned char uchar;

typedef struct{
  unsigned int rel_pos;
  uchar A:2;
  uchar C:2;
  uchar G:2;
  uchar T:2;    
}counts;


aMap readvcf(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(0);
  }
  int at=0;
  char buf[1024];
  aMap am;
  while(gzgets(gz,buf,1024)){
    char chr[1024];
    int pos;
    char al1,al2;
    if(4!=sscanf(buf,"%s\t%d\t%c\t%c\n",chr,&pos,&al1,&al2)){
      fprintf(stderr,"\t-> problem parsing line: %d which looks like: %s\n",at,buf);
    }
    char tmpnam[1024];
    sprintf(tmpnam,"%s %d",chr,pos);
    aMap::iterator it=am.find(tmpnam);
    if(it!=am.end()){
      fprintf(stderr,"\t-> Problem with duplicate positions: %s \n",tmpnam);
      //      exit(0);
    }
    am[strdup(tmpnam)]= at++;
  }
  fprintf(stderr,"\t-> read nsites from vcf:%lu\n",am.size());
  return am;
}




void abcScounts::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doScounts\t%d \n",doScounts);
  fprintf(argFile,"-vcfname\t%s \n",vcfname);
}

void abcScounts::getOptions(argStruct *arguments){
  doScounts=angsd::getArg("-doScounts",doScounts,arguments);
  vcfname = angsd::getArg("-vcfname",vcfname,arguments);
  
  if(doScounts==0){
    shouldRun[index]=0;
    return;
  }
}

abcScounts::abcScounts(const char *outfiles,argStruct *arguments,int inputtype){
  vcfname = NULL;
  doScounts = 0;
  outfile = NULL;

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doScounts")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  

  if(doScounts==0)
    return ;
  printArg(arguments->argumentFile);
  
  outfile = aio::openFileBG(outfiles,".scounts.gz");
  if(!vcfname){
    fprintf(stderr,"\t-> Must supply a file for screening\n");
    exit(0);
  }
  int docounts =0;
  docounts = angsd::getArg("-docounts",docounts,arguments);
  if(docounts==0){
    fprintf(stderr,"\t-> -doScounts requires -doCounts\n");
    exit(0);
  }
    
  am = readvcf(vcfname);
}

abcScounts::~abcScounts(){
  if(outfile!=NULL) 
    bgzf_close(outfile);
}

void abcScounts::clean(funkyPars *pars){

}

void abcScounts::print(funkyPars *pars){
  if(doScounts==0)
    return;
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]!=0){
      char tmpname[1024];
      sprintf(tmpname,"%s %d",header->target_name[pars->refId],pars->posi[s]+1);
      aMap::iterator it = am.find(tmpname);
      if(it==am.end()){
	fprintf(stderr,"\t-> problem finding site: %s\n",tmpname);
	continue;
      }
      counts cnts;
      if(pars->counts[s][0]>3||pars->counts[s][1]>3||pars->counts[s][2]>3||pars->counts[s][3]>3){
	fprintf(stderr,"\t-> skipping posi tmpname:%s du to depth>3\n",tmpname);
	continue;
      }
      if(pars->counts[s][0]+pars->counts[s][1]+pars->counts[s][2]+pars->counts[3]==0)
	continue;
      cnts.rel_pos = it->second;
      // fprintf(stderr,"realpos: %d rel_pos:%d\n",pars->posi[s]+1,it->second);
      cnts.A = pars->counts[s][0];
      cnts.C = pars->counts[s][1];
      cnts.G = pars->counts[s][2];
      cnts.T = pars->counts[s][3];
      assert(sizeof(counts)==bgzf_write(outfile,&cnts,sizeof(counts)*1));
      
    }
  }

}


void abcScounts::run(funkyPars *pars){

}


