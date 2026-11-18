#include "shared.h"
#include "analysisFunction.h"
#include "abcScounts.h"
#include "aio.h"

typedef unsigned char uchar;

typedef struct{
  unsigned int rel_pos;
  uchar A:2;
  uchar C:2;
  uchar G:2;
  uchar T:2;
}counts;


//read the sites file (chr pos allele1 allele2). The value of each site is its line number,
//which is what misc/scounts uses as rel_pos, so every line must be counted even when skipped.
aMap readvcf(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(1);
  }
  int at=0;
  char buf[1024];
  aMap am;
  while(gzgets(gz,buf,1024)){
    size_t blen = strlen(buf);
    if(blen==1023&&buf[blen-1]!='\n'){
      fprintf(stderr,"\t-> Line %d in sites file \'%s\' is longer than 1023 characters\n",at+1,fname);
      exit(1);
    }
    char chr[1024];
    int pos;
    char al1,al2;
    if(4!=sscanf(buf,"%1023s\t%d\t%c\t%c\n",chr,&pos,&al1,&al2)){
      fprintf(stderr,"\t-> Problem parsing line %d in sites file \'%s\', expected \'chr pos allele1 allele2\' but got: %s\n",at+1,fname,buf);
      exit(1);
    }
    char tmpnam[2048];
    snprintf(tmpnam,2048,"%s %d",chr,pos);
    aMap::iterator it=am.find(tmpnam);
    if(it!=am.end()){
      //keep the first occurrence, but still count the line so later indices match the sites file
      fprintf(stderr,"\t-> Duplicate position in sites file, using first occurrence: %s \n",tmpnam);
      at++;
      continue;
    }
    am[strdup(tmpnam)]= at++;
  }
  gzclose(gz);
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
  fprintf(stderr,"\t-> NOTE: -doScounts is deprecated and will be removed in a future version of ANGSD\n");
}

abcScounts::abcScounts(const char *outfiles,argStruct *arguments,int inputtype){
  vcfname = NULL;
  doScounts = 0;
  outfile = NULL;
  nSkippedDepth = 0;

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

  if(!vcfname){
    fprintf(stderr,"\t-> Must supply a file for screening (-vcfname)\n");
    exit(1);
  }
  int docounts =0;
  docounts = angsd::getArg("-docounts",docounts,arguments);
  if(docounts==0){
    fprintf(stderr,"\t-> -doScounts requires -doCounts\n");
    exit(1);
  }

  am = readvcf(vcfname);
  outfile = aio::openFileBG(outfiles,".scounts.gz");
}

abcScounts::~abcScounts(){
  if(nSkippedDepth>0)
    fprintf(stderr,"\t-> -doScounts skipped %lu sites where a base had depth>3\n",(unsigned long)nSkippedDepth);
  if(outfile!=NULL&&bgzf_close(outfile)!=0){
    fprintf(stderr,"\t-> Problem closing .scounts.gz output file\n");
    exit(1);
  }
  for(aMap::iterator it=am.begin();it!=am.end();++it)
    free((char*)it->first);
  free(vcfname);
}

void abcScounts::clean(funkyPars *pars){

}

void abcScounts::print(funkyPars *pars){
  if(doScounts==0)
    return;
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]!=0){
      char tmpname[1024];
      snprintf(tmpname,1024,"%s %d",header->target_name[pars->refId],pars->posi[s]+1);
      aMap::iterator it = am.find(tmpname);
      if(it==am.end()){
	fprintf(stderr,"\t-> problem finding site: %s\n",tmpname);
	continue;
      }
      counts cnts;
      //each count is stored in 2 bits, so sites with depth>3 for a base are skipped
      if(pars->counts[s][0]>3||pars->counts[s][1]>3||pars->counts[s][2]>3||pars->counts[s][3]>3){
	nSkippedDepth++;
	continue;
      }
      if(pars->counts[s][0]+pars->counts[s][1]+pars->counts[s][2]+pars->counts[s][3]==0)
	continue;
      cnts.rel_pos = it->second;
      // fprintf(stderr,"realpos: %d rel_pos:%d\n",pars->posi[s]+1,it->second);
      cnts.A = pars->counts[s][0];
      cnts.C = pars->counts[s][1];
      cnts.G = pars->counts[s][2];
      cnts.T = pars->counts[s][3];
      aio::doAssert(sizeof(counts)==bgzf_write(outfile,&cnts,sizeof(counts)*1),1,AT,"");

    }
  }

}


void abcScounts::run(funkyPars *pars){

}
