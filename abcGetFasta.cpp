/*
  thorfinn thorfinn@binf.ku.dk 18dec 2012
  assumes zero index positions

  should remove mutexs, and make function thread safe. dragon daemon

  part of angsd using faidx
 
  Actually no need to mutex, except for the first chunk. Lets fix it at some point.
*/

#include <cassert>

#include "analysisFunction.h"
#include "shared.h"

#include "pthread.h"
#include "abc.h"
#include "abcGetFasta.h"


//this will initialize our data
perFasta *init(const char *fname){

  fprintf(stderr,"\t-> Reading fasta: %s\n",fname);
  if(aio::fexists(fname)==0){
    fprintf(stderr,"\t-> fastafile: \'%s\' doesn't exists, will exit\n",fname);
    exit(0);
  }
  //check that fa hasn't been updated
  char *strtsk=NULL;
  strtsk = (char*)calloc(strlen(fname) + 5, 1);
  sprintf(strtsk, "%s.fai", fname);
  if(aio::isNewer(fname,strtsk)){
    fprintf(stderr,"\t-> fai index file: \'%s\' looks older than corresponding fastafile: \'%s\'.\n\t-> Please reindex fasta file\n",strtsk,fname);
    exit(0);
  }
  free(strtsk);
    
  perFasta *r= new perFasta;
  r->fastaname=NULL;
  r->fai = NULL;
  r->seqs = NULL;
  r->curChr = -1;
  if(pthread_mutex_init(&r->aMut,NULL)){
    fprintf(stderr,"[%s:%s] error initializing mutex \n",__FILE__,__FUNCTION__);
    exit(0);
  }

  if(NULL==(r->fai = fai_load(fname))){
    fprintf(stderr,"[%s:%s] error reading fai file:%s\n",__FILE__,__FUNCTION__,fname);
    exit(0);
  }
  r->fastaname=strdup(fname);
  return r;
}



void abcGetFasta::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n\n",__FILE__);
  fprintf(argFile,"\t-ref\t%s\t(afile.fasta)\n",refName);
  fprintf(argFile,"\t-anc\t%s\t(afile.fasta)\n",ancName);
  fprintf(argFile,"\tNB these fasta files should be indexed 'samtools faidx'\n");
  fprintf(argFile,"\n");

}


//this will destroy a perfasta structure
void destroy(perFasta *f){
  fai_destroy(f->fai);
  free(f->seqs);
  free(f->fastaname);
  delete f;
  f=NULL;
}

abcGetFasta::~abcGetFasta(){
  if(ref)
    destroy(ref);
  if(anc)
    destroy(anc);
  free(refName);
  free(ancName);
}

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));//<-superstupid hack to be able to parse crams. This should be fixed
void abcGetFasta::getOptions(argStruct *arguments){

  refName = NULL;
  ancName = NULL;
  //from command line
  ancName = angsd::getArg("-anc",ancName,arguments);
  refName = angsd::getArg("-ref",refName,arguments);
  if(refName==NULL)
    refName = angsd::getArg("-T",refName,arguments);
  //  fprintf(stderr,"refName=%s\tancName=%s\n",refName,ancName);

  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-ref")|| !strcasecmp(arguments->argv[1],"-anc") ){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    sprintf(ref, "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding2->specific,ref);
    free(ref);
  }
}

char *abcGetFasta::loadChr(perFasta *f, char*chrName,int chrId){
  //  fprintf(stderr,"[%s] \t->loading chr:%s from faidx=%p curchr=%d\n",__FUNCTION__,chrName,f,f->curChr);
  free(f->seqs);
  f->seqs=NULL;
  //fprintf(stderr,"f->curChr=%d chrId=%d\n",f->curChr,chrId);
  f->curChr = chrId;
  f->seqs = faidx_fetch_seq(f->fai, chrName, 0, 0x7fffffff, &f->chrLen);
  if(f->seqs==NULL){
    fprintf(stderr,"\n[%s] Error loading fasta info from chr:\'%s\' \n",__FUNCTION__,chrName);
    f->chrLen=0;
  }
  //  fprintf(stderr,"[%s] done\n",__FUNCTION__);
  if(f->chrLen!=header->target_len[chrId]){
    fprintf(stderr,"\t-> Problem with length of fastafile vs length of chr in BAM header\n");
    fprintf(stderr,"\t-> Chromosome name: \'%s\' length from BAM header:%d length from fai file:%d\n",chrName,header->target_len[chrId],f->chrLen);
    extern int SIG_COND;
    SIG_COND = 0;
  }
  return f->seqs;
}



char *abcGetFasta::magic(int refId,int *posi,int numSites,perFasta *f){
  pthread_mutex_lock(&f->aMut);
  assert(refId!=-1);
  //load new chr if different from last or if nothing has been loaded before
  if(f->curChr==-1||refId!=f->curChr){
    //    fprintf(stderr,"[%s] chaning to chr: %d\n",__FUNCTION__,refId);
    //fflush(stderr);
    loadChr(f,header->target_name[refId],refId);
  }
  //first check that last position is not to long after.
  if(posi[numSites-1] > f->chrLen+200){
    fprintf(stderr,"Trying to access fasta efter end of chromsome+200:%s/%s pos=%d ref_len=%d\n",header->target_name[f->curChr],header->target_name[refId],posi[numSites-1],f->chrLen);
  }
  //now loop over all positions. We allow +200 since a read might match to the end of a chr
  char *ret=new char[numSites];
    
  for(int i=0;i<numSites;i++){
    if(f->seqs&&posi[i]<f->chrLen)
      ret[i] = refToInt[f->seqs[posi[i]]];//offseting by one, to make in work like samtools
    else{
      ret[i] = 4;
    }
  }  

  pthread_mutex_unlock(&f->aMut);
  return ret;
}


 abcGetFasta::abcGetFasta(argStruct *arguments){
  refName = NULL;
  ancName = NULL;
  ref=NULL;
  anc=NULL;
   if(arguments->argc==2){
     if(!strcasecmp(arguments->argv[1],"-ref")||!strcasecmp(arguments->argv[1],"-anc")){
       printArg(stdout);
       exit(0);
    }else
       return;
   }

   getOptions(arguments);


   if(ref!=NULL || anc!=NULL)
     printArg(arguments->argumentFile);


   if(ancName!=NULL)
     anc = init(ancName);
  
   if(refName!=NULL)
    ref = init(refName);
   
   if(ref==NULL&&anc==NULL)
     shouldRun[index]=0;
}

void abcGetFasta::run(funkyPars *pars){
   if(pars->numSites==0)
    return;
  if(ref)
    pars->ref = magic(pars->refId,pars->posi,pars->numSites,ref);
  if(anc)
    pars->anc = magic(pars->refId,pars->posi,pars->numSites,anc);
  
}


void abcGetFasta::print(funkyPars *f){}
void abcGetFasta::clean(funkyPars *f){}
