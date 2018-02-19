#include <cmath>
#include <cstdlib>
#include <htslib/bgzf.h>
#include <assert.h>
#include <ctime>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abc.h"
#include "abcWriteFasta.h"

void abcWriteFasta::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doFasta\t%d\n",doFasta);
  fprintf(argFile,"\t1: use a random (non N) base (needs -doCounts 1)\n");
  fprintf(argFile,"\t2: use the most common (non N) base (needs -doCounts 1)\n");
  fprintf(argFile,"\t3: use the base with highest ebd (under development) \n");
  fprintf(argFile,"\t4: output iupac codes (under development) \n");
  fprintf(argFile,"\t-basesPerLine\t%d\t(Number of bases perline in output file)\n",NbasesPerLine);
  fprintf(argFile,"\t-explode\t%d\t print chromosome where we have no data (0:no,1:yes)\n",explode);
  fprintf(argFile,"\t-rmTrans\t%d\t remove transitions as different from -ref bases (0:no,1:yes)\n",rmTrans);
  fprintf(argFile,"\t-ref\t%s\t reference fasta, only used with -rmTrans 1\n",ref);
  fprintf(argFile,"\t-iupacRatio\t%.3f\t (Remove nucleotide below total depth ratio for IUPAC assignment)\n",iupacRatio);
  fprintf(argFile,"\t-seed\t%d\t use non random seed of value 1\n",seed);
  fprintf(argFile,"\n");
}

void abcWriteFasta::getOptions(argStruct *arguments){

  //from command line
  doFasta=angsd::getArg("-doFasta",doFasta,arguments);
  if(doFasta==0)
    return;  
  doCount=angsd::getArg("-doCounts",doCount,arguments);
  explode=angsd::getArg("-explode",explode,arguments);
  NbasesPerLine = angsd::getArg("-basesPerLine",NbasesPerLine,arguments);
  rmTrans=angsd::getArg("-rmTrans",rmTrans,arguments);
  ref=angsd::getArg("-ref",ref,arguments);
  iupacRatio=angsd::getArg("-iupacRatio",iupacRatio,arguments);
  seed=angsd::getArg("-seed",seed,arguments);


  if(doFasta){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doFasta \n");
      exit(0);
    }
    if(doFasta==3 && arguments->nInd!=1){
      fprintf(stderr,"Error: -doFasta 3 only works for a single individual\n");
      exit(0);
    }
    if((doFasta==2||doFasta==1) && doCount==0){
      fprintf(stderr,"Error: -doFasta 1 or 2 needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
  }
  if(rmTrans && ref==NULL){
    fprintf(stderr,"\t-> Must supply reference with -rmTrans 1\n");
    exit(0);
  }

}
void writeChr(kstring_t *bufstr,size_t len,char *nam,char*d,int nbpl){
  fprintf(stderr,"\t-> [%s] writing chr:%s\n",__FUNCTION__,nam);
  ksprintf(bufstr,">%s",nam);
  for(size_t i=0;i<len;i++){
    if(i % nbpl == 0)
      kputc('\n',bufstr);
    kputc(d!=NULL?d[i]:'N',bufstr);
  }
  kputc('\n',bufstr);

}


abcWriteFasta::abcWriteFasta(const char *outfiles,argStruct *arguments,int inputtype){
  explode = 0;
  myFasta = NULL;
  doFasta=0;
  doCount=0;
  currentChr=-1;
  NbasesPerLine=50;
  hasData =0;
  ref = NULL;
  iupacRatio = 0.0;
  rmTrans = 0;
  seed=0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doFasta")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);


  if(doFasta==0){
    shouldRun[index] = 0;
    return;
  }
  if(seed)
    srand(seed);
  else
    srand(time(0));

  
  printArg(arguments->argumentFile);


  for(int i=0;i<256;i++)
    lphred[i] =log(1.0-pow(i,-1.0*i/10.0));
  //make output files
  const char* postfix;
  postfix=".fa.gz";
  outfileZ = NULL;
  outfileZ = aio::openFileBG(outfiles,postfix);
  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;

}


abcWriteFasta::~abcWriteFasta(){
  if(ref)
    free(ref);
  if(doFasta==0)
    return;
  changeChr(-1);
  if(outfileZ!=NULL) bgzf_close(outfileZ); 
  if(bufstr.s!=NULL)
    free(bufstr.s);
}

void abcWriteFasta::changeChr(int refId) {
  if(doFasta==0)
    return;
 
  if(myFasta!=NULL){//proper case we have data
    if(explode||hasData){
      writeChr(&bufstr,header->target_len[currentChr],header->target_name[currentChr],myFasta,NbasesPerLine);
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    }
  }
  
  //ANDERS FILL IN MISSING CHRS IF YOU WANT HERE
  //ANDERS IS APPRANTLY LAZY SO NOW I'VE DONE IT FOR HIM
  //THORFINN THANK YOU FOR FIXING YOUR OWN MISTAKE


  if(refId!=-1){//-1 = destructor
    for(int i=currentChr+1;explode&&i<refId;i++){
      writeChr(&bufstr,header->target_len[i],header->target_name[i],NULL,NbasesPerLine);
      aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    }
    currentChr=refId;
    free(myFasta);
    myFasta=(char*)malloc(header->target_len[currentChr]);
    memset(myFasta,'N',header->target_len[currentChr]);
  }else{
    free(myFasta);
    for(int i=currentChr+1;explode&&i<header->n_targets;i++){
      writeChr(&bufstr,header->target_len[i],header->target_name[i],NULL,NbasesPerLine);
      aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    }
  }
  aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
}




void abcWriteFasta::run(funkyPars *pars){

  if(doFasta==0)
    return;


  hasData=1;
  if(doFasta==1){//random number read
    for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
      if(pars->keepSites[s]==0)
	continue;
      if(pars->nInd==1)
	myFasta[pars->posi[s]] = intToRef[ angsd::getRandomCount(pars->counts[s],0) ];
      else
	myFasta[pars->posi[s]] = intToRef[ angsd::getRandomCountTotal(pars->counts[s],pars->nInd) ];
    }     
  }
  else if(doFasta==2) {//most common
    for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
      if(pars->keepSites[s]==0)
	continue;
      if(pars->nInd==1)
	myFasta[pars->posi[s]] = intToRef[ angsd::getMaxCount(pars->counts[s],0) ];
      else
	myFasta[pars->posi[s]] = intToRef[ angsd::getMaxCountTotal(pars->counts[s],pars->nInd) ];
    }
  }else if(doFasta==3){
    for(int i=0;i<pars->nInd;i++){
      //      fprintf(stderr,"numSites: %d\n",pars->numSites);
      for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
	tNode *tn = pars->chk->nd[s][i];
	if(tn==NULL)
	  continue;
	double ebds[]= {0.0,0.0,0.0,0.0};
	for(int b=0;b<tn->l;b++){
	  int bof = refToInt[tn->seq[b]];
	  if(bof==4)
	    continue;
	  //	  fprintf(stderr,"pos:%d b:%d  bas:%c mapQ:%d qxcore:%c bof:%d qs:%d mapQ:%d %f %f \n",pars->posi[s],b,tn.seq[b],tn.mapQ[b],33+tn.qs[b],bof,tn.qs[b],tn.mapQ[b],lphred[tn.qs[b]],lphred[tn.mapQ[b]]);
	  ebds[bof] += exp(lphred[tn->qs[b]]+lphred[tn->mapQ[b]]);

	}
	for(int b=0;0&&b<4;b++)
	  fprintf(stderr,"b:%d %f\n",b,ebds[b]);
	int wh = angsd::whichMax(ebds,4);
	if(wh==-1) wh=4;//catch no information
	myFasta[pars->posi[s]] = intToRef[wh];
      }
    }
  }else if(doFasta==4){
    //supplied by kristian ullrich
    for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){
      if(pars->keepSites[s]==0)
  continue;
      if(pars->nInd==1)
        myFasta[pars->posi[s]] = intToIupac[ angsd::getIupacCount(pars->counts[s],0,iupacRatio) ];
      else
        myFasta[pars->posi[s]] = intToIupac[ angsd::getIupacCountTotal(pars->counts[s],pars->nInd,iupacRatio) ];
    }
  }
  //Do transitions removal
  if(rmTrans){
    assert(pars->ref!=NULL);
    for(int s=0;s<pars->numSites&&pars->posi[s]<header->target_len[pars->refId];s++){

      int ob = refToInt[myFasta[pars->posi[s]]];
      int rb = refToInt[pars->ref[pars->posi[s]]];
      //A <-> G, C <-> T
      if(ob==0&&rb==2)
	myFasta[pars->posi[s]]=intToRef[4];
      else if(ob==2&&rb==0)
	myFasta[pars->posi[s]]=intToRef[4];
      else if(ob==1&&rb==3)
	myFasta[pars->posi[s]]=intToRef[4];
      else if(ob==3&&rb==1)
	myFasta[pars->posi[s]]=intToRef[4];

    }

  }
}

