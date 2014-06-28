/*
  make a fasta file from a bam file. 
  thorfinn thorfinn@binf.ku.dk nov 28 2013   
  anders albrecht@binf.ku.dk made this.
  tsk tweaked it, to buffered gz dumping.
  tsk added ebd -doFasta 3
  part of angsd
*/



#include <cmath>
#include <cstdlib>
#include <zlib.h>
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include "kstring.h"//<-used for buffered output
#include "abc.h"
#include "abcWriteFasta.h"

typedef struct {
  char *seq;
  int start;
  int stop;
}funkyFasta;


void abcWriteFasta::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doFasta\t%d\n",doFasta);
  fprintf(argFile,"\t1: use a random base\n");
  fprintf(argFile,"\t2: use the most common base (needs -doCounts 1)\n");
  fprintf(argFile,"\t3: use the base with highest ebd (under development) \n");
  fprintf(argFile,"\t-basesPerLine\t%d\t(Number of bases perline in output file)\n",NbasesPerLine);
  fprintf(argFile,"\n");
}

void abcWriteFasta::getOptions(argStruct *arguments){

  //from command line
  doFasta=angsd::getArg("-doFasta",doFasta,arguments);
  doCount=angsd::getArg("-doCounts",doCount,arguments);
  NbasesPerLine = angsd::getArg("-basesPerLine",NbasesPerLine,arguments);
  
  if(doFasta){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doFasta \n");
      exit(0);
    }
    if(arguments->nInd!=1){
      fprintf(stderr,"Error: -doFasta only uses a single individual\n");
      exit(0);
    }
    if(doFasta==2 & doCount==0){
      fprintf(stderr,"Error: -doFasta 2 needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
  }


}

abcWriteFasta::abcWriteFasta(const char *outfiles,argStruct *arguments,int inputtype){
  doFasta=0;
  doCount=0;
  currentChr=-1;
  NbasesPerLine=50;
  lphred = NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doFasta")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doFasta==0){
    shouldRun[index] = 0;
    return;

  }
  lphred = new double[256];
  for(int i=0;i<256;i++)
    lphred[i] =log(1.0-pow(i,-1.0*i/10.0));
  //make output files
  const char* postfix;
  postfix=".fa.gz";
  outfileZ = NULL;
  outfileZ = aio::openFileGz(outfiles,postfix,GZOPT);
  //  const char* postfix2;
  //  postfix2=".fastaChr";
  //  outfile2 = openFile(outfiles,postfix2);

  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;

  currentPos=0;

}


abcWriteFasta::~abcWriteFasta(){

  if(doFasta==0)
    return;

  assert(bufstr.l==0);

  bufstr.l=0;
  
  extern int SIG_COND;

  while(currentChr<header->n_ref){
    //    fprintf(stderr,"currentpos %lu %d bufstr.l=%zu curChr:%s\n",currentPos,header->l_ref[currentChr],bufstr.l,header->name[currentChr]); 
    for(size_t p=currentPos;p<header->l_ref[currentChr];p++){
      if(p % NbasesPerLine == 0)
	kputc('\n',&bufstr);//fprintf(outfile,"\n");
      kputc('N',&bufstr);
      //      fprintf(outfile,"N");
    }
    if(bufstr.l!=0)
      gzwrite(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
  
    currentPos=0;
    currentChr++;
    if(currentChr<header->n_ref)
      ksprintf(&bufstr,"\n>%s",header->name[currentChr]);//fprintf(outfile,"\n>%s",header->name[currentChr]);
    
  }
  if(bufstr.l!=0)
    gzwrite(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
  

  if(outfileZ) gzclose(outfileZ); 
  if(bufstr.s!=NULL)
    free(bufstr.s);

}


void abcWriteFasta::clean(funkyPars *pars){

  if(doFasta==0)
    return;

  funkyFasta *fastaStruct =(funkyFasta *) pars->extras[index];
  delete[] fastaStruct->seq;
  delete fastaStruct;

}



void abcWriteFasta::print(funkyPars *pars){
  //   fprintf(stderr,"currentpos %lu %d\n",currentPos,header->l_ref[currentChr]);

  if(doFasta==0)
    return;

  funkyFasta *fastaStruct = (funkyFasta *) pars->extras[index];//new
  
  bufstr.l=0;
  
  if(currentChr==-1){
    currentPos=0;
    currentChr=0;
    //  fprintf(outfile,">%s",header->name[currentChr]);
    ksprintf(&bufstr,">%s",header->name[currentChr]);
  }

  while(currentChr!=pars->refId) {
    for(int p=currentPos;p<header->l_ref[currentChr];p++){
      if(p % NbasesPerLine == 0)
	kputc('\n',&bufstr);//fprintf(outfile,"\n");
      kputc('N',&bufstr);//fprintf(outfile,"N");
    }
    currentPos=0;
    currentChr++;
    //fprintf(outfile,"\n>%s",header->name[currentChr]);
    ksprintf(&bufstr,"\n>%s",header->name[currentChr]);
    if(bufstr.l!=0)
      gzwrite(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
  }

  for(int p=currentPos;p<fastaStruct->start;p++){
    if(p % NbasesPerLine == 0)
      kputc('\n',&bufstr);//fprintf(outfile,"\n");
    
    kputc('N',&bufstr);//fprintf(outfile,"N");
    currentPos++;
  
    
  }
  int c=0;
  for(int p=currentPos;p<fastaStruct->stop;p++){
    if(p==pars->posi[c]){
      if(p % NbasesPerLine == 0)
	kputc('\n',&bufstr);//fprintf(outfile,"\n");
      kputc(fastaStruct->seq[c],&bufstr);//fprintf(outfile,"%c",fastaStruct->seq[c]);
      //      fprintf(stdout,"%lu%d\t%c\n",currentPos+1,pars->posi[c],fastaStruct->seq[c]);
      c++;
    }
    else{
      if(p % NbasesPerLine == 0)
	kputc('\n',&bufstr);//fprintf(outfile,"\n");
      kputc('N',&bufstr);//fprintf(outfile,"N");
    }
    currentPos++;
    if(bufstr.l!=0)
      gzwrite(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
  }
    /*
  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    fprintf(stdout,"%c",fastaStruct->seq[s]);
  }
  fprintf(stdout,"\n");
  */
 
  if(bufstr.l>0)
    int ret =gzwrite(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    
}


void abcWriteFasta::run(funkyPars *pars){

  if(doFasta==0)
    return;

  funkyFasta *fastaStruct = new funkyFasta;//new

  char *seq;
  int start;
  int stop;
  start = pars->posi[0];
  stop = pars->posi[pars->numSites-1];
  //  seq = new char[stop-start+1];
  seq = new char[pars->numSites];
  for(int m=0;m<pars->numSites;m++)
    seq[m]='N';

 
  if(doFasta==1){//random number read
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      if(pars->chk->nd[s][0].l==0)
	continue;
      int j = std::rand() % pars->chk->nd[s][0].l;
      seq[s] = pars->chk->nd[s][0].seq[j];
      
    }
  }else if(doFasta==2) {//most common
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;

      int max=0;
      for( int b = 0; b < 4; b++ ){
	if( pars->counts[s][b] >= max){
	  max=pars->counts[s][b];
	}
      }
      if(max==0)
	continue;

      int nmax=0;
      int w=-1;
      for( int b = 0; b < 4; b++ ){
	if( pars->counts[s][b] == max ){
	  w=b;
	  nmax++;
	}
      }

      if(nmax==1)
	seq[s] = intToRef[w];
      else{
	int cumsum[4] = {-1,-1,-1,-1};
	int i=0;
	for( int b = 0; b < 4; b++ ){
	  if( pars->counts[s][b] == max ){
	    cumsum[b] = i;
	    i++;
	  }
 	}

	int j = std::rand() % i;

	for( int b = 0; b < 4; b++ ){
	  if( cumsum[b] ==j ){
	    seq[s] = intToRef[b];
	    break;
	  }
	}
      }
    }
  }else if(doFasta==3){
    for(int i=0;i<pars->nInd;i++){
      //      fprintf(stderr,"numSites: %d\n",pars->numSites);
      for(int s=0;s<pars->numSites;s++){
	tNode &tn = pars->chk->nd[s][i];
	double ebds[]= {0.0,0.0,0.0,0.0};
	for(int b=0;b<tn.l;b++){
	  int bof = refToInt[tn.seq[b]];
	  if(bof==4)
	    continue;
	  //	  fprintf(stderr,"pos:%d b:%d  bas:%c mapQ:%d qxcore:%c bof:%d qs:%d mapQ:%d %f %f \n",pars->posi[s],b,tn.seq[b],tn.mapQ[b],33+tn.qs[b],bof,tn.qs[b],tn.mapQ[b],lphred[tn.qs[b]],lphred[tn.mapQ[b]]);
	  ebds[bof] += exp(lphred[tn.qs[b]]+lphred[tn.mapQ[b]]);

	}
	for(int b=0;0&&b<4;b++)
	  fprintf(stderr,"b:%d %f\n",b,ebds[b]);
	int wh = angsd::whichMax(ebds,4);
	if(wh==-1) wh=4;//catch no information
	seq[s] = intToRef[wh];
	//exit(0);
      }
      //      exit(0);
    }
  }



  fastaStruct->seq=seq;
  fastaStruct->start=start;
  fastaStruct->stop=stop;
  pars->extras[index] = fastaStruct;
}

