/*
  make a fasta file from a bam file. 
  thorfinn thorfinn@binf.ku.dk dec17 2012   
  anders albrecht@binf.ku.dk made this.
  part of angsd
*/


#include <cmath>
#include <cstdlib>

#include "shared.h"
#include "analysisFunction.h"

#include <htslib/kstring.h>
#include "abc.h"

#include "abcDstat.h"



void abcDstat::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa\t\t%d\trun the abbababa analysis\n",doAbbababa);
  // fprintf(argFile,"\t1: use a random base\n");
  //  fprintf(argFile,"\t2: use the most common base (needs -doCounts 1)\n");
  fprintf(argFile,"\t-rmTrans\t\t%d\tremove transitions\n",rmTrans);
  fprintf(argFile,"\t-blockSize\t\t%d\tsize of each block in bases\n",blockSize);
  fprintf(argFile,"\t-enhance\t\t%d\toutgroup must have same base for all reads\n",enhance);
  fprintf(argFile,"\t-ans\t\t\t%s\tfasta file with outgroup\n",ancName);
  fprintf(argFile,"\t-useLast\t\t%d\tuse the last individuals as outgroup instead of -anc\n",useLast);
  fprintf(argFile,"\t-printEmpty\t\t%d\t print blocks without any ABBA or BABA sites \n",printEmpty);
  fprintf(argFile,"\t-seed\t%d\t use non random seed of value 1\n",seed);

  fprintf(argFile,"\n");
}

void abcDstat::getOptions(argStruct *arguments){
    //from command line
  doAbbababa=angsd::getArg("-doAbbababa",doAbbababa,arguments);
  if(doAbbababa==0)
    return;

  doCount=angsd::getArg("-doCounts",doCount,arguments);
  blockSize=angsd::getArg("-blockSize",blockSize,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);
  rmTrans = angsd::getArg("-rmTrans",rmTrans,arguments);
  Aanc = angsd::getArg("-Aanc",Aanc,arguments);
  useLast =  angsd::getArg("-useLast",useLast,arguments);
  enhance =  angsd::getArg("-enhance",enhance,arguments);
  printEmpty =  angsd::getArg("-printEmpty",printEmpty,arguments);
  seed=angsd::getArg("-seed",seed,arguments);
  if(useLast != 0)
    useLast = 1;


  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"Error: bam or soap input needed for -doAbbababa \n");
    exit(0);
  }
  if(arguments->nInd<3+useLast){
    fprintf(stderr,"Error: -doAbbababa needs atleast 3 individual + outgroup\n");
    exit(0);
  }
  if( useLast == 0 && enhance == 1){
    fprintf(stderr,"Error: -enhance only works with -useLast 1\n");
    exit(0);
  }
  if( doCount==0){
    fprintf(stderr,"Error: -doAbbababa needs allele counts (use -doCounts 1)\n");
    exit(0);
  }
  if(ancName==NULL & Aanc==0 & useLast == 0){
    fprintf(stderr,"Error: -doAbbababa needs an outgroup in fasta format (use -anc fastaFileName ) or (-useLast 1)\n");
    exit(0);
    }
  


}

abcDstat::abcDstat(const char *outfiles,argStruct *arguments,int inputtype){
  rmTrans = 0;
  outfile=NULL;
  ancName=NULL;
  doAbbababa=0;
  doCount=0;
  Aanc = 0;
  printEmpty = 1;
  useLast = 0;
  currentChr=-1;
  NbasesPerLine=50;
  blockSize=5000000;
  enhance=0;
  seed=0;
  block=0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAbbababa")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);


  if(doAbbababa==0){
    shouldRun[index] = 0;
    return;
  }
  printArg(arguments->argumentFile);


  if(seed)
    srand48(seed);
  //make output files
  const char* postfix;
  postfix=".abbababa";
  outfile = aio::openFile(outfiles,postfix);
  //  const char* postfix2;
  //  postfix2=".fastaChr";
  //  outfile2 = openFile(outfiles,postfix2);

  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;


  //the number of combinations
  nComb=(arguments->nInd -useLast)*(arguments->nInd -1 -useLast)*(arguments->nInd -2 -useLast)/2;
  ABBA = new int[nComb];
  BABA = new int[nComb];


  calcMatCat();
  //pre calc all posible combinations


}

int category(int x1, int x2, int x3, int x4){
  if(x1==4||x2==4||x3==4||x4==4)
    return 0; //missing
  if(x1==x2&&x3==x4&&x1==x3)
    return 1; //AAAA
  else if(x1==x2&&x4==x2&&x1!=x3)
    return 2; //AABA
  else if(x1==x3&&x3==x4&&x1!=x2)
    return 3; //ABAA
  else if(x1==x4&&x3==x2&&x1!=x3)
    return 4; //ABBA
  else if(x1==x4&&x3!=x2&&x1!=x2&&x3!=x1)
    return 5; //ABCA
  else if(x2==x4&&x3==x2&&x1!=x3)
    return 6; //BAAA
  else if(x1==x3&&x4==x2&&x1!=x2)
    return 7; //BABA
  else if(x2==x4&&x3!=x2&&x1!=x2&&x1!=x3)
    return 8; //BACA
  else if(x2==x1&&x3==x4&&x1!=x3)
    return 9; //BBAA
  else if(x2==x1&&x3==x2&&x1!=x4)
    return 10; //BBBA
  else if(x2==x1&&x3!=x4&&x1!=x3&&x1!=x4)
    return 11; //BBCA
  else if(x3==x4&&x3!=x2&&x1!=x2&&x1!=x3)
    return 12; //BCAA
  else if(x3==x1&&x3!=x4&&x1!=x2&&x2!=x4)
    return 13; //BCBA
  else if(x3==x2&&x3!=x4&&x3!=x1&&x1!=x4)
    return 14; //CBBA
  else if(x2!=x1&&x1!=x3&&x1!=x4&&x2!=x3&&x2!=x4&&x3!=x4)
    return 15; //BCDA
  else{
    return -999;
    fprintf(stdout,"Warning: Not in any category!");
  }
}


void abcDstat::calcMatCat(){

  for(int b1=0; b1<5; b1++){
    for(int b2=0; b2<5; b2++){
      for(int b3=0; b3<5; b3++){
	for(int b4=0; b4<5; b4++){
	  matcat[b1][b2][b3][b4]=category(b1,b2,b3,b4);
	}    
      }    
    }    
  }
}



abcDstat::~abcDstat(){


  free(ancName);
  if(doAbbababa==0)
    return;


  if(currentChr > -1)//if first chunk
    printAndEmpty();
  if(outfile) fclose(outfile);
  //  fclose(outfile2);
  if(bufstr.s!=NULL)
    free(bufstr.s);
  delete[] ABBA;
  delete[] BABA;

}


void abcDstat::clean(funkyPars *pars){

  if(doAbbababa==0)
    return;

  funkyAbbababa *abbababaStruct =(funkyAbbababa *) pars->extras[index];

  for(int i=0; i < abbababaStruct->nBlocks ; i++)
    for(int c=0; c<nComb ; c++)
      delete[] abbababaStruct->ABBABABAblocks[i][c];
  for(int i=0; i < abbababaStruct->nBlocks ; i++)
    delete[] abbababaStruct->ABBABABAblocks[i];
  delete[] abbababaStruct->ABBABABAblocks;
  
  delete[] abbababaStruct->blockPos;
  
  delete abbababaStruct;

}

void abcDstat::printAndEmpty(){


  
  if(printEmpty == 0 && 1){
    int total=0;
    for(int j=0;j<nComb;j++)
      total += ABBA[j] + BABA[j];      
    if(total==0)
      return;
    
  }
 
  fprintf(outfile,"%s\t%d\t%d",header->target_name[currentChr],block*blockSize,block*blockSize+blockSize-1);
  for(int j=0;j<nComb;j++){
    fprintf(outfile,"\t%d\t%d",ABBA[j],BABA[j]);
    ABBA[j]=0;
    BABA[j]=0;
  }
  fprintf(outfile,"\n");

  fflush(outfile);
}


void abcDstat::getBlockNum(int pos){
  block=(int)((pos+1)/blockSize);
}


// returns number of blocks in chunk
int abcDstat::getNumBlocks(funkyPars *pars){

  int nBlocks = 1;
  int blockHere = -1;// (int)((pars->posi[0]+1)/blockSize);
  for(int s = 0 ; s < pars->numSites; s++){
    if( pars->keepSites[s]==0 )
      continue;

    if( pars->posi[s]+1 >= blockHere*blockSize+blockSize ){
      nBlocks++;
      blockHere =  (int)((pars->posi[s]+1)/blockSize);
    }
  }
  return(nBlocks);
}


void abcDstat::print(funkyPars *pars){

  //fprintf(stderr,"currentChr %d %d\t refid %s \t pos[0] %d-%d\tnumSites %d\n",currentChr,header->target_name[currentChr],header->target_name[pars->refId],pars->posi[0]+1,pars->posi[pars->numSites-1]+1,pars->numSites);

  // for(int i=0;i<pars->numSites;i++)
  //if(pars->keepSites[i])
  //  fprintf(stderr,"pos %d\n",pars->posi[i]+1);
 
  
  if(doAbbababa==0)
    return;

  

  funkyAbbababa *abbababaStruct = (funkyAbbababa *) pars->extras[index];//new
  

  
  if(currentChr==-1){//if first chunk
    //start new block
    for(int j=0;j<nComb;j++){
      ABBA[j]=0;
      BABA[j]=0;
    }
    getBlockNum(pars->posi[0]);
    if(currentChr > pars->refId)
      fprintf(stdout,"Warning. Your regions are not sorted. Becarefull...\n");
    currentChr=pars->refId;
  }

  while(currentChr!=pars->refId){//if new chr (not first)
    //start new block
    printAndEmpty();
    //    currentChr++;
    currentChr = pars->refId;
    block = abbababaStruct->blockPos[0];
  }

  // fprintf(stdout,"nBlocks=%d\t",abbababaStruct->nBlocks);
  for(int b=0; b < abbababaStruct->nBlocks;b++){
    //  fprintf(stdout,"block %d-%d\t",abbababaStruct->blockPos[b],abbababaStruct->ABBABABAblocks[b][0][4]);
    
    if(abbababaStruct->blockPos[b] > block){
      printAndEmpty();
      block = abbababaStruct->blockPos[b];
    }
    //   fprintf(stdout,"insert: bPos-%d-%d ABBA-%d\n",abbababaStruct->blockPos[0],block,abbababaStruct->ABBABABAblocks[0][0][4]);
    for(int comp=0; comp<nComb;comp++){

	  ABBA[comp] += abbababaStruct->ABBABABAblocks[b][comp][4] ;
	  BABA[comp] += abbababaStruct->ABBABABAblocks[b][comp][7] ;
    }   
  }
  
  //fprintf(stdout,"\n");
}


void abcDstat::run(funkyPars *pars){

  if(doAbbababa==0)
    return;

  funkyAbbababa *abbababaStruct = new funkyAbbababa;//new

  abbababaStruct->nBlocks = getNumBlocks(pars);
  abbababaStruct->ABBABABAblocks = new int**[abbababaStruct->nBlocks];
  
  for(int i=0; i < abbababaStruct->nBlocks ; i++){
    abbababaStruct->ABBABABAblocks[i] = new int*[nComb];
    for(int c=0; c<nComb ; c++){
      abbababaStruct->ABBABABAblocks[i][c] = new int[256];
      for(int type=0;type<256;type++)
	abbababaStruct->ABBABABAblocks[i][c][type] = 0;

    }
  }
  abbababaStruct->blockPos = new int[abbababaStruct->nBlocks];
  for(int b=0;b<abbababaStruct->nBlocks;b++)
    abbababaStruct->blockPos[b] = 0 ;
  
  
  

  int **pattern;
  pattern = new int*[pars->numSites];
  for(int s=0;s<pars->numSites;s++){
    pattern[s]=new int[pars->nInd];
    for(int i=0;i<pars->nInd;i++)
      pattern[s][i]=4;


  }


 

  // sample an allele for each individuals
  if( doAbbababa!=0 ){//random number read
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      
      for(int i=0;i<pars->nInd;i++){
	int dep=0;
	for( int b = 0; b < 4; b++ ){
	  dep+=pars->counts[s][i*4+b];
	}
	if(dep==0)
	  continue;

	if( enhance == 1 && useLast == 1 && i == pars->nInd -1 ){ //only use outgroup if all bases are the same
	  pattern[s][i] = 4;
	  
	  for( int b = 0; b < 4; b++ )
	    if(dep == pars->counts[s][i*4+b]){
	      pattern[s][i] = b;
	      //fprintf(stdout,"s \n");
	      break;
	    }
	  
	 
	}
	else{
	  int j = std::rand() % dep;
	  int cumSum=0;
	  for( int b = 0; b < 4; b++ ){
	    cumSum+=pars->counts[s][i*4+b];
	    if( cumSum > j ){
	      pattern[s][i] = b;
	      break;
	    }
	  }
	}
      }     
    }
  }


  //get pattern
  int comp=0;
  for(int h3=0; h3<(pars->nInd-useLast); h3++){
    for(int h2=0; h2<(pars->nInd-useLast); h2++){
      if(h2==h3)
	continue;
      for(int h1=0; h1<h2; h1++){
	if(h1==h3)
	  continue;
	if(rmTrans){
	  int blockHere = (int)((pars->posi[0]+1)/blockSize);
	  int nBlocks = 0;
	  abbababaStruct->blockPos[nBlocks] = blockHere;
	  for(int s=0;s<pars->numSites;s++){
	    int theAnc=4;
	    if(Aanc==0 & useLast==0)
	      theAnc =pars->anc[s];
	    else if(useLast == 1){
	      theAnc=pattern[s][pars->nInd-1];
	    }

	    if(pars->keepSites[s]==0)
	      continue;

	    if( pars->posi[s]+1 >= blockHere*blockSize+blockSize ){
	      nBlocks++;
	      blockHere =  (int)((pars->posi[s]+1)/blockSize);
	      abbababaStruct->blockPos[nBlocks] = blockHere;
	    }

	    int ABCDsite=0;
	    if(theAnc==0 && pattern[s][h3]==2)
	      ABCDsite=0;
	    else if(theAnc==2 && pattern[s][h3]==0)
	      ABCDsite=0;
	    else if(theAnc==1 && pattern[s][h3]==3)
	      ABCDsite=0;
	    else if(theAnc==3 && pattern[s][h3]==1)
	      ABCDsite=0;
	    else 
	      ABCDsite = matcat[pattern[s][h1]] [pattern[s][h2]] [pattern[s][h3]] [theAnc];
	    abbababaStruct->ABBABABAblocks[nBlocks][comp][ABCDsite]++;
	    
	  }  
	}
	else{
	  int blockHere = (int)((pars->posi[0]+1)/blockSize);
	  int nBlocks = 0;
	  
	  abbababaStruct->blockPos[nBlocks] = blockHere;

	  for(int s=0;s<pars->numSites;s++){
	      int theAnc=4;
	    if(Aanc==0 & useLast==0)
	      theAnc =pars->anc[s];
	    else if(useLast){
	      theAnc=pattern[s][pars->nInd-1];
	    }

	    if(pars->keepSites[s]==0)
	      continue;

	    if( pars->posi[s]+1 >= blockHere*blockSize+blockSize ){
	      nBlocks++;
	      blockHere =  (int)((pars->posi[s]+1)/blockSize);
	      abbababaStruct->blockPos[nBlocks] = blockHere;
	    }

  	    int ABCDsite = matcat[pattern[s][h1]] [pattern[s][h2]] [pattern[s][h3]] [theAnc];
	    abbababaStruct->ABBABABAblocks[nBlocks][comp][ABCDsite]++;
	    //	    if(comp==0 && ABCDsite==4)
	    // fprintf(stdout,"ABCD=%d-%d-%d\t",ABCDsite,nBlocks,abbababaStruct->ABBABABAblocks[nBlocks][comp][ABCDsite]);

	    //  	    ABCD[s][comp] = matcat[2] [3] [2] [3];
	    //	    if(ABCD[s][comp]!=1)
	    //	      fprintf(stdout,"%d \n",ABCD[s][comp]);
	  }
	}
	comp++;
      }
    }
  }



  // fprintf(stdout,"run: bPos-%d ABBA-%d\n",abbababaStruct->blockPos[0],abbababaStruct->ABBABABAblocks[0][0][4]);
 
  
  for(int s=0;s<pars->numSites;s++)
    delete[] pattern[s];
  delete[] pattern;


  pars->extras[index] = abbababaStruct;
}

