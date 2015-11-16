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
typedef struct {
  int **ABCD;
}funkyAbbababa;


void abcDstat::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa\t\t%d\trun the abbababa analysis\n",doAbbababa);
  // fprintf(argFile,"\t1: use a random base\n");
  //  fprintf(argFile,"\t2: use the most common base (needs -doCounts 1)\n");
  fprintf(argFile,"\t-rmTrans\t\t%d\tremove transitions\n",rmTrans);
  fprintf(argFile,"\t-blockSize\t\t%d\tsize of each block in bases\n",blockSize);
  fprintf(argFile,"\t-ans\t\t\t%s\tfasta file with outgroup\n",ancName);
  fprintf(argFile,"\n");
}

void abcDstat::getOptions(argStruct *arguments){
    //from command line
  doAbbababa=angsd::getArg("-doAbbababa",doAbbababa,arguments);
 
  doCount=angsd::getArg("-doCounts",doCount,arguments);
  blockSize=angsd::getArg("-blockSize",blockSize,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);
  rmTrans = angsd::getArg("-rmTrans",rmTrans,arguments);
    
  if(doAbbababa){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doAbbababa \n");
      exit(0);
    }
    if(arguments->nInd<3){
      fprintf(stderr,"Error: -doAbbababa needs atleast 3 individual\n");
      exit(0);
    }
    if( doCount==0){
      fprintf(stderr,"Error: -doAbbababa needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
    if(ancName==NULL){
      fprintf(stderr,"Error: -doAbbababa needs an outgroup in fasta format (use -anc fastaFileName )\n");
      exit(0);
    }
  }


}

abcDstat::abcDstat(const char *outfiles,argStruct *arguments,int inputtype){
  rmTrans = 0;
  outfile=NULL;
  ancName=NULL;
  doAbbababa=0;
  doCount=0;
  
  currentChr=-1;
  NbasesPerLine=50;
  blockSize=5000000;
  block=0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAbbababa")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAbbababa==0){
    shouldRun[index] = 0;
    return;
}
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
  nComb=arguments->nInd*(arguments->nInd-1)*(arguments->nInd-2)/2;
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

  for(int s=0;s<pars->numSites;s++){
    delete[] abbababaStruct->ABCD[s];
  }
  delete[] abbababaStruct->ABCD;
  delete abbababaStruct;

}

void abcDstat::printAndEmpty(){

  fprintf(outfile,"%s\t%d\t%d",header->target_name[currentChr],block*blockSize+1,block*blockSize+blockSize);
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


void abcDstat::print(funkyPars *pars){
  //   fprintf(stderr,"currentpos %lu %d\n",currentPos,header->l_ref[currentChr]);

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
    currentChr=0;
  }

  while(currentChr!=pars->refId){//if new chr (not first)
    //start new block
    printAndEmpty();
    currentChr++;
    getBlockNum(pars->posi[0]);
  }





  for(int s=0;s<pars->numSites;s++){
    int comp=0;

    if(pars->posi[s]>=block*blockSize+blockSize){
      printAndEmpty();
      getBlockNum(pars->posi[s]);
    }

    if(pars->keepSites[s]==0)
	continue;
    for(int h3=0; h3<(pars->nInd); h3++){
      for(int h2=0; h2<(pars->nInd); h2++){
	if(h2==h3)
	  continue;
	for(int h1=0; h1<h2; h1++){
	  if(h1==h3)
	    continue;
	  if(abbababaStruct->ABCD[s][comp]==4)
	    ABBA[comp]++;
	  else if(abbababaStruct->ABCD[s][comp]==7){
	    BABA[comp]++;
	  }
	    comp++;
	}
      }
    }
  }

}


void abcDstat::run(funkyPars *pars){

  if(doAbbababa==0)
    return;

  funkyAbbababa *abbababaStruct = new funkyAbbababa;//new

  int **ABCD;
  int **pattern;
  ABCD = new int*[pars->numSites];
  pattern = new int*[pars->numSites];
  for(int s=0;s<pars->numSites;s++){
    pattern[s]=new int[pars->nInd];
    for(int i=0;i<pars->nInd;i++)
      pattern[s][i]=4;

    ABCD[s] = new int[nComb];
  }


 

  // sample an allele for each individuals
  if(doAbbababa!=0){//random number read
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


  //get pattern
  int comp=0;
  for(int h3=0; h3<(pars->nInd); h3++){
    for(int h2=0; h2<(pars->nInd); h2++){
      if(h2==h3)
	continue;
      for(int h1=0; h1<h2; h1++){
	if(h1==h3)
	  continue;
	if(rmTrans){
	  for(int s=0;s<pars->numSites;s++){
	    if(pars->keepSites[s]==0)
	      continue;
 
	    if(pars->anc[s]==0 && pattern[s][h3]==2)
	      ABCD[s][comp]=0;
	    else if(pars->anc[s]==2 && pattern[s][h3]==0)
	      ABCD[s][comp]=0;
	    else if(pars->anc[s]==1 && pattern[s][h3]==3)
	      ABCD[s][comp]=0;
	    else if(pars->anc[s]==3 && pattern[s][h3]==1)
	      ABCD[s][comp]=0;
	    else 
	      ABCD[s][comp] = matcat[pattern[s][h1]] [pattern[s][h2]] [pattern[s][h3]] [pars->anc[s]];
	  }  
	}
	else
	  for(int s=0;s<pars->numSites;s++){
	    if(pars->keepSites[s]==0)
	      continue;
  	    ABCD[s][comp] = matcat[pattern[s][h1]] [pattern[s][h2]] [pattern[s][h3]] [pars->anc[s]];
	    //  	    ABCD[s][comp] = matcat[2] [3] [2] [3];
	    //	    if(ABCD[s][comp]!=1)
	    //	      fprintf(stdout,"%d \n",ABCD[s][comp]);
	  }
	comp++;
      }
    }
  }


  for(int s=0;s<pars->numSites;s++)
    delete[] pattern[s];
  delete[] pattern;

  abbababaStruct->ABCD=ABCD;
  pars->extras[index] = abbababaStruct;
}

