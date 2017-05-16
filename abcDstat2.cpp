#include <cmath>
#include <cstdlib>
#include <ctime>
#include "shared.h"
#include "analysisFunction.h"
#include <htslib/kstring.h>
#include "abc.h"
#include "abcDstat2.h"

typedef struct {
  double **NUM;
  double **DEN;
  double ***COMB;
  int *BLOCKNUM;
  int **NSITE;
  int NUMBLOCK;
}funkyAbbababa2;

// shows up when you run ./angsd -doAbbababa2
void abcDstat2::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa2\t\t%d\trun the abbababa analysis\n",doAbbababa2);
  fprintf(argFile,"\t-rmTrans\t\t%d\tremove transitions\n",rmTrans);
  fprintf(argFile,"\t-blockSize\t\t%d\tsize of each block in bases\n",blockSize);
  fprintf(argFile,"\t-anc\t\t\t%s\tfasta file with outgroup\n",ancName);
  fprintf(argFile,"\t-sample\t\t\t%d\tsample a single base\n",sample);
  fprintf(argFile,"\t-maxDepth\t\t%d\tmax depth of each site allowed\n",maxDepth);
  fprintf(argFile,"\t-sizeFile\t\t%s\tfile with size of populations\n",sizeFile);
  fprintf(argFile,"\t-enhance\t\t%d\tonly analyze sites where outgroup H4 is non poly\n",enhance);
  fprintf(argFile,"\t-Aanc\t\t\t%d\tset H4 outgroup allele as A in each site\n",Aanc);
  fprintf(argFile,"\t-useLast\t\t%d\tuse last group of bam files as outgroup in the D-stat\n",useLast);
  fprintf(argFile,"\n");
}

// get you arguments
void abcDstat2::getOptions(argStruct *arguments){
  //from command line
  // 0: ignore this class, non zero: run this class
  doAbbababa2 = angsd::getArg("-doAbbababa2",doAbbababa2,arguments);

  if(doAbbababa2==0)
    return;

  doCount = angsd::getArg("-doCounts",doCount,arguments);
  blockSize = angsd::getArg("-blockSize",blockSize,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);
  rmTrans = angsd::getArg("-rmTrans",rmTrans,arguments);
  sample = angsd::getArg("-sample",sample,arguments);
  maxDepth = angsd::getArg("-maxDepth",maxDepth,arguments);
  enhance = angsd::getArg("-enhance",enhance,arguments);
  Aanc = angsd::getArg("-Aanc",Aanc,arguments);
  sizeFile = angsd::getArg("-sizeFile",sizeFile,arguments);
  useLast = angsd::getArg("-useLast",useLast,arguments);

  if(doAbbababa2){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doAbbababa2 \n");
      exit(0);
    }
    if(arguments->nInd==3 && ancName==NULL){
      fprintf(stderr,"Error: -doAbbababa2 needs at least 4 individual\n");
      exit(0);
    }
    if(doCount==0){
      fprintf(stderr,"Error: -doAbbababa2 needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
  }
} //---end of abcDstat2::getOptions(argStruct *arguments)

// Construction
abcDstat2::abcDstat2(const char *outfiles, argStruct *arguments,int inputtype){
  //default values
  rmTrans = 0;
  sample = 0;
  outfile=NULL;
  ancName=NULL;
  sizeFile=NULL;
  doAbbababa2=0;
  doCount=0;
  blockSize=5000000;
  block=0;
  enhance = 0;
  maxDepth=1000;
  useLast = 0;

 
  // you are starting before chromosome 0
  currentChr=-1;
  NbasesPerLine=50;
  Aanc=0;

  //if you dont use the -doAbbababa2 argument then return
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doAbbababa2")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  //get all options
  getOptions(arguments);

  //ignore
  if(doAbbababa2==0){
    shouldRun[index] = 0;
    return;
  }
  printArg(arguments->argumentFile);

  //correctly counting number of individuals
  nIndFasta = arguments -> nInd;
  //int numPop;
  if(sizeFile != NULL){
    sizeMat = angsd::getMatrixInt(sizeFile,nIndFasta);
    numPop = sizeMat.x;
  }
  else{
    numPop = nIndFasta;
  }
  
  if(ancName != NULL && useLast == 0){
    nIndFasta += 1;
    if(sizeFile == NULL)
      numPop = nIndFasta;
  }
  

  
  /*------ get population sizes and indexes for combinations of populations ------*/
  /*------------------------------------------------------------------------------*/
  
  // populations sizes from text file
  if( (sizeFile != NULL && ancName == NULL) || (sizeFile != NULL && ancName != NULL && useLast == 1) ){
    POPSIZE = new int[numPop];
    for(int a=0; a<numPop; a++){
      POPSIZE[a] = sizeMat.matrix[a][0];
    }
  }
  else if( (sizeFile != NULL && ancName != NULL && useLast == 0) ){
    POPSIZE = new int[numPop+1];
    for(int a=0; a<numPop; a++)
      POPSIZE[a] = sizeMat.matrix[a][0];
    POPSIZE[numPop] = 1;     
  }
  else if(sizeFile != NULL && useLast == 1){
    POPSIZE = new int[numPop];
    for(int a=0; a<numPop-1; a++)
      POPSIZE[a] = sizeMat.matrix[a][0];
    POPSIZE[numPop-1] = 1; 
  }
  else{ //default with no sizeFile: every individual is a population
    POPSIZE = new int[numPop];
    for(int i=0; i<numPop; i++)
      POPSIZE[i] = 1;
  }

  //check if the sizes of populations sum to the number of data files
  int sumCheck=0;
  for(int i=0; i<numPop; i++)
    sumCheck += POPSIZE[i];
  if(sumCheck != nIndFasta){
    fprintf(stderr,"\t-> Error: Num of individuals in file \"%s\" is %d and different from num %d of data files\n",sizeFile,sumCheck,nIndFasta);
    exit(0);
  }

  // cumulative populatios sizes (to get indexes for reading data in ABCD)
  CUMPOPSIZE = new int[numPop+1];
  CUMPOPSIZE[0] = 0;
  int cumCont = 0;
  for(int a=0; a<numPop; a++){
    cumCont += POPSIZE[a];
    CUMPOPSIZE[a+1] = cumCont;
  }

  //combinations of populations
  numComb =  (numPop-1)*(numPop-2)*(numPop-3)/2;
 
  //print some information
  if((useLast==1) || (useLast==0 && ancName==NULL))
    fprintf(stderr,"\t-> %d Populations | %lu trees | %d Individuals\n", numPop, numComb, nIndFasta);
  else if(useLast==0 && ancName!=NULL) 
    fprintf(stderr,"\t-> %d Populations | %lu trees | %d Individuals | %s Outgroup\n", numPop, numComb, nIndFasta, ancName);
  
  SIZEABCD = new int*[numComb]; //indexes for reading ABCD2 data  
  int cont=0;
    for(int i=0; i<numPop-1; i++){
      for(int j=0; j<numPop-1; j++){
	for(int k=0; k<numPop-1; k++){
	  if(j>i && j!=k && i!=k){
	    SIZEABCD[cont] = new int[3];   
	    SIZEABCD[cont][0] = 4*i;
	    SIZEABCD[cont][1] = 4*j;
	    SIZEABCD[cont][2] = 4*k;
	    cont += 1;
	  }
	}
      }
    }
  
  /*ENDEND get population sizes and indexes for combinations of populations ------*/
  /*------------------------------------------------------------------------------*/

  //make output files
  const char* postfix;
  postfix=".abbababa2";
  outfile = aio::openFile(outfiles,postfix);

  // store large amounts of data for printing fast
  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;
  
  //alleles pattern variable for printing
  COMBprint = new double*[numComb];
  for(int m=0; m<numComb; m++){
    COMBprint[m] = new double[256];
    for(int i=0; i<256; i++)
      COMBprint[m][i] = 0;
  }


  //numerator and denominator for printing
  NUMprint = new double[numComb];
  DENprint = new double[numComb];
  NSITEprint = new double[numComb];
  for(int m=0; m<numComb; m++){
    NUMprint[m]=0;
    DENprint[m]=0;
    NSITEprint[m]=0;
  }

  
  //print header
  fprintf(outfile,"CHR\tBLOCKstart\tBLOCKend\tNumer\tDenom\tnumSites");
  for(int a=0;a<4;a++)
    for(int b=0;b<4;b++)
      for(int c=0;c<4;c++)
	for(int d=0;d<4;d++)
	  fprintf(outfile,"\t%d%d%d%d",a,b,c,d);
  fprintf(outfile,"\n");

}//---end of abcDstat2::abcDstat2(const char *outfiles, argStruct *arguments,int inputtype)

//destructor
abcDstat2::~abcDstat2(){
  free(ancName);
  if(doAbbababa2==0)
    return;
  
  printAndEmpty(block*blockSize+1,currentChr);
  
  if(outfile) fclose(outfile);
  if(bufstr.s!=NULL)
    free(bufstr.s);
  if(sizeFile != NULL)
    angsd::deleteMatrixInt(sizeMat);
  
  for(int m=0; m<numComb; m++)
    delete[] SIZEABCD[m];
  delete[] SIZEABCD;
 
  delete[] POPSIZE;
  delete[] CUMPOPSIZE;
  delete[] DENprint;
  delete[] NUMprint;
  delete[] NSITEprint;

  for(int m=0; m<numComb; m++)
    delete[] COMBprint[m];
  delete[] COMBprint;
}

void abcDstat2::getBlockNum(int pos){
  block=(int)((pos+1)/blockSize);
}


// returns number of blocks in chunk
int abcDstat2::getNumBlocks(funkyPars *pars){

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




//clean after yourself
void abcDstat2::clean(funkyPars *pars){
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct =(funkyAbbababa2 *) pars->extras[index];

  int totBlocks = abbababaStruct->NUMBLOCK;
  for(int m=0; m<numComb; m++)
    delete[] abbababaStruct->NUM[m];
  delete[] abbababaStruct->NUM;

  for(int m=0; m<numComb; m++)
    delete[] abbababaStruct->DEN[m];
  delete[] abbababaStruct->DEN;

  for(int m=0; m<numComb; m++){
    for(int j=0; j<totBlocks; j++)
      delete[] abbababaStruct->COMB[m][j];
    delete[] abbababaStruct->COMB[m];
  }

  for(int m=0; m<numComb; m++)
    delete[] abbababaStruct->NSITE[m];
  delete[] abbababaStruct->NSITE;
  
  delete[] abbababaStruct->COMB;
  delete[] abbababaStruct->BLOCKNUM;

  delete abbababaStruct;
}//---end of abcDstat2::clean(funkyPars *pars)
//print and eventually reset counters

void abcDstat2::printAndEmpty(int blockAddress,int theChr){
  double denCont=0;
  for(int m=0; m<numComb; m++)
    denCont += DENprint[m];
  if(denCont != 0){
    for(int m=0; m<numComb; m++){
      fprintf(outfile,"%s\t%d\t%d\t%f\t%f\t%f",header->target_name[theChr],blockAddress-1,blockAddress+blockSize-2,NUMprint[m],DENprint[m],NSITEprint[m]);
      for(int i=0;i<256;i++)
	fprintf(outfile,"\t%f",COMBprint[m][i]);
      fprintf(outfile,"\n");
      }
    }
  
  for(int m=0; m<numComb; m++){
    DENprint[m]=0;
    NUMprint[m]=0;
    NSITEprint[m]=0;
  }

  for(int m=0; m<numComb; m++)
    for(int i=0; i<256; i++)
      COMBprint[m][i]=0;

  
  fflush(outfile);
  
}//---end of abcDstat2::printAndEmpty(int blockStart,int theChr)

void abcDstat2::print(funkyPars *pars){
  
  if(doAbbababa2==0)
    return;
  funkyAbbababa2 *abbababaStruct = (funkyAbbababa2 *) pars->extras[index];//new
  
  if(currentChr==-1){//if first chunk
    for(int m=0; m<numComb; m++){
      DENprint[m] = 0; //    
      NUMprint[m] = 0; //numerator for current block
      NSITEprint[m] = 0;
    }
    for(int m=0; m<numComb; m++){
      for(int i=0;i<256;i++)
    	COMBprint[m][i]=0;
    }
    getBlockNum(pars->posi[0]);
    if(currentChr > pars->refId)
      fprintf(stdout,"Warning. Your regions are not sorted. Becarefull...\n");
    currentChr=pars->refId;
    //start new block
    
    
    //block = abbababaStruct->BLOCKNUM[0];
  }

  while(currentChr!=pars->refId){ //if new chr (not first)
    //start new block
    printAndEmpty(block*blockSize+1,currentChr);
    currentChr=pars->refId;
    block = abbababaStruct->BLOCKNUM[0];
    //block = (int)((pars->posi[s]+1)/blockSize);
  }
  
  for(int b=0;b<abbababaStruct->NUMBLOCK;b++){
    
    
    if(abbababaStruct->BLOCKNUM[b] > block){
      printAndEmpty(block*blockSize+1,pars->refId);
      block = abbababaStruct->BLOCKNUM[b];
    }

    for(int c=0;c<numComb;c++)
      for(int pat=0;pat<256;pat++)
	COMBprint[c][pat] += abbababaStruct->COMB[c][b][pat];

    for(int c=0;c<numComb;c++){
      NUMprint[c] += abbababaStruct->NUM[c][b];
      DENprint[c] += abbababaStruct->DEN[c][b];
      NSITEprint[c] += abbababaStruct->NSITE[c][b];
    }


  }//---end of for(int s=0;s<pars->numSites;s++)
  
}//---end of abcDstat2::print(funkyPars *pars)


void abcDstat2::run(funkyPars *pars){
  if(doAbbababa2==0)
    return;

  //number of blocks
  int totBlocks = 0;
  int blockHere = -1;
  for(int s = 0 ; s < pars->numSites; s++){
    if( pars->keepSites[s]==0 )
      continue;
    
    if( pars->posi[s]+1 >= blockHere*blockSize+blockSize ){
      totBlocks++;
      blockHere =  (int)((pars->posi[s]+1)/blockSize);
    }
  } 
 
  funkyAbbababa2 *abbababaStruct = new funkyAbbababa2; //new structure
  //fprintf(stderr,"nindfasta %d\n",nIndFasta);
  //abbababaStruct->NUMBLOCK = totBlocks;
  abbababaStruct->NUMBLOCK = totBlocks;
  //double ABCD[4*nIndFasta];
  //for(int i = 0; i < 4*nIndFasta; i++)
  //  ABCD[i] = 0;
  
  //double ABCD2[numPop*4];
  //for(int i=0;i<numPop*4;i++)
  //  ABCD2[i*4] = 0;

  double ***COMB;
  COMB = new double**[numComb];
  for(int m=0; m<numComb; m++){
    COMB[m] = new double*[totBlocks];
    for(int i=0; i<totBlocks; i++){
      COMB[m][i] = new double[256];
      for(int j=0; j<256; j++)
	COMB[m][i][j] = 0;
    }
  }
    
  double **NUM;
  NUM = new double*[numComb];
  for(int m=0; m<numComb; m++){
    NUM[m] = new double[totBlocks];
    for(int i=0; i<totBlocks; i++)
      NUM[m][i] = 0;
  }


  double **DEN;
  DEN= new double*[numComb];
  for(int m=0; m<numComb; m++){
    DEN[m] = new double[totBlocks];
    for(int i=0; i<totBlocks; i++)
      DEN[m][i] = 0;
  }

  int *BLOCKNUM;
  BLOCKNUM = new int[totBlocks];
  for(int i=0; i<totBlocks; i++)
    BLOCKNUM[i] = 0;
  
  int **NSITE;
  NSITE = new int*[numComb];
  for(int m=0; m<numComb; m++){
    NSITE[m] = new int[totBlocks];
    for(int i=0; i<totBlocks; i++)
      NSITE[m][i] = 0;
  }
  

  //normalizing constants
  double somma;
  double normc;

 
  int blockIdx = -1;
  blockHere = -1;
  
  for(int s=0;s<pars->numSites;s++){
    //fprintf(stderr,"nindFasta2 %d\n",nIndFasta);
    double ABCD[4*nIndFasta];
    double ABCD2[4*numPop];
    for(int i = 0; i < 4*nIndFasta; i++)
      ABCD[i] = 0;
    for(int i=0; i<numPop*4; i++)
      ABCD2[i] = 0;
    
    if(pars->keepSites[s]==0)
      continue;
    
    if( pars->posi[s]+1 >= blockHere*blockSize + blockSize ){
      blockIdx++;
      blockHere =  (int)((pars->posi[s]+1)/blockSize);
      BLOCKNUM[blockIdx] = blockHere;
    }
    for(int i=0;i<pars->nInd;i++){
      //read the data
      if(pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3] == 0)
	continue;//no data
      if(pars->counts[s][i*4]<maxDepth && pars->counts[s][i*4+1]<maxDepth && pars->counts[s][i*4+2]<maxDepth && pars->counts[s][i*4+3]<maxDepth){
	
	if(sample==1){
	  int dep=0;
	  for( int b = 0; b < 4; b++ )
	    dep+=pars->counts[s][i*4+b];	  	  
	  int j;
	  j = std::rand()%dep;
	  int cumSum=0;
	  for( int b = 0; b < 4; b++ ){
	    cumSum+=pars->counts[s][i*4+b];
	    if(cumSum > j){
	      ABCD[i*4+b] = 1;
	      break;
	    }
	  }	    
	}
	else{
	  for( int b = 0; b < 4; b++ )//{ //bases
	    ABCD[i*4+b]=pars->counts[s][i*4+b];
	}
      }
    }
    
    if(ancName != NULL && useLast == 0){
      if(pars->anc[s] < 4)
	ABCD[pars->nInd * 4 + pars->anc[s]] = 1;
      else
	continue;
    }
    
    //------do all the populationss------
    double w[nIndFasta];//weights
    for(int p=0; p<numPop; p++){
      
      //---------------building weighted individual 1. written in ABCD2.     
      somma = 0;
      normc = 0;
      for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++){
	w[i]=0;
	somma = ABCD[i*4+0]+ABCD[i*4+1]+ABCD[i*4+2]+ABCD[i*4+3]; 
	w[i] = (2*somma)/(somma + 1);
	normc += w[i];
      }
      if(normc!=0){
	for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++)
	  w[i] = w[i]/normc;
      }
      for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++)
	//---------------building ABCD2 - weighted (pseudo) individuals
	for(int al=0;al<4;al++){
	  for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++){
	    ABCD2[p*4 + al] += w[i]*ABCD[i*4+al];
	  }
	}
      //normalize
      somma = 0;
      for(int al=0;al<4;al++)
	somma += ABCD2[p*4+al];
      if(somma!=0){
	for(int al=0;al<4;al++)
	  ABCD2[p*4+al] /= somma;}     
    }
    
    /*-ENDENDEND--count WEIGHTED normalized allele combinations -----------------------*/
    /*-------------------------------------------------------------------------------- */
        
    if(enhance==1){
      int enh=0;
      for(int j=0;j<4;j++)
	if(ABCD2[ (numPop-1)*4 + j]==0)
	  enh++;
      if(enh!=3)
	continue;	  
    }
	  
    for(int m=0; m<numComb; m++){

      double abba=0;
      double baba=0;
      double h1=0;
      double h12=0;
      double h123=0;
      double h4=0;
      double h1234=0;
      int pattern=0;
      double siteCont = 0;
      
      for(int i=0;i<4;i++){
	h1 = ABCD2[ SIZEABCD[m][0] + i];
	if(h1==0){
	  pattern+=64; continue;}
	for(int j=0;j<4;j++){
	  h12 = h1 * ABCD2[ SIZEABCD[m][1] + j];
	  if(h12==0){
	    pattern+=16; continue;}
	  for(int k=0;k<4;k++){
	    h123 = h12 * ABCD2[ SIZEABCD[m][2] + k];
	    if(h123==0){
	      pattern+=4; continue;}
	    for(int l=0;l<4;l++){
	      if( rmTrans==1 && (pattern==40 || pattern==130 || pattern==34 || pattern==136 ||pattern==125 || pattern==215 || pattern==119 || pattern==221) ){
		pattern++; continue;}
	      
	      h4=ABCD2[ (numPop-1)*4 + l];
	      if(h4==0){
		pattern++;
	      }
	      else{
		h1234 = h123 * h4;
		
		COMB[m][blockIdx][pattern] += h1234;
		siteCont += h1234;
		
		if(i==l && j==k && i!=j)
		  abba += h1234;
		
		
		if(i==k && j==l && i!=j)
		  baba += h1234;
		
		pattern++;
	      }
	    }
	  }
	}
      }
      NUM[m][blockIdx] += abba - baba;
      DEN[m][blockIdx] += abba + baba;
      NSITE[m][blockIdx] += siteCont;
    }
    

    
  }//---end for(int s=0;s<pars->numSites;s++)
  
  abbababaStruct -> NUM=NUM;
  abbababaStruct -> DEN=DEN;
  abbababaStruct -> COMB=COMB;
  abbababaStruct -> BLOCKNUM = BLOCKNUM;
  abbababaStruct -> NSITE = NSITE;
  pars -> extras[index] = abbababaStruct;      
}
  



