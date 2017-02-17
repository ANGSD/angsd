#include <cmath>
#include <cstdlib>
#include <ctime>
#include "shared.h"
#include "analysisFunction.h"
#include <htslib/kstring.h>
#include "abc.h"
#include "abcDstat2.h"

typedef struct {
  double **ABCD; //counts
  double **ABCD2;
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
  fprintf(argFile,"\t-useLast\t\t%d\tuse fasta file as outgroup in the D-stat\n",useLast);
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
  NSITEprint=0;
  enhance = 0;
  maxDepth=100;
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
  
  if(ancName != NULL && useLast == 1)
    nIndFasta += 1;

  
  /*------ get population sizes and indexes for combinations of populations ------*/
  /*------------------------------------------------------------------------------*/
  
  // populations sizes from text file
  if( (sizeFile != NULL && ancName == NULL) || (sizeFile != NULL && ancName != NULL && useLast == 0) ){
    POPSIZE = new int[numPop];
    for(int a=0; a<numPop; a++)
      POPSIZE[a] = sizeMat.matrix[a][0];
  }
  else if(sizeFile != NULL && ancName != NULL && useLast == 1){
    numPop += 1;
    POPSIZE = new int[numPop];
    for(int a=0; a<numPop-1; a++)
      POPSIZE[a] = sizeMat.matrix[a][0];
    POPSIZE[numPop-1] = 1; 
  }
  else{ //default with no sizeFile: every individual is a population
    POPSIZE = new int[nIndFasta];
    for(int i=0; i<nIndFasta; i++)
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
  if(useLast==0)
    fprintf(stderr,"\t-> %d Populations | %d trees | %d Individuals\n", numPop, numComb, nIndFasta);
  else if(useLast==1 && ancName!=NULL) 
    fprintf(stderr,"\t-> %d Populations | %d trees | %d Individuals | %s Outgroup\n", numPop, numComb, nIndFasta, ancName);
  
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
  for(int m=0; m<numComb; m++){
    NUMprint[m] = 0;
    DENprint[m] = 0;
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
  
  printAndEmpty(block*blockSize,currentChr);
  
  if(outfile) fclose(outfile);
  if(bufstr.s!=NULL)
    free(bufstr.s);

  angsd::deleteMatrixInt(sizeMat);

  for(int m=0; m<numComb; m++)
    delete[] SIZEABCD[m];
  delete[] SIZEABCD;

  for(int m=0; m<numComb; m++)
    delete[] COMBprint[m];
  delete[] COMBprint;

  delete[] NUMprint;
  delete[] DENprint;
  delete[] POPSIZE;
  delete[] CUMPOPSIZE;
}

//clean after yourself
void abcDstat2::clean(funkyPars *pars){
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct =(funkyAbbababa2 *) pars->extras[index];
  
  for(int s=0;s<pars->numSites;s++)
    delete[] abbababaStruct->ABCD[s];
  delete[] abbababaStruct->ABCD;

  for(int s=0;s<pars->numSites;s++)
    delete[] abbababaStruct->ABCD2[s];
  delete[] abbababaStruct->ABCD2;

  delete abbababaStruct;
}//---end of abcDstat2::clean(funkyPars *pars)

//print and eventually reset counters
void abcDstat2::printAndEmpty(int blockStart,int theChr){

  for(int m=0; m<numComb; m++){
    if(NUMprint[m] != 0){ //avoid to print 0s in the first line of the output file
    fprintf(outfile,"%s\t%d\t%d\t%f\t%f\t%d",header->target_name[theChr],blockStart,blockStart+blockSize,NUMprint[m],DENprint[m],NSITEprint);
    for(int i=0;i<256;i++)
      fprintf(outfile,"\t%f",COMBprint[m][i]);
    fprintf(outfile,"\n");
    }
  }
  
  for(int m=0; m<numComb; m++){
    DENprint[m]=0;
    NUMprint[m]=0;
  }
  
  NSITEprint = 0; 
  for(int m=0; m<numComb; m++){
    for(int i=0; i<256; i++)
      COMBprint[m][i]=0;
  }
  fflush(outfile);
  
}//---end of abcDstat2::printAndEmpty(int blockStart,int theChr)


void abcDstat2::getBlockNum(int pos){
  block=(int)((pos+1)/blockSize);
}

void abcDstat2::print(funkyPars *pars){
  
  if(doAbbababa2==0)
    return;
  funkyAbbababa2 *abbababaStruct = (funkyAbbababa2 *) pars->extras[index];//new
  
  if(currentChr==-1){//if first chunk
    for(int m=0; m<numComb; m++){
      DENprint[m] = 0; //    
      NUMprint[m] = 0; //numerator for current block
    }
    for(int m=0; m<numComb; m++){
      for(int i=0;i<256;i++)
	COMBprint[m][i]=0;
    }
    //start new block
    getBlockNum(pars->posi[0]);
    currentChr=0;
  }

  while(currentChr!=pars->refId){ //if new chr (not first)
    //start new block
    printAndEmpty(block*blockSize,currentChr);
    currentChr=pars->refId;
    getBlockNum(pars->posi[0]);
  }
  
  for(int s=0;s<pars->numSites;s++){
    
    if(pars->posi[s]>=block*blockSize+blockSize){
      printAndEmpty(block*blockSize,pars->refId);
      getBlockNum(pars->posi[s]);
    }
    
    if(pars->keepSites[s]==0)
      continue;

    
    for(int m=0; m<numComb; m++){
      double abba=0;
      double baba=0;
      int pattern=0;
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  for(int k=0;k<4;k++){
	    for(int l=0;l<4;l++){
	      
	      if(i==l && j==k && i!=j)
	        abba += abbababaStruct->ABCD2[s][ SIZEABCD[m][0] + i] * abbababaStruct->ABCD2[s][ SIZEABCD[m][1] + j] * abbababaStruct->ABCD2[s][ SIZEABCD[m][2] + k] * abbababaStruct->ABCD2[s][ (numPop-1)*4 + l];
	      
	      if(i==k && j==l && i!=j)
	        baba += abbababaStruct->ABCD2[s][ SIZEABCD[m][0] + i] * abbababaStruct->ABCD2[s][ SIZEABCD[m][1] + j] * abbababaStruct->ABCD2[s][ SIZEABCD[m][2] + k] * abbababaStruct->ABCD2[s][ (numPop-1)*4 + l];
	      
	      COMBprint[m][pattern] += abbababaStruct->ABCD2[s][ SIZEABCD[m][0] + i] * abbababaStruct->ABCD2[s][ SIZEABCD[m][1] + j] * abbababaStruct->ABCD2[s][ SIZEABCD[m][2] + k] * abbababaStruct->ABCD2[s][ (numPop-1)*4 + l];
	      pattern++;
	    }
	  }
	}
      }
      //NUMprint[m] += COMBprint[m][20] + COMBprint[m][40] + COMBprint[m][60] +COMBprint[m][65] +COMBprint[m][105] +COMBprint[m][125] +COMBprint[m][130] +COMBprint[m][150] +COMBprint[m][190] +COMBprint[m][195] +COMBprint[m][215] +COMBprint[m][235] - COMBprint[m][17] - COMBprint[m][34] - COMBprint[m][51] - COMBprint[m][68] - COMBprint[m][102] - COMBprint[m][119] - COMBprint[m][136] - COMBprint[m][153] - COMBprint[m][187] - COMBprint[m][204] - COMBprint[m][221] - COMBprint[m][238];
      //DENprint[m] += COMBprint[m][20] + COMBprint[m][40] + COMBprint[m][60] +COMBprint[m][65] +COMBprint[m][105] +COMBprint[m][125] +COMBprint[m][130] +COMBprint[m][150] +COMBprint[m][190] +COMBprint[m][195] +COMBprint[m][215] +COMBprint[m][235] + COMBprint[m][17] + COMBprint[m][34] + COMBprint[m][51] + COMBprint[m][68] + COMBprint[m][102] + COMBprint[m][119] + COMBprint[m][136] + COMBprint[m][153] + COMBprint[m][187] + COMBprint[m][204] + COMBprint[m][221] + COMBprint[m][238];
      NUMprint[m] += abba - baba;
      DENprint[m] += abba + baba;
    }
    



    

    NSITEprint++;
  }//---end of for(int s=0;s<pars->numSites;s++)
  
}//---end of abcDstat2::print(funkyPars *pars)


void abcDstat2::run(funkyPars *pars){
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct = new funkyAbbababa2; //new structure

  //allele counts from the data .bam files
  double **ABCD; 
  ABCD = new double*[pars->numSites]; 

  for(int s=0;s<pars->numSites;s++){//set to zero
    ABCD[s] = new double[4*nIndFasta];
    for(int b = 0; b < 4*nIndFasta; b++)
      ABCD[s][b] = 0;   
  }
  
  double **ABCD2;
  ABCD2 = new double*[pars->numSites];
  for(int s = 0; s < pars->numSites; s++){
    ABCD2[s] = new double[numPop * 4];
    for(int i=0;i<numPop * 4;i++)
       ABCD2[s][i] = 0;
  }

  //normalizing constants
  double somma;
  double normc;
  
  if(doAbbababa2==1){
    
    for(int s=0;s<pars->numSites;s++){
    
      if(pars->keepSites[s]==0)
	continue;
      
      for(int i=0;i<pars->nInd;i++){
	//read the data
	if(pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3] == 0) //if no data at site s
	  continue;
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
		ABCD[s][i*4+b] = 1;
		break;
	      }
	    }	    
	  }
	  else{
	    for( int b = 0; b < 4; b++ )//{ //bases
	      ABCD[s][i*4+b]=pars->counts[s][i*4+b];
	  }
	}
      }

      if(ancName != NULL && useLast == 1)
	if(pars->anc[s] < 4)
	  ABCD[s][pars->nInd * 4 + pars->anc[s]] = 1;
      
      //------do all the populationss------
      double w[nIndFasta];//weights
      for(int p=0; p<numPop; p++){
      
      //---------------building weighted individual 1. written in ABCD2.     
      somma = 0;
      normc = 0;
      for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++){
	w[i]=0;
	somma = ABCD[s][i*4+0]+ABCD[s][i*4+1]+ABCD[s][i*4+2]+ABCD[s][i*4+3]; 
	w[i] = (2*somma)/(somma + 1);
	normc += w[i];
      }
      if(normc!=0){
	for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++)
	  w[i] = w[i]/normc;
      }
      
      //---------------building ABCD2 - weighted (pseudo) individuals
      for(int al=0;al<4;al++){
	for(int i=CUMPOPSIZE[p];i<CUMPOPSIZE[p+1];i++){
	  ABCD2[s][p*4 + al] += w[i]*ABCD[s][i*4+al];
	}
      }
      //normalize
      somma = 0;
      for(int al=0;al<4;al++)
	somma += ABCD2[s][p*4+al];
      if(somma!=0){
	for(int al=0;al<4;al++)
	  ABCD2[s][p*4+al] /= somma;}
      }

      
      /*-ENDENDEND--count WEIGHTED normalized allele combinations -----------------------*/
      /*-------------------------------------------------------------------------------- */


      /*---'-enhance' option for analyzing only non-polymorphic sites of the outgroup----*/ 
      /*-------------------------------------------------------------------------------- */
      /*
      if(enhance==1){
	int enh=0;
	for(int j=0;j<4;j++)
	  if(ABCD2OUT[j]==0)
	    enh++;
	if(enh!=3){
	  DEN[comb][s]=0;
	  NUM[comb][s]=0;
	}	  
      }
      */
      /*END-'-enhance' option for analyzing only non-polymorphic sites of the outgroup---*/ 
      /*-------------------------------------------------------------------------------- */  

      }//---end for(combinations of populations)

      abbababaStruct -> ABCD=ABCD;
      abbababaStruct -> ABCD2=ABCD2;
      pars -> extras[index] = abbababaStruct;
      
  }//---end for(int s=0;s<pars->numSites;s++)
    //}//---end if(doAbbababa2==1)

  //abbababaStruct -> ABCD=ABCD;
  //abbababaStruct -> DEN=DEN;
  //abbababaStruct -> NUM=NUM;
  //abbababaStruct -> INNERPOP=ABCD2;
  //abbababaStruct -> OUTPOP=ABCD2OUT;
  //abbababaStruct -> COMB=ALLELES;

  //pars -> extras[index] = abbababaStruct;
}


