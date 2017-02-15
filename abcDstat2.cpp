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
  double **NUM;
  double **DEN;
  double **COMB;  //alleles combinations in weighted individuals
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
  else{ //default with no sizeFile
    numPop = 4;
    POPSIZE = new int[4];
    POPSIZE[0]=1;
    POPSIZE[1]=1;
    POPSIZE[2]=1;
    POPSIZE[3]=1;
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
  numComb = 3;
  if(numPop > 4){
    numComb =  (numPop-1)*(numPop-2)*(numPop-3)/2;
  }

  //print some information
  if(useLast==0)
    fprintf(stderr,"\t-> %d Populations | %d trees | %d Individuals\n", numPop, numComb, nIndFasta);
  else if(useLast==1 && ancName!=NULL) 
    fprintf(stderr,"\t-> %d Populations | %d trees | %d Individuals | %s Outgroup\n", numPop, numComb, nIndFasta, ancName);
  
  SIZEIDX = new int*[numComb]; //populations' sizes for each possible combination
  SIZEABCD = new int*[numComb]; //indexes for reading ABCD data
  
  int cont=0;
    for(int i=0; i<numPop-1; i++){
      for(int j=0; j<numPop-1; j++){
	for(int k=0; k<numPop-1; k++){
	  if(j>i && j!=k && i!=k){
	    SIZEIDX[cont] = new int[3];   
	    SIZEIDX[cont][0] = POPSIZE[i];
	    SIZEIDX[cont][1] = POPSIZE[j];
	    SIZEIDX[cont][2] = POPSIZE[k];
	    SIZEABCD[cont] = new int[3];   
	    SIZEABCD[cont][0] = 4*CUMPOPSIZE[i];
	    SIZEABCD[cont][1] = 4*CUMPOPSIZE[j];
	    SIZEABCD[cont][2] = 4*CUMPOPSIZE[k];
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
    delete[] SIZEIDX[m];
  delete[] SIZEIDX;

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

  for(int m=0; m<numComb; m++){
    delete[] abbababaStruct->NUM[m];
    delete[] abbababaStruct->DEN[m];
  }

  for(int m=0; m<numComb*256; m++)
    delete[] abbababaStruct->COMB[m];
  
  delete[] abbababaStruct->NUM;
  delete[] abbababaStruct->DEN;
  delete[] abbababaStruct->COMB;

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
      NUMprint[m] += abbababaStruct->NUM[m][s];
      DENprint[m] += abbababaStruct->DEN[m][s];
      if(abbababaStruct->NUM[m][s]>abbababaStruct->DEN[m][s])
	Eprint += 1;
    }

    for(int m=0; m<numComb; m++){
      for(int i=0; i<256; i++)
	COMBprint[m][i] += abbababaStruct->COMB[m*256 + i][s];
    }
    NSITEprint++;
  }
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
  
  double **NUM;//numerator of D-statistic
  NUM = new double*[numComb];
  for(int m=0; m<numComb; m++){
    NUM[m] = new double[pars->numSites];
    for(int b = 0; b < pars->numSites; b++)
      NUM[m][b] = 0;   
  }

  double **DEN;//denominator of D-statistic
  DEN = new double*[numComb];
  for(int m=0; m<numComb; m++){
    DEN[m] = new double[pars->numSites];
    for(int b = 0; b < pars->numSites; b++)
      DEN[m][b] = 0;   
  }
  
  double **ALLELES;//allele patterns
  ALLELES = new double*[numComb*256];
  for(int m=0; m<numComb*256; m++){
    ALLELES[m] = new double[pars->numSites];{
      for(int s = 0; s < pars->numSites; s++)
	ALLELES[m][s] = 0;
    }
  }

  double ABCD2[numComb * 12]; //three inner populations weighted allele counts
  double ABCD2OUT[4]; //outgroup weighted allele counts

  for(int i=0;i<4;i++)
    ABCD2OUT[i]=0;

  //indexes for data reading and populations sizes
  int size1[numComb];
  int size2[numComb];
  int size3[numComb];
  int id1[numComb];
  int id2[numComb];
  int id3[numComb];
  
  for(int comb=0; comb<numComb; comb++){  
    size1[comb] = SIZEIDX[comb][0];
    size2[comb] = SIZEIDX[comb][1];
    size3[comb] = SIZEIDX[comb][2];
    id1[comb] = SIZEABCD[comb][0];
    id2[comb] = SIZEABCD[comb][1];
    id3[comb] = SIZEABCD[comb][2];
    }

  //normalizing constants
  double somma;
  double normc;
  
  if(doAbbababa2==1){
    
    for(int s=0;s<pars->numSites;s++){
    
      if(pars->keepSites[s]==0)
	continue;
      
      if(ancName != NULL && useLast == 1){ //if fasta file has no data
	if(pars->anc[s] == 4)
	  continue;
      }
      
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
	  ABCD[s][pars->nInd * 4 + pars->anc[s]] = 1;

      //------fix weights and alleles for the outgroup------     
      double w4[POPSIZE[numPop-1]];
      
      if(POPSIZE[numPop-1] == 1)
	w4[0] = 1;
      else{
	somma = 0;
	normc = 0;
	for(int i=0;i<POPSIZE[numPop-1];i++){
	  somma = ABCD[s][4*CUMPOPSIZE[numPop-1]+i*4]+ABCD[s][4*CUMPOPSIZE[numPop-1]+i*4+1]+ABCD[s][4*CUMPOPSIZE[numPop-1]+i*4+2]+ABCD[s][4*CUMPOPSIZE[numPop-1]+i*4+3]; 
	  w4[i] = (2*somma)/(somma + 1);
	  normc += w4[i];
	}
	if(normc!=0){
	  for(int i=0;i<POPSIZE[numPop-1];i++)
	    w4[i] = w4[i]/normc;
	}
      }
      
      somma=0;
      
      for(int j=0;j<4;j++){
	for(int i=0; i<POPSIZE[numPop-1]; i++){
	  ABCD2OUT[j] += w4[i]*ABCD[s][4*CUMPOPSIZE[numPop-1]+i*4+j];
	}
	somma += ABCD2OUT[j];
      }

      if(somma!=0){
	for(int j=0;j<4;j++)
	  ABCD2OUT[j] = ABCD2OUT[j]/somma;
      }
    
      //------do all the combinations------
      for(int comb=0; comb<numComb; comb++){
	double w1[size1[comb]];//
	double w2[size2[comb]];//
	double w3[size3[comb]];//
	double sum1[size1[comb]];//
	double sum2[size2[comb]];//
	double sum3[size3[comb]];//
	
      //---------------building weighted individual 1. written in ABCD2.     
      somma = 0;
      normc = 0;      
      for(int i=0;i<size1[comb];i++){
	w1[i]=0;
	somma = ABCD[s][id1[comb]+i*4+0]+ABCD[s][id1[comb]+i*4+1]+ABCD[s][id1[comb]+i*4+2]+ABCD[s][id1[comb]+i*4+3]; 
	w1[i] = (2*somma)/(somma + 1);
	normc += w1[i];
      }
      if(normc!=0){
	for(int i=0;i<size1[comb];i++)
	  w1[i] = w1[i]/normc;
      }
      //---------------building weighted individual 2. written in ABCD2.
      somma = 0;
      normc = 0;
      for(int i=0;i<size2[comb];i++){
	somma = ABCD[s][id2[comb]+i*4]+ABCD[s][id2[comb]+i*4+1]+ABCD[s][id2[comb]+i*4+2]+ABCD[s][id2[comb]+i*4+3]; 
	w2[i] = (2 * somma)/(somma + 1);
	normc += w2[i];
      }
      if(normc!=0){
	for(int i=0;i<size2[comb];i++)
	  w2[i] = w2[i]/normc;
      }
      //---------------building weighted individual 3. written in ABCD2.

      somma = 0;
      normc = 0;
      for(int i=0;i<size3[comb];i++){
	somma = ABCD[s][id3[comb]+i*4]+ABCD[s][id3[comb]+i*4+1]+ABCD[s][id3[comb]+i*4+2]+ABCD[s][id3[comb]+i*4+3]; 
	w3[i] = (2*somma)/(somma + 1);
	normc += w3[i];
      }
      if(normc!=0){
	for(int i=0;i<size3[comb];i++)
	  w3[i] = w3[i]/normc;
      }

      //---------------building ABCD2 - weighted (pseudo) 3 inner individuals
      for(int j=0;j<12;j++)
	ABCD2[12 * comb + j] = 0;

      for(int j=0;j<4;j++){// first pseudo individual
	for(int i=0; i<size1[comb]; i++){
	  ABCD2[12*comb + j] += w1[i]*ABCD[s][id1[comb]+i*4+j];
	}
	for(int i=0; i<size2[comb]; i++){
	  ABCD2[12*comb + 4 + j] += w2[i]*ABCD[s][id2[comb]+i*4+j];
	}
	for(int i=0; i<size3[comb]; i++){
	  ABCD2[12*comb + 8 + j] += w3[i]*ABCD[s][id3[comb]+i*4+j];
	}
      }
    
      //---------------count alleles combination of 4 pseudoindividials
      somma = 0;
      //normalizing observation vector
      for(int i=0;i<3;i++){
	somma = ABCD2[12*comb + i*4] + ABCD2[12*comb + i*4 + 1] + ABCD2[12*comb + i*4 + 2] + ABCD2[12*comb + i*4 + 3];
	if(somma != 0){
	  ABCD2[12*comb + i*4] = ABCD2[12*comb + i*4]  / somma;
	  ABCD2[12*comb + i*4 + 1] = ABCD2[12*comb + i*4 + 1]  / somma;
	  ABCD2[12*comb + i*4 + 2] = ABCD2[12*comb + i*4 + 2]  / somma;
	  ABCD2[12*comb + i*4 + 3] = ABCD2[12*comb + i*4 + 3]  / somma;
	}	 
      }

      /*-------------------------------------------------------------------------------- */
      /*------------count WEIGHTED normalized allele combinations,---------------------- */ 
      /*---------------------numerator and denominator ----------------------------------*/
      int posiz = 0;
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  for(int k=0;k<4;k++){
	    for(int l=0;l<4;l++){
	      ALLELES[256*comb + posiz][s] += ABCD2[comb*12 + i] * ABCD2[comb*12 + 4 + j] * ABCD2[comb*12 + 8 + k] * ABCD2OUT[l];
	      posiz++;
	    }
	  }
	}
      }

      NUM[comb][s] = NUM[comb][s] + ALLELES[256*comb+20][s] + ALLELES[256*comb+40][s] + ALLELES[256*comb+60][s] +ALLELES[256*comb+65][s] +ALLELES[256*comb+105][s] +ALLELES[256*comb+125][s] +ALLELES[256*comb+130][s] +ALLELES[256*comb+150][s] +ALLELES[256*comb+190][s] +ALLELES[256*comb+195][s] +ALLELES[256*comb+215][s] +ALLELES[256*comb+235][s] - ALLELES[256*comb+17][s] - ALLELES[256*comb+34][s] - ALLELES[256*comb+51][s] - ALLELES[256*comb+68][s] - ALLELES[256*comb+102][s] - ALLELES[256*comb+119][s] - ALLELES[256*comb+136][s] - ALLELES[256*comb+153][s] - ALLELES[256*comb+187][s] - ALLELES[256*comb+204][s] - ALLELES[256*comb+221][s] - ALLELES[256*comb+238][s];
      DEN[comb][s] = DEN[comb][s] + ALLELES[256*comb+20][s] + ALLELES[256*comb+40][s] + ALLELES[256*comb+60][s] +ALLELES[256*comb+65][s] +ALLELES[256*comb+105][s] +ALLELES[256*comb+125][s] +ALLELES[256*comb+130][s] +ALLELES[256*comb+150][s] +ALLELES[256*comb+190][s] +ALLELES[256*comb+195][s] +ALLELES[256*comb+215][s] +ALLELES[256*comb+235][s] + ALLELES[256*comb+17][s] + ALLELES[256*comb+34][s] + ALLELES[256*comb+51][s] + ALLELES[256*comb+68][s] + ALLELES[256*comb+102][s] + ALLELES[256*comb+119][s] + ALLELES[256*comb+136][s] + ALLELES[256*comb+153][s] + ALLELES[256*comb+187][s] + ALLELES[256*comb+204][s] + ALLELES[256*comb+221][s] + ALLELES[256*comb+238][s];
      
      /*-ENDENDEND--count WEIGHTED normalized allele combinations -----------------------*/
      /*-------------------------------------------------------------------------------- */


      /*---'-enhance' option for analyzing only non-polymorphic sites of the outgroup----*/ 
      /*-------------------------------------------------------------------------------- */    
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
      /*END-'-enhance' option for analyzing only non-polymorphic sites of the outgroup---*/ 
      /*-------------------------------------------------------------------------------- */  

      }//---end for(combinations of populations)
    }//---end for(int s=0;s<pars->numSites;s++)
  }//---end if(doAbbababa2==1)

  abbababaStruct -> ABCD=ABCD;
  abbababaStruct -> DEN=DEN;
  abbababaStruct -> NUM=NUM;
  abbababaStruct -> COMB=ALLELES;

  pars -> extras[index] = abbababaStruct;
}


