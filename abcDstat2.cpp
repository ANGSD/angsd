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
  double *NUM;
  double *DEN;
  double *COMB;  //alleles combinations
 }funkyAbbababa2;

// shows up when you run ./angsd -doAbbababa2
void abcDstat2::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa2\t\t%d\trun the abbababa analysis\n",doAbbababa2);
  fprintf(argFile,"\t-rmTrans\t\t%d\tremove transitions\n",rmTrans);
  fprintf(argFile,"\t-blockSize\t\t%d\tsize of each block in bases\n",blockSize);
  fprintf(argFile,"\t-ans\t\t\t%s\tfasta file with outgroup\n",ancName);
  fprintf(argFile,"\t-sample\t\t\t%d\tsample a single base\n",sample);
  fprintf(argFile,"\t-maxDepth\t\t\t%d\tsample a single base\n",maxDepth);
  fprintf(argFile,"\t-sizeH1\t\t\t%d\tnum of individuals in group H1",sizeH1);
  fprintf(argFile,"\t-sizeH2\t\t\t%d\tnum of individuals in group H2",sizeH2);
  fprintf(argFile,"\t-sizeH3\t\t\t%d\tnum of individuals in group H3",sizeH3);
  fprintf(argFile,"\t-sizeH4\t\t\t%d\tnum of individuals in group H4",sizeH4);
  fprintf(argFile,"\t-enhance\t\t\t%d\tonly analyze sites where outgroup H4 is non poly. Put 2 if the order of populations is reverted.",enhance);
  fprintf(argFile,"\t-Aanc\t\t\t%d\tset H4 outgroup allele as A in each site",Aanc);
  fprintf(argFile,"\n");
}

// get you arguments
void abcDstat2::getOptions(argStruct *arguments){
  //from command line
  // 0: ignore this class, non zero: run this class
  doAbbababa2=angsd::getArg("-doAbbababa2",doAbbababa2,arguments); 
  doCount=angsd::getArg("-doCounts",doCount,arguments);
  blockSize=angsd::getArg("-blockSize",blockSize,arguments);
  ancName = angsd::getArg("-anc",ancName,arguments);
  rmTrans = angsd::getArg("-rmTrans",rmTrans,arguments);
  sample = angsd::getArg("-sample",sample,arguments);
  maxDepth = angsd::getArg("-maxDepth",maxDepth,arguments);
  sizeH1 = angsd::getArg("-sizeH1",sizeH1,arguments);
  sizeH2 = angsd::getArg("-sizeH2",sizeH2,arguments);
  sizeH3 = angsd::getArg("-sizeH3",sizeH3,arguments);
  sizeH4 = angsd::getArg("-sizeH4",sizeH4,arguments);
  enhance = angsd::getArg("-enhance",enhance,arguments);
  Aanc = angsd::getArg("-Aanc",Aanc,arguments);
  
  if(doAbbababa2){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doAbbababa2 \n");
      exit(0);
    }
    if(arguments->nInd==3){
      fprintf(stderr,"Error: -doAbbababa2 needs 4 individual\n");
      exit(0);
    }
    if(doCount==0){
      fprintf(stderr,"Error: -doAbbababa2 needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
    if(arguments->nInd != sizeH1 + sizeH2 + sizeH3 + sizeH4){
      fprintf(stderr,"Error: the number of individuals in sizeH* does not match the number of individuals in bam files\n");
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
  doAbbababa2=0;
  doCount=0;
  blockSize=5000000;
  block=0;
  NSITEprint=0;
  sizeH1 = 1;
  sizeH2 = 1;
  sizeH3 = 1;
  sizeH4 = 1;
  enhance = 1;
  DENprint=0;
  NUMprint=0;
  maxDepth=100;
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
  printArg(arguments->argumentFile);
  //ignore
  if(doAbbababa2==0){
    shouldRun[index] = 0;
    return;
  }

  //make output files
  const char* postfix;
  postfix=".abbababa2";
  outfile = aio::openFile(outfiles,postfix);

  //const char* postfix2; //new
  //postfix2=".abbababa"; //new
  //outfile2 = aio::openFile(outfiles,postfix2); //new

  // store large amounts of data for printing fast
  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;
  
  //the number of combinations
  //nComb=arguments->nInd*(arguments->nInd-1)*(arguments->nInd-2)/2; 
  //ABBACOMB = new int[nComb]; //new
  //BABACOMB = new int[nComb]; //new

  //calcMatCat(); //new //pre calc all possible combinations

  COMBprint = new double[256];
  
  fprintf(outfile,"CHR\tBLOCKstart\tBLOCKend\tNumer\tDenom\tnumSites");
  for(int a=0;a<4;a++)
    for(int b=0;b<4;b++)
      for(int c=0;c<4;c++)
	for(int d=0;d<4;d++)
	  fprintf(outfile,"\t%d%d%d%d",a,b,c,d);
  fprintf(outfile,"\n");
}//---end of abcDstat2::abcDstat2(const char *outfiles, argStruct *arguments,int inputtype)


abcDstat2::~abcDstat2(){
  free(ancName);
  if(doAbbababa2==0)
    return;

  printAndEmpty(block*blockSize,currentChr);

  if(outfile) fclose(outfile);
  if(bufstr.s!=NULL)
    free(bufstr.s);
}

void abcDstat2::clean(funkyPars *pars){
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct =(funkyAbbababa2 *) pars->extras[index];
  delete[] abbababaStruct->NUM;
  delete[] abbababaStruct->DEN;
  delete[] abbababaStruct->COMB;
  for(int s=0;s<pars->numSites;s++)
    delete[] abbababaStruct->ABCD[s];
  delete[] abbababaStruct->ABCD;
  delete abbababaStruct;
}//---end of abcDstat2::clean(funkyPars *pars)

void abcDstat2::printAndEmpty(int blockStart,int theChr){

  if(NUMprint!=0){ //avoid to print 0s in the first line of the output file
    fprintf(outfile,"%s\t%d\t%d\t%f\t%f\t%d",header->target_name[theChr],blockStart,blockStart+blockSize,NUMprint,DENprint,NSITEprint);
    for(int i=0;i<256;i++)
      fprintf(outfile,"\t%f",COMBprint[i]);
    fprintf(outfile,"\n");
  }

  DENprint=0;
  NUMprint=0;
  NSITEprint = 0;
  Eprint = 0; //overflow errors (for debugging - not printed anymore in the final version)
  
  for(int i=0;i<256;i++)
    COMBprint[i]=0;
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
    DENprint=0; //    
    NUMprint=0; //numerator for current block
    Eprint = 0;
    for(int i=0;i<256;i++)
      COMBprint[i]=0; 
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
 
  for(int i=0;i<256;i++)
    COMBprint[i] += abbababaStruct->COMB[i]; //?????
  
  for(int s=0;s<pars->numSites;s++){
    //int comp=0; //new
    if(pars->posi[s]>=block*blockSize+blockSize){
      printAndEmpty(block*blockSize,pars->refId);
      getBlockNum(pars->posi[s]);
    }    
    if(pars->keepSites[s]==0)
      continue;
    NUMprint+=abbababaStruct->NUM[s];
    DENprint+=abbababaStruct->DEN[s];
    NSITEprint++;
    if(abbababaStruct->NUM[s]>abbababaStruct->DEN[s])
      Eprint += 1;    
  }
}//---end of abcDstat2::print(funkyPars *pars)


void abcDstat2::run(funkyPars *pars){
  // pars->nInd
  // pars->numSites
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct = new funkyAbbababa2; //new structure

  double **ABCD; //pointer to nSites counts of allele
  ABCD = new double*[pars->numSites]; 

  for(int s=0;s<pars->numSites;s++){
    ABCD[s] = new double[4*pars->nInd];
    for(int b = 0; b < 4*pars->nInd; b++)
      ABCD[s][b] = 0;   
  }

  double *NUM = new double[pars->numSites];
  double *DEN = new double[pars->numSites];
  double *ALLELES = new double[256];//----observation of four alleles *1*2*3*4* counter
  double *ABCD2 = new double[16];  //----weighted site of the 4 populations
  double *w1 = new double[sizeH1]; // Weights
  double *w2 = new double[sizeH2]; // of
  double *w3 = new double[sizeH3]; // the 4
  double *w4 = new double[sizeH4]; // populations

  double somma;
  double normc;
  int contNeg;

  for(int i=0;i<256;i++)
    ALLELES[i] = 0;
  for(int i=0;i<16;i++)
    ABCD2[i] = 0;

  if(doAbbababa2==1){
    
    for(int s=0;s<pars->numSites;s++){
      contNeg = 0;
      if(pars->keepSites[s]==0)
	continue;

      for(int i=0;i<pars->nInd;i++){ 
	if(pars->counts[s][i*4]<0 || pars->counts[s][i*4+1]<0 || pars->counts[s][i*4+2]<0 || pars->counts[s][i*4+3]<0)
	  contNeg += 1;                //count negative occurrences in site s
	if(pars->counts[s][i*4] + pars->counts[s][i*4+1] + pars->counts[s][i*4+2] + pars->counts[s][i*4+3] == 0)
	  continue;
	if(pars->counts[s][i*4]<maxDepth && pars->counts[s][i*4+1]<maxDepth && pars->counts[s][i*4+2]<maxDepth && pars->counts[s][i*4+3]<maxDepth){
	 
	  if(sample==1){
	    int dep=0;
	    for( int b = 0; b < 4; b++ ){
	      dep+=pars->counts[s][i*4+b];
	    }	  
	    srand(time(0));
	    int j;
	    j = std::rand()%dep;
	    int cumSum=0;
	    for( int b = 0; b < 4; b++ ){
	      cumSum+=pars->counts[s][i*4+b];
	      if( cumSum > j){
		ABCD[s][i*4+b] = 1;
	        break;
	      }
	    }	    
	  }
	  else{
	    for( int b = 0; b < 4; b++ ) //bases
	      ABCD[s][i*4+b]=pars->counts[s][i*4+b];
	  }
	}
      }

      if(contNeg>0)//throw out the site if there are negative values in it
	continue;
      
      //---------------building weighted individual 1. written in ABCD2.
      
      somma = 0;
      normc = 0;
       
      for(int i=0;i<sizeH1;i++){
	w1[i]=0;
	somma = ABCD[s][i*4+0]+ABCD[s][i*4+1]+ABCD[s][i*4+2]+ABCD[s][i*4+3]; 
	w1[i] = (2*somma)/(somma + 1);
	normc += w1[i];
      }
      for(int i=0;i<sizeH1;i++){
	if(normc!=0)
	  w1[i] = w1[i]/normc;
      }
          
      //---------------building weighted individual 2. written in ABCD2.

      somma = 0;
      normc = 0;

      for(int i=0;i<sizeH2;i++){
	somma = ABCD[s][sizeH1*4+i*4]+ABCD[s][sizeH1*4+i*4+1]+ABCD[s][sizeH1*4+i*4+2]+ABCD[s][sizeH1*4+i*4+3]; 
	w2[i] = (2*somma)/(somma + 1);
	normc += w2[i];
      }
      for(int i=0;i<sizeH2;i++){
	if(normc!=0)
	  w2[i] = w2[i]/normc;
      }

      //---------------building weighted individual 3. written in ABCD2.

      somma = 0;
      normc = 0;

      for(int i=0;i<sizeH3;i++){
	somma = ABCD[s][sizeH1*4+sizeH2*4+i*4]+ABCD[s][sizeH1*4+sizeH2*4+i*4+1]+ABCD[s][sizeH1*4+sizeH2*4+i*4+2]+ABCD[s][sizeH1*4+sizeH2*4+i*4+3]; 
	w3[i] = (2*somma)/(somma + 1);
	normc += w3[i];
      }
      for(int i=0;i<sizeH3;i++){
	if(normc!=0)
	  w3[i] = w3[i]/normc;
      }

      //---------------building weighted individual 4. written in ABCD2.

      somma = 0;
      normc = 0;
 
      for(int i=0;i<sizeH4;i++){
	somma = ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+i*4]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+i*4+1]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+i*4+2]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+i*4+3]; 
	w4[i] = (2*somma)/(somma + 1);
	normc += w4[i];
      }
      for(int i=0;i<sizeH4;i++){
	if(normc!=0)
	  w4[i] = w4[i]/normc;
      }
      
      //---------------building ABCD2 - weighted (pseudo) 4 individuals

      for(int j=0;j<16;j++)
	ABCD2[j] = 0;

      for(int j=0;j<4;j++){// first pseudo individual
	for(int i=0; i<sizeH1; i++){
	  ABCD2[j] += w1[i]*ABCD[s][i*4+j];
	}
	for(int i=0; i<sizeH2; i++){
	  ABCD2[4+j] += w2[i]*ABCD[s][sizeH1*4+i*4+j];
	}
	for(int i=0; i<sizeH3; i++){
	  ABCD2[8+j] += w3[i]*ABCD[s][sizeH1*4+sizeH2*4+i*4+j];
	}
	for(int i=0; i<sizeH4; i++){
	  ABCD2[12+j] += w4[i]*ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+i*4+j];
	}
      }
      
      //---------------count alleles combination of 4 pseudoindividials


      int posiz = 0;
      
      somma = 0;

      //normalizing observation vector
      for(int i=0;i<4;i++){
	somma = ABCD2[i*4]+ABCD2[i*4+1]+ABCD2[i*4+2]+ABCD2[i*4+3];
	if(somma != 0){
	  ABCD2[i*4]=ABCD2[i*4]/somma;
	  ABCD2[i*4+1]=ABCD2[i*4+1]/somma;
	  ABCD2[i*4+2]=ABCD2[i*4+2]/somma;
	  ABCD2[i*4+3]=ABCD2[i*4+3]/somma;
	}	 
      }
   

      // allele counts
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  for(int k=0;k<4;k++){
	    for(int l=0;l<4;l++){
	      ALLELES[posiz] += ABCD2[i] * ABCD2[4+j] * ABCD2[8+k] * ABCD2[12+l];
	      posiz++;
	    }
	  }
	}
      }

      DEN[s]=0;
      NUM[s]=0;
      
      // numerator and denominator for the error-free estimation of the D-statistic
      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  if(i!=j){
	    NUM[s] += ABCD2[i]*ABCD2[4+j]*ABCD2[8+j]*ABCD2[12+i] - ABCD2[j]*ABCD2[4+i]*ABCD2[8+j]*ABCD2[12+i];
	    DEN[s] += ABCD2[i]*ABCD2[4+j]*ABCD2[8+j]*ABCD2[12+i] + ABCD2[j]*ABCD2[4+i]*ABCD2[8+j]*ABCD2[12+i];
	    if(NUM[s] > 1 || NUM[s] < -1){fprintf(stderr,"error NUM>1\n%f\t\n",NUM[s]);}
	    if(DEN[s] > 1 || DEN[s] < 0){fprintf(stderr,"error DEN>1\n%f\t\n",DEN[s]);}
	  }
	}
      }
      
      //'-enhance' option for analyzing only non-polymorphic sites of the outgroup
      if(enhance==1){
	int enh=0;
	for(int j=0;j<4;j++)
	  if(ABCD2[12+j]==0)
	    enh++;
	if(enh!=3){
	  DEN[s]=0;
	  NUM[s]=0;
	}	  
      }

    }//---end for(int s=0;s<pars->numSites;s++)
  }//---end if(doAbbababa2==1)

  abbababaStruct -> ABCD=ABCD;
  abbababaStruct -> DEN=DEN;
  abbababaStruct -> NUM=NUM;
  abbababaStruct -> COMB=ALLELES;
  delete[] w1;
  delete[] w2;
  delete[] w3;
  delete[] w4;
  delete[] ABCD2;
  pars -> extras[index] = abbababaStruct;
}

