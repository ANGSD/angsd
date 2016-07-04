
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
  double *COMB;  //alleles combinations in weighted individuals
  double *ALLCOMB;  //alleles combinations in unweighted individuals
 }funkyAbbababa2;

// shows up when you run ./angsd -doAbbababa2
void abcDstat2::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa2\t\t%d\trun the abbababa analysis\n",doAbbababa2);
  fprintf(argFile,"\t-rmTrans\t\t%d\tremove transitions\n",rmTrans);
  fprintf(argFile,"\t-blockSize\t\t%d\tsize of each block in bases\n",blockSize);
  fprintf(argFile,"\t-ans\t\t\t%s\tfasta file with outgroup\n",ancName);
  fprintf(argFile,"\t-sample\t\t\t%d\tsample a single base\n",sample);
  fprintf(argFile,"\t-maxDepth\t\t\t%d\tmax depth of each site allowed\n",maxDepth);
  fprintf(argFile,"\t-sizeH1\t\t\t%d\tnum of individuals in group H1",sizeH1);
  fprintf(argFile,"\t-sizeH2\t\t\t%d\tnum of individuals in group H2",sizeH2);
  fprintf(argFile,"\t-sizeH3\t\t\t%d\tnum of individuals in group H3",sizeH3);
  fprintf(argFile,"\t-sizeH4\t\t\t%d\tnum of individuals in group H4",sizeH4);
  fprintf(argFile,"\t-enhance\t\t\t%d\tonly analyze sites where outgroup H4 is non poly.",enhance);
  fprintf(argFile,"\t-Aanc\t\t\t%d\tset H4 outgroup allele as A in each site",Aanc);
  fprintf(argFile,"\t-combFile\t\t\t%d\tcreate a .abbababa2counts file where are printed the numbers of alleles combinations without having weighted the individuals",combFile);
  fprintf(argFile,"\n");
}

// get you arguments
void abcDstat2::getOptions(argStruct *arguments){
  //from command line
  // 0: ignore this class, non zero: run this class
  doAbbababa2 = angsd::getArg("-doAbbababa2",doAbbababa2,arguments); 
  doCount = angsd::getArg("-doCounts",doCount,arguments);
  blockSize = angsd::getArg("-blockSize",blockSize,arguments);
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
  combFile = angsd::getArg("-combFile",combFile,arguments);
  
  if(doAbbababa2){
    if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
      fprintf(stderr,"Error: bam or soap input needed for -doAbbababa2 \n");
      exit(0);
    }
    if(arguments->nInd==3 && ancName==NULL){
      fprintf(stderr,"Error: -doAbbababa2 needs 4 individual\n");
      exit(0);
    }
    if(doCount==0){
      fprintf(stderr,"Error: -doAbbababa2 needs allele counts (use -doCounts 1)\n");
      exit(0);
    }
    if((arguments->nInd != sizeH1 + sizeH2 + sizeH3 + sizeH4) && ancName == NULL){
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
  combFile = 0;
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

  //update number of individuals if using fasta outgroup
  //with the option -anc ancName
  nIndFasta = arguments->nInd;
  if(ancName != NULL)
    nIndFasta++;

  //make output files
  const char* postfix;
  postfix=".abbababa2";
  outfile = aio::openFile(outfiles,postfix);
  
  
  //if(combFile == 1){// unweighted combinations count file
  const char* postfix2; 
  postfix2=".abbababa2counts"; 
  outfile2 = aio::openFile(outfiles,postfix2); //}

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
  ALLCOMBprint = new double[256];
  
  fprintf(outfile,"CHR\tBLOCKstart\tBLOCKend\tNumer\tDenom\tnumSites");
  for(int a=0;a<4;a++)
    for(int b=0;b<4;b++)
      for(int c=0;c<4;c++)
	for(int d=0;d<4;d++)
	  fprintf(outfile,"\t%d%d%d%d",a,b,c,d);
  fprintf(outfile,"\n");

  //if(combFile == 1){
    fprintf(outfile2,"CHR\tBLOCKstart\tBLOCKend\tnumSites");
    for(int a=0;a<4;a++)
      for(int b=0;b<4;b++)
	for(int c=0;c<4;c++)
	  for(int d=0;d<4;d++)
	    fprintf(outfile2,"\t%d%d%d%d",a,b,c,d);
    fprintf(outfile2,"\n");
    //}
}//---end of abcDstat2::abcDstat2(const char *outfiles, argStruct *arguments,int inputtype)


abcDstat2::~abcDstat2(){
  free(ancName);
  if(doAbbababa2==0)
    return;
  
  printAndEmpty(block*blockSize,currentChr);
  
  if(outfile) fclose(outfile);
  if(outfile2) fclose(outfile2);
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
  delete[] abbababaStruct->ALLCOMB;
  for(int s=0;s<pars->numSites;s++)
    delete[] abbababaStruct->ABCD[s];
  delete[] abbababaStruct->ABCD;
  delete abbababaStruct;
}//---end of abcDstat2::clean(funkyPars *pars)

void abcDstat2::printAndEmpty(int blockStart,int theChr){

  if(NUMprint!=0){ //avoid to print 0s in the first line of the output file
    fprintf(outfile,"%s\t%d\t%d\t%f\t%f\t%d",header->target_name[theChr],blockStart,blockStart+blockSize,NUMprint,DENprint,NSITEprint);
    fprintf(outfile2,"%s\t%d\t%d\t%d",header->target_name[theChr],blockStart,blockStart+blockSize,NSITEprint);
    for(int i=0;i<256;i++){
      fprintf(outfile,"\t%f",COMBprint[i]);
      //if(combFile==1)
	fprintf(outfile2,"\t%f",ALLCOMBprint[i]);
	}
    fprintf(outfile,"\n");
    //if(combFile == 1)
    fprintf(outfile2,"\n");
      //}
  }

  DENprint=0;
  NUMprint=0;
  NSITEprint = 0;
  Eprint = 0; //overflow errors (for debugging - not printed anymore in the final version) 
  for(int i=0;i<256;i++)
    COMBprint[i]=0;
  fflush(outfile);
  //if(combFile == 1){
    for(int i=0;i<256;i++)
      ALLCOMBprint[i]=0;
    fflush(outfile2);
    //}
  
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
    //if(combFile == 1){
    for(int i=0;i<256;i++)
      ALLCOMBprint[i]=0;
      //}
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
    COMBprint[i] += abbababaStruct->COMB[i]; 
  //if(combFile == 1){
  for(int i=0;i<256;i++)
    ALLCOMBprint[i] += abbababaStruct->ALLCOMB[i];
    //}
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
  if(doAbbababa2==0)
    return;

  funkyAbbababa2 *abbababaStruct = new funkyAbbababa2; //new structure

  double **ABCD; //pointer to nSites counts of allele
  ABCD = new double*[pars->numSites]; 

  for(int s=0;s<pars->numSites;s++){
    ABCD[s] = new double[4*nIndFasta];
    for(int b = 0; b < 4*nIndFasta; b++)
      ABCD[s][b] = 0;   
  }

  
  double *ALLCOMB = new double[256];
  double *NUM = new double[pars->numSites];
  double *DEN = new double[pars->numSites];
  double *ALLELES = new double[256];//----observation of four alleles *1*2*3*4* counter
  double *ABCD2 = new double[16];  //----weighted site of the 4 populations
  double *w1 = new double[sizeH1]; // Weights
  double *w2 = new double[sizeH2]; // of
  double *w3 = new double[sizeH3]; // the 4
  double *w4 = new double[sizeH4]; // populations
  double *sum1 = new double[sizeH1];
  double *sum2= new double[sizeH2];
  double *sum3 = new double[sizeH3];
  double *sum4 = new double[sizeH4];
      
  double somma;
  double normc;
  int contNeg;

  for(int i=0;i<256;i++)
    ALLELES[i] = 0;
  for(int i=0;i<16;i++)
    ABCD2[i] = 0;
  for(int j=0;j<256;j++)
    ALLCOMB[j] = 0;
    
  if(doAbbababa2==1){
    
    for(int s=0;s<pars->numSites;s++){
      contNeg = 0;
      if(pars->keepSites[s]==0)
	continue;

      for(int i=0;i<nIndFasta;i++){ 
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
      /*
      fprintf(stderr,"Weigths ");
     for(int i=0;i<sizeH1;i++) 
       fprintf(stderr,"%.2e ",w1[i]);
     for(int i=0;i<sizeH2;i++)
       fprintf(stderr,"%.2e ",w2[i]);
     for(int i=0;i<sizeH3;i++)
       fprintf(stderr,"%.2e ",w3[i]);
     for(int i=0;i<sizeH4;i++)
       fprintf(stderr,"%.2e ",w4[i]);
       fprintf(stderr,"\n");*/

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
      somma = 0;
      //normalizing observation vector
      for(int i=0;i<4;i++){
	somma = ABCD2[i*4]+ABCD2[i*4+1]+ABCD2[i*4+2]+ABCD2[i*4+3];
	if(somma != 0){
	  /*fprintf(stderr,"Define Weighted ABCD alleles\n ");
	  for(int mm=0;mm<16;mm++)
	    fprintf(stderr,"%.2e ",ABCD2[mm]);
	    fprintf(stderr,"\n");*/
	  ABCD2[i*4] = ABCD2[i*4]  / somma;
	  ABCD2[i*4+1] = ABCD2[i*4+1]  / somma;
	  ABCD2[i*4+2] = ABCD2[i*4+2]  / somma;
	  ABCD2[i*4+3] = ABCD2[i*4+3]  / somma;
	}	 
      }

      /*-------------------------------------------------------------------------------- */
      /*---------------count normalized allele combinations without weighting individuals*/
      int posiz = 0;

      for(int h1=0;h1<sizeH1;h1++)
	sum1[h1] = ABCD[s][h1*4]+ABCD[s][h1*4+1]+ABCD[s][h1*4+2]+ABCD[s][h1*4+3];
      for(int h2=0;h2<sizeH2;h2++)
	sum2[h2] = ABCD[s][sizeH1*4+h2*4]+ABCD[s][sizeH1*4+h2*4+1]+ABCD[s][sizeH1*4+h2*4+2]+ABCD[s][sizeH1*4+h2*4+3];
      for(int h3=0;h3<sizeH3;h3++)
	sum3[h3] = ABCD[s][sizeH1*4+sizeH2*4+h3*4]+ABCD[s][sizeH1*4+sizeH2*4+h3*4+1]+ABCD[s][sizeH1*4+sizeH2*4+h3*4+2]+ABCD[s][sizeH1*4+sizeH2*4+h3*4+3];
      for(int h4=0;h4<sizeH4;h4++)
	sum4[h4] = ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+h4*4]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+h4*4+1]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+h4*4+2]+ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+h4*4+3];


      if(combFile == 1){
	double allele1 = 0;
	double allele2 = 0;
	double allele3 = 0;
	double allele4 = 0;
	for(int h1=0;h1<sizeH1;h1++){
	  posiz = 0;
	  for(int h2=0;h2<sizeH2;h2++){
	    posiz = 0;
	    for(int h3=0;h3<sizeH3;h3++){
	      posiz = 0;
	      for(int h4=0;h4<sizeH4;h4++){
		posiz = 0;
		for(int i=0;i<4;i++){
		  for(int j=0;j<4;j++){
		    for(int k=0;k<4;k++){
		      for(int l=0;l<4;l++){
			allele1 = ABCD[s][h1*4+i];
			allele2 = ABCD[s][sizeH1*4+h2*4+j];
			allele3 = ABCD[s][sizeH1*4+sizeH2*4+h3*4+k];
			allele4 = ABCD[s][sizeH1*4+sizeH2*4+sizeH3*4+h4*4+l];
			if(sum1[h1]*sum2[h2]*sum3[h3]*sum4[h4] != 0)
			  ALLCOMB[posiz] += (allele1/sum1[h1]) * (allele2/sum2[h2]) * (allele3/sum3[h3]) * (allele4/sum4[h4]);
			allele1 = 0;
			allele2 = 0;
			allele3 = 0;
			allele4 = 0;
			posiz++;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      /*-ENDENDEND-----count normalized allele combinations without weighting individuals*/
      /*-------------------------------------------------------------------------------- */





      /*-------------------------------------------------------------------------------- */
      /*------------count WEIGHTED normalized allele combinations -----------------------*/
      posiz = 0;
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
      /*-ENDENDEND--count WEIGHTED normalized allele combinations -----------------------*/
      /*-------------------------------------------------------------------------------- */


      /*-------------------------------------------------------------------------------- */
      /*------------numerator and denominator for the D-statistic -----------------------*/
      /*-------HERE IS WHERE I PRINT OUT THINGS ORIGINATING THE PROBLEM------------------*/
      DEN[s]=0;
      NUM[s]=0;
      double addNum;
      double addDen;

      for(int i=0;i<4;i++){
	for(int j=0;j<4;j++){
	  if(i!=j){
	    addNum = ABCD2[i]*ABCD2[4+j]*ABCD2[8+j]*ABCD2[12+i] - ABCD2[j]*ABCD2[4+i]*ABCD2[8+j]*ABCD2[12+i];
	    addDen = ABCD2[i]*ABCD2[4+j]*ABCD2[8+j]*ABCD2[12+i] + ABCD2[j]*ABCD2[4+i]*ABCD2[8+j]*ABCD2[12+i];
	    //check if there are problems in the numerator
	    if(addNum > 1 || addNum < -1){
	      fprintf(stderr,"error NUM>1\n%f\t\n",NUM[s]);
	      fprintf(stderr,"NON-weighted ABCD alleles\n ");
	      for(int mm=0;mm<16;mm++)
		fprintf(stderr,"%.2e ",ABCD[s][mm]);
	      fprintf(stderr,"\n");
	      fprintf(stderr,"weighted ABCD alleles\n ");
	      for(int mm=0;mm<16;mm++)
		fprintf(stderr,"%.2e ",ABCD2[mm]);
	      fprintf(stderr,"\n");
	    }
	    if(addDen > 1 || addDen < 0){
	      fprintf(stderr,"error DEN>1\n%f\t\n",DEN[s]);
	      fprintf(stderr,"NON-weighted ABCD alleles\n ");
	      for(int mm=0;mm<16;mm++)
		fprintf(stderr,"%.2e ",ABCD[s][mm]);
	      fprintf(stderr,"\n");
	      fprintf(stderr,"weighted ABCD alleles\n ");
	      for(int mm=0;mm<16;mm++)
		fprintf(stderr,"%.2e ",ABCD2[mm]);
	      fprintf(stderr,"\n");
	    }
	    NUM[s] += addNum;
	    DEN[s] += addDen;
	    //if(NUM[s] > 1 || NUM[s] < -1){fprintf(stderr,"error NUM>1\n%f\t\n",NUM[s]); continue;}
	    //if(DEN[s] > 1 || DEN[s] < 0){fprintf(stderr,"error DEN>1\n%f\t\n",DEN[s]); continue;}
	  }
	}
      }
      /*---ENDEND---numerator and denominator for the D-statistic -----------------------*/ 
      /*-------------------------------------------------------------------------------- */


      /*---'-enhance' option for analyzing only non-polymorphic sites of the outgroup----*/ 
      /*-------------------------------------------------------------------------------- */    
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
      /*END-'-enhance' option for analyzing only non-polymorphic sites of the outgroup---*/ 
      /*-------------------------------------------------------------------------------- */  


    }//---end for(int s=0;s<pars->numSites;s++)
  }//---end if(doAbbababa2==1)

  abbababaStruct -> ABCD=ABCD;
  abbababaStruct -> DEN=DEN;
  abbababaStruct -> NUM=NUM;
  abbababaStruct -> COMB=ALLELES;
  abbababaStruct -> ALLCOMB=ALLCOMB;
  delete[] w1;
  delete[] w2;
  delete[] w3;
  delete[] w4;
  delete[] ABCD2;
  delete[] sum1;
  delete[] sum2;
  delete[] sum3;
  delete[] sum4;
  //delete[] ALLCOMB;
  //delete[] ALLELES;
  pars -> extras[index] = abbababaStruct;
}


