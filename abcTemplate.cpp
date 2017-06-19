/*
  thorfinn 26june 2014 thorfinn@binf.ku.dk, part of angsd

  generally angsd will send the data for a genomic region through all analysis classes (abc*).
  the data is encapsulated in a 'funkyPars' struct (found in shared.h)
  
  For each abcClass it will do
  1) ::run
  2) ::print
  3) ::clean
  
  below are some toy example that will show
  a) how the class system works ()
  b) how to access the internal datastructure

  -doTemplate 1:
  Will simply count the number of A, C, G, T, N's for both strand

  We will put all analysis in the ::print function. therefore it will NOT be threaded
  the raw analysis should be put into ::run (this is threaded in angsd)
  but we wont do this in case of simplicity

  Run command examples:
  ---------------------
  ./angsd -i YanaRef.bam -doTemplate 2 -r 3:100000-50000000 -ref ../ProbabilisticAncientDNA/HumanGenomeReference/hs.build37.1.fa  

 */

#include <ctype.h> //<-used for isupper/islower
#include <cmath> //in order to construct large array/objects
#include <stdlib.h>

#include "shared.h" //<-contains the struct defintino for funkypars
#include "analysisFunction.h" //<-contains some utility functions, usefull for parsing args
#include "abcTemplate.h"//contains the analysis class defition associated with this file

//#include "phys_genolike_calc.h"


// Define pointer s.t. it can:
//phys_genolike_calc *like_calc;

#define ReadLengthMax 110
#define DepthMax 200
#define QscoreMax 61
#define MapScoreMax 61
#define Nallele 4
#define UndamagedDepth 15

// Basic distributions:
int Lread[ReadLengthMax];
int DepthDist[DepthMax];
long int QscoreDist[QscoreMax];
long int MapScoreDist[MapScoreMax];

// Higher dimensional distributions:
int BaseDistInRead[ReadLengthMax][Nallele][2][2];      // We consider the allele frequency as a function of distance from ends.
int QscoreDistInRead[ReadLengthMax][QscoreMax][2][2];  // We consider Qscores as a function of distance from ends.
long int QscoreVsErrorFreq[QscoreMax][Nallele][2];          // For each Qscore and base type, we consider the number of errors (1) vs. total number of bases (0).
int StrandDist[DepthMax][DepthMax];                    // Distribution (binomial?) of strands
long int ReadsVsRef[2];                                // Count the number of times a "sure" base read matches the reference (or not).
long int GenotypeFreq[11];                             // Frequency of genotypes (11: Undetermined from our requirements, i.e. "efficiency")
int CorWro[ReadLengthMax][2][2][QscoreMax][Nallele][Nallele];
int SigVsReadL[ReadLengthMax][ReadLengthMax][2][2][Nallele][Nallele];    // First "ReadLength" is the length of the read, the second is the actual position!
int HeteroZygDist[46][46];                             // The distribution of alleles in heterozygous positions (of depth 30-45).

// Cross check summations:
double AveEE[Nallele][Nallele][3];
double ActTT[Nallele][Nallele][3];
double ActEE[Nallele][Nallele][3];
double CheckQscore[Nallele][4];
double CheckPQ[Nallele][Nallele][3];


double min(double a, double b) {if (a<b) return a; else return b;}
double max(double a, double b) {if (a>b) return a; else return b;}



//this msg is shown on screen if you type
// ./angsd -doTemplate
//otherwise this is printet to the file .arg
void abcTemplate::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"-doTemplate\t%d (Which analysis should we perform?)\n",doTemplate);
  fprintf(argFile,"\t\t1: Count and print basetypes in combination with strand\n");
  if (doTemplate == 2) fprintf(argFile,"-ref\t%s\n",refName);
}

//this is the function that parses the parameters used for this analysis class
void abcTemplate::getOptions(argStruct *arguments){
  //from command line
  doTemplate=angsd::getArg("-doTemplate",doTemplate,arguments);

  
  if (doTemplate == 2) refName = angsd::getArg("-ref", refName, arguments);

  if ((doTemplate == 2) && (refName==NULL)) {
    fprintf(stderr, "\t-> Must supply -ref \n");
    printArg(stderr);
    exit(0);
  }

  if(doTemplate==0){
    /*
      if this class shouldnt do any analysis,
      then setting this to zero will make sure nothing is run (apart from destructor)
    
      this could also have been accomplished by
      if(doTemplate==0) return
      in ::run ::clean ::print
    */
    
    shouldRun[index]=0;
    return;
  }
}

//constructor
abcTemplate::abcTemplate(const char *outfiles,argStruct *arguments,int inputtype){
  doTemplate = 0; //defaults= dont do analysis
  refName = NULL;
  outfile = NULL;

  char *empty;
  //like_calc = new phys_genolike_calc( empty );

  //first a hook for the interactive help:
  //  ./angsd -doTemplate
  if (arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doTemplate")){
      printArg(stdout);
      exit(0);
    } else
      return;
  }
  
  //now parse the arguments
  getOptions(arguments);
  //now print the arguments

  if(doTemplate==0)
    return ;

  printArg(arguments->argumentFile);
  //initalize outputfile
  outfile = aio::openFile(outfiles,".results");
  // fprintf(outfile,"Chromo\tPosition\t+A\t+C\t+G\t+T\t-A\t-C\t-G\t-T\n");


  // Initialize arrays:
  for (int i=0; i<ReadLengthMax; i++) Lread[i] = 0;
  for (int i=0; i<DepthMax; i++) DepthDist[i] = 0;
  for (int i=0; i<QscoreMax; i++) QscoreDist[i] = 0;
  for (int i=0; i<MapScoreMax; i++) MapScoreDist[i] = 0;

  for (int i=0; i<ReadLengthMax; i++) {
    for (int j=0; i<Nallele; i++) {
      BaseDistInRead[i][j][0][0] = 0;
      BaseDistInRead[i][j][1][0] = 0;
      BaseDistInRead[i][j][0][1] = 0;
      BaseDistInRead[i][j][1][1] = 0;
    }
  }
  for (int i=0; i<ReadLengthMax; i++) {
    for (int j=0; i<QscoreMax; i++) {
      QscoreDistInRead[i][j][0][0] = 0;
      QscoreDistInRead[i][j][1][0] = 0;
      QscoreDistInRead[i][j][0][1] = 0;
      QscoreDistInRead[i][j][1][1] = 0;
    }
  }
  for (int i=0; i<QscoreMax; i++) {
    for (int j=0; i<Nallele; i++) {
      QscoreVsErrorFreq[i][j][0] = 0;
      QscoreVsErrorFreq[i][j][1] = 0;
    }
  }

  for (int i=0; i<DepthMax; i++)
    for (int j=0; j<DepthMax; j++)
      StrandDist[i][j] = 0;

  ReadsVsRef[0] = 0;
  ReadsVsRef[1] = 0;
  for (int i=0; i<11; i++) {
    GenotypeFreq[i] = 0;
  }

  for (int i=0; i<ReadLengthMax; i++) {
    for (int j=0; j<QscoreMax; j++) {
      for (int k=0; k<4; k++){
	for (int l=0;l<4; l++) {
	  CorWro[i][0][0][j][k][l] = 0;
	  CorWro[i][0][1][j][k][l] = 0;
	  CorWro[i][1][0][j][k][l] = 0;
	  CorWro[i][1][1][j][k][l] = 0;
        }
      }
    }
  }

  for (int i=0; i<ReadLengthMax; i++) {
    for (int j=0; j<ReadLengthMax; j++) {
      for (int k=0; k<4; k++){
	for (int l=0;l<4; l++) {
	  SigVsReadL[i][j][0][0][k][l] = 0;
	  SigVsReadL[i][j][0][1][k][l] = 0;
	  SigVsReadL[i][j][1][0][k][l] = 0;
	  SigVsReadL[i][j][1][1][k][l] = 0;
        }
      }
    }
  }

  for (int i=0; i<46; i++) {
    for (int j=0; j<46; j++) {
      HeteroZygDist[i][j] = 0;
    }
  }

  for (int i=0; i<Nallele; i++) {
    for (int j=0; j<Nallele; j++) {
      for (int k=0; k<3; k++){
	AveEE[i][j][k] = 0.0;
	ActTT[i][j][k] = 0.0;
	ActEE[i][j][k] = 0.0;
      }
    }
  }

  // Check of Q and PQ scores:
  for (int i=0; i<Nallele; i++) {
    for (int k=0; k<4; k++){
      CheckQscore[i][k] = 0.0;
    }
    for (int j=0; j<Nallele; j++) {
      for (int k=0; k<3; k++){
	CheckPQ[i][j][k] = 0.0;
      }
    }
  }

}



//destructor
abcTemplate::~abcTemplate(){
  if(doTemplate==0)
    return;
  fprintf(outfile,"\n\nRead length distribution:\n");
  for (int i=0; i<ReadLengthMax; i++) {
    fprintf(outfile,"\t%d \t%d \n", i, Lread[i]);
  }

  fprintf(outfile,"\n\nDepth distribution:\n");
  for (int i=0; i<DepthMax; i++) {
    fprintf(outfile,"\t%d \t%d \n", i, DepthDist[i]);
  }

  fprintf(outfile,"\n\nQscore distribution:\n");
  for (int i=0; i<QscoreMax; i++) {
    fprintf(outfile,"\t%d \t%ld \n", i, QscoreDist[i]);
  }

  fprintf(outfile,"\n\nMapScore distribution:\n");
  for (int i=0; i<MapScoreMax; i++) {
    fprintf(outfile,"\t%d \t%ld \n", i, MapScoreDist[i]);
  }

  fprintf(outfile,"\n\nAllele Frequency vs. Position in Read distribution:\n");
  fprintf(outfile,"  From start of Read:\n");
  for (int i=0; i<ReadLengthMax; i++) {
    fprintf(outfile,"\t  Posi: %2d    \tA: %d %d   \tC: %d %d   \tG: %d %d   \tT: %d %d \n", i,
	    BaseDistInRead[i][0][0][0], BaseDistInRead[i][0][0][1], BaseDistInRead[i][1][0][0], BaseDistInRead[i][1][0][1],
	    BaseDistInRead[i][2][0][0], BaseDistInRead[i][2][0][1], BaseDistInRead[i][3][0][0], BaseDistInRead[i][3][0][1]);
  }
  fprintf(outfile,"  From end of Read:\n");
  for (int i=0; i<ReadLengthMax; i++) {
    fprintf(outfile,"\t  Posi: %2d    \tA: %d %d   \tC: %d %d   \tG: %d %d   \tT: %d %d \n", i,
	    BaseDistInRead[i][0][1][0], BaseDistInRead[i][0][1][1], BaseDistInRead[i][1][1][0], BaseDistInRead[i][1][1][1],
	    BaseDistInRead[i][2][1][0], BaseDistInRead[i][2][1][1], BaseDistInRead[i][3][1][0], BaseDistInRead[i][3][1][1]);
  }

  // Q-scores vs. Read position:
  fprintf(outfile,"\n\nQscore Frequency vs. Position in Read distribution:\n");
  fprintf(outfile,"  From start/end of Read:\n");
  for (int i=0; i<ReadLengthMax; i++) {
    for (int j=0; j<QscoreMax; j++) {
      fprintf(outfile,"\t  Position from end: %3d   Q-score: %2d     \tStart (U/L strand): %d %d  \tEnd (U/L strand): %d %d \n", i, j,
	      QscoreDistInRead[i][j][0][0], QscoreDistInRead[i][j][0][1], QscoreDistInRead[i][j][1][0], QscoreDistInRead[i][j][1][1]);
    }
  }


  // Q-scores vs. Error frequencies:
  fprintf(outfile,"\n\nQscore vs. error frequency:\n");
  for (int i=0; i<QscoreMax; i++) {

    // Only print if there is significant data:
    int Nsignif = 100;
    if (QscoreVsErrorFreq[i][0][0] > Nsignif && QscoreVsErrorFreq[i][1][0] > Nsignif &&
	QscoreVsErrorFreq[i][2][0] > Nsignif && QscoreVsErrorFreq[i][3][0] > Nsignif) {
      double ErrFreq[4], eErrFreq[4];
      for (int j=0; j<4; j++) {
	ErrFreq[j]  = double(QscoreVsErrorFreq[i][j][1]) / double(QscoreVsErrorFreq[i][j][0]);
	eErrFreq[j] = sqrt(ErrFreq[j] * (1.0-ErrFreq[j]) / double(QscoreVsErrorFreq[i][j][0]));
      }
      fprintf(outfile,"\t  Qscore: %2d    \tA: %7.5f (%5d/%7d)   \tC: %7.5f (%5d/%7d)   \tG: %7.5f (%5d/%7d)   \tT: %7.5f (%5d/%7d) \n", i,
	      ErrFreq[0],(int) QscoreVsErrorFreq[i][0][1],(int) QscoreVsErrorFreq[i][0][0],
	      ErrFreq[1],(int) QscoreVsErrorFreq[i][1][1],(int) QscoreVsErrorFreq[i][1][0],
	      ErrFreq[2],(int) QscoreVsErrorFreq[i][2][1],(int) QscoreVsErrorFreq[i][2][0],
	      ErrFreq[3],(int) QscoreVsErrorFreq[i][3][1],(int) QscoreVsErrorFreq[i][3][0]);
    }
  }

  // Strand distribution:
  fprintf(outfile,"\n\nStrand distribution:\n");
  for (int i=1; i<DepthMax; i++) {
    fprintf(outfile,"  Depth: %3d   ", i);
    for (int j=0; j<i+1; j++) fprintf(outfile,"%d ", StrandDist[i][j]);
    fprintf(outfile,"\n");
  }

  // Number of times, that the read combination matches the reference base:
  fprintf(outfile,"\n\n Read correctness:\n");
  double fracWrong = double(ReadsVsRef[1]) / double(ReadsVsRef[0]);
  fprintf(outfile,"  Nwrong: %d   Ntotal: %d    frac = %8.6f \n",(int) ReadsVsRef[1],(int) ReadsVsRef[0], fracWrong);

  // Number of times, that the read combination matches the reference base:
  fprintf(outfile,"\n\n Genotype frequencies (and undetermined sites):\n");
  char* GTnames[10] = {(char*)"AA",(char*) "AC",(char*) "AG",(char*) "AT",(char*) "CC",(char*) "CG",(char*) "CT",(char*) "GG",(char*) "GT",(char*) "TT"};
  for (int i=0; i<10; i++) {
    fprintf(outfile, "  Genotype %s:  %8d \n", GTnames[i], (int)GenotypeFreq[i]);
  }
  fprintf(outfile, "  Genotype not determined:  %8d \n",(int) GenotypeFreq[10]);



  // The 85 * 2 * 2 * 40 * 4 * 4 = 217600 table of values for:    p_base = f_correct(posi/isop, q-score, base)
  fprintf(outfile,"\n\n TransitionsBaseQscore (i.e. matrix of error rates from which model is obtained):\n");
  for (int k=0; k<2; k++) {
    for (int l=0; l<2; l++) {
      for (int i=0; i<ReadLengthMax; i++) {
	for (int j=0; j<QscoreMax; j++) {
	  fprintf(outfile,"\t%d \t%d \t%d \t%d \t %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d \n", k, l, i , j,
		  CorWro[i][k][l][j][0][0], CorWro[i][k][l][j][0][1], CorWro[i][k][l][j][0][2], CorWro[i][k][l][j][0][3],
		  CorWro[i][k][l][j][1][0], CorWro[i][k][l][j][1][1], CorWro[i][k][l][j][1][2], CorWro[i][k][l][j][1][3],
		  CorWro[i][k][l][j][2][0], CorWro[i][k][l][j][2][1], CorWro[i][k][l][j][2][2], CorWro[i][k][l][j][2][3],
		  CorWro[i][k][l][j][3][0], CorWro[i][k][l][j][3][1], CorWro[i][k][l][j][3][2], CorWro[i][k][l][j][3][3]);
	}
      }
    }
  }


  // The 85 * 85 * 2 * 2 * 4 * 4 = 400.000 values for detecting aDNA signal strength vs. ReadLength:
  fprintf(outfile,"\n\n SignalStrenghVsReadLength:\n");
  for (int i=0; i<ReadLengthMax; i++) {
    for (int k=0; k<2; k++) {
      for (int l=0; l<2; l++) {
	for (int j=0; j<ReadLengthMax; j++) {
	  fprintf(outfile,"\t%d \t%d \t%d \t%d \t %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d %8d \n", i, k, l, j,
		  SigVsReadL[i][j][k][l][0][0], SigVsReadL[i][j][k][l][0][1], SigVsReadL[i][j][k][l][0][2], SigVsReadL[i][j][k][l][0][3],
		  SigVsReadL[i][j][k][l][1][0], SigVsReadL[i][j][k][l][1][1], SigVsReadL[i][j][k][l][1][2], SigVsReadL[i][j][k][l][1][3],
		  SigVsReadL[i][j][k][l][2][0], SigVsReadL[i][j][k][l][2][1], SigVsReadL[i][j][k][l][2][2], SigVsReadL[i][j][k][l][2][3],
		  SigVsReadL[i][j][k][l][3][0], SigVsReadL[i][j][k][l][3][1], SigVsReadL[i][j][k][l][3][2], SigVsReadL[i][j][k][l][3][3]);
	}
      }
    }
  }

  // Check of allele distribution in heterozygous cases:
  fprintf(outfile,"\n\n HeterozygousAlleleDistribution:\n");
  for (int i=0; i<46; i++) {
    fprintf(outfile, "  Depth: %2d   ", i);
    for (int j=0; j<46; j++) fprintf(outfile, " %d", HeteroZygDist[i][j]);
    fprintf(outfile, "\n");
  }


  // Simple check of Qscores:
  fprintf(outfile,"\n\n Quick check if probabilities from corrected Qscores are reasonable:\n");
  for (int i=0; i<4; i++) {
    fprintf(outfile, "  Base %1d:     SumProbQscore = %7.3f    SumProbQscoreCorr = %7.3f    NbaseWrong = %7.1f     NbasesChecked = %10.1f \n",
	    i, CheckQscore[i][0], CheckQscore[i][1], CheckQscore[i][2], CheckQscore[i][3]);
  }


  // Larger check of Qscores and transition probabilities together:
  fprintf(outfile,"\n\n Check of Q-scores + transition probability model validity:\n");
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      fprintf(outfile, "  True %1d:   Obs %1d:     SpQscore = %12.3f     SpModel = %12.3f     Nobs = %17.1f \n",
	      i, j, CheckPQ[i][j][0], CheckPQ[i][j][1], CheckPQ[i][j][2]);
    }
    fprintf(outfile, "  True %1d:                       Ntotal = %17.1f \n",
	    i, CheckPQ[i][0][2]+CheckPQ[i][1][2]+CheckPQ[i][2][2]+CheckPQ[i][3][2]);
  }


  //delete like_calc;

  if(outfile!=NULL) 
    fclose(outfile);
}



//this function is run, after ::run and ::print
void abcTemplate::clean(funkyPars *pars){
  //we havent done any allocation so we dont need to cleanup
}



void abcTemplate::print(funkyPars *pars){
  
  if(doTemplate==1){
    //count bases by strand

    //rawseqdata is in chunkyT struct (bambi_interface.h)
    chunkyT *chk = pars->chk;

    //loop over sites;
    for(int s=0;s<pars->numSites;s++){
      int bases[2][5] = {{0,0,0,0,0},{0,0,0,0,0}};      
      
      //loop over samples
      for(int i=0;i<pars->nInd;i++){
	//all seqdata associated with single bamfile is in a tNode
	tNode *nd = chk->nd[s][i];
	//loop over the individual bases
	for(int l=0;l<nd->l;l++){
	  char c = nd->seq[l]; //this is the base
	  char q = nd->qs[l]; //this is the associated qscore, fancy shit
	  int strand = isupper(nd->seq[l])==0; //strand is defined as either small/big letters
	  
	  //there is a lookuptable called refToInt which maps
	  //a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->5
	  bases[strand][refToInt[c]]++;
	}

	//print chr and position
	fprintf(outfile,"%s\t%d",header->target_name[pars->refId],pars->posi[s]+1);//position is zero index internally
	
	//print the basecount
	for(int i=0;i<2;i++)
	  for(int j=0;j<5;j++)
	    fprintf(outfile,"\t%d",bases[i][j]);
	fprintf(outfile,"\n");
	

      }
      
    }
  }
  
}



// -------------------------------------------------------------------------------------------------- //
void abcTemplate::run(funkyPars *pars){
// -------------------------------------------------------------------------------------------------- //

  if (doTemplate==2) {
    int count = 0;
    chunkyT *chk = pars->chk;
    
    // Point the likelihood calculator to the correct chunkyT
    //like_calc->update_chunkyT( chk );
    
    for(int s=0;s<pars->numSites;s++){
      for(int i=0;i<pars->nInd;i++){
	tNode *nd = chk->nd[s][i];

	// Reference base:
	int refB = -1;
	if (pars->ref == 0) fprintf(outfile,"  Warning: Ref not defined! %s \n", pars->ref);
	else refB = refToInt[pars->ref[s]];
	
	
	// Basic counting for distributions:
	// ---------------------------------
	// Depth:
	int Depth = min(DepthMax, max(0, nd->l));
	DepthDist[Depth] += 1;

	// Loop over the individual bases to count these:
	int Nbases[5] = {0, 0, 0, 0, 0};
	for (int l=0; l<nd->l; l++) {
	  char c = nd->seq[l];          // This is the base (a char)
	  int base = refToInt[c];      // This is the base (an int)
	  if (base > -1 && base < 5) Nbases[base]++;
	}


	// Check if position is of high depth, and if it can be used for control analysis:
	// -------------------------------------------------------------------------------
	int correct_base = -1;               // 0: AA, 1: CC, 2: GG, 3: TT
	int correct_genotype = -1;           // 0: AA, 1: AC, 2: AG, 3: AT, 4: CC, 5: CG, 6: CT, 7: GG, 8: GT, 9: TT
	int b2gt_homo[4] = {0, 4, 7, 9};
	int b2gt_hetero[4][4] = {{0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
	bool UsePosCorr = false;

	// Decide that base is correct, if 30 < depth < 45 and 90% are of one type:
	if (nd->l >= 30 && nd->l <= 45) {
	  for (int i=0; i<4; i++) {
	    if (float(Nbases[i]) > 0.9*float(nd->l))
	      correct_base = i;   // This could be done smarter (i.e. taking specific bases into account!!!)
	  }
	  
	  if (correct_base > -1 && refB > -1 && correct_base == refB)
	    UsePosCorr = true;

	  // Low-interlectual-budget code for GenoType:
	  int major = -1;
	  int Nmajor = 0;
	  int minor = -1;
	  int Nminor = 0;
	  for (int i=0; i<4; i++) {
	    if (Nbases[i] > Nmajor) {
	      Nminor = Nmajor;
	      minor = major;
	      Nmajor = Nbases[i];
	      major = i;
	    } else if (Nbases[i] > Nminor) {
	      Nminor = Nbases[i];
	      minor = i;
	    }
	  }

	  if (major > -1 && float(Nmajor) > 0.9*float(nd->l)) {
	    correct_genotype = b2gt_homo[major];
	    GenotypeFreq[correct_genotype]++;
	  } else if (major > -1 && minor > -1 && float(Nmajor+Nminor) > 0.9*float(nd->l) && Nminor > 6) {
	    correct_genotype = b2gt_hetero[major][minor];
	    GenotypeFreq[correct_genotype]++;
	    HeteroZygDist[nd->l][Nmajor]++;
	  } else {
	    GenotypeFreq[10]++;
	  }

	  // Smarter way for GenoType (for smarter guys than us!):
	  // std::vector<std::pair<K,V>> site_bases;
	  // std::sort(site_bases.begin(), site_bases.end(), value_comparer);   ...something something!
	}




	// ------------------------------------------------------------------------- //
	// Loop over the individual bases:
	// ------------------------------------------------------------------------- //
	int UpperStrands = 0;	// Count number of upper strands.
	//like_calc->update_tNode(nd);
	
	for (int l=0; l<nd->l; l++) {
	  char c = nd->seq[l];          // This is the base (a char)
	  int base = refToInt[c];      // This is the base (an int)
	  int qscore = nd->qs[l];       // This is the associated qscore (NOTE: Implicit type casting!!!)
	  int strand = (isupper(nd->seq[l]) == 0) ? 0 : 1;   // Strand is defined as either small/big letters
	  int mscore = nd->mapQ[l];              // Map score - to be understood in the same way as a Phred score (Yana: cut at 33)
	  int posi_here = nd->posi[l];           // NOTE: Implicit type casting!!!
	  int isop_here = nd->isop[l];           // NOTE: Implicit type casting!!!
	  
	  if (qscore < 1) continue;
	  if (base < 0 || base > 3) continue;
	  
	  // Read Length:
	  int ReadLength = min(ReadLengthMax, max(0, nd->posi[l]+nd->isop[l]) + 1);
	  if (nd->posi[l] == 0) {
	    Lread[ReadLength]++;
	  }

	  // Map scores:
	  int MapScore = min(MapScoreMax, max(0, mscore));
	  MapScoreDist[MapScore]++;

	  // Q scores:
	  int Qscore = min(QscoreMax, max(0, qscore));
	  QscoreDist[Qscore]++;
	  
	  // TODO: Get this to have length 61 (to include modern samples with Q-scores up to 60):
	  double qscore_corr[42] = { 1.000, 1.000, 1.000 , 1.000 , 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 ,
				     1.000, 1.000, 1.000 , 1.821 , 0.999, 0.923, 0.927, 7.469, 8.233, 12.997,
				     1.966, 8.325, 13.066, 13.422, 5.006, 1.653, 4.479, 1.827, 2.802, 1.266 ,
				     14.635, 1.927, 12.169, 2.192, 2.229, 2.891, 2.757, 2.328, 4.187, 7.902 , 5.391, 4.820};
	  double prob_Qscore      = pow( 10, -0.1 * qscore );
	  double prob_Qscore_corr = prob_Qscore * qscore_corr[qscore] / 1.0;

	  // Q-score vs. distance from read ends:
	  QscoreDistInRead[posi_here][qscore][0][strand]++;
	  QscoreDistInRead[isop_here][qscore][1][strand]++;

	  // Allele/base frequency vs. distance from read ends:
	  if (isop_here > UndamagedDepth) BaseDistInRead[posi_here][base][0][strand]++;
	  if (posi_here > UndamagedDepth) BaseDistInRead[isop_here][base][1][strand]++;

	  // Qscore calibration:
	  if (correct_base > -1 && posi_here > UndamagedDepth && isop_here > UndamagedDepth) {
	    QscoreVsErrorFreq[Qscore][base][0]++;
	    if (base != correct_base) QscoreVsErrorFreq[Qscore][base][1]++;
	  }

	  // Correctness/matching of reads vs. reference base:
	  if (correct_base > -1 && refB > -1) {
	    ReadsVsRef[0]++;
	    if (correct_base != refB) ReadsVsRef[1]++;
	  }

	  if (strand == 0) UpperStrands++;
	  
	  
	  // The major table for mapping error rates as a function of:
	  //   Position, Read-direction, Strand, Qscore, reference base, observed base
	  if (correct_base > -1 && refB > -1 && correct_base == refB && posi_here < ReadLengthMax && isop_here > UndamagedDepth && Qscore < QscoreMax) {
	    CorWro[posi_here][0][strand][Qscore][refB][base]++;
	    SigVsReadL[ReadLength][posi_here][0][strand][refB][base]++;
	  }
	  if (correct_base > -1 && refB > -1 && correct_base == refB && isop_here < ReadLengthMax && posi_here > UndamagedDepth && Qscore < QscoreMax) {
	    CorWro[isop_here][1][strand][Qscore][refB][base]++;
	    SigVsReadL[ReadLength][isop_here][1][strand][refB][base]++;
	  }



	  // Check if there is any correlation between adjacent (+2,+1,0,-1,-2 positions) allele misreads:
	  // ---------------------------------------------------------------------------------------------
	  // If position can be used by correlation analysis:
	  if (UsePosCorr) {
	    
	    // fprintf(outfile, "  %1d  %3d %3d  %2d  %1d  %3d \n", base, posi_here, isop_here, Qscore, strand, l);
	    
	    
	    // If this base can be used by correlation analysis (not too close to edge of read):
	    if (posi_here > UndamagedDepth && UndamagedDepth < isop_here) {
	      
	      // Check Q-scores, both raw and corrected:
	      CheckQscore[correct_base][0] += prob_Qscore;
	      CheckQscore[correct_base][1] += prob_Qscore_corr;
	      if (base != correct_base) CheckQscore[correct_base][2] += 1.0;
	      CheckQscore[correct_base][3] += 1.0;
	    }
	    

	    // Check model of Q-scores and transition probabilities (should work for all positions!):
	    // --------------------------------------------------------------------------------------
	    int notbase[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
	    CheckPQ[base][base][0] += 1.0 - prob_Qscore;
	    CheckPQ[notbase[base][0]][base][0] += prob_Qscore / 3.0;
	    CheckPQ[notbase[base][1]][base][0] += prob_Qscore / 3.0;
	    CheckPQ[notbase[base][2]][base][0] += prob_Qscore / 3.0;
#if 0
	    like_calc->update_pbase(l);
	    CheckPQ[0][base][1] += like_calc->base_prob[0];
	    CheckPQ[1][base][1] += like_calc->base_prob[1];
	    CheckPQ[2][base][1] += like_calc->base_prob[2];
	    CheckPQ[3][base][1] += like_calc->base_prob[3];
#endif	    
	    // CheckPQ[correct_base][base][1] += like_calc->base_prob[correct_base];
	    CheckPQ[correct_base][base][2] += 1.0;
	  }
	  
	  // There is a lookuptable called refToInt which maps
	  //a->0,A->0,c->1,C->1,g->2,G->2,t->3,T=>3,n->4,N->5
	  
	  // Random number:
	  // double rand_uniform = distribution(generator);
	  // double rand_uniform = double((rand() % 1000000) + 0.5) / 1000000.0;

	} // End of loop over depth.



	// Strand distribution:
	StrandDist[Depth][UpperStrands]++;
	
	if (UsePosCorr) {
	  // fprintf(outfile, " \n");
	}


      }    // End of loop over individuals
    }      // End of loop over sites

  }        // End of "doTemplate 2"

}


