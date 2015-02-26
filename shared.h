
extern int SIG_COND;




#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <map>
#include <zlib.h>
#include "bambi_interface.h"
#include "argStruct.h"

#ifndef _types_t
#define _types_t
typedef unsigned int suint;

#define LENS 10000
#define GZOPT "w6h"


enum{INPUT_BAM,INPUT_GLF,INPUT_GLF3,INPUT_BEAGLE,INPUT_PILEUP,INPUT_VCF_GL,INPUT_VCF_GP};

typedef struct {

  //primitives
  int chunkNumber; //4
  int numSites; //8
  int nInd; //12


  int refId;//reference number from header //16
  int *posi;//position is 0 indexed //24

  suint **counts;// contains counts of A,C,G,T //32 
  double **likes;// contains logscaled likelihoodratios, 10 per sample //40
  double **post;// contains genotype posteriors, 3 persample DRAGON //48
  
  char *major;//major //
  char *minor;//minor
  char *ref;//reference
  char *anc;//ancestral
 
  //sites removed in analysis and print
  int *keepSites;
 
  
  //stuff needed for bamreader
  fcb *for_callback;
  chunkyT *chk;
  int killSig;
  
  //extra stuff associated with each analysis module
  void **extras;
  
}funkyPars;

funkyPars *allocFunkyPars();
void deallocFunkyPars(funkyPars *p);


void init(argStruct *arguments);//intialize all needed objects
void destroy();//destroy all initialized objects
void selector(funkyPars *p);

#define LOGMAX 20000   // pre-computed logfactorial 

#endif
