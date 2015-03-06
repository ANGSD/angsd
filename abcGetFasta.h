#pragma once

#include <htslib/faidx.h>
#include "abc.h"//ONLY USED for reference stuff in mpileup output.
#include "shared.h"//ONLY USED for reference stuff in mpileup output.
//generic struct for holding the faidx, the mutex and the actual data.
typedef struct{
  char *fastaname;
  pthread_mutex_t aMut;//mutex will be locked unlocked whenever we do something
  faidx_t *fai;//contains the faidx structure
  char *seqs;//contains the reference for the current chr;
  int curChr;//the exact chromosome name for the seqs above
  int chrLen;//length of chromosome
}perFasta;

class abcGetFasta : public abc{
  char *ancName;
  char *refName;
  char *getRef(int refId,int *posi,int lens,faidx_t *fai);
  char *getAnc(int refId,int *posi,int lens,faidx_t *fai);
  faidx_t *faiRef;
  faidx_t *faiAnc;
public:
  char *loadChr(perFasta *f, char*chrName,int chrId);
  char *magic(int refId,int *posi,int numSites,perFasta *f);
  perFasta *ref;
  perFasta *anc;  
  abcGetFasta(argStruct *arguments);
  ~abcGetFasta();
  void printArg(FILE *argFile);
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);//plugs in reference and ancestral
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
};
