#pragma once

#include <pthread.h>
#include <htslib/faidx.h>
#include "abc.h"
#include "shared.h"
//generic struct for holding the faidx, the mutex and the actual data.
typedef struct{
  char *fastaname;
  pthread_mutex_t aMut;//mutex will be locked unlocked whenever we do something
  faidx_t *fai;//contains the faidx structure
  char *seqs;//contains the reference for the current chr;
  int curChr;//the chromosome id (index in the BAM header) of seqs above, -1 if none
  int chrLen;//length of chromosome
}perFasta;

class abcGetFasta : public abc{
  char *ancName;
  char *refName;
  int isVcf;//input is -vcf-gl/-vcf-gp, its ##contig lengths are not checked as strictly as a BAM header
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
