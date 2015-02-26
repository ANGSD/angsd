#pragma once
//DEFAULT TRIMMING
#define MINQ 13
#define EM_START 0.001
#define EM_NITER 100

#include "argStruct.h"
#include "analysisFunction.h"
class abc{
public:
  static char *shouldRun;
  const static bam_hdr_t *header;//contains the header of a single bam;
  const static aMap *revMap;
  int index;
  static int tot_index;
  //  virtuel general()
  virtual void run(funkyPars *f)=0;
  virtual void print( funkyPars *f)=0;
  virtual void printArg(FILE *fp)=0;
  virtual void clean(funkyPars *f)=0;
  abc(){index=tot_index++;};
  virtual ~abc(){};
};


abc **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments);
