
#include "abc.h"

//these are the major builtin analysis that angsd can perform
#include "abcFilter.h"
#include "abcMajorMinor.h"
#include "abcFreq.h"
#include "abcError.h"
#include "abcGL.h"
#include "abcAsso.h"
#include "abcHWE.h"
#include "abcHWE_F.h"
#include "abcAncError.h"
#include "abcHaploCall.h"
#include "abcDstat.h"
#include "abcDstat2.h"
#include "abcWriteFasta.h"
#include "abcCallGenotypes.h"
#include "abcGetFasta.h"//for reading fasta; ancestral and refernce
#include "abcCounts.h" //generate counts from reads
#include "abcSaf.h" //original
#include "abcCovar.h" //calculate covar
#include "abcMismatch.h" //mismatch matrix used for GL project
#include "abcFilterSNP.h" //some snp filters, not finished yet
#include "abcSnpTools.h" //<-implementation of some stuff from snptools, not finished yet 
#include "abcHetPlas.h" //<-implementation of hetero plasmic
#include "abcWritePlink.h" //<- dump plink files.
#include "abcSmartCounts.h"
#include "abcTemplate.h"
#include "abcAncestry.h"
#include "abcIBS.h"
#include "abcWriteVcf.h" //<- dump plink files.
#include "abcPSMC.h"
#include "abcScounts.h"
#include "abcWriteBcf.h"
//below we set some variables that are shared between all the analysis classes
#define MAX_CLASS 30

int abc::tot_index =0;
const bam_hdr_t *abc::header = NULL;
const aMap *abc::revMap = NULL;
char *abc::shouldRun = new char[MAX_CLASS]; 

//here we instiantiate all the analysis classes, and return them as an array of pointers to the <abc> class. First parameter will contain the number of analysis classes.
abc **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments){
  memset(abc::shouldRun,'1',MAX_CLASS);
  
  int nit=0;
  
  abc **tskStuff =new abc*[MAX_CLASS];
  tskStuff[nit++] = new abcFilter(arguments);//0
  tskStuff[nit++] = new abcGetFasta(arguments);//1
  tskStuff[nit++] = new abcCounts(outfiles,arguments,inputtype);//2
  tskStuff[nit++] = new abcError(outfiles,arguments,inputtype);//3
  tskStuff[nit++] = new abcGL(outfiles,arguments,inputtype);//4
  tskStuff[nit++] = new abcMajorMinor(outfiles,arguments,inputtype);//5
  tskStuff[nit++] = new abcFreq(outfiles,arguments,inputtype);//6
  tskStuff[nit++] = new abcAsso(outfiles,arguments,inputtype);//7
  tskStuff[nit++] = new abcHWE(outfiles,arguments,inputtype); // 8
  tskStuff[nit++] = new abcAncError(outfiles,arguments,inputtype);//9
  tskStuff[nit++] = new abcCallGenotypes(outfiles,arguments,inputtype);//10
  tskStuff[nit++] = new abcSaf(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcCovar(outfiles,arguments);
  tskStuff[nit++] = new abcTsk(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcFilterSNP(outfiles,arguments,inputtype);//14
  tskStuff[nit++] = new abcSnpTools(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcHetPlas(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcWritePlink(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcDstat(outfiles,arguments,inputtype);//18
  tskStuff[nit++] = new abcWriteFasta(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcSmartCounts(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcTemplate(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcWriteVcf(outfiles,arguments,inputtype);//22
  tskStuff[nit++] = new abcHaploCall(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcDstat2(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcHWE_F(outfiles,arguments,inputtype); // 25
  tskStuff[nit++] = new abcIBS(outfiles,arguments,inputtype); // 
  tskStuff[nit++] = new abcPSMC(outfiles,arguments,inputtype); //
  tskStuff[nit++] = new abcScounts(outfiles,arguments,inputtype); //
  tskStuff[nit++] = new abcWriteBcf(outfiles,arguments,inputtype); // 29


  nItem = nit;
  return tskStuff;
}
