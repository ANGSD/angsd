
#include "abc.h"

//these are the major builtin analysis that angsd can perform
#include "abcFilter.h"
#include "abcMajorMinor.h"
#include "abcFreq.h"
#include "abcError.h"
#include "abcGL.h"
#include "abcAsso.h"
#include "abcHWE.h"
#include "abcAncError.h"
#include "abcDstat.h"
#include "abcWriteFasta.h"
#include "abcCallGenotypes.h"
#include "abcGetFasta.h"//for reading fasta; ancestral and refernce
#include "abcCounts.h" //generate counts from reads
#include "abcSaf.h" //original
#include "abcSaf2.h" //<- EJ version, banded dynamic programming recursion (cite paper, add more infor etc)
#include "abcCovar.h" //calculate covar
#include "abcTsk.h" //mismatch matrix used for GL project
#include "abcFilterSNP.h" //some snp filters, not finished yet
#include "abcSnpTools.h" //<-implementation of some stuff from snptools, not finished yet 
#include "abcHetPlas.h" //<-implementation of hetero plasmic
#include "abcWritePlink.h" //<- dump plink files.
#include "abcSmartCounts.h"
#include "abcTemplate.h"

//below we set some variables that are shared between all the analysis classes
#define MAX_CLASS 25
int abc::tot_index =0;
const aHead *abc::header = NULL;
const aMap *abc::revMap = NULL;
char *abc::shouldRun = new char[MAX_CLASS]; 

//here we instiantiate all the analysis classes, and return them as an array of pointers to the <abc> class. First parameter will contain the number of analysis classes.
abc **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments){
  memset(abc::shouldRun,'1',MAX_CLASS);
  
  int nit=0;
  
  abc **tskStuff =new abc*[MAX_CLASS];
  tskStuff[nit++] = new abcFilter(arguments);
  tskStuff[nit++] = new abcGetFasta(arguments);
  tskStuff[nit++] = new abcCounts(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcError(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcGL(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcMajorMinor(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcFreq(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcAsso(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcHWE(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcAncError(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcCallGenotypes(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcSaf(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcCovar(outfiles,arguments);
  tskStuff[nit++] = new abcTsk(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcFilterSNP(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcSnpTools(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcHetPlas(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcWritePlink(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcDstat(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcWriteFasta(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcSmartCounts(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcSaf2(outfiles,arguments,inputtype);
  tskStuff[nit++] = new abcTemplate(outfiles,arguments,inputtype);

  //don't touch below
  nItem = nit;
  return tskStuff;
}
