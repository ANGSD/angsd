#include "soap_likes.h"
#include "gatk_likes.h"
#include "bam_likes.h"
#include "simple_likes.h"
#include "ancestral_likes.h"

class abcGL:public abc{
private:
  BGZF *outfileSAF;
  FILE *outfileSAFIDX;
  BGZF *outfileSAFPOS;

  const char * postfix;
  const char * beaglepostfix;

  //below is for SYK
  double **errors;
  char * errorFname;
  double ****errorProbs;
  //below is used for output
  kstring_t bufstr;
  //soap method is implemented as a class
  soap_likes soap;
  // ancestral method is implemented as a class
  anc_likes ancestral_lik;
  //a tempdir used by SOAPsnp for the recalibration matrices
  char *angsd_tmpdir;

  //below used for filtering qscore,mapQ and trimming
  int trim;

  int64_t offs[2];
  BGZF *gzoutfile;
  BGZF *gzoutfile2;
  int GL;
  int doGlf;
  int minInd;
  char *tmpChr;
  int nnnSites;
  void getLikesFullError10Genotypes(int numSites,int nInd,suint **counts,double ****errorProbs,int *keepSites,double **loglikes);
  void printLike(funkyPars *pars);

public:
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  void changeChr(int refId);
  abcGL(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcGL();

};
