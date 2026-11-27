#include "soap_likes.h"
#include "gatk_likes.h"
#include "bam_likes.h"
#include "simple_likes.h"
#include "ancestral_likes.h"

class abcGL:public abc{
private:
  const char * postfix;
  const char * beaglepostfix;

  //below is for SYK
  double **errors;
  char * errorFname;
  double ****errorProbs;
  //below is used for output
  kstring_t bufstr;
  //bam, gatk, simple and soap methods are implemented as classes
  bam_likes *bam;//only allocated with -GL 1, its tables are 32MB
  gatk_likes gatk;
  simple_likes simple;
  soap_likes soap;
  // ancestral method is implemented as a class
  anc_likes ancestral_lik;
  //a tempdir used by SOAPsnp for the recalibration matrices
  char *angsd_tmpdir;

  //below used for filtering qscore,mapQ and trimming
  int trim;

  BGZF *gzoutfile;
  BGZF *gzoutfile2;
  int GL;
  int doGlf;
  int minInd;
  int seed;//for -GL 1 (sampling of bases above depth 255) and -GL 6
  void getLikesFullError10Genotypes(int numSites,int nInd,suint **counts,double ****errorProbs,int *keepSites,double **loglikes);
  void printLike(funkyPars *pars);

public:
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcGL(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcGL();

};
