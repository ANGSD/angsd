
#include "soap_likes.h"
#include "gatk_likes.h"
#include "bam_likes.h"

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
  //soap method is implemented as a class
  soap_likes soap;
  //a tempdir used by SOAPsnp for the recalibration matrices
  char *angsd_tmpdir;

  //below used for filtering qscore,mapQ and trimming
  int trim;

  
  gzFile gzoutfile;
  int GL;
  int doGlf;
  int minInd;
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
