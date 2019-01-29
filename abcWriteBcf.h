#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

class abcWriteBcf:public abc{
private:
  htsFile *fp;
  bcf_hdr_t *hdr;
  bcf1_t *rec; 
  char* refName;
  char* ancName;
  int doPost;
  int doMajorMinor;
  int gl;
  int doMaf;
  int doCounts;
  int doGeno;
  argStruct *args;
  const char *outfiles;
public:
  int doBcf;
  abcWriteBcf(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWriteBcf();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
