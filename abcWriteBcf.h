#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

class abcWriteBcf:public abc{
private:
  htsFile *fp;
  kstring_t *kstr;
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
public:
  int doVcf;
  abcWriteBcf(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWriteBcf();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
