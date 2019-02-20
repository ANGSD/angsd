#include <htslib/vcf.h>
#include "argStruct.h"
#include "analysisFunction.h"

#define PHREDMAX 256

typedef struct{
  char *fname;
  char *seek;
  htsFile *hts_file;
  bcf_hdr_t *hdr;
  hts_idx_t *idx;
  hts_itr_t *iter;
  int nsamples;
}htsstuff;


void htsstuff_seek(htsstuff *hs,char *seek);
htsstuff *htsstuff_init(char *fname,char *seek);
void htsstuff_destroy(htsstuff *hs);
class vcfReader{
private:
  bcf1_t *acpy;
  int32_t *pl;
float pl2ln[PHREDMAX];
  int curChr;
  float pl2ln_f(int32_t & val){
    if(val>=PHREDMAX){
      return log(pow(10.0,-0.1*val));
    } else {
      return pl2ln[val];
    }
  }
  int parseline(bcf1_t *rec,htsstuff *hs,funkyPars *r,int &balcon);

public:
  htsstuff *hs;
  funkyPars *fetch(int chunkSize);
  void seek(char *seek){htsstuff_seek(hs,seek);}
  vcfReader(char *fname,char *seek){
    hs=htsstuff_init(fname,seek);
    acpy=NULL;
    for(int i=0;i<PHREDMAX;i++)
      pl2ln[i] = log(pow(10.0,-0.1*i));
    curChr=-1;
    pl=NULL;
  }
  ~vcfReader(){htsstuff_destroy(hs);free(pl);}
};
