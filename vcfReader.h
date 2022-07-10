#include <htslib/vcf.h>
#include<cstdio>
#include <cstring>

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

bam_hdr_t *bcf_hdr_2_bam_hdr_t (htsstuff *hs);
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
  //type=0 -> PL
  //type=1 -> GL
  int parseline(bcf1_t *rec,htsstuff *hs,funkyPars *r,int &balcon,int type);
  int pl_gl_gp;
  int ln_gl_m;
  float *ln_gl = NULL;
  float *farr = NULL;
  int32_t *iarr = NULL;
  int mfarr;
  int miarr;
  std::vector<regs> *regions;//this is a pointer. I should really be more consistent, but i dont want to copy construct
  kstring_t itrname;
public:
  htsstuff *hs;
  bam_hdr_t *bamhdr; //wrapper for converting bcf/vcf style header to bam header
  funkyPars *fetch(int chunkSize);
  void seek(char *seek){htsstuff_seek(hs,seek);}
  int vcfReaderwrap_reader(htsstuff *hts,bcf1_t *rec);
  vcfReader(char *fname,char *seek,int pl_or_gl_a,std::vector<regs> *regions_a);
  ~vcfReader(){
    htsstuff_destroy(hs);
    if(pl)
      free(pl);
    if(ln_gl!=NULL)
      free(ln_gl);
    if(iarr)
      free(iarr);
    if(farr)
      free(farr);
    free(itrname.s);}
};
