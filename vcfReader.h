#pragma once
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <htslib/vcf.h>
#include <htslib/tbx.h>

#include "argStruct.h"
#include "analysisFunction.h"

#define PHREDMAX 256

typedef struct{
  char *fname;
  char *seek;
  htsFile *hts_file;
  bcf_hdr_t *hdr;
  hts_idx_t *idx;//index for bcf files
  tbx_t *tbx;//index for bgzipped vcf files
  kstring_t line;//text line read with tbx
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
      //not log(pow(10,-0.1*val)), the pow is 0 for PL above ~3080 and the log is then -inf
      return -0.1*M_LN10*val;
    } else {
      return pl2ln[val];
    }
  }
  //type=0 -> PL
  //type=1 -> GL
  //type=2 -> GP
  int parseline(bcf1_t *rec,htsstuff *hs,funkyPars *r,int &balcon,int type);
  int isindel(const bcf1_t *v);
  int pl_gl_gp;
  int ln_gl_m;
  float *ln_gl,*farr;
  int32_t *iarr;
  int mfarr;
  int miarr;
  std::vector<regs> *regions;//this is a pointer. I should really be more consistent, but i dont want to copy construct
  kstring_t itrname;
  int whichRegion;//region currently being read, when regions are given
  int hasSeeked;//has the iterator been set for whichRegion
  int indelTag;//id of the INFO/INDEL tag in the header, -1 if absent
  int nWarn[5];//number of times each kind of warning has been printed
  int warnOnce(int which);//returns 1 if a warning of this kind should still be printed
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
    free(itrname.s);
  }
};
