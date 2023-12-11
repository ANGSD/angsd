#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <htslib/vcf.h>

#include "argStruct.h"
#include "analysisFunction.h"
#include "aio.h"

#define PHREDMAX 256

#define BCF_TAG_FMT_PL (1<<0)
#define BCF_TAG_FMT_GL (1<<1)
#define BCF_TAG_FMT_GP (1<<2)
#define BCF_TAG_FMT_AD (1<<3)
#define BCF_TAG_FMT_DP (1<<4)

#define BASE_NONREF 4



typedef struct {
  char* fname;
  char* seek;
  htsFile* hts_file;
  bcf_hdr_t* hdr;
  hts_idx_t* idx;
  hts_itr_t* iter;
  int nsamples;
}htsFileStruct;

bam_hdr_t* bcf_hdr_2_bam_hdr_t(htsFileStruct* hs);
void htsFileStruct_seek(htsFileStruct* hs, char* seek);
htsFileStruct* htsFileStruct_init(char* fname, char* seek);
void htsFileStruct_destroy(htsFileStruct* hs);
class bcfReader {
private:
  bcf1_t* acpy;
  int32_t* pl;
  float pl2ln[PHREDMAX];
  int curChr;
  float pl2ln_f(int32_t& val) {
    if (val >= PHREDMAX) {
      return log(pow(10.0, -0.1 * val));
    } else {
      return pl2ln[val];
    }
  }
  //type=0 -> PL
  //type=1 -> GL
  int read_bcf_record(bcf1_t* rec, htsFileStruct* hs, funkyPars* r, int& balcon);
  uint32_t bcf_tags = 0;
  int ln_gl_m;
  float* ln_gl, * farr;
  int32_t* iarr;
  int32_t* ad_arr;
  int mfarr;
  int miarr;
  int32_t m_ad_arr;
  std::vector<regs>* regions;//this is a pointer. I should really be more consistent, but i dont want to copy construct
  kstring_t itrname;
public:
  htsFileStruct* hs;
  bam_hdr_t* bamhdr; //wrapper for converting bcf/vcf style header to bam header
  funkyPars* fetch(int chunkSize);

  // for doFasta
  suint* count_bases(funkyPars* pars, int32_t* ad_arr, const int n, bcf1_t* rec);

  void seek(char* seek) { htsFileStruct_seek(hs, seek); }
  int bcfReaderwrap_reader(htsFileStruct* hts, bcf1_t* rec);
  bcfReader(char* fname, char* seek, uint32_t bcf_tags_to_use, std::vector<regs>* regions_a);
  ~bcfReader() {
    htsFileStruct_destroy(hs);
    if (pl)
      free(pl);
    if (ln_gl != NULL)
      free(ln_gl);
    if (iarr)
      free(iarr);
    if (ad_arr)
      free(ad_arr);
    if (farr)
      free(farr);
    free(itrname.s);
  }
};
