#include <htslib/vcf.h>
#include "argStruct.h"
#include "analysisFunction.h"

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
  htsstuff *hs;
  int curChr;
public:
  funkyPars *fetch(int chunkSize);
  vcfReader(char *fname,char *seek){ hs=htsstuff_init(fname,seek);}
  ~vcfReader(){htsstuff_destroy(hs);}
};
