#include <htslib/hts.h>
#include <htslib/sam.h>
int cigstat_calc(bam1_t *rd);
int cigstat_init(const char *fname);
int cigstat_close();
