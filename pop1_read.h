#pragma once
#include <htslib/hts.h>
#include <htslib/sam.h>

typedef struct{
  hts_idx_t *hts_idx;
  hts_itr_t *hts_itr;
}iter_t;


int bam_iter_read2(htsFile *fp, iter_t *iter,bam1_t *b,bam_hdr_t *hdr);
int restuff(bam1_t *b);

