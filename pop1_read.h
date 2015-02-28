#pragma once
#include <htslib/hts.h>
#include <htslib/sam.h>

//function to fetch one read from fp,
int pop1_read(htsFile *fp, hts_itr_t *itr,bam1_t *b,bam_hdr_t *hdr);

//sometimes we need to recheck that a read is fine. This is what restuff does;
int restuff(bam1_t *b);
