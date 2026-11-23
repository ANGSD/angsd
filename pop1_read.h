#pragma once
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <cstdio>

//fetch the next read that passes the read filters from fp (or from itr if not NULL)
//returns >=0 for a read, -1 at the end of the file/region. Exits on read errors.
int pop1_read(htsFile *fp, hts_itr_t *itr,bam1_t *b,bam_hdr_t *hdr);

//redo -baq and -C for a read that was buffered while the reference was on another chromosome
//returns 1 if the read should be kept, 0 if it is now below -minMapQ or removed by -C
int restuff(bam1_t *b);

//print how many reads were read, used and removed by each read filter (BAM/CRAM input only)
void printReadCounts(FILE *fp);
