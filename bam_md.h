#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  int bam_prob_realn_core(bam1_t *b, const char *ref,int flag);
  int bam_cap_mapQ(bam1_t *b, char *ref, int thres);


#ifdef __cplusplus
}
#endif
