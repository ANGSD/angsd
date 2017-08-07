/*
  These files contains modified versions from SAMtools (1.5-1-g27b628e (using htslib 1.5-2-g2739558))
  thorfinn@binf.ku.dk 7aug 2017 cambridge
 */
#pragma once
#include <stdlib.h>
#include <htslib/khash_str2int.h>
#ifdef __cplusplus
extern "C" {
#endif
  
  void* add_read_group_single(char *name);
  void* add_read_groups_file(char *fn);
#ifdef __cplusplus
}
#endif
