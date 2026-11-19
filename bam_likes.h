#pragma once
#include "bambi_interface.h"

struct errmod_t;

//-GL 1, the SAMtools model from htslib's errmod
class bam_likes{
  //lookup tables from errmod_init, only read by run so the threads can share it
  errmod_t *mod;
public:
  bam_likes();
  ~bam_likes();
  void run(chunkyT *chk,double **lk,int trim,int *keepSites,int seed) const;
};
