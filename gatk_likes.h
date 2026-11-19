#pragma once
#include "bambi_interface.h"

//-GL 2, the model from the original GATK paper
class gatk_likes{
  //log probability of a base with a given quality score under a genotype
  //[0] homozygous for the base, [1] heterozygous with it, [2] without it
  double probs[3][256];
public:
  gatk_likes();
  void run(chunkyT *chk,double **lk,int trim,int *keepSites) const;
};
