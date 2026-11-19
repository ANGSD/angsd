// Dear emacs, this is -*- c++ -*-
/*
  Genotype likelihoods for -GL 5 (phys), with a per sample model of
  sequencing errors and ancient DNA damage. DEPRECATED: will be removed
  in a future version.

  Each bam file needs a coefficient file <bam>.phys with two lines:
    QscoreCorrections c0 c1 ... (1 to 61 values, one per quality score)
    ModelParameters p0 ... p17  (18 values)
  p0-p11 are the off diagonal relative error rates (true base -> observed
  base) for reads on the forward strand, reverse strand reads use the
  complemented matrix. p12-p14 and p15-p17 are a*exp(-b*dist)+c for the
  C->T and G->A transitions as a function of the distance to the nearest
  read end.
*/
#ifndef PHYS_GENOLIKE_CALC_H
#define PHYS_GENOLIKE_CALC_H

#include <vector>
#include "bambi_interface.h"

class phys_genolike_calc {
 private:
  struct physPars{
    float qscore_corr[61];//correction factor per quality score
    float parlist[18];
    float pararray[2][4][4];//[strand][true base][observed base]
  };
  std::vector<physPars> pars;//one per sample

  void update_pbase(int depth,const physPars &p,tNode *nd,float *base_prob);

 public:
  void read_coef(const char *parspath);
  void get_genolikes(int site,int sample,chunkyT *chk,int trim,double *return_likes);
};

#endif
