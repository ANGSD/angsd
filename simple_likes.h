#pragma once
#include "bambi_interface.h"

//-GL 6, samples one base per individual and site, for testing
class simple_likes{
  //log likelihood of the 10 genotypes given the sampled base A,C,G,T
  double probs[4][10];
public:
  simple_likes();
  void run(suint **counts,int *keepSites,double **lk,int nSites,int nSamples) const;
};
