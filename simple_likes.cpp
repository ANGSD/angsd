/*
  -GL 6, a very simple model meant for testing.

  For every individual and site one base is sampled from the counts
  (-doCounts 1), and the genotype likelihoods are 1 for the homozygote of
  that base, 0.5 for the three heterozygotes with it and 0 for the other
  six genotypes. The sampling uses drand48, which abcGL seeds from -seed
  or the clock. With -nThreads above 1 the threads share drand48, so the
  same seed only gives the same result with -nThreads 1.
 */

#include <cmath>
#include "bambi_interface.h"
#include "simple_likes.h"
#include "analysisFunction.h"

// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
//AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
simple_likes::simple_likes(){
  for(int b=0;b<4;b++){
    for(int i=0;i<4;i++)
      for(int ii=i;ii<4;ii++)
	probs[b][angsd::majorminor[i][ii]] = -INFINITY;//genotypes without b
    for(int i=0;i<4;i++)
      probs[b][angsd::majorminor[b][i]] = log(0.5);//heterozygotes with b
    probs[b][angsd::majorminor[b][b]] = 0.0;//homozygote for b, log(1)
  }
}

void simple_likes::run(suint **counts,int *keepSites,double **lk,int nSites,int nSamples) const{
  for(int s=0;s<nSites;s++){
    lk[s] = new double[10*nSamples];
    for(int ii=0;ii<10*nSamples;ii++)
      lk[s][ii] = -0.0;//set default values
    if(keepSites[s]==0)
      continue;

    for(int i=0;i<nSamples;i++){
      int r = angsd::getRandomCount(counts[s],i,-1);
      if(r==4)//no data for this individual
	continue;
      for(int j=0;j<10;j++)
	lk[s][i*10+j]=probs[r][j];
    }
  }
}
