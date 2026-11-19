/*
  The GLs directly from the qscore, assuming diploid.

  This is the model from the original GATK paper, and is most likely very
  different from the newest GATK. A base with error probability e has
  probability 1-e under its own allele and e/3 under each of the others.
 */

#include <cmath>
#include <cstdio>
#include "bambi_interface.h"
#include "gatk_likes.h"

extern int refToInt[256];

//e above this would make a base more likely under the other alleles than
//its own, which happens for quality scores 0 and 1. At 0.75 the base has
//probability 1/4 under every allele, so it does not change the GLs.
static const double MAX_ERROR = 0.75;

//for each observed base, the genotypes homozygous for it, the three
//heterozygous with it, and the six without it
static const int offsets[4][10]={
  {0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
  {4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
  {7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
  {9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};

gatk_likes::gatk_likes(){
  for(int i=0;i<256;i++){
    double e = pow(10,-i/10.0);
    if(e>MAX_ERROR)
      e = MAX_ERROR;
    double p = 1-e;
    double p3 = (1-p)/3.0;
    probs[0][i] = log(p*0.5+p*0.5);
    probs[1][i] = log(p*0.5+p3*0.5);
    probs[2][i] = log(0.5*p3+0.5*p3);
  }
}

void gatk_likes::run(chunkyT *chk,double **lk,int trim,int *keepSites) const{
  for(int s=0;s<chk->nSites;s++){
    lk[s] = new double[10*chk->nSamples];
    for(int ii=0;ii<10*chk->nSamples;ii++)
      lk[s][ii] = -0.0;//set default values
    if(keepSites[s]==0)
      continue;

    for(int i=0;i<chk->nSamples;i++){
      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;
      //calc like persample
      double *likes1 = lk[s]+10*i;
      for(int j=0;j<nd->l;j++){
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
	//skip trimmed bases and n/N, quality and mapQ are filtered before this
	if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
	  continue;
	}

	likes1[offsets[allele][0]] += probs[0][qs]; //'homo'zygotic hit
	for(int o=1;o<4;o++){//'heterozygotic' hit
	  likes1[offsets[allele][o]] += probs[1][qs];
	}
	for(int o=4;o<10;o++){//non hit
	  likes1[offsets[allele][o]] += probs[2][qs];
	}
      }
    }
  }
}
