
#include <cmath>
#include <cstdio>
#include "bambi_interface.h"

extern int refToInt[256];
double **probs = NULL;


/*
  allocate the prob matrix with
  prob[1][qs] = homohit
  prob[2][qs] = hethit
  prob[3][qs] = homoNonhit

  I'm perfecly aware that these are not genotypes, but what should these be called
 */

double **genLikes(int len){
  double **ret = new double*[3];
  for(int i=0;i<3;i++)
    ret[i] = new double[256];
  
  for(int i=0;i<len;i++){
    double p = 1-pow(10,-i/10.0);
    double p3 = (1-p)/3.0;
    ret[0][i] = log(p*0.5+p*0.5);
    ret[1][i] = log(p*0.5+p3*0.5);
    ret[2][i] = log(0.5*p3+0.5*p3);
  }
  return ret;
}

int offsets[4][10]={
  {0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
  {4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
  {7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
  {9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};

void gatk_init(){
  probs = genLikes(256);
}

void gatk_destroy(){
  for(int i=0;i<3;i++)
    delete [] probs[i];
  delete [] probs;
  probs=NULL;
}

void call_gatk(chunkyT *chk,double **lk,int trim){
  
  for(int s=0;s<chk->nSites;s++){
    for(int i=0;i<chk->nSamples;i++){
      //allocate for all samples, for first sample
      if(i==0){
	lk[s] = new double[10*chk->nSamples];
	for(int ii=0;ii<10*chk->nSamples;ii++)
	  lk[s][ii] = -0.0;//set default values
      }
      
      tNode &nd = chk->nd[s][i];
      
      //calc like persample
      double *likes1 = lk[s]+10*i;
      for(int j=0;j<nd.l;j++){
	int allele = refToInt[nd.seq[j]];
	int qs = nd.qs[j];
	//filter qscore, mapQ,trimming, and always skip n/N
	if(nd.posi[j]<trim||nd.isop[j]<trim||allele==4){
	  continue;
	}

	likes1[offsets[allele][0]] += probs[0][qs]; //'homo'zygotic hit
	for(int o=1;o<4;o++){//'heterozygotic' hit{
	  likes1[offsets[allele][o]] += probs[1][qs];
	}
	for(int o=4;o<10;o++){//non hit
	  likes1[offsets[allele][o]] += probs[2][qs];
	}
	
      }
     
    }
  }
}
