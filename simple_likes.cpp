/*
 
 */


#include <cmath>
#include <cstdio>
#include <cassert>
#include "bambi_interface.h"
#include "simple_likes.h"
#include "analysisFunction.h"

double **probs_simple = NULL;

// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
//AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
double **genLikes_simple(){
  double **ret = new double*[4];
  for(int b=0;b<4;b++){
    ret[b]=new double[10];
    ret[b][angsd::majorminor[b][b]] = log(1.0);
    //    fprintf(stderr,"offss:%d\n",angsd::majorminor[b][b]);
    for(int i=0;i<4;i++){
      if(b!=i){
	ret[b][angsd::majorminor[b][i]] = log(0.5);
	//	fprintf(stderr,"offss2:%d\n",angsd::majorminor[b][i]);
      }
    }
    for(int i=0;i<4;i++)
      for(int ii=i;ii<4;ii++){
	if(i!=b&&ii!=b){
	  ret[b][angsd::majorminor[i][ii]] = log(0.0);
	  //fprintf(stderr,"offss3:%d\n",angsd::majorminor[i][ii]);
	}
      }
    //exit(0);
  }
  
  return ret;
}


void simple_init(){
  probs_simple = genLikes_simple();
  for(int b=0;0&&b<4;b++){
    fprintf(stderr,"%d) ",b);
    for(int i=0;i<10;i++)
      fprintf(stderr,"%f ",probs_simple[b][i]);
    fprintf(stderr,"\n");
  }
}

void simple_destroy(){
  for(int i=0;i<3;i++)
    delete [] probs_simple[i];
  delete [] probs_simple;
  probs_simple=NULL;
}


void call_simple(suint **counts,int *keepSites,double **lk,int nSites,int nSamples){
  
  for(int s=0;s<nSites;s++){
    for(int i=0;i<nSamples;i++){
      //allocate for all samples, for first sample
      if(i==0){
	lk[s] = new double[10*nSamples];
	for(int ii=0;ii<10*nSamples;ii++)
	  lk[s][ii] = -0.0;//set default values
      }
      if(keepSites[s]==0)
	continue;

      int r = angsd::getRandomCount(counts[s],i,-1);
      //      fprintf(stdout,"ris\t%d\n",r);
      assert(r>-1);
      if(r==4)
	continue;
      for(int j=0;j<10;j++){
	lk[s][i*10+j]=probs_simple[r][j];
      }
    }
  }
}
