#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <htslib/bgzf.h>

#include <cstring>
#include "safstat.h"
#include "realSFS_args.h"
//#include "realSFS.h"
//lazy binomial
int choose(int n,int m){
  if(n==2&&m==2)
    return 1;
  else if(n==3&&m==2)
    return 3;
  else{
    fprintf(stderr,"\t-> Never here\n");
  }
}


void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat){
  fprintf(stderr,"\t-> [%s] sfs1:%d sfs2:%d dimspace:%d \n",__FUNCTION__,sfs1,sfs2,(sfs1+1)*(sfs2+1));
  *aMat = new double[(sfs1+1)*(sfs2+1)];
  *bMat = new double[(sfs1+1)*(sfs2+1)];
  int at=0;
  for(int a1=0;a1<=sfs1;a1++)
    for(int a2=0;a2<=sfs2;a2++){
      double p1 = 1.0 * a1/(1.0*sfs1);
      double p2 = 1.0 * a2/(1.0*sfs2);
      double q1 = 1 - p1;
      double q2 = 1 - p2;
      double alpha1 = 1 - (p1*p1 + q1*q1);
      double alpha2 = 1 - (p2*p2 + q2*q2);
      
      double al =  0.5 * ( pow(p1-p2,2.0) + pow(q1-q2,2)) - (sfs1+sfs2) *  (sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      double bal = 0.5 * ( pow(p1-p2,2) + pow(q1-q2,2)) + (4*sfs1*sfs2-sfs1-sfs2)*(sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      (*aMat)[at] = al;
      (*bMat)[at] = bal;
      //      fprintf(stderr,"p1:%f p2:%f q1:%f q2:%f alhpa1:%f alpha:%f al:%f bal:%f\n",p1,p2,q1,q2,alpha1,alpha2,al,bal);
      at++;
    }
}

void block_coef(Matrix<float > *gl1,Matrix<float> *gl2,double *prior,double *a1,double *b1,std::vector<double> &ares,std::vector<double> &bres){
  double tre[3]={0,0,0};//a/b,sum(a),sum(0)
  for(int s=0;s<gl1->x;s++){
    int inc =0 ;
    double tmp[(gl1->y+1)*(gl2->y+1)];
    
    for(int i=0;i<gl1->y;i++)
      for(int j=0;j<gl2->y;j++)
	tmp[inc++] = prior[inc]* gl1->mat[s][i] *gl2->mat[s][j];
    
    double as=0;
    double bs=0;
    void normalize(double *,int);
    normalize(tmp,inc);
    for(int i=0;i<inc;i++){
      as += a1[i]*tmp[i];
      bs += b1[i]*tmp[i];
    }
    tre[0] += as/bs;
    tre[1] += as;
    tre[2] += bs;
    ares.push_back(as);
    bres.push_back(bs);
      
    //    fprintf(stdout,"%f %f\n",as,bs);
  }
  fprintf(stderr,"bres.size:%lu\n",bres.size());
  fprintf(stderr,"u:%f w:%f\n",tre[0]/(1.0*gl1->x),tre[1]/tre[2]);
}


int fst_print(int argc,char **argv){
  return 0;
}
