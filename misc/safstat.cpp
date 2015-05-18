//input are the dimension of the sfs for the two pops.
//for diploid humands sfs can be 2,4,6 etc. but we will ofcause adjoint with the zero category
#include <cmath>
void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat){
  *aMat = new double[(sfs1+1)*(sfs2+1)];
  *bMat = new double[(sfs1+1)*(sfs2+1)];
  for(int a1=0;a1<=sfs1;a1++)
    for(int a2=0;a2<=sfs1;a2++){
      double p1 = a1/sfs1;
      double p2 = a2/sfs2;
      double q1 = 1 - p1;
      double q2 = 1 - p2;
      double alpha1 = 1 - (p1*p1 + q1*q1);
      double alpha2 = 1 - (p2*p2 + q2*q2);
      
      double al =  1/2 * ( pow(p1-p2,2.0) + pow(q1-q2,2)) - (sfs1+sfs2) *  (sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      double bal = 1/2 * ( pow(p1-p2,2) + pow(q1-q2,2)) + (4*sfs1*sfs2-sfs1-sfs2)*(sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      *aMat[a1*sfs1+a2] = al;
      *bMat[a2*sfs2+a1] = bal;
    }
  
}
