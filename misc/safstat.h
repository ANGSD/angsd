#include <vector>
#include <cmath>
#include "Matrix.hpp"

//estimate coefficient need for fsts
void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat,int whichFst);
//use coefficients to estimate persite fsts along with alpha and beta used for multilocus weighted estimator
void block_coef(Matrix<float > *gl1,Matrix<float> *gl2,double *prior,double *a1,double *b1,std::vector<double> &ares,std::vector<double> &bres,int *remap,int *adjust);
int choose(int n,int m);
int choose(size_t n,int m);
int fst_print(int argc,char **argv);
int fst_stat(int argc,char **argv);
int fst_stat2(int argc,char **argv);
