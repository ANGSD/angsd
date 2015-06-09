#include <vector>
#include <cmath>
#include "Matrix.hpp"

void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat);
void block_coef(Matrix<float > *gl1,Matrix<float> *gl2,double *prior,double *a1,double *b1,std::vector<double> &ares,std::vector<double> &bres);
int choose(int n,int m);
int fst_print(int argc,char **argv);
int fst_stat(int argc,char **argv);
