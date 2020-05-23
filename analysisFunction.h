#include <zlib.h>
#ifndef _angsd_
#define _angsd_


#define FOR(i,n) for(i=0; i<n; i++)
#include <cmath>
#include <vector>
#include "shared.h"
#include "chisquare.h" 
#include <htslib/bgzf.h>

extern int refToInt[256];
extern char intToRef[5];
extern char intToIupac[15];




namespace angsd {

  double lbico(double n, double k);
  double myComb2(int k,int r, int j);
  
  //assuming samtools ordering
 
  static int majorminor[4][4] = 
    {
      {0,1,2,3},
      {1,4,5,6},
      {2,5,7,8},
      {3,6,8,9}
    };
 

  template<typename T>
  struct Matrix {
    int x;
    int y;
    T **matrix;
  };
  void norm(double *d,size_t len);
  double getMax(double a,double b, double c);
  double addProtect2(double a,double b);
  double addProtect3(double a,double b, double c);
  double addProtectN(double a[],int N);
  Matrix<double> getMatrix(const char *name,int doBinary,int lens);
  Matrix<int> getMatrixInt(const char *name,int lens);
  int fexists(const char* str);
  double sigm(double x);
  double **get3likes(funkyPars *pars);
  double **get3likesRescale(funkyPars *pars);
  double **get3likes(funkyPars *pars,int *keepInd);
  double **get3likesRMlow(funkyPars *pars,int *keepInd);
  double **getlikes(funkyPars *pars,int *keepInd);
  void swapDouble (double& first, double& second);
  int matinv( double x[], int n, int m, double space[]);
  void deleteMatrix(Matrix<double> mat);
  void deleteMatrixInt(Matrix<int> mat);
  void printMatrix(Matrix<double> mat,FILE *file);
  void logrescale(double *ary,int len);
  int svd_inverse(double mat[],int xLen, int yLen);
  double dnorm(double x,double mean,double sd,int ifLog);
  double bernoulli(int k, double p, int ifLog);
  double sd(double* phe, int size );
  double poisson(double k,  double lambda, int ifLog);
  double to_pval(Chisqdist *chisq,double f);
  char *strpop(char **str,char split);
  int getRandomCount(suint *d, int i,  int depth = -1);
  int getMaxCount(suint *d, int i, int depth = -1);
  int getIupacCount(suint *d, int i, double iRatio, int depth = -1);
  int getRandomCountTotal(suint *d, int nInd);
  int getMaxCountTotal(suint *d, int nInd);
  int getIupacCountTotal(suint *d, int nInd, double iRatio);
  double estFreq(double *loglike,int numInds);


  template <typename T>
  T * allocArray(size_t len,T defval){
    T *ret= new T[len];
    for(size_t i=0;i<len;i++)
      ret[i]=defval;
    return ret;
    

  }
  template <typename T>
  T * allocArray(size_t len){
    T *ret= new T[len];
    return ret;
    
  }

  template <typename T>
  T sum(const T *ary,size_t len){
  T s =0;
  for(size_t i=0;i<len ; i++)
    s+=ary[i];
  //  printf("sum:%f\n",s);
  return s;
  }

  void print_array(FILE *fp,double *ary,int len);
  void print_array(FILE *fp,int *ary,int len);

  double *readDouble(const char*fname,int hint);
  int whichMax(double *d,int len);
}


#endif


double phi(double x);


