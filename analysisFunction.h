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
  int getArg(const char* argName,int type,argStruct *arguments);
  float getArg(const char* argName,float type,argStruct *arguments);
  char* getArg(const char* argName,char* type,argStruct *arguments);
  char* getArg(const char* argName,const char* type,argStruct *arguments);
  double getArg(const char* argName,double type,argStruct *arguments);
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
  std::vector<char*> getFilenames(const char * name,int nInd);
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

//angsd io
namespace aio{
  size_t fsize(const char* fname);
  int fexists(const char* str);//{///@param str Filename given as a string.
  FILE *openFile(const char* a,const char* b);
  FILE *getFILE(const char*fname,const char* mode);
  BGZF *openFileBG(const char* a,const char* b);
  htsFile *openFileHts(const char * a, const char*b);
  int isNewer(const char *newer,const char *older);
  ssize_t bgzf_write(BGZF *fp, const void *data, size_t length);
  int tgets(gzFile gz,char**buf,int *l);
}
#endif


double phi(double x);
