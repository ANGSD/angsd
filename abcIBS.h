#include "abc.h"
typedef struct{
  int **dat;
  int *major;
  int *minor;
  int *ibsMatrix; 
  int *nonMis;// N x N matrix of non missing comparisons for ibs matrix
  double *covMatrix;
  int *covMis; // N x N matrix of non missing comparisons for covariance matrix

}IBSstruct;

class abcIBS:public abc{
private:

  //non optional arguments
  int doIBS;
  int majorminor;
  int doCount;
  int output01;
  int doCov;
  int intToMajorMinorAA[5];
  //optional arguments
  int maxMis;
  int minMinor;
  double minFreq;
  int makeMatrix;
  int *ibsMatrixAll;
  int *nonMisAll;
  double *covMatrixAll;
  int *nonMisCov;
  //out file
  BGZF* outfileZ;
  FILE* outfileMat;
  FILE* outfileCov;
  int nInd;

  //functions
  void printHaplo(funkyPars *pars);
  void getHaplo(funkyPars *pars);
  void makeIBSmatrix(funkyPars *pars);
  void makeCovMatrix(funkyPars *pars);

  //print buffer
  kstring_t bufstr;

public:
  abcIBS(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcIBS();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};
