#include "abc.h"

typedef struct{

  double **stat;
  int **keepInd;
  int *highHe;
  int *highWt;
  int *highHo;
}assoStruct;



class abcAsso:public abc{
private:
  kstring_t bufstr;
public:
  //none optional stuff
  BGZF **multiOutfile;
  int doPrint;
  int minCov; //not for users
  int doMaf;
  int dynCov;//not for users
  int doAsso;
  int doPost;
  int GL;
  int sitePerm;  //not for users
  int isBinary;
  int minHigh;
  int minCount;
  int adjust;  //not for users
  int model;
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);

  abcAsso(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcAsso();
  //other stuff
  char *covfile;
  char *yfile;
  
  angsd::Matrix<double> ymat;
  angsd::Matrix<double> covmat;
  void scoreAsso(funkyPars  *pars,assoStruct *assoc);
  void frequencyAsso(funkyPars  *pars,assoStruct *assoc);
  double doAssociation(funkyPars *pars,double *post,double *y,int keepInd,int *keepList,double freq,int s,assoStruct *assoc);
  int getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv);
  int getFitBin(double *res,double *Y,double *covMatrix,int nInd,int nEnv);
  double normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s);
  double binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,assoStruct *assoc,int s);
  void printDoAsso(funkyPars *pars);
};
