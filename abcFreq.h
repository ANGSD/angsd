//DRAGON fix degree of freedom in chisq unknown
#include "abc.h"

typedef struct{
  double *freq_EM; 
  double *freq_EM_unknown; 
  double *lrt_EM; 
  double *lrt_EM_unknown;
  double *freq;//can be copy of the two above, or the pp
  double *lrt;//can be a copy of hte two above
  double *phat; //counts based freq
}freqStruct;



class abcFreq:public abc{
private:
  double to_pval(Chisqdist *str,double f);
  Chisqdist *chisq1;
  Chisqdist *chisq3;
  int nInd;
  char *refName;
  char *ancName;
  char *indFname;
  static double EM_start;//static due to called from analysisEstErorr,analysisAsso
  static int emIter;//static due to called from analysisEstErorr,analysisAsso
  double eps;
  
  void prepPrint(funkyPars *pars);

  gzFile outfileZ;
  gzFile outfileZ2;

  int doMaf;
 
  int doPost;
  int doSNP;
  int doMajorMinor;
  int GL;

  double minMaf;
  double SNP_pval;
   int beagleProb;

  int inputIsBeagle;
public:
  static double *indF;
  void run(funkyPars  *pars);
  void likeFreq(funkyPars *pars,freqStruct *freq);
  void postFreq(funkyPars *pars,freqStruct *freq);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void openfile(const char *outfiles);
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);
  abcFreq(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcFreq();
  static double emFrequencyNoFixed_F(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int posi);
  static double emFrequencyNoFixed(double *loglike,int numInds, int iter,double start,int *keep,int keepInd,int major,int posi);
  static double emFrequencyNoFixed_ext(double *loglike,int numInds,int *keep,int keepInd,int major,int posi);
  static double likeNoFixedMinor(double p,double *logLikes,int numInds,int major);
  static double emFrequency_F(double *loglike,int numInds, int iter,double start,int *keep,int keepInd);
  static double emFrequency(double *loglike,int numInds, int iter,double start,int *keep,int keepInd);
  static double emFrequency_ext(double *loglike,int numInds,int *keep,int keepInd);
  static double likeFixedMinor(double p,double *logLikes,int numInds);
};


