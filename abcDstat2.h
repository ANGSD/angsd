class abcDstat2:public abc{
private:
  kstring_t bufstr;
  int currentChr;
  int NbasesPerLine;
  int nBlocks;
  int block;
  int blockSize;
  double *DENprint;
  double *NUMprint;
  double *NSITEprint;
  int Eprint;
  char *ancName;
  double **COMBprint;
  
public:
  int doAbbababa2;
  FILE *outfile;
  int sample;
  int doCount;
  int useLast;
  int maxDepth;
  int enhance;
  int nIndFasta;
  int rmTrans;
  int Aanc;
  int *POPSIZE;
  int *CUMPOPSIZE;
  int **SIZEABCD;
  char *sizeFile;
  long int numComb;
  int numPop;


  angsd::Matrix<int> sizeMat;
  
  abcDstat2(const char *outfiles, argStruct *arguments,int inputtype);
  ~abcDstat2();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);  //not protected
  void print(funkyPars *pars); // protect (MUTEX)
  void clean(funkyPars *pars); //
  void printArg(FILE *argFile);
  void printAndEmpty(int blockAddress,int theChr);
  void getBlockNum(int pos);
  int getNumBlocks(funkyPars *pars);
};
