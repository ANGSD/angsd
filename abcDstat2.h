class abcDstat2:public abc{
private:
  kstring_t bufstr;
  int currentChr;
  int NbasesPerLine;
  int nBlocks;
  int block; //the current block number
  int blockSize;
  double DENprint;
  double NUMprint;
  int NSITEprint;
  int Eprint;
  char *ancName;
  double *COMBprint;
  double *ALLCOMBprint;
public:
  int doAbbababa2;
  FILE *outfile;
  FILE *outfile2;
  int sample;
  int doCount;
  int maxDepth;
  int sizeH1;
  int sizeH2;
  int sizeH3;
  int sizeH4;
  int enhance;
  //int nComb;
  int rmTrans;
  int Aanc;
  int combFile;
  abcDstat2(const char *outfiles, argStruct *arguments,int inputtype);
  ~abcDstat2();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);  //not protected
  void print(funkyPars *pars); // protect (MUTEX)
  void clean(funkyPars *pars); //
  void printArg(FILE *argFile);
  void printAndEmpty(int blockStart,int theChr);
  void getBlockNum(int pos);
};
