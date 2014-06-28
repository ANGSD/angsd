
class abcDstat:public abc{
private:
  kstring_t bufstr;
  int currentChr;
  int NbasesPerLine;
  int nBlocks;
  int block; //the current block number
  int blockSize;
  int matcat[5][5][5][5];
  char *ancName;
  int *ABBA;
  int *BABA;
public:
  int doAbbababa;
  FILE *outfile;
  int doCount;
  int nComb;
  int rmTrans;
  abcDstat(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcDstat();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);
  void calcMatCat();
  void printAndEmpty();
  void getBlockNum(int pos);
};
