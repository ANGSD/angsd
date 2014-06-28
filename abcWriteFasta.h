
class abcWriteFasta:public abc{
private:
  kstring_t bufstr;
  size_t currentPos;
  int currentChr;
  int NbasesPerLine;
  double *lphred;//these are log phread scores log(10^(1:255)/(-10))
public:
  int doFasta;
  gzFile outfileZ;
  int doCount;

  abcWriteFasta(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWriteFasta();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};
