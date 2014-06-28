class abcWritePlink:public abc{
private:
  FILE *fp1;//tfam/fam/bim
  FILE *fp2;//tped/bed
public:
  int doPlink;
  abcWritePlink(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWritePlink();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
