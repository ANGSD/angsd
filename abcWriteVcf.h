class abcWriteVcf:public abc{
private:
  gzFile fp;//tfam/fam/bim
  kstring_t kstr;
public:
  int doVcf;
  abcWriteVcf(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcWriteVcf();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
