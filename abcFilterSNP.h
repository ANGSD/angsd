

class abcFilterSNP:public abc{
  BGZF* outfileZ;
public:
  int doSnpStat;
  
  abcFilterSNP(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcFilterSNP();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
