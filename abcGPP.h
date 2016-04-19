
class abcGPP:public abc{
private:
public:
  int doGPP;
  
  abcGPP(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcGPP();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
};
