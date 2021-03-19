
class abcMcall:public abc{
private:
public:
  int domcall;
  
  abcMcall(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcMcall();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
