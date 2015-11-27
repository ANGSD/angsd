
class abcSnpTools:public abc{
private:
  uint16_t *ebd;
  int curChr;
  double *phred;
public:
  int doSnpTools;
  
  abcSnpTools(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcSnpTools();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
