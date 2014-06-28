
class abcHetPlas:public abc{
private:
  double **probs;
  FILE *outputFile;
  int maxIter;
  double minLRT;
public:

  int doHetPlas;
  int makellhs(tNode*,double**,int*);
  void doNew(funkyPars *pars);
  abcHetPlas(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcHetPlas();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
