
typedef struct{
  char *hasAlloced; 
  double **lh3; 
}lh3struct;




class abcMajorMinor:public abc{
  int doMajorMinor;
  int doSaf;
  char *pest;
  int skipTriallelic;
  void majorMinorGL(funkyPars *pars,int doMajorMinor);
  int doVcf;
  int doGlf;
  int rmTrans;
public:
  void printArg(FILE *argFile);
  void run(funkyPars *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  abcMajorMinor(const char *outfiles,argStruct *arguments,int inputtype);
  void getOptions(argStruct *arguments);
  ~abcMajorMinor();
};
