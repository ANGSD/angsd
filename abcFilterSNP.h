

class abcFilterSNP:public abc{
  BGZF* outfileZ;
  double edge_pval;
  double mapQ_pval;
  double sb_pval;
  double hwe_pval;
  double qscore_pval;
  double hetbias_pval;
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
