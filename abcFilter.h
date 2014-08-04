
typedef struct{
  int64_t offs;
  int len;
}asdf_dats;



typedef struct{
  BGZF *bg;
  FILE *fp;
  std::map<int,asdf_dats> offs;
  char *keeps;
  char *major;
  char *minor;
  int hasMajMin;
  int curLen;//<-length of the char arrays;
}filt;




class abcFilter : public abc{
  filt *fl;
  int doMajorMinor;
  char *fname;
  FILE *fp;

  char *keepsChr;//<- a char array of length=reflength, used for .keeps
  int minInd;
  //  int nInd;
  int curChr;
  int isBed;
public:
  void readSites(int refId);
  //none optional stuff
  FILE *outfile;
  abcFilter(argStruct *arguments);
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  ~abcFilter();
};
