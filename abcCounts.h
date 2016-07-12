
class abcCounts : public abc{
private:
  const char* postfix1; //.pos
  const char* postfix2;//.counts;
  const char* postfix3;//.qs
  const char* postfix4;//.depthSample
  const char* postfix5;//.depthGlobal
  
  char *qfileFname;
  char *ffileFname;
  char *oFiles;
  int dumpCounts;
  int iCounts;
  int doCounts;
  int doQsDist;//0=nothing 1=overall 
  int setMaxObs;
  char *minQfile;
  angsd::Matrix<double> minQmat;


  kstring_t bpos;
  kstring_t bbin;
  kstring_t bufstr;
  BGZF* oFileCountsBin;
  BGZF* oFileCountsPos;
  BGZF* oFileIcounts;

  size_t *qsDist;
  int nInd;
  int minInd;
  int minQ;
  int setMaxDepth;
  int setMinDepth;
  int setMaxDepthInd;
  int setMinDepthInd;

  //from depth class
  int doDepth;
  int maxDepth;

  size_t **depthCount;
  size_t *globCount;
  unsigned char lookup[256];
  int qCutoff[9];
  double fCutoff[9];
public:

  
  //none optional stuff
  suint **countNucs(const chunkyT *chk,int *keepSites,int mmin,int mmax);  
  abcCounts(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcCounts();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};
