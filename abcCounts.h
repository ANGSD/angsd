
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

  char *minQfile;
  angsd::Matrix<double> minQmat;


  kstring_t bpos;
  kstring_t bbin;

  gzFile oFileCountsBin;
  gzFile oFileCountsPos;
  gzFile oFileIcounts;

  size_t *qsDist;
  int nInd;
  int minInd;

  int setMaxDepth;
  int setMinDepth;
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
  suint **countNucs(const chunkyT *chk,int *keepSites);  
  abcCounts(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcCounts();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};
