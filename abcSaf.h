
typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  float **pLikes;
}realRes;


class abcSaf : public abc{
  int doSaf;
  BGZF *outfileGprobs;
  BGZF *outfileSAF;
  FILE *outfileSAFIDX;
  BGZF *outfileSAFPOS;
  int underFlowProtect;
  int isSim;
  int noTrans;
  char *anc;
  char *pest;
  int doPost;
  int tsktsktsk;
  int isHap;
  double *filipeIndF;
  int ishap;
  int newDim;
  int64_t offs[2];
  int nnnSites;
  char *tmpChr;
  void algoJointPost(double **post,int nSites,int nInd,int *keepSites,realRes *r);
  void algoGeno(int refId,double **liks,char *major,char *minor,int nsites,int numInds,kstring_t *sfsfile,int underFlowProtect, int *posi,int *keepSites,double *pest);
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int *keepSites,realRes *r,int noTrans);
  void algoJointHap(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int *keepSites,realRes *r,int noTrans);
  void algoJointMajorMinor(double **liks,int nsites,int numInds, int *keepSites,realRes *r, char *major, char *minor);

  void writeAll();
  int mynchr;
public:
  static double *prior; //<- outputfile form pest;
  static double *lbicoTab; //dim = [2*numInds+1]
  static double **myComb2Tab;
  //none optional stuff
  FILE *outfile;
  abcSaf(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcSaf();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  void changeChr(int refId);
  
};
