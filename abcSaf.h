
typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  float **pLikes;
}realRes;


class abcSaf : public abc{
  std::vector<float> *theta_res;
  std::vector<int> theta_pos;
  int doSaf;
  BGZF *outfileGprobs;
  BGZF *outfileSAF;
  FILE *outfileSAFIDX;
  BGZF *outfileSAFPOS;
  BGZF *theta_dat;
  FILE *theta_idx;
  int underFlowProtect;
  int fold;
  int isSim;
  int noTrans;
  char *anc;
  char *pest;
  int doPost;
  int doThetas;
  void calcThetas(funkyPars *p,int index,double *prior,std::vector<float> *vecs,std::vector<int> &myposi,int newdim);

  double aConst;
  double aConst2;
  double aConst3;
  double *scalings;
  int tsktsktsk;
  int isHap;
  double *filipeIndF;
  int ishap;
  int newDim;
  int64_t offs[2];
  int64_t offs_thetas;
  int nnnSites;
  char *tmpChr;
  void algoJointPost(double **post,int nSites,int nInd,int *keepSites,realRes *r,int doFold);
  void algoGeno(int refId,double **liks,char *major,char *minor,int nsites,int numInds,kstring_t *sfsfile,int underFlowProtect, int *posi,int *keepSites,double *pest);
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans);
  void algoJointHap(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans);
  void algoJointMajorMinor(double **liks,int nsites,int numInds, int *keepSites,realRes *r,int fold,char *major, char *minor);

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
