
typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  double **pLikes;
}realRes;


class abcSaf : public abc{
  int doSaf;
  gzFile outfileGprobs;
  FILE *outfileSFS;
  gzFile outfileSFSPOS;
  gzFile theta_fp;
  int underFlowProtect;
  int fold;
  int isSim;
  int noTrans;
  char *anc;
  char *pest;
  double *prior; //<- outputfile form pest;
  int doThetas;
  void calcThetas(funkyPars *p,int index,double *prior,gzFile fpgz);

  double aConst;
  double aConst2;
  double aConst3;
  double *scalings;


  double *filipeIndF;

  int newDim;
  void algoJointPost(double **post,int nSites,int nInd,int *keepSites,realRes *r,int doFold);
  void algoGeno(int refId,double **liks,char *major,char *minor,int nsites,int numInds,kstring_t *sfsfile,int underFlowProtect, int *posi,int *keepSites,double *pest);
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes *r,int noTrans);
  double *lbicoTab; //dim = [2*numInds+1]
  double **myComb2Tab;
public:
  //none optional stuff
  FILE *outfile;
  abcSaf(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcSaf();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

  
};
