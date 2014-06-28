
typedef struct{
  char *oklist;//<- {0,1,2}, length=numSites, 0 don't keep 1 do keep 2 error
  double **pLikes;
}realRes2;


class abcSaf2 : public abc{
  int doSaf2;
  FILE *outfileSFS;
  gzFile outfileSFSPOS;
  int underFlowProtect;
  int fold;
  int isSim;
  int noTrans;
  char *anc;

  int mode; //[EJ]
  int newDim;
  void algoJoint(double **liks,char *anc,int nsites,int numInds,int underFlowProtect, int fold,int *keepSites,realRes2 *r,int noTrans);
public:
  //none optional stuff
  FILE *outfile;
  abcSaf2(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcSaf2();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

  
};

