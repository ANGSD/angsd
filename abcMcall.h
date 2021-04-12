typedef struct{
  float *quals;//phred of hyp(alt) hyp(nul) (2xnumsites)
  float *QS;//sum of phreds for each allele 5*numsites
  int **gcdat;
  char *isvar;
  /* contains -1,0,1,2.
     -1 no call
     0 AA
     1 Aa
     2 aa
  */
  char *als;
}angsd_mcall;


class abcMcall:public abc{
private:
public:
  int domcall;
  
  abcMcall(const char *outfiles,argStruct *arguments,int inputtype);
  ~abcMcall();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
