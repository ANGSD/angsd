/*
  Class that contains the perChr hmm

 */
double qkFunction(unsigned i, double pix, unsigned numWind,double **nP,double **PP);
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti);
void setEPSize(double *ary,int l,double *from_infile);
#define PSMC_T_INF 1000.0

struct wins{
  int from;//inclusive
  int to;//inclusive
};
void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double mu,double **emis);

class fastPSMC {
  double pix;
  int tk_l;
  double max_t;
  //  unsigned maxTime; //number of time intervals
  double rho;
  double mu;
  // double *tk;//tk is tk_l long
  //  double *epsize;//tk_l long
  double **P;
  double **PP;
  double *stationary,*R1,*R2;//tk_l long
  double **fw;//tk_l x nWindows+1
  double **bw;//tk_l x nWindows+1
  double **pp;//tk_l x nWindows+1
  double **emis;//tk_l x nWindows+1
  std::vector<wins> windows;
public:
  fastPSMC(){
    pix = -666;
    max_t = 15;
    rho = 0.207;
    mu = 0.0001;
    fprintf(stderr,"\t-> rho:%f mu:%f\n",rho,mu);
  }
  
  void setWindows(perpsmc *perc,char *chooseChr,int start,int stop,int block);
  void printWindows(FILE *fp){
    //print indices for endpoint of windows
    for(int w=0;w<windows.size();w++)
      fprintf(stdout,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);
  //print out data:
#if 0
  for(int w=0;0&&w<windows.size();w++)
    for(int p=windows[w].from;p<windows[w].to;p++)
      fprintf(stdout,"%d\t%d\t%f\t%f\n",p,w,gls[2*p],gls[2*p+1]);
#endif
  }
  void allocate(int tk_l);
  void ComputePii(unsigned numWind,int tk_l,double **P,double **PP,double **fw,double **bw,double *stationary);
  void calculate_stationary(double *tk,int tk_l,double *lambda,double *results,double *P2);
  void calculate_FW_BW_PP_Probs();
  void make_hmm(double *tk,int tk_l,double *gls,double *epsize);
  
private:
  void ComputeR2(int v,double **mat){
    //    fprintf(stderr,"tk_l:%d\n",tk_l);
    double tmp = 0;
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = tmp*P[2][i]+mat[i][v]*P[6][i]+R1[i]*P[7][i];
      tmp = R2[i];
    }
  }
  
  void ComputeR1(int v,double **mat){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + mat[i+1][v];
  }
  
  void ComputeRs(int v,double **mat){
    ComputeR1(v,mat);
    ComputeR2(v,mat);
  }

};
