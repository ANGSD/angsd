#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"
#include "compute.c"
#define PSMC_T_INF 1000.0

typedef struct{
  double **P;
  double **PP;
}oPars;



double qkFunction(unsigned K, double pix, unsigned numWind,double **P,double **PP){

  /*
  for (unsigned l = 1; l < numWind + 1; l++) //This block is needed if eimission probabilities depend on estimated parameters, e.g. on time disctretisation 
    qi += log(emis[K][l])*fw[K][l]*bw[K][l];
    qi /= pix;
  */
  double qi = 0;
  //qi += nP1[i]*PP[1][i] + nP2[i]*PP[2][i] + nP3[i]*PP[3][i] + nP4[i]*PP[4][i] + nP5[i]*PP[5][i] + nP6[i]*PP[6][i] + nP7[i]*PP[7][i];
  return qi;
}


/*
  objective function. Function to be optimized


double qFunction(double *params , void *d){
  oPars *data = (oPars*) d;

  double Q = 0;
  //Compute new values of probabilities P_i with parameters params
  cout << "COMPUTE THEM HERE" << enld;
  //Compute P_ii
  ComputePii(numWind);//pull before optimization
  
  for (unsigned i = 0; i < maxTime; i++)
    Q += qkFunction(i, pix);
  return Q;
}

*/
struct wins{
  int from;//inclusive
  int to;//inclusive
};


double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;
  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


void printarray(FILE *fp,double *ary,int l){
  for(int i=0;i<l;i++)
    fprintf(fp,"%d)\t%f\n",i,ary[i]);

}
void printmatrix(FILE *fp,double **mat,int x,int y){
  fprintf(fp,"#printmatrix with x:%d y:%d\n",x,y);
  for(int i=0;i<y;i++){
    for(int j=0;j<x-1;j++)
      fprintf(fp,"%f\t",mat[j][i]);
    fprintf(fp,"%f\n",mat[x-1][i]);
  }

}


void ComputeSW(int maxTime,double W[],double sW[]){
  double tmp = 0;
  for (int i = maxTime; i >=0 ; i--){
    tmp += W[i];
    sW[i] = tmp;
  }
}

void UpdateEPSize(int maxTime, double W[],double sW[],double epSize[],double T[]){
  for (unsigned i = 1; i < maxTime; i++){
    double tmp;
    tmp = W[i]/sW[i];
    epSize[i] = -log(1 - tmp)/(T[i+1]-T[i]);
  }
}


/*
  Calculate stationary distrubution
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  stationary distribution will be put in results array, also of length tk
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */
void calculate_stationary(double *tk,int tk_l,double *lambda,double *results,double *P2){
  results[0] = 1;//fix this
  for(int i=1;i<tk_l;i++){
    double tmp =0;
    for(int j=0;j<i-1;j++)
      tmp += (tk[j+1]-tk[j])/lambda[i];
    results[i] = tmp*P2[i];
  }

}

/* 
  Calculate emission probabilityes
  tk array of length tk_l
  lambda array effective population sizes
  both has length tk_l
  
  emission probablities will be put in the **emis
  
  stationary(i) = exp(-sum_{j=0}^{i-1}{tau_j/lambda_j}*P2[i])
 */
void calculate_emissions(double *tk,int tk_l,double *gls,std::vector<wins> &windows,double mu,double **emis){
  fprintf(stderr,"[Calculating emissions with tk_l:%d and windows.size():%lu:%s ] start\n",tk_l,windows.size(),__TIME__);
  //initialize the first:
  for(int j=0;j<tk_l;j++)
    emis[j][0] = 0;
 

  for(int v=0;v<windows.size();v++){
    for(int j=0;j<tk_l;j++){
      emis[j][v+1] = 0;
      double inner = exp(-2.0*tk[j]*mu);
      for(int i=windows[v].from;i<windows[v].to;i++)
	emis[j][v+1] += log(exp(gls[i*2])*inner + exp(gls[2*i+1])*(1-inner));

    }
  }
  fprintf(stderr,"[Calculating emissions with tk_l:%d and windows.size():%lu:%s ] stop\n",tk_l,windows.size(),__TIME__);
}


/*
  This functions either set the tk, NOT the intervals.
  n, is the true length of tk. First entry zero, last entry INF
 */
void setTk(int n, double *t, double max_t, double alpha, double *inp_ti){
  //  fprintf(stderr,"[%s] (n,tk,max_t,alpha,inp_ti)=(%d,%p,%f,%f,%p)\n",__FUNCTION__,n,t,max_t,alpha,inp_ti);
  int k;
  if (inp_ti == 0) {
    double beta;
    beta = log(1.0 + max_t / alpha) / n; // beta controls the sizes of intervals
    for (k = 0; k < n; ++k)
      t[k] = alpha * (exp(beta * k) - 1);
    t[n-1] = max_t;
    t[n] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
  } else {
    memcpy(t, inp_ti, (n-1) * sizeof(double));
    t[n-1] = PSMC_T_INF;
  }
}


void setEPSize(double *ary,int l,double *from_infile){
  if(!from_infile)
    for (int i = 0; i <l; i++)
      ary[l]=1;
  else
    memcpy(ary,from_infile,(l-1)*sizeof(double));
	
}
/*
  pix = 0;
  for (unsigned i = 0; i < maxTime; i++)
  pix += fw[i][last_window]; 
  

 */


class fastPSMC {
  double pix;
  int tk_l;
  double max_t;
  //  unsigned maxTime; //number of time intervals
  double rho;
  double mu;
  double *tk;//tk is tk_l long
  double *epsize;//tk_l long
  double **P;
  double **PP;
  double *stationary;//tk_l long
  //  std::vector<double> P1, P2, P3, P4, P5, P6, P7;//each has length of timePoints.size()
  std::vector<double> R1, R2; ////each has length of timePoints.size()
  //  std::vector<double> epSize; //effective population size, S_k in latex
  //  std::vector<double> fbProb; //length = timePoints.size()
  double **fw;//tk_l x nWindows+1
  double **bw;//tk_l x nWindows+1
  double **pp;//tk_l x nWindows+1
  double **emis;//tk_l x nWindows+1
public:
  fastPSMC(){
    pix = -666;
    max_t = 15;
    rho = 0.207;
    mu = 0.0001;
    fprintf(stderr,"\t-> rho:%f mu:%f\n",rho,mu);
  }

  void init(int numWin,psmc_par *p);
 
  void ComputePii(unsigned numWind){
    ComputeP11(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP22(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP33(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP44(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP55(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP66(numWind,tk_l,P,PP,fw,bw,stationary);
    ComputeP77(numWind,tk_l,P,PP,fw,bw,stationary);
  }

  void ComputeFBProbs(double *gls,std::vector<wins> &windows,int n){
    fprintf(stderr,"[%s] start\n",__FUNCTION__ );

    //calculate emissions
    calculate_emissions(tk,tk_l,gls,windows,mu,emis);exit(0);
    //    printmatrix(stdout,emis,tk_l,(int)windows.size());exit(0);
    //we first set the initial fwprobs to stationary distribution
    double stationary[tk_l];
    calculate_stationary(tk,tk_l,epsize,stationary,P[2]);
    for(int i=0;i<tk_l;i++)
      fw[i][0] =stationary[i];

  

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      ComputeRs(v,fw);//<-prepare R1,R2
      fw[0][v+1] = (fw[0][v]*P[1][0] + R1[0]*P[3][0] + fw[0][v]*P[4][0])*emis[0][v+1] ;
      for (unsigned i = 1; i < tk_l; i++)
	fw[i][v+1]= (fw[i][v+1]*P[1][i] + R2[i-1]*P[2][i-1] + R1[i]*P[3][i] + fw[i][v+1]*P[4][i])*emis[i][v+1];
    }
    double tmp =0;
    for(int i=0;i<tk_l;i++){
      fprintf(stderr,"fw[%d][%lu]:%f\n",i,windows.size(),fw[i][windows.size()]);
      tmp += log(fw[i][windows.size()]);
    }
    fprintf(stderr,"forward llh:%f\n",tmp);
    exit(0);

    pix = 0;
    for (unsigned i = 0; i < tk_l; i++)
      pix += fw[i][windows.size()-1]; 
 


    //now do backward algorithm
    for(int i=0;i<n;i++)
      bw[i][windows.size()] = stationary[i];

    //we plug in values at v-1, therefore we break at v==1
    for(int v=windows.size();v>0;v--){
      ComputeRs(v,bw);//<-prepare R1,R2
      bw[0][v-1] = (bw[0][v]*P[1][0] + R1[0]*P[3][0] + bw[0][v]*P[4][0])*emis[0][v] ;
      for (unsigned i = 1; i < tk_l; i++)
	bw[i][v-1] = (stationary[i]*bw[i][v]*emis[i][v]*P[1][i] + R2[i-1]*P[2][i-1] + R1[i]*P[3][i] + stationary[i]*bw[i][v]*emis[i][v]*P[4][i])/stationary[i];
    }
    
    tmp =0;
    for(int i=0;i<tk_l;i++)
      tmp += log(bw[i][windows.size()]);
    fprintf(stderr,"backward llh:%f\n",tmp);


    //calculate post prob per window per state
    for(int v=1;v<windows.size();v++){
      for(int j=0;j<tk_l;j++)
	pp[j][v] = fw[j][v]*bw[j][v];
	
    }
    
    fprintf(stderr,"[%s] stop\n",__FUNCTION__ );
  }  
private:
  void ComputeP1(){ //TODO: tk_l is in fact the number of time intervals
    for (unsigned i = 0; i < tk_l-1; i++){
      P[1][i] = 1.0/(1.0+epsize[i]*2.0*rho);
      P[1][i] *= exp( -rho*2.0*tk[i] ) - exp(-rho*2.0*tk[i+1]-(tk[i+1]-tk[i])/epsize[i]);
      P[1][i] /= 1.0 - exp( -(tk[i+1]-tk[i])/epsize[i] );
    }
    //Last interval ends with +infinity
    unsigned i = tk_l - 1;
    P[1][i] = 1.0/(1.0+epsize[i]*2.0*rho)* exp( -rho*2.0*tk[i] );
    
  }
  void ComputeP2(){
    for (unsigned i = 0; i < tk_l; i++)
      P[2][i] = 1.0 - P[5][i];
  }
  
  void ComputeP3(){
    for (unsigned i = 0; i < tk_l - 1; i++){
      P[3][i] = exp(-tk[i]*2.0*rho);
      P[3][i] += epsize[i]*2.0*rho/(1.0 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
      P[3][i] -= 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i+1]*2.0*rho);
    }
    unsigned i = tk_l - 1;
    P[3][i] = exp(-tk[i]*2.0*rho);
  }
	
  void ComputeP4(){
    for (unsigned i = 0; i < tk_l-1; i++){
      P[4][i] = 1.0/(1.0 - exp(-(tk[i+1]-tk[i])/epsize[i]) );
      double tmp = 2.0*rho/(1.0 + 2*rho*epsize[i])*exp(-2*rho*tk[i]);
      tmp -= 2.0*exp(-(tk[i+1] - tk[i])/epsize[i] - 2.0*rho*tk[i] );
      tmp -= 2.0*rho*epsize[i]/(1.0 - epsize[i]*2.0*rho)*exp(-2.0*rho*tk[i]-2.0*(tk[i+1]-tk[i])/epsize[i]);
      tmp += 2.0/(1.0-epsize[i]*2.0*rho)/(1.0 + 2.0*rho)*exp(-rho*tk[i+1]-(tk[i+1]-tk[i])/epsize[i]);
      P[4][i] *= tmp;
    }
    unsigned i = tk_l - 1;
    P[4][i] = 2.0*rho/(1.0 + 2.0*rho*epsize[i])*exp(-2.0*rho*tk[i]);
  }
	
  void ComputeP5(){
    for (unsigned i = 0; i < tk_l-1; i++)
      P[5][i] = exp( -(tk[i+1] - tk[i])/epsize[i] );
    P[5][tk_l-1] = 0.0;
  }
	
  void ComputeP6(){
    for (unsigned i = 0; i < tk_l-1; i++){
      P[6][i] = 1/(1-exp(-(tk[i+1]-tk[i])/epsize[i]));
      P[6][i] *= exp(-(tk[i+1]-tk[i])/epsize[i]);
      double tmp = exp(-2*rho*tk[i]);
      tmp -= 1/(1-2*rho*epsize[i])*exp(-2*rho*tk[i+1]);
      tmp += 2*rho*epsize[i]/(1 - 2*rho*epsize[i])*exp(-2*rho*tk[i]-(tk[i+1]-tk[i])/epsize[i]);
      P[6][i] *= tmp;
    }
    P[6][tk_l - 1] = 0.0;
  }
	
  void ComputeP7(){
    for (unsigned i = 0; i < tk_l - 1; i++){
      P[7][i] = 1.0 - exp(-(tk[i+1]-tk[i])*2.0*rho) - exp(-tk[i]*2.0*rho);
      P[7][i] -= epsize[i]*2*rho/(1 - epsize[i]*2.0*rho)*exp(-(tk[i+1]-tk[i])/epsize[i]-tk[i]*2.0*rho);
      P[7][i] += 1.0/(1.0 - epsize[i]*2.0*rho)*exp(-tk[i]*2.0*rho);
    }
    unsigned i = tk_l - 1;
    P[7][i] = 1.0 - exp(-2.0*rho*tk[i]);
  }
  
  void ComputeP0(){ //probability P(T > i)
    P[0][0] = P[5][0];
    for (unsigned i = 1; i < tk_l; i++)
      P[0][i] = P[0][i-1]*P[5][i];
  }

  void ComputeR1(int v,double **mat){
    R1[tk_l - 1] = 0;
    for (int i = tk_l - 2; i >= 0 ; i--)
      R1[i] = R1[i+1] + mat[i+1][v];

  }
	
  void ComputeR2(int v,double **mat){
    //    fprintf(stderr,"tk_l:%d\n",tk_l);
    double tmp = 0;
    for (unsigned i = 0; i < tk_l ; i++){
      R2[i] = tmp*P[2][i]+mat[i][v]*P[6][i]+R1[i]*P[7][i];
      tmp = R2[i];
    }
  }
	
  void ComputeRs(int v,double **mat){
    //    fprintf(stderr,"v:%d\n",v);
    ComputeR1(v,mat);
    ComputeR2(v,mat);
  }
	
  void ComputeGlobalProbabilities(){
    ComputeP1();
    ComputeP5();
    ComputeP6();
    ComputeP2();
    ComputeP3();
    ComputeP4();
    ComputeP7();
    //    ComputePigi();
  }

};




void fastPSMC::init(int numWindows,psmc_par *p){
  fprintf(stderr,"[%s]: pars->n:%d will allocate tk with length n+2:%d\n",__FUNCTION__,p->n,p->n+2);
  tk_l = p->n+2;
  tk = new double[tk_l];
  setTk(tk_l,tk,max_t,0.01,p->times);//<- last position will be infinity
  //  printarray(stderr,tk,tk_l);
  epsize = new double[tk_l];
  setEPSize(epsize,tk_l,p->params);
  epsize[tk_l-1] = 100;//<- set the last element; DRAGON what should last epsize be?
  //  printarray(stderr,epsize,tk_l);exit(0);
  //  exit(0);
  fw = new double *[tk_l];
  bw = new double *[tk_l];
  pp = new double *[tk_l];
  emis = new double *[tk_l];
  for(int i=0;i<tk_l;i++){
    emis[i] = new double[numWindows+1];
    fw[i] = new double[numWindows+1];
    bw[i] = new double[numWindows+1];
    pp[i] = new double[numWindows+1];

  }
  fprintf(stderr,"emission allocated with [%d][%d]\n",tk_l,numWindows+1);
  P = new double *[8];
  PP= new double *[8];
  for(int i=0;i<8;i++){
    P[i] = new double[tk_l];
    PP[i]= new double[tk_l];
  }
  ComputeGlobalProbabilities();//only the P* ones
}

//function to print the data we need
int main_analysis(double *gls,std::vector<wins> &windows,psmc_par *p){

  fastPSMC obj;
  //prepare datastructures
  obj.init(windows.size(),p);

#if 0
  //print indices for endpoint of windows
  for(int w=0;w<windows.size();w++)
    fprintf(stdout,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);
  exit(0);
  //print out data:
  for(int w=0;0&&w<windows.size();w++)
    for(int p=windows[w].from;p<windows[w].to;p++)
      fprintf(stdout,"%d\t%d\t%f\t%f\n",p,w,gls[2*p],gls[2*p+1]);
#endif
  obj.ComputeFBProbs(gls,windows,p->n);
  return 1;
}


int psmc_wrapper(args *pars,int block) {

  //  fprintf(stderr,"[%s]\n",__FUNCTION__);
#if 0 //print pars
  psmc_par *p=pars->par;
  fprintf(stderr,"par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;i<pars->par->n+1;i++)
    fprintf(stderr,"%i)\t%f\t%f\n",i,pars->par->times[i],pars->par->params[i]);
  
#endif
  
 

  //loop over chrs;

  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    //set perchr iterator, if pars->chooseChr, then we have only use a single chr
    std::vector<wins> windows;

    it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);

    int beginIndex =0;
    int endIndex=0;
    int beginPos = 1;
    int endPos = beginPos+block-1;

    while(1){
      wins w;
      if(endPos>pars->perc->pos[pars->perc->last-1])
	break;

      while(pars->perc->pos[beginIndex]<beginPos)
	beginIndex++;
      while(pars->perc->pos[endIndex]<endPos)
	endIndex++;
      //      fprintf(stdout,"endpiadsf:%d\n",pars->perc->pos[endIndex]);
      //      fprintf(stdout,"winsize:%d bp:%d,ep:%d bi:%d ei:%d ei-bi:%d\n",pars->winSize,beginPos,endPos,beginIndex,endIndex,endIndex-beginIndex);
      w.from = beginIndex;
      w.to = endIndex+1;
      windows.push_back(w);
      beginPos+=block;
      endPos+=block;

    }
    fprintf(stderr,"\t->[%s] number of windows:%lu\n",__FUNCTION__,windows.size());
    main_analysis(pars->perc->gls,windows,pars->par);
    break;

    for(size_t s=pars->perc->first;0&&s<pars->perc->last;s++)
      fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->pos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
    
    if(pars->chooseChr!=NULL)
      break;
  }
  return 1;
}


