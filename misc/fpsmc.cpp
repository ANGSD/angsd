#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"

void ComputeSW(int maxTime,double W[],double sW[]){
  double tmp = 0;
  for (unsigned i = maxTime; i >=0 ; i--){
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




using namespace std;

struct Site{
  double g0, g1;
};

struct wins{
  int from;//inclusive
  int to;//inclusive
};

class fastPSMC {
  unsigned maxTime; //number of time intervals
  double rho;
  double mu;
  std::vector<double> timePoints;//coalescence times. Called T_k in document
  std::vector<double> P1, P2, P3, P4, P5, P6, P7;//each has length of timePoints.size()
  std::vector<double> R1, R2; ////each has length of timePoints.size()
  std::vector<double> epSize; //effective population size, S_k in latex
  //  std::vector<double> fbProb; //length = timePoints.size()
  double **fw;//maxTime x nWindows+1
  double **bw;//maxTime x nWindows+1
  double **pp;
public:
  fastPSMC(){
    rho = 0.207;
    mu = 0.0001;
    maxTime=100;
    fprintf(stderr,"\t-> rho:%f maxTime:%d mu:%f\n",rho,maxTime,mu);
  }
  fastPSMC(double rho_a,double mu_a,int maxTime_a){
    rho=rho_a;
    mu=mu_a;
    maxTime=maxTime_a;
    fprintf(stderr,"\t-> rho:%f maxTime:%d mu:%f\n",rho,maxTime,mu);
  }
  void init(size_t numWin);
  double llh(std::vector<Site> &data);
	//Initialisation
  void SetTimeIntervals(){
    for (unsigned i = 0; i < maxTime+1; i++){
      timePoints.push_back( i );
      //prob1.push_back( 0 );
    }
  }
  
  void SetEPSize(){
    epSize.clear();
    for (unsigned i = 0; i < maxTime; i++)
      epSize.push_back( 1 );
  }

  /*
    initial distribution:
    exp(-sum_{<=k-1} lambda_i*(T_i-T_{i-1})))(1-exp(lambda_k(T_{k-1}-T_k)))
  */
  
  
  void ComputeFBProbs(std::vector<Site> &data,std::vector<wins> &windows){
    fprintf(stderr,"[%s] start\n",__FUNCTION__ );
    //we first set the initial fwprobs to uniform
    for(int i=timePoints.size();i>0;i++){
      double tmp =0;
      for(int t=timePoints.size()-1;t>=0;t--)
	tmp += 1;//epsize[t];
      fw[i][i] = exp(-tmp)*(1-exp(epSize[i]));
    }

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      
      //calculate emissions
      ComputeRs(v,fw);//<-prepare R1,R2

      //Calcute emission probs, looping over all sites from: windows[v].from to: windows[v].to
      //we are calculating emission probs for each T_k?!
      double logemis[maxTime];

      for(int j=0;j<maxTime;j++){
	logemis[j] = 0;
	for(int i=windows[v].from;i<windows[v].to;i++){
	  logemis[j] += log(exp(data[i].g0)*exp(-timePoints[j]*mu) + exp(data[i].g1)*(1-exp(-timePoints[j]*mu) ));
	}
      }
      fw[v+1][0] = (fw[v][0]*P1[0] + R1[0]*P3[0] + fw[v][0]*P4[0])*logemis[0] ;
      for (unsigned i = 1; i < maxTime; i++)
	fw[v+1][i] = (fw[v+1][i]*P1[i] + R2[i-1]*P2[i-1] + R1[i]*P3[i] + fw[v+1][i]*P4[i])*logemis[i];

      
      
      
    }
    

    //now do backward algorithm
    for(int i=0;i<timePoints.size();bw++)
      bw[i][windows.size()] = 1.0;

    //we plug in values at v-1, therefore we break at v==1
    for(int v=windows.size();v>0;v--){
      
      //calculate emissions
      
      
      ComputeRs(v,bw);//<-prepare R1,R2

      //Calcute emission probs, looping over all sites from: windows[v].from to: windows[v].to
      //we are calculating emission probs for each T_k?!
      double logemis[maxTime];

      for(int j=0;j<maxTime;j++){
	logemis[j] = 0;
	for(int i=windows[v].from;i<windows[v].to;i++){
	  logemis[j] += log(exp(data[i].g0)*exp(-timePoints[j]*mu) + exp(data[i].g1)*(1-exp(-timePoints[j]*mu) ));
	}
      }
      
      bw[v-1][0] = (bw[v][0]*P1[0] + R1[0]*P3[0] + bw[v][0]*P4[0])*logemis[0] ;
      for (unsigned i = 1; i < maxTime; i++)
	bw[v-1][i] = (bw[v][i]*P1[i] + R2[i-1]*P2[i-1] + R1[i]*P3[i] + bw[v][i]*P4[i])*logemis[i];
    }

    //calculate post prob per site.
    for(int v=1;v<windows.size();v++){
      for(int j=0;j<maxTime;j++)
	pp[v][j] = fw[v][j]*bw[v][j];
	
    }
    
    fprintf(stderr,"[%s] stop\n",__FUNCTION__ );
  }  
private:
  double timePoint(){

  }
  
  void ComputeP1(){ //TODO: maxTime is in fact the number of time intervals
    for (unsigned i = 0; i < maxTime-1; i++){
      P1[i] = 1.0/(1.0+epSize[i]*2.0*rho);
      P1[i] *= exp( -rho*2.0*timePoints[i] ) - exp(-rho*2.0*timePoints[i+1]-(timePoints[i+1]-timePoints[i])/epSize[i]);
      P1[i] /= 1.0 - exp( -(timePoints[i+1]-timePoints[i])/epSize[i] );
    }
    //Last interval ends with +infinity
    unsigned i = maxTime - 1;
    P1[i] = 1.0/(1.0+epSize[i]*2.0*rho)* exp( -rho*2.0*timePoints[i] );
    
  }
  void ComputeP2(){
    for (unsigned i = 0; i < maxTime; i++)
      P2[i] = 1.0 - P5[i];
  }
  
  void ComputeP3(){
    for (unsigned i = 0; i < maxTime - 1; i++){
      P3[i] = exp(-timePoints[i]*2.0*rho);
      P3[i] += epSize[i]*2.0*rho/(1.0 - epSize[i]*2.0*rho)*exp(-(timePoints[i+1]-timePoints[i])/epSize[i]-timePoints[i]*2.0*rho);
      P3[i] -= 1.0/(1.0 - epSize[i]*2.0*rho)*exp(-timePoints[i+1]*2.0*rho);
    }
    unsigned i = maxTime - 1;
    P3[i] = exp(-timePoints[i]*2.0*rho);
  }
	
  void ComputeP4(){
    for (unsigned i = 0; i < maxTime-1; i++){
      P4[i] = 1.0/(1.0 - exp(-(timePoints[i+1]-timePoints[i])/epSize[i]) );
      double tmp = 2.0*rho/(1.0 + 2*rho*epSize[i])*exp(-2*rho*timePoints[i]);
      tmp -= 2.0*exp(-(timePoints[i+1] - timePoints[i])/epSize[i] - 2.0*rho*timePoints[i] );
      tmp -= 2.0*rho*epSize[i]/(1.0 - epSize[i]*2.0*rho)*exp(-2.0*rho*timePoints[i]-2.0*(timePoints[i+1]-timePoints[i])/epSize[i]);
      tmp += 2.0/(1.0-epSize[i]*2.0*rho)/(1.0 + 2.0*rho)*exp(-rho*timePoints[i+1]-(timePoints[i+1]-timePoints[i])/epSize[i]);
      P4[i] *= tmp;
    }
    unsigned i = maxTime - 1;
    P4[i] = 2.0*rho/(1.0 + 2.0*rho*epSize[i])*exp(-2.0*rho*timePoints[i]);
  }
	
  void ComputeP5(){
    for (unsigned i = 0; i < maxTime-1; i++)
      P5[i] = exp( -(timePoints[i+1] - timePoints[i])/epSize[i] );
    P5[maxTime-1] = 0.0;
  }
	
  void ComputeP6(){
    for (unsigned i = 0; i < maxTime-1; i++){
      P6[i] = 1/(1-exp(-(timePoints[i+1]-timePoints[i])/epSize[i]));
      P6[i] *= exp(-(timePoints[i+1]-timePoints[i])/epSize[i]);
      double tmp = exp(-2*rho*timePoints[i]);
      tmp -= 1/(1-2*rho*epSize[i])*exp(-2*rho*timePoints[i+1]);
      tmp += 2*rho*epSize[i]/(1 - 2*rho*epSize[i])*exp(-2*rho*timePoints[i]-(timePoints[i+1]-timePoints[i])/epSize[i]);
      P6[i] *= tmp;
    }
    P6[maxTime - 1] = 0.0;
  }
	
  void ComputeP7(){
    for (unsigned i = 0; i < maxTime - 1; i++){
      P7[i] = 1.0 - exp(-(timePoints[i+1]-timePoints[i])*2.0*rho) - exp(-timePoints[i]*2.0*rho);
      P7[i] -= epSize[i]*2*rho/(1 - epSize[i]*2.0*rho)*exp(-(timePoints[i+1]-timePoints[i])/epSize[i]-timePoints[i]*2.0*rho);
      P7[i] += 1.0/(1.0 - epSize[i]*2.0*rho)*exp(-timePoints[i]*2.0*rho);
    }
    unsigned i = maxTime - 1;
    P7[i] = 1.0 - exp(-2.0*rho*timePoints[i]);
  }

  void ComputeR1(int v,double **mat){
    double tmp = 0;
    for (unsigned i = maxTime - 1; i > 0 ; i--){
      R1[i] = tmp + mat[i][v];
      tmp = R1[i];
    }
  }
	
  void ComputeR2(int v,double **mat){
    double tmp = 0;
    for (unsigned i = 0; i < maxTime ; i++){
      R2[i] = tmp*P2[i]+mat[v][i]*P6[i]+R1[i]*P7[i];
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
  }

};

void fastPSMC::init(size_t numWindows){
  SetTimeIntervals();
  SetEPSize();

  fw = new double *[timePoints.size()];
  bw = new double *[timePoints.size()];
  pp = new double *[timePoints.size()];
  for(int i=0;i<timePoints.size();i++){
    pp[i] = new double[numWindows+1];
    pp[i] = new double[numWindows+1];
  }
  maxTime++;
  P1.reserve(maxTime);
  P2.reserve(maxTime);
  P2.reserve(maxTime);
  P3.reserve(maxTime);
  P4.reserve(maxTime);
  P5.reserve(maxTime);
  P6.reserve(maxTime);
  P7.reserve(maxTime);
  R1.reserve(maxTime);
  R2.reserve(maxTime);
  maxTime--;
  ComputeGlobalProbabilities();
}

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

//function to print the data we need
int main_analysis(std::vector<Site> &data,std::vector<wins> &windows){
  fastPSMC obj;
  //prepare datastructures
  obj.init(windows.size());
  
  //print indices for endpoint of windows
  for(int w=0;0&w<windows.size();w++)
    fprintf(stdout,"win[%d]=(%d,%d)\n",w,windows[w].from,windows[w].to);

  //print out data:
  for(int w=0;0&w<windows.size();w++)
    for(int p=windows[w].from;p<windows[w].to;p++)
      fprintf(stdout,"%d\t%d\t%f\t%f\n",p,w,data[p].g0,data[p].g1);


  obj.ComputeFBProbs(data,windows);
}


int psmc_wrapper(args *pars){
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);



  //loop over chrs;

  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    //set perchr iterator, if pars->chooseChr, then we have only use a single chr
    std::vector<Site> data;//really stupid, but lets go with it for now.
    std::vector<wins> windows;
    it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);

    //generate window data
    int beginIndex = pars->perc->first;
    while(1){
      int beginPos=pars->perc->pos[beginIndex]+1;//add one due to positions being offset with one.
      int endPos = beginPos+pars->winSize;

      int at=beginIndex;
      while(at<pars->perc->last&&pars->perc->pos[at]<endPos){
      Site d;
	d.g0 = pars->perc->gls[at];
	d.g1 = pars->perc->gls[at+1];
	data.push_back(d);
	at++;
      }
      if(at>=pars->perc->last)
	break;
      wins w;
      w.from=beginIndex;
      w.to=at;
      windows.push_back(w);
      beginIndex = at;
      //      fprintf(stderr,"(%d,%d)\n",beginPos,endPos);
    }
    main_analysis(data,windows);
    break;
    //print_psmc_print_windows(data);
    /*
    for(size_t s=pars->perc->first;0&&s<pars->perc->last;s++){
      fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->pos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
    }
    */
    
    if(pars->chooseChr!=NULL)
      break;
  }

}


