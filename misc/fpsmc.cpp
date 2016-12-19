#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"

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
  std::vector<double> epSize; //effective population size
  //  std::vector<double> fbProb; //length = timePoints.size()
  double **fw;//maxTime x nWindows+1
  double **bw;//maxTime x nWindows+1
	
public:
  fastPSMC(){
    rho = 0.207;
    mu = 0.0001;
    maxTime=100;
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

  
  
  void ComputeFBProbs(std::vector<Site> &data,std::vector<wins> &windows){
    
    //we first set the initial fwprobs to uniform
    for(int i=0;i<timePoints.size();i++)
      fw[i][0] = 1.0/timePoints.size();

    //we now loop over windows.
    //v=0 is above and is the initial distribution, we therefore plug in at v+1
    for(int v=0;v<windows.size();v++){
      //calculate emissions
      
      
      ComputeRs(v);//<-prepare R1,R2

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
      for (unsigned i = 1; i < maxTime; i++){
	fw[v+1][i] = (fw[v+1][i]*P1[i] + R2[i-1]*P2[i-1] + R1[i]*P3[i] + fw[v+1][i]*P4[i])*logemis[i];
      }
    }
  }  
private:
  void ComputeP1(){
    for (unsigned i = 0; i < maxTime; i++){
      P1[i] = exp( -rho*1*(timePoints[i+1] - timePoints[i]) );
    }
  }
  void ComputeP2(){
    for (unsigned i = 0; i < maxTime; i++)
      P2[i] = 1 - P5[i];
  }
  
  void ComputeP3(){
    for (unsigned i = 0; i < maxTime; i++)
      P3[i] = exp(-2*rho*timePoints[i]) - exp(-2*rho*timePoints[i+1]) - P6[i];
  }
	
  void ComputeP4(){
    for (unsigned i = 0; i < maxTime; i++)
      P4[i] = P3[i];
  }
	
  void ComputeP5(){
    for (unsigned i = 0; i < maxTime; i++)
      P5[i] = exp( -epSize[i]*(timePoints[i+1] - timePoints[i]) );
  }
	
  void ComputeP6(){
    for (unsigned i = 0; i < maxTime; i++){
      double tmp = 0;
      tmp = exp((epSize[i]-2*rho)*timePoints[i]) - exp((epSize[i]-2*rho)*timePoints[i+1]);
      tmp = tmp*2*rho*exp(-epSize[i]*timePoints[i+1])/(epSize[i] - 2*rho);
      P6[i] = tmp;
    }
  }
	
  void ComputeP7(){
    for (unsigned i = 0; i < maxTime; i++)
      P7[i] = P6[i];
  }

  void ComputeR1(int v){
    double tmp = 0;
    for (unsigned i = maxTime - 1; i > 0 ; i--){
      //      fprintf(stderr,"i:%u\n",i);
      //      fprintf(stderr,"i:%d cap:%lu val:(%d,%d)\n",i,R1.capacity(),i,v);
      R1[i] = tmp + fw[i][v];
      tmp = R1[i];
    }
  }
	
  void ComputeR2(int v){
    double tmp = 0;
    for (unsigned i = 0; i < maxTime ; i++){
      R2[i] = tmp*P2[i]+fw[v][i]*P6[i]+R1[i]*P7[i];
      tmp = R2[i];
    }
  }
	
  void ComputeRs(int v){
    //    fprintf(stderr,"v:%d\n",v);
    ComputeR1(v);
    ComputeR2(v);
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
  for(int i=0;i<timePoints.size();i++){
    fw[i] = new double[numWindows+1];
    bw[i] = new double[numWindows+1];
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


  //  obj.ComputeFBProbs(data,windows);
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


