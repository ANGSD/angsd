#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"

using namespace std;

struct Site{
  double g0, g1;
  size_t siteId;//might not be needed
  char *name;//chr_firstIdx:lastIdx_winStart:winStop //not needed
};

class fastPSMC {
  unsigned maxTime; //number of time intervals
  double rho;
  std::vector<double> timePoints;//coalescence times.
  std::vector<double> P1, P2, P3, P4, P5, P6, P7;//each has length of timePoints.size()
  std::vector<double> R1, R2; ////each has length of timePoints.size()
  std::vector<double> epSize; //effective population size
  std::vector<double> fbProb;
	
public:
  void init();
  void llh(std::vector<Site> &data);
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
  
private:
  void ComputeP1(){
    for (unsigned i = 0; i < maxTime; i++)
      P1[i] = exp( -rho*1*(timePoints[i+1] - timePoints[i]) );
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

  void ComputeR1(){
    double tmp = 0;
    for (unsigned i = maxTime - 1; i >= 0 ; i++){
      R1[i] = tmp + fbProb[i];
      tmp = R1[i];
    }
  }
	
  void ComputeR2(){
    double tmp = 0;
    for (unsigned i = 0; i < maxTime ; i++){
      R2[i] = tmp*P2[i]+fbProb[i]*P6[i]+R1[i]*P7[i];
      tmp = R2[i];
    }
  }
	
  void ComputeRs(){
    ComputeR1();
    ComputeR2();
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
	
  void ComputeFBProbs(){
    fbProb[0] = fbProb[0]*P1[0] + R1[0]*P3[0] + fbProb[0]*P4[0];
    for (unsigned i = 1; i < maxTime; i++){
      fbProb[i] = (fbProb[i]*P1[i] + R2[i-1]*P2[i-1] + R1[i]*P3[i] + fbProb[i]*P4[i]);//removed last factor otherwise woulndt compile
      //fbProb[i] = (fbProb[i]*P1[i] + R2[i-1]*P2[i-1] + R1[i]*P3[i] + fbProb[i]*P4[i])*emProb[i];
    }
  }
};

void fastPSMC::init(){
  ComputeGlobalProbabilities();
  ComputeRs();
}

double fastPSMC::llh(std::vector<Site> &data){
  double llh=0;
  for(int v=0;v<data.length();v++){//loop over windows
    double g0=exp(data[v].g0);//homo gl
    double g1=exp(data[v].g1);//het gl
    double tmp =0;
    for (unsigned t = 1; t < maxTime; t++)
      tmp += log((1-exp(t))*g0+exp(t)*g1);
	
    llh+=tmp;
  }

  return llh;
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
int print_psmc_print_windows(std::vector<Site> &data){
  for(int i=0;i<data.size();i++)
    fprintf(stdout,"%lu\t%s\t%f\t%f\n",data[i].siteId,data[i].name,data[i].g0,data[i].g1);

}

//function to print the data we need
int main_analysis(std::vector<Site> &data){
  for(int i=0;0&&i<data.size();i++)
    fprintf(stdout,"%lu\t%s\t%f\t%f\n",data[i].siteId,data[i].name,data[i].g0,data[i].g1);
  fastPSMC obj;
  obj.init();
}


int psmc_wrapper(args *pars){
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);
  //loop over chrs;
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    //set perchr iterator, if pars->chooseChr, then we have only use a single chr
    it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);

    std::vector<Site> data;//this should be moved up so the end data will contain data for all chrs and not just a single chrs

    //generate window data
    int beginIndex = pars->perc->first;
    while(1){
      int beginPos=pars->perc->pos[beginIndex]+1;//add one due to positions being offset with one.
      int endPos = beginPos+pars->winSize;
      Site d;
      d.g0=d.g1=log(.0);//-Inf
      d.name=NULL;
      int at=beginIndex;
      while(at<pars->perc->last&&pars->perc->pos[at]<endPos){
	d.g0 = addProtect2(d.g0,pars->perc->gls[at]);
	d.g1 = addProtect2(d.g1,pars->perc->gls[at+1]);
	at++;
      }
      if(at>=pars->perc->last)
	break;
      d.name =(char*) calloc(1024,sizeof(char));
      snprintf(d.name,1024,"%s_%d:%d_%d:%d",it->first,beginIndex,at,beginPos,endPos);
      d.siteId=data.size();//very redundant. but lets keep for now.
      data.push_back(d);
      beginIndex = at;
    }
    main_analysis(data);
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


