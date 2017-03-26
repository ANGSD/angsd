#include <vector>
#include <cassert>
#include <cmath>
#include "psmcreader.h"
#include "main_psmc.h"
#include "hmm_psmc.h"
#include "../bfgs.h"

typedef struct{
  double **nP;
  double **PP;
  double *tk;
  int tk_l;
  double pix;
  int numWind;

  double rho;
  double *epsize;
}oPars;

/*
  objective function. Function to be optimized
*/

double qFunction(const double *params ,const void *d){
  oPars *data = (oPars*) d;
  void ComputeGlobalProbabilities(double *tk,int tk_l,double **P,double *epsize,double rho);
  ComputeGlobalProbabilities(data->tk,data->tk_l,data->nP,data->epsize,data->rho);
  double Q = 0;
  for (unsigned i = 0; i < data->tk_l; i++)
    Q += qkFunction(i, data->pix,data->numWind,data->nP,data->PP);
  return Q;
}


void runoptim(double *tk,int tk_l,double *epsize,double rho,double **PP,double pix,int numWin){
  double **nP = new double *[8];
  for(int i=0;i<8;i++)
    nP[i] = new double[tk_l];

  //get start
  double pars[tk_l];
  for(int i=0;i<tk_l;i++)
    pars[i] = drand48();
  //set bounds
  int nbd[tk_l];
  double lbd[tk_l];
  double ubd[tk_l];
  for(int i=0;i<tk_l;i++){
    nbd[i]=1;
    lbd[i]=0.000001;
    ubd[i]=PSMC_T_INF;
  }

  oPars data;
  data.nP = nP;
  data.PP = PP;
  data.tk = tk;
  data.tk_l = tk_l;
  data.pix = pix;
  data.numWind=numWin;
  data.rho= rho;
  data.epsize=epsize;
  
  double max_llh = findmax_bfgs(tk_l,pars,(void *)&data,qFunction,NULL,lbd,ubd,nbd,1);
}


int psmc_wrapper(args *pars,int block) {
#if 1 //print pars
  psmc_par *p=pars->par;
  fprintf(stderr,"par->n:%d\tpar->n_free:%d\tpar_map:%p\tpar->pattern:%s\tpar->times:%p\tpar->params:%p\n",p->n,p->n_free,p->par_map,p->pattern,p->times,p->params);
  for(int i=0;i<pars->par->n+1;i++)
    fprintf(stderr,"%i)\t%f\t%f\n",i,pars->par->times[i],pars->par->params[i]);
  
#endif
  int tk_l = pars->par->n+2;
  double *tk = new double [tk_l];
  double *epsize = new double [tk_l];
  setEPSize(epsize,tk_l,p->params);
  //(nelems,array,max_t,alpha,array with values from file, can be NULL)
  setTk(tk_l,tk,15,0.01,p->times);//<- last position will be infinity
#if 1
  for(int i=0;i<tk_l;i++)
    fprintf(stderr,"[%d]:(%f,%f)\n",i,tk[i],epsize[i]);
#endif
  
  //initialize all hmm (one for each chr), for now just a single
  fastPSMC obj;
  myMap::const_iterator it = pars->perc->mm.begin();
  //calculate window end points
  obj.setWindows(pars->perc,it->first,pars->start,pars->stop,pars->block);
  //  obj.printWindows(stdout);
  //allocate internal structures needed
  obj.allocate(tk_l);
  //make an hmm

  obj.make_hmm(tk,tk_l,pars->perc->gls,epsize);

  return 1;
}

