/*
  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its 10 sep 2016, it is hot outside. this code is a hack from the realSFS.
  
  
*/

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>
#include <unistd.h>
#include <zlib.h>
#include "safstat.h"
#include <algorithm>
#include <libgen.h>
#include "psmcreader.h"
#include "header.h"
#include "main_psmc.h"

int SIG_COND =1;
double ttol = 1e-16; 

int really_kill =3;
int VERBOSE = 1;
void handler(int s) {
  if(s==13)//this is sigpipe
    exit(0);
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want \'ngsPSMC\' to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;

}




int print_header(int argc,char **argv){

  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx \n");
    return 0; 
  }
  
  args *pars = getArgs(argc,argv);
  if(!pars)
    return 1;
  
  writepsmc_header(stdout,pars->perc);
  
  destroy_args(pars);
  return 0;
}

void print(int argc,char **argv){
  
  if(argc<1){
    fprintf(stderr,"\t-> Must supply afile.saf.idx files \n");
    fprintf(stderr,"\t-> Examples \n");
    fprintf(stderr,"\t-> ./ngsPSMC print pop1.saf.idx \n");
    fprintf(stderr,"\t-> ./ngsPSMC print pop1.saf.idx -r chr1:10000000-12000000\n");
    return; 
  }
  
  args *pars = getArgs(argc,argv);
  writepsmc_header(stderr,pars->perc);
  
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    
    if(pars->chooseChr!=NULL)
      it = iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop);
    else
      it = iter_init(pars->perc,it->first,pars->start,pars->stop);

    for(size_t s=pars->perc->first;s<pars->perc->last;s++)
      fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->pos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
    
    if(pars->chooseChr!=NULL)
      break;
  }
  
  destroy_args(pars);
}
double em(double &x,double *gls,int nSites,double tol){
  fprintf(stderr,"[%s] x:%f gls:%p nSites:%d\n",__FUNCTION__,x,gls,nSites);
  double llh = 0;
  double est = 0;
  double start = x;
  double lastllh =0;
  for(int iter=0;iter<20;iter++){
    fprintf(stderr,"\r %d/%d     ",iter,20);fflush(stderr);
    for(int i=0;i<nSites;i++) {
      //  fprintf(stderr,"gls=(%f,%f)\n",gls[2*i],gls[2*i+1]);
      double tmp[2];
      tmp[0] = exp(gls[2*i])*(1-start);
      tmp[1] = exp(gls[2*i+1])*(start);
      est += tmp[1]/(tmp[0]+tmp[1]);
      llh -= log(tmp[0]+tmp[1]);
    }
    est = est/((double) nSites);
    fprintf(stderr,"iter:%d est: %f llh: %f diffInLlh:%e diffInPars:%e\n",iter,est,llh,llh-lastllh,est-start);
    if(llh<lastllh){
      fprintf(stderr,"\t-> Problem with EM newllh is larger than lastllh, will break\n");
      break;
    }
    if(fabs(est-start)<tol){
      fprintf(stderr,"\t-> Difference in estimated pars is smaller than tol, convergence achieved\n");
      start=est;
      lastllh=llh;
      break;
    }
    start=est;
    lastllh=llh;
  }
  x=start;
  return lastllh;
}
void calcpost(double &x,double *gls,int nSites,double *pp){
  double start =x;
  for(int i=0;i<nSites;i++) {
    //  fprintf(stderr,"gls=(%f,%f)\n",gls[2*i],gls[2*i+1]);
    double tmp[2];
    tmp[0] = exp(gls[2*i])*(1-start);
    tmp[1] = exp(gls[2*i+1])*(start);
    double est = tmp[1]/(tmp[0]+tmp[1]);
    pp[i] = est;
    //    fprintf(stdout,"%f\n",est);
  }

}

void writefa(FILE *fp,double *pp,int ppl, int ll,double cutoff,char *chrname){
  fprintf(stdout,">%s\n",chrname);
  int i;
  for(i=0;i<ppl;i++){
    if(i>0&&((i %ll )==0)){
      // fprintf(stderr,"i:%d\n",i);
      fprintf(fp,"\n");
    }
    if(pp[i]<cutoff)
      fprintf(fp,"A");
    else
      fprintf(fp,"M");
  }
  if(i>0&&((i %ll )==0))
    return;
  else
    fprintf(fp,"\n");
  return;
}


int makeold(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"\t-> output is a vcf2fq style file \n");
    return 0; 
  }
  args *pars = getArgs(argc,argv);
  writepsmc_header(stderr,pars->perc);
  
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    if(pars->chooseChr!=NULL)
      it = iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop);
    else
      it = iter_init(pars->perc,it->first,pars->start,pars->stop);
    double opt = 0.01;
    //    double llh = em(opt,pars->perc->gls,pars->perc->last,1e-8);
    double *pp = new double[pars->perc->last];
    calcpost(opt,pars->perc->gls,pars->perc->last,pp);
    writefa(stdout,pp,pars->perc->last,50,0.9,it->first);
    break;
  }
  return 0;
}

int main(int argc,char **argv){
  //start of signal handling
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags =  0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  if(argc==1){
    fprintf(stderr, "\t-> ---./ngsPSMC\n");
    fprintf(stderr,"\t-> ./ngsPSMC [print print_header makeold] afile.psmc.idx \n");
    return 0;
  }

  ++argv;
  --argc;
  if(!strcasecmp(*argv,"print"))
    print(--argc,++argv);
  else if(!strcasecmp(*argv,"print_header"))
    print_header(--argc,++argv);
  else if(!strcasecmp(*argv,"makeold"))
    makeold(--argc,++argv);
  else {
    if(isatty(fileno(stdout))){
      fprintf(stderr,"\t-> You are printing results to the terminal consider dumping into a file\n");
      fprintf(stderr,"\t-> E.g.: \'./realSFS");
      for(int i=0;i<argc;i++)
	fprintf(stderr," %s",argv[i]);
      fprintf(stderr," >sfs.ml.txt\'\n");   
    }
    main_psmc(argc,argv);
  }
  return 0;
}
