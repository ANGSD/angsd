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


int main(int argc,char **argv){
  //start of signal handling
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  if(argc==1){
    fprintf(stderr, "\t-> ---./ngsPSMC\n");
    fprintf(stderr,"\t-> ./ngsPSMC afile.psmc.idx \n");
    return 0;
  }
  ++argv;
  --argc;
  if(!strcasecmp(*argv,"print"))
    print(--argc,++argv);
  else if(!strcasecmp(*argv,"print_header"))
    print_header(--argc,++argv);
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
