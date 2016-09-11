#include "main_psmc.h"

args * getArgs(int argc,char **argv){
  args *p = new args;
  p->chooseChr=NULL;
  p->start=p->stop=-1;
  p->maxIter=1e2;
  p->tole=1e-6;
  p->nSites =0;
  p->fname = NULL;
  p->onlyOnce = 0;
  p->seed =0;
  if(argc==0)
    return p;

  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-seed"))
      p->seed = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = get_region(*(++argv),p->start,p->stop);
      if(!p->chooseChr)
	return NULL;
    }
    else{
      p->perc = perpsmc_init(*argv);
      p->fname = *argv;
    }
    argv++;
  }
  if(p->seed==0)
    p->seed = time(NULL);
  srand48(p->seed);
  fprintf(stderr,"\t-> args: tole:%f maxiter:%d chr:%s start:%d stop:%d fname:%s seed:%ld \n",p->tole,p->maxIter,p->chooseChr,p->start,p->stop,p->fname,p->seed);

  return p;
}

//made a seperate function for this. Im assuming our args will contain allocated data at some point.
void destroy_args(args *p){
  perpsmc_destroy(p->perc);
  delete p;
}


//simple function 
int main_psmc(int argc, char **argv){
  fprintf(stderr,"\t-> we are in file: %s function: %s line:%d\n",__FILE__,__FUNCTION__,__LINE__);

  //we loop over the single chromosomes
  args *pars = getArgs(argc,argv);
  if(!pars)
    return 0;
  //this will printout the header
  writesaf_header(stderr,pars->perc);
  
  for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
    //set perchr iterator, if pars->chooseChr, then we have only use a single chr
    it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);

    //print out the chromosome position and the two gls
    for(size_t s=pars->perc->first;s<pars->perc->last;s++)
      fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->ppos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
    
    if(pars->chooseChr!=NULL)
      break;
  }
  destroy_args(pars);
  return 0;
}
