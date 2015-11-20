#include "realSFS_args.h"
#include <ctime>
char * get_region(char *extra,int &start,int &stop) {
  if(!extra){
    fprintf(stderr,"Must supply parameter for -r option\n");
    return NULL;
  }
  if(strrchr(extra,':')==NULL){//only chromosomename
    char *ref = extra;
    start = stop = -1;;
    return ref;
  }
  char *tok=NULL;
  tok = strtok(extra,":");

  char *ref = tok;

  start =stop=-1;

  tok = extra+strlen(tok)+1;//tok now contains the rest of the string
 
  if(strlen(tok)==0)//not start and/or stop ex: chr21:
    return ref;
  

  if(tok[0]=='-'){//only contains stop ex: chr21:-stop
    tok =strtok(tok,"-");
    stop = atoi(tok);
  }else{
    //catch single point
    int isProper =0;
    for(size_t i=0;i<strlen(tok);i++)
      if(tok[i]=='-'){
	isProper=1;
	 break;
      }
    //fprintf(stderr,"isProper=%d\n",isProper);
    if(isProper){
      tok =strtok(tok,"-");
      start = atoi(tok)-1;//this is important for the zero offset
      tok = strtok(NULL,"-");
      if(tok!=NULL)
	stop = atoi(tok);
    }else{
      //single point
      stop = atoi(tok);
      start =stop -1;
      
    }
    
  }
  if(stop!=-1&&stop<start){
    fprintf(stderr,"endpoint:%d is larger than startpoint:%d\n",start,stop);
    exit(0);
    
  }
  if(0){
    fprintf(stderr,"[%s] ref=%s,start=%d,stop=%d\n",__FUNCTION__,ref,start,stop);
    exit(0);
  }
  return ref;
}



args * getArgs(int argc,char **argv){
  args *p = new args;
  p->bootstrap = 0;
  p->chooseChr=NULL;
  p->start=p->stop=-1;
  p->maxIter=1e2;
  p->tole=1e-6;
  p->nThreads=4;
  p->nSites =0;
  p->fname = NULL;
  p->onlyOnce = 0;
  p->emAccl =1;
  p->fstout = NULL;
  p->start=p->stop=-1;
  p->win=p->step=-1;
  p->type =0;
  p->oldout =0;
  p->seed =0;
  if(argc==0)
    return p;

  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-P"))
      p->nThreads = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-win"))
      p->win = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-type"))
      p->type = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-step"))
      p->step = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-bootstrap"))
      p->bootstrap = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-oldout"))
      p->oldout = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-m"))
      p->emAccl = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-seed"))
      p->seed = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-onlyOnce"))
      p->onlyOnce = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-r")){
      p->chooseChr = get_region(*(++argv),p->start,p->stop);
      if(!p->chooseChr)
	return NULL;
    }
    else  if(!strcasecmp(*argv,"-sfs")){
      p->sfsfname.push_back(strdup(*(++argv)));
    }
    else  if(!strcasecmp(*argv,"-fstout")){
      p->fstout = strdup(*(++argv));
    }else{
      p->saf.push_back(persaf_init<float>(*argv));
      p->fname = *argv;
      //   fprintf(stderr,"toKeep:%p\n",p->saf[p->saf.size()-1]->toKeep);
    }
    argv++;
  }
  if(p->seed==0)
    p->seed = time(NULL);
  srand48(p->seed);
  for(int i=0;(p->saf.size()>1)&&(i<p->saf.size());i++)
    p->saf[i]->kind =2;
  fprintf(stderr,"\t-> args: tole:%f nthreads:%d maxiter:%d nsites:%lu start:%s chr:%s start:%d stop:%d fname:%s fstout:%s oldout:%d seed:%ld bootstrap:%d\n",p->tole,p->nThreads,p->maxIter,p->nSites,p->sfsfname.size()!=0?p->sfsfname[0]:NULL,p->chooseChr,p->start,p->stop,p->fname,p->fstout,p->oldout,p->seed,p->bootstrap);

  if((p->win==-1 &&p->step!=-1) || (p->win!=-1&&p->step==-1)){
    fprintf(stderr,"\t-> Both -win and -step must be supplied for sliding window analysis\n");
    exit(0);

  }
  return p;
}


void destroy_safvec(std::vector<persaf *> &saf){
  //fprintf(stderr,"destroy &saf\n");
  for(int i=0;i<saf.size();i++)
    persaf_destroy(saf[i]);
}


void destroy_args(args *p){
  //  fprintf(stderr,"destroy args\n");
  destroy_safvec(p->saf);
  delete p;
}
