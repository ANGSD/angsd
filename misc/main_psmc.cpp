#include "main_psmc.h"
#include "fpsmc.h"
#include <ctime>
#include <ctype.h>
#define DEFAULT_PATTERN "4+5*3+4"
// parse a pattern like "4+5*3+4"
// the returned array holds which parameters are linked together
// number of parameters and number of free parameters will be also returned

int *psmc_parse_pattern(const char *pattern, int *n_free, int *n_pars)
{
  fprintf(stderr,"parsing pattern :\"%s\"\n",pattern);
	char *q, *p, *tmp;
	int top = 0, *stack = (int*)malloc(sizeof(int) * 0x100);
	int *pars_map, k, l, i;
	p = q = tmp = strdup(pattern);
	k = 1;
	while (1) {
		assert(isdigit(*p) || *p == '*' || *p == '+' || *p == '\0'); // allowed characters
		if (*p == '+' || *p == '\0') {
			int is_end = (*p == 0)? 1 : 0;
			*p++ = '\0';
			l = atoi(q); q = p;
			for (i = 0; i < k; ++i) {
				stack[top++] = l;
				assert(top <= 0xff);
			}
			k = 1;
			if (is_end) break;
		} else if (*p == '*') {
			*p = '\0';
			k = atoi(q); // number of repeats
			*p++ = '*'; q = p;
		} else ++p;
	}
	for (k = l = 0; k != top; ++k) l += stack[k];
	*n_pars = l - 1; *n_free = top;
	pars_map = (int*)malloc(sizeof(int) * (*n_pars + 1));
	for (k = i = 0; k != top; ++k)
		for (l = 0; l < stack[k]; ++l)
			pars_map[i++] = k;
	free(tmp); free(stack);
	return pars_map;
}


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
  p->winSize = 100;//default 100bp
  p->par =(psmc_par*) calloc(1,sizeof(psmc_par));
  if(argc==0)
    return p;

  while(*argv){
    //    fprintf(stderr,"%s\n",*argv);
    if(!strcasecmp(*argv,"-tole"))
      p->tole = atof(*(++argv));
    else  if(!strcasecmp(*argv,"-maxIter"))
      p->maxIter = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-winSize"))
      p->winSize = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-p"))
      p->par->pattern =  strdup(*(++argv));
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
  fprintf(stderr,"\t-> args: tole:%f maxiter:%d chr:%s start:%d stop:%d fname:%s seed:%ld winsize:%d\n",p->tole,p->maxIter,p->chooseChr,p->start,p->stop,p->fname,p->seed,p->winSize);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  if(p->par->pattern==NULL)
    p->par->pattern = strdup(DEFAULT_PATTERN);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  p->par->par_map = psmc_parse_pattern(p->par->pattern, &p->par->n_free, &p->par->n);
  
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
  writepsmc_header(stderr,pars->perc);

  if(1){
    psmc_wrapper(pars);
  }  else{
    //below is old printout, keeping for reference
    for(myMap::iterator it=pars->perc->mm.begin();it!=pars->perc->mm.end();++it){
      //set perchr iterator, if pars->chooseChr, then we have only use a single chr
      it = pars->chooseChr?iter_init(pars->perc,pars->chooseChr,pars->start,pars->stop):iter_init(pars->perc,it->first,pars->start,pars->stop);
      
      //print out the chromosome position and the two gls
      for(size_t s=pars->perc->first;0&&s<pars->perc->last;s++)
	fprintf(stdout,"%s\t%d\t%e\t%e\n",it->first,pars->perc->pos[s]+1,pars->perc->gls[2*s],pars->perc->gls[2*s+1]);
      
      if(pars->chooseChr!=NULL)
	break;
    }
  }
  destroy_args(pars);
  return 0;
}
