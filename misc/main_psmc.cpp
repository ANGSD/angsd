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



void setpars( char *fname,psmc_par *pp) {
  fprintf(stderr,"[%s]:%s\n",__FUNCTION__,fname);
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(!fp){
    fprintf(stderr,"\t-> Problem opening file:%s\n",fname);
    exit(0);
  }
  char *buf = new char[fsize(fname)+10];
  memset(buf,0,fsize(fname)+10);
  assert(fsize(fname)==fread(buf,sizeof(char),fsize(fname),fp));
  fclose(fp);
  char *slashslash[100];
  
  //stupid loop below....
  int n=0;

  for(int i=0;i<strlen(buf)-1;i++){//offset with one so we dont get the last empty output from PSCMC
    if(strncmp(buf+i,"\n//\n",4)==0)
      slashslash[n++] = buf+i;
  }

  char *last= slashslash[n-2];
  char *line = NULL;
  strtok(last,"\n");

  line=strtok(NULL,"\n");
  int IT=-1;
  sscanf(line,"IT\t%d",&IT);
  int RD=-1;
  line=strtok(NULL,"\n"); sscanf(line,"RD\t%d",&RD);
  double LK=-1;
  line=strtok(NULL,"\n"); sscanf(line,"LK\t%lf",&LK);
  double QD[2]={-1,-1};
  line=strtok(NULL,"\n"); sscanf(line,"QD\t%lf -> %lf",&QD[0],&QD[1]);
  double RI=-1;
  line=strtok(NULL,"\n"); sscanf(line,"RI\t%lf",&RI);
  double TR[2]={-1,-1};
  line=strtok(NULL,"\n"); sscanf(line,"TR\t%lf\t%lf",&TR[0],&TR[1]);
  double MT={-1};
  line=strtok(NULL,"\n"); sscanf(line,"MT\t%lf",&MT);
  double C_pi=-1;
  double n_recomb=-1;
  line=strtok(NULL,"\n"); sscanf(line,"MM\tC_pi: %lf, n_recomb: %lf",&C_pi,&n_recomb);
  std::vector<char *> RS;
  fprintf(stderr,"IT:%d RD:%d lk:%f qd[0]:%f qd[1]:%f ri:%f tr[0]:%f tr[1]:%f mt:%f c_pi:%f n_rebomc:%f\n",IT,RD,LK,QD[0],QD[1],RI,TR[0],TR[1],MT,C_pi,n_recomb);
  while(((line=strtok(NULL,"\n")))){
    if(line[0]=='R'&&line[1]=='S')
      RS.push_back(line);
    else{
      break;
    }
  }
  //  fprintf(stderr,"number of lines with RS:%lu\n",RS.size());exit(0);
  char *nline = strdup(line);
  char *tok = strtok(nline,"\n\t ");
  tok = strtok(NULL,"\n\t ");
  pp->pattern=strdup(tok);

  pp->par_map= psmc_parse_pattern(pp->pattern,&pp->n_free,&pp->n);
  assert(RS.size()-1==pp->n);
  pp->params = new double[RS.size()];
  pp->times = new double[RS.size()];

  for(int i=0;i<RS.size();i++){
    int val;
    sscanf(RS[i],"RS\t%d\t%lf\t%lf\t",&val,&pp->times[i],&pp->params[i]);
    assert(val==i);
  }
  fprintf(stderr,"\t-> Done reading parameters from file: \'%s\'\n",fname);
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
  p->block = 100;//default 100bp
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
      p->block = atoi(*(++argv));
    else  if(!strcasecmp(*argv,"-p"))
      p->par->pattern =  strdup(*(++argv));
    else  if(!strcasecmp(*argv,"-nSites"))
      p->nSites = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-seed"))
      p->seed = atol(*(++argv));
    else  if(!strcasecmp(*argv,"-infile"))
      setpars(*++argv,p->par);
    
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
  fprintf(stderr,"\t-> args: tole:%f maxiter:%d chr:%s start:%d stop:%d fname:%s seed:%ld winsize:%d\n",p->tole,p->maxIter,p->chooseChr,p->start,p->stop,p->fname,p->seed,p->block);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  if(p->par->pattern==NULL)
    p->par->pattern = strdup(DEFAULT_PATTERN);
  //  fprintf(stderr,"par:%p par->pattern:%p DEFAULT_PATTERN:%s\n",p->par,p->par->pattern,DEFAULT_PATTERN);
  if(p->par->pattern!=NULL&&p->par->params==NULL)
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
    psmc_wrapper(pars,100);
  }else{
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
