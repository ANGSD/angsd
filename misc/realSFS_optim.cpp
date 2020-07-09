#include <pthread.h>
#include <algorithm>
#include <numeric>
#include "realSFS_args.h"
#include "realSFS_optim.h"
extern int howOften;
#include "multisafreader.hpp"
void readSFS(const char *fname, size_t hint, double *ret);

size_t * foldremapper = NULL;
int * foldfactors = NULL;
int * foldkeep = NULL;

template<typename T>
struct emPars
{
  int threadId;
  std::vector <Matrix<T> *> gls;
  size_t from;
  size_t to;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  double *inner;//allocated for every thread
  int dim;
};

pthread_t *thd = NULL;
size_t *bootstrap = NULL;
size_t bootnSites = 0;//nspope; number of sites changes if bootstrapping chromosomes
emPars<float> *emp = NULL;

double ttol = 1e-16; 

bool multi_saf_loop (int *bin, int *nchr, int i);
size_t sfs_linear_index (int *bin, int *nchr, int npop, bool rev);

// use a recursive nested loop so we don't need to rewrite this for each # of populations
template <typename T>
double lik(double *sfs, std::vector< Matrix<T> *> &gls, size_t from, size_t to)
{
  assert(from >= 0 && to >=0);
  double r = 0;
  size_t ss, lind;

  int chr[gls.size()];
  int rbin[gls.size()];
  int bin[gls.size()];
  int stt[gls.size()];
  int len[gls.size()];

  for (int p = 0; p<gls.size(); ++p)
    chr[p] = gls[p]->y-1;

  for (size_t s = from; SIG_COND && s < to; ++s)
  {
    ss = !bootstrap ? s : bootstrap[s];

    for (int p = 0; p<gls.size(); ++p)
    {
      bin[p] = 0;
      stt[p] = gls[p]->mat[ss][0];
      len[p] = gls[p]->mat[ss][1]-1;
    }

    double tmp = 0;
    do 
    {
      for (int p = 0; p<gls.size(); ++p)
        rbin[p] = stt[p] + bin[p];
      lind = sfs_linear_index(rbin, chr, gls.size(), false);

      double tmp2 = sfs[foldremapper[lind]] * foldfactors[lind];
      for (int p=0; p<gls.size(); ++p)
        tmp2 *= gls[p]->mat[ss][2 + bin[p]];
      tmp += tmp2;
    }
    while(multi_saf_loop(bin, len, gls.size()-1));

    r += log(tmp);
  }

  return r;
}

// use a recursive nested loop so we don't need to rewrite this for each # of populations
template <typename T>
void emStep(double *pre, std::vector<Matrix<T> *> &gls, double *post, size_t start, size_t stop, int dim, double *inner)
{
  for(int x=0; x<dim; x++)
    post[x] = 0.0;
  size_t ss, lind;
   
  int chr[gls.size()];
  int rbin[gls.size()];
  int bin[gls.size()];
  int stt[gls.size()];
  int len[gls.size()];

  for (int p = 0; p<gls.size(); ++p)
    chr[p] = gls[p]->y-1;

  for(size_t s=start; SIG_COND && s<stop; s++)
  {
    ss = !bootstrap ? s : bootstrap[s];

    for (int i = 0; i<gls.size(); ++i)
    {
      bin[i] = 0;
      stt[i] = gls[i]->mat[ss][0];
      len[i] = gls[i]->mat[ss][1] - 1;
    }

    for(int i=0; i<dim; i++)
      inner[i] = 0;

    do {
      for (int p = 0; p<gls.size(); ++p)
        rbin[p] = stt[p] + bin[p];
      lind = sfs_linear_index(rbin, chr, gls.size(), false);

      double tmp2 = pre[foldremapper[lind]] * foldfactors[lind];
      for (int i=0; i<gls.size(); ++i)
        tmp2 *= gls[i]->mat[ss][2 + bin[i]];
      inner[foldremapper[lind]] += tmp2;
    }
    while(multi_saf_loop(bin, len, gls.size()-1));

    normalize(inner, dim);
    for(int x=0; x<dim; x++)
      post[x] += inner[x];
  }
}

template <typename T>
void *like_slave(void *p)
{
  emPars<T> &pars = emp[(size_t) p];
  pars.lik = lik(pars.sfs, pars.gls, pars.from, pars.to);
  pthread_exit(NULL);
}

template <typename T>
double like_master(int nThreads)
{
  normalize(emp[0].sfs,emp[0].dim);
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,like_slave<T>,(void*) i);
    if(rc)
      fprintf(stderr,"Error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  double res=0;
  for(int i=0;i<nThreads;i++){
    //    fprintf(stderr,"lik=%f\n",emp[i].lik);
    res += emp[i].lik;
  
  }
  return res;
}

template <typename T>
void *emStep_slave(void *p)
{
  emPars<T> &pars = emp[(size_t) p];
  emStep<T>(pars.sfs, pars.gls, pars.post, pars.from, pars.to, pars.dim, pars.inner);
  pthread_exit(NULL);
}

template<typename T>
void emStep_master(double *post,int nThreads){
  for(size_t i=0;i<nThreads;i++){
    int rc = pthread_create(&thd[i],NULL,emStep_slave<T>,(void*) i);
    if(rc)
      fprintf(stderr,"Error creating thread\n");
    
  }
  for(int i=0;i<nThreads;i++)
    pthread_join(thd[i], NULL);
    
  memcpy(post,emp[0].post,emp[0].dim*sizeof(double));
  for(int i=1;i<nThreads;i++){
    for(int j=0;j<emp[0].dim;j++)
      post[j] += emp[i].post[j];
  }
  
  for(int j=0;j<emp[0].dim;j++)
    if(bootnSites)
      post[j] /= double(bootnSites);
    else
      post[j] /= double(emp[0].gls[0]->x);//nspope; rescale *after* threads have merged
}

template<typename T>
emPars<T> *setThreadPars(std::vector<Matrix<T> * > &gls,double *sfs,int nThreads,int dim,size_t nSites)
{
  //  fprintf(stderr,"nSites:%d\n",nSites);
  emPars<T> *temp = new emPars<T>[nThreads];
  size_t blockSize = nSites/nThreads;
  for(int i=0;i<nThreads;i++){
    temp[i].threadId = i;
    temp[i].gls=gls;
    temp[i].from =0;
    temp[i].to=blockSize;
    temp[i].sfs = sfs;
    temp[i].post=new double[dim];
    temp[i].inner=new double[dim];
    temp[i].dim = dim;
  }
  //redo the from,to
  for(int i=1;i<nThreads;i++){
    temp[i].from = temp[i-1].to;
    temp[i].to = temp[i].from+blockSize;
  }
  //fix last end point
  temp[nThreads-1].to=nSites;
  
  thd= new pthread_t[nThreads];
  return temp;
}

template<typename T>
void destroy(emPars<T> *a,int nThreads)
{
  for(int i=0;i<nThreads;i++){
    delete [] a[i].inner;
    delete [] a[i].post;
  }
  delete [] a;
  delete [] thd;
}

template <typename T>
double em(double *sfs,double tole,int maxIter,int nThreads,int dim,std::vector<Matrix<T> *> &gls,int verbose){
  if(bootnSites)
    emp = setThreadPars<T>(gls,sfs,nThreads,dim,bootnSites);
  else
    emp = setThreadPars<T>(gls,sfs,nThreads,dim,gls[0]->x);
  fprintf(stderr,"------------\n");
  double oldLik,lik;
  oldLik = like_master<T>(nThreads);
  
  fprintf(stderr,"startlik=%f\n",oldLik);fflush(stderr);
  if(verbose){
    fprintf(stdout,"%d\t%f",-1,oldLik);
    for(int i=0;i<dim;i++)
      fprintf(stdout,"\t%f",sfs[i]);
    fprintf(stdout,"\n");
  }

  double tmp[dim];
  double expected[dim];
  for(int i=0;i<dim;i++)
    expected[i] = sfs[i]; 

  for(int it=0;SIG_COND&&it<maxIter;it++) {
    emStep_master<T>(tmp,nThreads);
    double sr2 = 0;
    double exp_diff = 0;
    int exp_which_diff = -1;
    for(int i=0;i<dim;i++){
      sr2 += (sfs[i]-tmp[i])*(sfs[i]-tmp[i]);
      sfs[i]= tmp[i];
      double tmpdiff = fabs(tmp[i]*gls[0]->x-expected[i]);
      if(tmpdiff>exp_diff){
        exp_diff = tmpdiff;
        exp_which_diff = i;
      }
      expected[i] = gls[0]->x*tmp[i];
    }

    lik = like_master<T>(nThreads);

    fprintf(stderr,"[%d] lik=%f diff=%e sr:%e nsites_difference[%d]: %e\n",it,lik,fabs(lik-oldLik),sr2,exp_which_diff,exp_diff);
    if(verbose){
      fprintf(stdout,"%d\t%f",it,lik);
      for(int i=0;i<dim;i++)
        fprintf(stdout,"\t%f",sfs[i]);
      fprintf(stdout,"\n");
    }
    if((fabs(lik-oldLik)<tole)||(0&&(sqrt(sr2)<tole))){//should update simfiles...
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  destroy<T>(emp,nThreads);
  return oldLik;
}


template <typename T>
double emAccl(double *p,double tole,int maxIter,int nThreads,int dim,std::vector<Matrix<T> *> &gls,int useSq,int verbose){
  if(bootnSites)
    emp = setThreadPars<T>(gls,p,nThreads,dim,bootnSites);
  else
    emp = setThreadPars<T>(gls,p,nThreads,dim,gls[0]->x);
  //double *p contains the 'start' that is used for em step and needs to be updated
  fprintf(stderr,"------------\n");
  double oldLik,lik;
  oldLik = like_master<T>(nThreads);

  fprintf(stderr,"startlik=%f\n",oldLik);fflush(stderr);
  if(verbose){
    fprintf(stdout,"%d\t%f",-1,oldLik);
    for(int i=0;i<dim;i++)
      fprintf(stdout,"\t%f",p[i]);
    fprintf(stdout,"\n");
  }

  double *p1 = new double[dim];// first guess
  double *p2 = new double[dim];// second guess
  double *q1 = new double[dim];
  double *q2 = new double[dim];
  double *pnew = new double[dim];//final guess
  double *oldp = new double[dim];
  double *tmp = new double[dim];
  int stepMax = 1;
  int mstep = 4;
  int stepMin = 1;

  int iter =0;
  double expected[dim];
  for(int i=0;i<dim;i++)
    expected[i] = p[i]; 

  while(SIG_COND&&iter<maxIter){
    emStep_master<T>(p1,nThreads);
    iter++;
    double sr2 =0;
    for(size_t i=0;i<dim;i++){
      q1[i] = p1[i]-p[i];
      sr2 += q1[i]*q1[i];
    }
    // oldp,p1,p2,pnew
    memcpy(oldp,p,sizeof(double)*dim);
    memcpy(p,p1,sizeof(double)*dim);  
    if(sqrt(sr2)<tole){
      fprintf(stderr,"\t-> Breaking EM(sr2) at iter:%d, sqrt(sr2):%e\n",iter,sqrt(sr2));
      break;
    }
    emStep_master<T>(p2,nThreads);
    iter++;


    double sq2 =0;
    for(size_t i=0;i<dim;i++){
      q2[i] = p2[i]-p1[i];
      sq2 += q2[i]*q2[i];
    }

    if(sqrt(sq2)<tole){
      fprintf(stderr,"\t-> Breaking EM(sq2) at iter:%d, sqrt(sq2):%e\n",iter,sqrt(sq2));
      break;
    }
    double sv2=0;
    for(size_t i=0;i<dim;i++)
      sv2 += (q2[i]-q1[i])*(q2[i]-q1[i]);

    double alpha = sqrt(sr2/sv2);
    alpha =std::max(stepMin*1.0,std::min(1.0*stepMax,alpha));

    //the magical step
    double ts =0;
    for(size_t j=0;j<dim;j++){
      pnew[j] = oldp[j] + 2.0 * alpha * q1[j] + alpha*alpha * (q2[j] - q1[j]);
      ts += pnew[j];
    }
    if(fabs(1-ts)>1e-8)
      fprintf(stderr,"\t-> Sum of parameters is not one? This could be due to loss of precision value:%e 1-value:%e \n",ts,1-ts);

#if 1 //fix for going out of bound
    int printOutOfBounds = 0;
    for(size_t j=0;j<dim;j++){
      if(pnew[j]<ttol){
        printOutOfBounds =1;
        pnew[j]=ttol;
      }if(pnew[j]>1-ttol){
        printOutOfBounds =1;
        pnew[j]=1-ttol;
      }
    }
    if(printOutOfBounds)
      fprintf(stderr,"\t-> Instability detected, accelerated guess is too close to bound or outside will fallback to regular EM for this step\n");
#endif

    int extra =0;
    if(printOutOfBounds==0&&fabs(alpha-1) >0.01){
      //this is clearly to stabilize
      //      double tmp[dim];
      memcpy(p,pnew,sizeof(double)*dim);
      emStep_master<T>(tmp,nThreads);
      memcpy(pnew,tmp,sizeof(double)*dim);
      extra =1;
      iter++;
    }
    memcpy(p,pnew,sizeof(double)*dim);

    //like_master is using sfs[] to calculate like
    lik = like_master<T>(nThreads);
    if(lik<oldLik||printOutOfBounds==1){
      if(printOutOfBounds==0)
        fprintf(stderr,"\t-> New like is worse? new:%e old:%e diff:%e will skip accelerated EM in this round\n",lik,oldLik,lik-oldLik);
      memcpy(p,p2,sizeof(double)*dim);
      lik = like_master<T>(nThreads);
      stepMax =1;
      mstep=4;
      stepMin =1;
      if(extra)//if have been doing the stabilizing step
        iter--;
    }
    double exp_diff = 0;
    int exp_which_diff = -1;
    for(int i=0;i<dim;i++){
      double tmpdiff = fabs(p[i]*gls[0]->x-expected[i]);
      if(tmpdiff>exp_diff){
        exp_diff = tmpdiff;
        exp_which_diff = i;
      }
      expected[i] = gls[0]->x*p[i];
    }    
    fprintf(stderr,"lik[%d]=%f diff=%e alpha:%f sr2:%e nsites_difference[%d]: %e\n",iter,lik,fabs(lik-oldLik),alpha,sr2,exp_which_diff,exp_diff);
    if(verbose){
      fprintf(stdout,"%d\t%f",iter,lik);
      for(int i=0;i<dim;i++)
        fprintf(stdout,"\t%f",p[i]);
      fprintf(stdout,"\n");
    }
    if(std::isnan(lik)) {
      //this should not happen
      fprintf(stderr,"\t-> Observed NaN in accelerated EM, will use last reliable value. Consider using as input for ordinary em\n");
      fprintf(stderr,"\t-> E.g ./realSFS -sfs current.output -m 0 >new.output\n");//thanks morten rasmussen
      memcpy(p,p2,sizeof(double)*dim);
      break;
    }

#if 1
    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
#endif
    if (alpha == stepMax) 
      stepMax = mstep * stepMax;
    if(stepMin<0 &&alpha==stepMin)
      stepMin = mstep*stepMin;

  }
  destroy<T>(emp,nThreads);
  delete [] p1;
  delete [] p2;
  delete [] q1;
  delete [] q2;
  delete [] pnew;
  delete [] oldp;
  delete [] tmp;
  return oldLik;
}

void make_folder (size_t*& mapping, int*& weight, int*& keep, int fold, int* chr, int npop);

template <typename T>
int main_opt(args *arg){

  int chr[arg->saf.size()];
  for (int i=0; i<arg->saf.size(); ++i)
    chr[i] = arg->saf[i]->nChr;

  make_folder(foldremapper, foldfactors, foldkeep, arg->fold==1, chr, arg->saf.size());
   
  std::vector<persaf *> &saf =arg->saf;
  for(int i=0;i<saf.size();i++)
    assert(saf[i]->pos!=NULL&&saf[i]->saf!=NULL);
  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=calc_nsites(saf,arg);
  }

  // for v4 memory check doesn't apply... TODO: could guesstimate based on file size??
  if(saf[0]->version==2 && fsizes<T>(saf,nSites)>getTotalSystemMemory())
    fprintf(stderr,"\t-> Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
  fprintf(stderr,"\t-> nSites: %lu\n",nSites);
  if(saf[0]->version==2)
  {
    float bytes_req_megs =(float) fsizes<T>(saf,nSites)/1024/1024;
    float mem_avail_megs =(float) getTotalSystemMemory()/1024/1024;//in percentile
    //fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
    fprintf(stderr,"\t-> The choice of -nSites will require atleast: %f megabyte memory, that is at least: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);
  }

  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim=(int) parspace(saf);
  double *sfs=new double[ndim];

  //nspope; saf input is sorted by chromosome, and read in that order
  //so only need to track starting pos of each chrom in gls[0]
  //if ever saf input is allowed to be unsorted, -resample_chr will break
  std::vector<size_t> chrStart(1,0);

  //temp used for checking pos are in sync
  setGloc(saf,nSites);
  while(1) {
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop,NULL,NULL,arg->fl,1);//read nsites from data
    int b=0;  

    //fprintf(stderr,"\t\tRET:%d gls->x:%lu\n",ret,gls[0]->x);
    if(ret==-2&&gls[0]->x==0)//no more data in files or in chr, eith way we break;
      break;

    if(ret==-2&&arg->bootstrap&&arg->resample_chr)//nspope; start pos for new chrom
      chrStart.push_back(gls[0]->x);

    {
      if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
        //fprintf(stderr,"continue continue\n");
        continue;
      }
    }

    if(gls[0]->x==0)
      continue;
    
    fprintf(stderr,"\t-> Will run optimization on nSites: %lu\n",gls[0]->x);

    //nspope; calculate number of sites for each chr
    std::vector<size_t> chrSize (chrStart.size()-1);
    if(arg->bootstrap&&arg->resample_chr&&chrStart.size()>1)
      std::adjacent_difference(++chrStart.begin(), chrStart.end(), chrSize.begin());
    chrSize.push_back(gls[0]->x-chrStart.back());
    std::vector<size_t> bootChr (chrStart.size(), 0);
  
  neverusegoto:
    if(arg->bootstrap)
      fprintf(stderr,"Will do bootstrap replicate %d/%d\n",b+1,arg->bootstrap);
    if(arg->sfsfname.size()!=0){
      fprintf(stderr,"\t-> Reading SFS from file: \'%s\'\n",arg->sfsfname[0]);
      readSFS(arg->sfsfname[0],ndim,sfs);
    } else{
      if(arg->seed==-1){
	for(int i=0;i<ndim;i++)
	  //	  sfs[i] = ((double)(ndim))/(i+1);
	  sfs[i] = 1;
      }else{
	for(int i=0;i<ndim;i++){
	  double r=drand48();
	  while(r==0.0)
	    r = drand48();
	  sfs[i] = r;
	}
      }
    }
    if(foldremapper) { //implies 1 or 2 pops
#if 0
      fprintf(stdout,"SFSIN");
      for(int i=0;i<ndim;i++)
	fprintf(stdout,"%f\t",sfs[i]);
      fprintf(stdout,"\n");
#endif
      
      //this block is debug block for validating that the input sfs gets folded correctly
      double newtmp[ndim];
      for(int i=0;i<ndim;i++)
	newtmp[i] = 0;
      for(int i=0;i<ndim;i++){
	newtmp[foldremapper[i]] += sfs[i];
	//	fprintf(stderr,"%d %f\n",i,newtmp[i]);
      }
      for(int i=0;i<ndim;i++){
	//	fprintf(stdout,"%f ",newtmp[i]);
	sfs[i] = newtmp[i];
      }
    }
    normalize(sfs,ndim);
      
      if(bootstrap==NULL &&arg->bootstrap)
	bootstrap = new size_t[gls[0]->x];
      
      if(bootstrap){
        if(arg->resample_chr){//nspope; resample chromosomes, delete/reallocate bootstrap[], fill in sites
          size_t newNsites = 0;
          for(size_t i=0; i<bootChr.size(); i++){
            bootChr.at(i) = lrand48()%bootChr.size();
            newNsites += chrSize.at(bootChr.at(i));
          }
  	  fprintf(stderr,"\t-> Resampled %lu chromosomes to get %lu total sites\n",bootChr.size(),newNsites);
          delete [] bootstrap; bootstrap = NULL;//no leaks here
          bootstrap = new size_t[newNsites];
          bootnSites = newNsites;//doesn't equal gls[0]->x thus is defined globally and used by setThreadPars 
          size_t k=0;
          for(size_t i=0;i<bootChr.size();i++)
            for(size_t j=0;j<chrSize.at(bootChr.at(i));j++)
              bootstrap[k++] = chrStart.at(bootChr.at(i))+j;
          std::sort(bootstrap,bootstrap+newNsites);
        }else{
	  for(size_t i=0;i<gls[0]->x;i++)
	    bootstrap[i] = lrand48() % gls[0]->x;
	  std::sort(bootstrap,bootstrap+gls[0]->x);
        }
      }
      double lik;

      if(arg->emAccl==0)
        lik = em<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim,gls,arg->verbose);
      else
        lik = emAccl<float>(sfs,arg->tole,arg->maxIter,arg->nThreads,ndim,gls,arg->emAccl,arg->verbose);

      fprintf(stderr,"likelihood: %f\n",lik);
      fprintf(stderr,"------------\n");

      // TODO: nsp check this. When bootstrapping contigs the number of sites can change.
      for(int x=0;x<ndim;x++)
        if (bootnSites)
          fprintf(stdout, "%f ", foldkeep[x]*bootnSites*sfs[x]);
        else
          fprintf(stdout, "%f ", foldkeep[x]*gls[0]->x*sfs[x]);
      fprintf(stdout,"\n");
      fflush(stdout);

      if(++b<arg->bootstrap)
        goto neverusegoto;

      for(int i=0;i<gls.size();i++)
        gls[i]->x =0;
    
      if(ret==-2&&arg->chooseChr!=NULL)
        break;
      if(arg->onlyOnce)
        break;
  }
  delGloc(saf,nSites);
  destroy(gls,nSites);
  destroy_args(arg);
  delete [] sfs;
  
  fprintf(stderr,"\n\t-> NB output is no longer log probs of the frequency spectrum!\n");
  fprintf(stderr,"\t-> Output is now simply the expected values! \n");
  fprintf(stderr,"\t-> You can convert to the old format simply with log(norm(x))\n");
  fprintf(stderr,"\n\t-> Please check that it has converged!\n\t-> You can add start new optimization by supplying -sfs FILE, where is >FILE from this run\n\t-> -maxIter INT -tole FLOAT\n");

  if(foldremapper)
    delete [] foldremapper;
  if(foldfactors)
    delete [] foldfactors;
  if(foldkeep)
    delete [] foldkeep;
  return 0;
}

template int main_opt<float>(args *arg);
