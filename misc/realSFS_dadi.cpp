#include <pthread.h>
#include <algorithm>
#include "realSFS_args.h"
#include "realSFS_dadi.h"
#include "realSFS_shared.h"
#include <htslib/faidx.h>
extern int howOften;
#include "multisafreader.hpp"

faidx_t *fai_ref=NULL;
faidx_t *fai_anc=NULL;
char *ref=NULL;
char *anc=NULL;
int ref_l=0;
int anc_l=0;

int whichmax(double *ary,int len){
  assert(len>0);
  int mmax = 0;
  for(int i=1;i<len;i++)
    if(ary[i]>ary[mmax])
      mmax=i;
  return mmax;
}


void printFastaFlank(FILE *fp,char *fasta,int center){
  if(fasta)
    fprintf(fp,"%c%c%c\t",fasta[center-1],fasta[center],fasta[center+1]);
  else
    fprintf(fp,"NNN\t");
}

template <typename T>
void print_dadi(std::vector<double *> &priors,std::vector<Matrix<T> *> &gls,int *posiToPrint,char *curChr){
  fprintf(stderr,"[%s]:\n",__func__);
  int nsites = gls[0]->x;
  if(nsites==0)
    return;

  int npop = priors.size();
  if(ref_l>0&&posiToPrint[gls[0]->x-1]>ref_l-1){
    fprintf(stderr,"\t-> Looks like position ofseqdata is larger than reference fasta\n");
    exit(0);
  }
  if(anc_l>0&&posiToPrint[gls[0]->x-1]>anc_l-1){
    fprintf(stderr,"\t-> Looks like position ofseqdata is larger than reference fasta\n");
    exit(0);
  }

  for(int s=0;s<nsites;s++){//sites
    int isvar=0;
    int counts[npop];
    double prop[npop];
    for(int p=0;p<npop;p++){//populations
      int ndim = gls[p]->y;
      double tmp[ndim];
      double *prior = priors[p];
      for(int i=0;i<ndim;i++) //categories
	tmp[i] = gls[p]->mat[s][i]*prior[i];
      normalize(tmp,ndim);
      counts[p] = whichmax(tmp,ndim);
      isvar += (counts[p]>1&&counts[p]<ndim-1)?1:0;
      prop[p] = tmp[counts[p]];
      isvar += counts[p];
    }
    if(isvar){
      printFastaFlank(stdout,ref,posiToPrint[s]);
      printFastaFlank(stdout,anc,posiToPrint[s]);
      fprintf(stdout,"%s\t%d",curChr,posiToPrint[s]+1);
      for(int i=0;i<npop;i++)
	fprintf(stdout,"\t%d\t%f",counts[i],prop[i]);
      fprintf(stdout,"\n");
    }
  }


}


template <typename T>
int main_dadi(int argc, char** argv){
  if(argc==0){
    fprintf(stderr,"\t-> Example: realSFS dadi p1.saf.idx p2.saf.idx -sfs p1.saf.idx.ml  -sfs p2.saf.idx.ml -ref hg19NoChr.fa -anc hg19ancNoChr.fa.gz \n");
    fprintf(stderr,"\t-> Possible options: []\n");
    return 0;
  }
  args *arg = getArgs(argc,argv);
  std::vector<persaf *> &saf =arg->saf;
  if(saf.size()==1)
    saf[0]->kind = 2;
  for(int i=0;i<saf.size();i++)
    assert(saf[i]->pos!=NULL&&saf[i]->saf!=NULL);
  if(saf.size()!=arg->sfsfname.size()){
    fprintf(stderr,"\t-> Must supply a perpopulation prior for each population\n");
    return 0;
  }
  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=calc_nsites(saf,arg);
  }
  if(fsizes<T>(saf,nSites)>getTotalSystemMemory())
    fprintf(stderr,"\t-> Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument\n"); 
    
  fprintf(stderr,"\t-> nSites: %lu\n",nSites);
  float bytes_req_megs =(float) fsizes<T>(saf,nSites)/1024/1024;
  float mem_avail_megs =(float) getTotalSystemMemory()/1024/1024;//in percentile
  //fprintf(stderr,"en:%zu to:%f\n",bytes_req_megs,mem_avail_megs);
  fprintf(stderr,"\t-> The choice of -nSites will require atleast: %f megabyte memory, that is at least: %.2f%% of total memory\n",bytes_req_megs,bytes_req_megs*100/mem_avail_megs);

  std::vector<Matrix<T> *> gls;
  std::vector<double *> priors;
  for(int i=0;i<saf.size();i++){
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));
    double *prior = new double[saf[i]->nChr+1];
    readSFS(arg->sfsfname[i],saf[i]->nChr+1,prior);
    priors.push_back(prior);
  }
  if(arg->ref)
    fai_ref = fai_load(arg->ref);
  if(arg->anc)
    fai_anc = fai_load(arg->anc);

  
  int *posiToPrint = new int[nSites];
  //temp used for checking pos are in sync
  setGloc(saf,nSites);
  char *lastchr=NULL;
  while(1) {
    static char *curChr=NULL;
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop,posiToPrint,&curChr,arg->fl,1);//read nsites from data
    //    int b=0;  

    char *thisChr=arg->chooseChr?arg->chooseChr:curChr;
    //    fprintf(stderr,"lastchr:%p curChr:%p\n",lastchr,thisChr);
    if(lastchr==NULL||strcmp(lastchr,thisChr)!=0){
      //cleanup first
      free(ref);ref=NULL;ref_l=0;
      free(anc);ref=NULL;anc_l=0;
      //then load;
      if(fai_ref){
	ref = faidx_fetch_seq(fai_ref, thisChr, 0, 0x7fffffff, &ref_l);
      }if(fai_anc)
	anc = faidx_fetch_seq(fai_anc, thisChr, 0, 0x7fffffff, &anc_l);
      free(lastchr);lastchr=NULL;
      lastchr=strdup(thisChr);
    }
    if(gls[0]->x>0)
      print_dadi(priors,gls,posiToPrint,thisChr);
    if(ret==-3&&gls[0]->x==0){//no more data in files or in chr, eith way we break;g
      //fprintf(stderr,"breaking\n");
      break;
    }
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
    if(arg->onlyOnce)
      break;
  }
  delete [] posiToPrint;
  free(lastchr);lastchr=NULL;
  delGloc(saf,nSites);
  destroy(gls,nSites);
  destroy_args(arg);
  for(int i=0;i<priors.size();i++)
    delete [] priors[i];
  fprintf(stderr,"\n\t-> NB NB output is no longer log probs of the frequency spectrum!\n");
  fprintf(stderr,"\t-> Output is now simply the expected values! \n");
  fprintf(stderr,"\t-> You can convert to the old format simply with log(norm(x))\n");
  return 0;
}

template int main_dadi<float>(int ,char**);

