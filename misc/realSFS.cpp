/*

  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

  april 13, safv3 added, safv2 removed for know. Will be reintroduced later.
  april 20, removed 2dsfs as special scenario
  april 20, split out the safreader into seperate cpp/h
  may 5, seems to work well now
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
#include <libgen.h>
#include <algorithm>
#include "realSFS_args.h"
#include "header.h"
#include "safcat.h"
#include "realSFS_optim.h"
#include "realSFS_dadi.h"
int SIG_COND =1;
int howOften =5e6;//how often should we print out (just to make sure something is happening)
#include "multisafreader.hpp"

int really_kill =3;
int VERBOSE = 1;
/*
  when doing folded spectra we reuse the unfolded datastructure containing sample allele frequencies (SAF).
  and we loop over the entire unfolded matrix. But (here comes the trick).
  categories of the spectra with an allelecount >0.5 should be 'merged' with the corresponding category. that is
  2dsfs[al1,al2] =  2dsfs[al1,al2]+2dsfs[dimpop1-al1,dimpop2-al2]
  so we populate a remap/lookup vector that maps the categories of the unfolded 2dsfs with the parameter(sfs) for the unfolded. Such that some of the entries in the unfoldedsfs gets reused.

example 6alleles in pop1, 5 alleles inpop2: 5+6 total
so entry 4,4 =freq 8/11 >0.5 gets should be added to (6-4,5-4)=(2,1)
  
> m
  0  1  2  3  4  5
0 1  8 15 22 29 36
1 2  9 16 23 30 37
2 3 10 17 24 31 38
3 4 11 18 25 32 39
4 5 12 19 26 33 40
5 6 13 20 27 34 41
6 7 14 21 28 35 42

> fold2d(m)
   0  1  2  3  4  5
0 43 43 43 43 43 43
1 43 43 43 43 43 NA
2 43 43 43 43 NA NA
3 43 43 43 NA NA NA
4 43 43 NA NA NA NA
5 43 NA NA NA NA NA
6 NA NA NA NA NA NA

the remap vector will therefore substitute position 33 with 10
   0  1  2  3  4  5
0  1  8 15 22 29 36
1  2  9 16 23 30  6
2  3 10 17 24 12  5
3  4 11 18 18 11  4
4  5 12 24 17 10  3
5  6 30 23 16  9  2
6 36 29 22 15  8  1

abit mindfuck, the above matrices are represented as flattened versions.
also some subtle things with row vs col in the implementation(since c is rowwise, R colwise)

*/

//

int *makefoldadjust(int *ary,int len){
  //  fprintf(stderr,"makefoldadjust:%d\n",len);
  int *ret = new int[len];
  for(int i=0;i<len;i++)
    ret[i] = 1;
  int shouldexit =1;
  for(int i=0;i<len;i++)
    if(ary[i]!=i)
      shouldexit =0;
  if(shouldexit){
    fprintf(stderr,"will breakk\n");
    return ret;
  }
  for(int i=0;i<len;i++){
    int nobs=0;
    for(int j=0;j<len;j++){
      //      fprintf(stderr,"i:%d j:%d\n",i,j);
      if(ary[i]==ary[j])
	nobs++;
    }
    if(nobs==1)
      ret[i] = 2;
    else if(nobs>2){
      fprintf(stderr,"nobs>2:%f should only be 1 or 2\n");
      exit(0);
    }
  }
  for(int i=0;1&&i<len;i++)
    fprintf(stderr,"%d %d\n",i, ret[i]);
  return ret;
}


int *makefoldremapper(args *arg,int pop1,int pop2){
  fprintf(stderr,"\t-> generating offset remapper lookup for flux cascade curves\n");
  int *mapper=NULL;
  if(pop1==pop2 &&arg->saf.size()==1){
    int ndim = arg->saf[pop1]->nChr+1;
    mapper = new int [ndim];
    for(int i=0;i<ndim;i++)
      mapper[i] = i;
    if(arg->fold==1){
      for(int i=0;i<ndim/2;i++)
	mapper[ndim-i-1] =i;
    }
    for(int i=0;0&&i<ndim;i++)
      fprintf(stderr,"%d:%d\n",i,mapper[i]);
  }else{
    int dimpop1 = arg->saf[pop1]->nChr+1;
    int dimpop2 = arg->saf[pop2]->nChr+1;
    int ndim = dimpop1*dimpop2;
    //fprintf(stderr,"ndim:%d (%lu,%lu)\n",ndim,dimpop1,dimpop2);
    
    int map[dimpop1][dimpop2];
    double tot=dimpop1+dimpop2-2;
    //fprintf(stderr,"tot:%.0f\n",tot);
    int inc=0;
    
    for(size_t x=0;x<dimpop1;x++)
      for(size_t y=0;y<dimpop2;y++)
	map[x][y] = inc++;
    mapper = new int [ndim];
    inc=0;
    
    for(size_t x=0;x<dimpop1;x++){
      for(size_t y=0;y<dimpop2;y++){
	double af= (x+y)/tot;
	//      fprintf(stderr,"x:%lu y:%lu af:%f ",x,y,af);
	if((arg->fold==1)&&(af>0.5)){
	  mapper[inc] = map[dimpop1-x-1][dimpop2-y-1];
	  //fprintf(stderr,"(%lu,%lu)->%d->%d\n",dimpop1-x-1,dimpop2-y-1,inc,mapper[inc]);
	}else{
	  mapper[inc] = inc;
	  //fprintf(stderr,"(%lu,%lu)->%d->%d\n",dimpop1-x-1,dimpop2-y-1,inc,mapper[inc]);
	}
	inc++;
      }
    }
#if 0
    for(int i=0;i<ndim;i++)
      fprintf(stdout,"%d %d\n",i,mapper[i]);
    // exit(0);
#endif
  }
  return mapper;
}


extern std::vector <char *> dumpedFiles;

void handler(int s) {
  if(s==13)//this is sigpipe
    exit(0);
  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want \'realSFS\' to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
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
  if(pars->saf.size()!=1){
    fprintf(stderr,"print_header only implemeted for single safs\n");
    exit(0);
  }
  writesaf_header(stdout,pars->saf[0]);
  
  destroy_args(pars);
  return 0;
}
template <typename T>
int printMulti(args *arg){
  //fprintf(stderr,"[%s]\n",__FUNCTION__);
  std::vector<persaf *> &saf =arg->saf;
  for(int i=0;i<saf.size();i++)
    assert(saf[i]->pos!=NULL&&saf[i]->saf!=NULL);

  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites=calc_nsites(saf,arg);
  }
  std::vector<Matrix<T> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<T>(nSites,saf[i]->nChr+1));

  int ndim=(int) parspace(saf);
  double *sfs=new double[ndim];
  
  //temp used for checking pos are in sync
  setGloc(saf,nSites);
  int *posiToPrint = new int[nSites];
  //used for printout old format
  FILE **oldfp = NULL;
  gzFile oldpos = Z_NULL;
  if(arg->oldout){
    oldfp = new FILE*[saf.size()];
    for(int i=0;i<saf.size();i++){
      size_t newlen = strlen(saf[i]->fname)+100;
      char *tmp =(char*) calloc(newlen,sizeof(char));
      tmp = strncpy(tmp,saf[i]->fname,strlen(saf[i]->fname)-4);
      fprintf(stderr,"\t-> Generating outputfile: %s\n",tmp);
      oldfp[i] = fopen(tmp,"wb");
      free(tmp);
    }
    size_t newlen = strlen(dirname(saf[0]->fname))+100;
    char *tmp = (char*) calloc(newlen,sizeof(char));
    snprintf(tmp,newlen,"%s/shared.pos.gz",dirname(saf[0]->fname));
    fprintf(stderr,"\t-> Generating outputfile: %s\n",tmp);
    oldpos = gzopen(tmp,"wb");
    free(tmp);
  }
  while(1) {
    static char *curChr=NULL;
    int ret=readdata(saf,gls,nSites,arg->chooseChr,arg->start,arg->stop,posiToPrint,&curChr,arg->fl,1);//read nsites from data
    if(arg->oldout==0){
      for(int s=0;s<gls[0]->x;s++){
	if(arg->chooseChr==NULL)
	  fprintf(stdout,"%s\t%d",curChr,posiToPrint[s]+1);
	else
	  fprintf(stdout,"%s\t%d",arg->chooseChr,posiToPrint[s]+1);
	for(int i=0;i<saf.size();i++)
	  for(int ii=0;ii<gls[i]->y;ii++)
	    fprintf(stdout,"\t%f",log(gls[i]->mat[s][ii]));
	fprintf(stdout,"\n");
      }
    }else{
      for(int s=0;s<gls[0]->x;s++){
	if(arg->chooseChr==NULL)
	  gzprintf(oldpos,"%s\t%d\n",curChr,posiToPrint[s]+1);
	else
	  gzprintf(oldpos,"%s\t%d\n",arg->chooseChr,posiToPrint[s]+1);
	for(int i=0;i<saf.size();i++){
	  double mytmp[gls[i]->y];
	  for(int ii=0;ii<gls[i]->y;ii++)
	    mytmp[ii] = log(gls[i]->mat[s][ii]);
	  fwrite(mytmp,sizeof(double),gls[i]->y,oldfp[i]);
	}
      }
    }
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
  delGloc(saf,nSites);
  destroy(gls,nSites);

  delete [] sfs;
  delete [] posiToPrint;

  if(arg->oldout==1){
    for(int i=0;i<saf.size();i++)
      fclose(oldfp[i]);
    delete [] oldfp;
    gzclose(oldpos);
  }
  destroy_args(arg);
  fprintf(stderr,"\t-> Run completed\n");
  return 0;
}


template<typename T>
void print(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"\t-> Must supply afile.saf.idx files \n");
    fprintf(stderr,"\t-> Examples \n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx \n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx -r chr1:10000000-12000000\n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx pop2.saf.idx -r chr2:10000000-12000000\n");
    fprintf(stderr,"\t-> You can generate the oldformat by appending the -oldout 1 to the print command like\n");
    fprintf(stderr,"\t-> ./realSFS print pop1.saf.idx pop2.saf.idx -oldout 1\n");
    return; 
  }
  
  args *pars = getArgs(argc,argv);
  for(int i=0;i<pars->saf.size();i++)
    pars->saf[0]->kind = 2;
  if(1||pars->saf.size()!=1){
    fprintf(stderr,"\t-> Will jump to multisaf printer and will only print intersecting sites between populations\n");
    printMulti<T>(pars);
    return;
  }

  writesaf_header(stderr,pars->saf[0]);
  
  float *flt = new float[pars->saf[0]->nChr+1];
  for(myMap::iterator it=pars->saf[0]->mm.begin();it!=pars->saf[0]->mm.end();++it){

    if(pars->chooseChr!=NULL)
      it = iter_init(pars->saf[0],pars->chooseChr,pars->start,pars->stop);
    else
      it = iter_init(pars->saf[0],it->first,pars->start,pars->stop);
 
    size_t ret;
    int pos;

    while((ret=iter_read(pars->saf[0],flt,sizeof(float)*(pars->saf[0]->nChr+1),&pos))){
      fprintf(stdout,"%s\t%d",it->first,pos+1);
      for(int is=0;is<pars->saf[0]->nChr+1;is++)
	fprintf(stdout,"\t%f",flt[is]);
      fprintf(stdout,"\n");
    }
 
    if(pars->chooseChr!=NULL)
      break;
  }
  
  delete [] flt;
  destroy_args(pars);
}


/*
  return value 
  -3 indicates that we are doing multi sfs and that we are totally and should flush

 */

int fst_index(int argc,char **argv){
  if(argc<1){
    fprintf(stderr,"Must supply afile.saf.idx [chrname, write more info]\n");
    return 0; 
  }
  args *arg = getArgs(argc,argv);
  if(!arg->outname){
    fprintf(stderr,"\t-> Must supply -fstout for doing fstindex\n");
    return 0;
  }

  std::vector<persaf *> &saf =arg->saf;
  //assert(saf.size()==2);
  size_t nSites = arg->nSites;
  if(nSites == 0){//if no -nSites is specified
    nSites = 100000;//<- set default to 100k sites, no need to load everything...
    // nSites=nsites(saf,arg);
  }
  fprintf(stderr,"\t-> nSites: %lu\n",nSites);
  std::vector<Matrix<float> *> gls;
  for(int i=0;i<saf.size();i++)
    gls.push_back(alloc<float>(nSites,saf[i]->nChr+1));

  //  int ndim= parspace(saf);
  if(arg->sfsfname.size()!=choose(saf.size(),2)){
    fprintf(stderr,"\t-> You have supplied: %lu populations, that is %d pairs\n",saf.size(),choose(saf.size(),2));
    fprintf(stderr,"\t-> You therefore need to supply %d 2dsfs priors instead of:%lu\n",choose(saf.size(),2),arg->sfsfname.size());
    exit(0);
  }
  fprintf(stderr,"\t-> IMPORTANT: please make sure that your saf files hasnt been folded with -fold 1 in -doSaf in angsd\n");
  
  double **a1,**b1;
  int **remaps,**remaps_scaling;
  a1=new double*[choose(saf.size(),2)];
  b1=new double*[choose(saf.size(),2)];
  remaps=new int*[choose(saf.size(),2)];
  remaps_scaling=new int*[choose(saf.size(),2)];
  int inc=0;
  for(int i=0;i<saf.size();i++)
    for(int j=i+1;j<saf.size();j++){
      calcCoef((int)saf[i]->nChr,(int)saf[j]->nChr,&a1[inc],&b1[inc],arg->whichFst);
      remaps[inc] = makefoldremapper(arg,i,j);
      remaps_scaling[inc] = makefoldadjust(remaps[inc],(arg->saf[i]->nChr+1)*(arg->saf[j]->nChr+1));
      //      fprintf(stderr,"a1[%d]:%p b1[%d]:%p\n",inc,&a1[inc][0],inc,&b1[inc][0]);
      inc++;
    }


  inc =0;
  std::vector<double *> sfs;
  for(int i=0;i<saf.size();i++)
    for(int j=i+1;j<saf.size();j++){
      size_t pairdim = (saf[i]->nChr+1)*(saf[j]->nChr+1);
      double *ddd=new double[pairdim];
      readSFS(arg->sfsfname[inc],pairdim,ddd);
      if(arg->fold==1){
#if 0
	fprintf(stdout,"SFSIN");
	  for(int i=0;i<ndim;i++)
	    fprintf(stdout,"%f\t",ddd[i]);
	  fprintf(stdout,"\n");
#endif
	  //this block is debug block for validating that the input sfs gets folded correctly
	  double newtmp[pairdim];
	  for(int i=0;i<pairdim;i++)
	    newtmp[i] = 0;
	  for(int i=0;i<pairdim;i++){
	    newtmp[remaps[inc][i]] += ddd[i];
	    //	fprintf(stderr,"%d %f\n",i,newtmp[i]);
	  }
	  for(int i=0;i<pairdim;i++){
	    //	fprintf(stdout,"%f ",newtmp[i]);
	    ddd[i] = newtmp[i];
	  }
      }
      
      normalize(ddd,pairdim);
      sfs.push_back(ddd);
      inc++;
    }

  BGZF *fstbg = openFileBG(arg->outname,".fst.gz");
  FILE *fstfp = openFile(arg->outname,".fst.idx");
  char buf[8]="fstv1";
  if(bgzf_write(fstbg,buf,8)!=8){
    fprintf(stderr,"\t-> Problem writing 8bytes into output file\n");
    exit(0);
  }
  fwrite(buf,1,8,fstfp);
#if 0
  for(int i=0;i<ndim;i++)
    fprintf(stdout,"%f %f\n",a1[i],b1[i]);
  exit(0);
#endif
#if 1
  size_t nsafs=saf.size();
  fwrite(&nsafs,sizeof(size_t),1,fstfp);
  for(int i=0;i<nsafs;i++){
    size_t clen= strlen(saf[i]->fname);
    fwrite(&clen,sizeof(size_t),1,fstfp);
    fwrite(saf[i]->fname,1,clen,fstfp);
  }
#endif
  int asdf = choose(saf.size(),2);
  std::vector<double> *ares = new std::vector<double> [choose(saf.size(),2)];
  std::vector<double> *bres = new std::vector<double> [choose(saf.size(),2)];
  //  for(int i=0;i<3;i++)
    //    fprintf(stderr,"ares.size():%lu bres.size():%lu sfs:%p\n",ares[i].size(),bres[i].size(),&sfs[i][0]);
  std::vector<int> posi;
  setGloc(saf,nSites);
  int *posiToPrint = new int[nSites];
  for(myMap::iterator it = saf[0]->mm.begin();it!=saf[0]->mm.end();++it) {
    //    fprintf(stderr,"doing chr:%s\n",it->first);
    if(arg->chooseChr!=NULL){
      it = saf[0]->mm.find(arg->chooseChr);
      if(it==saf[0]->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",arg->chooseChr);
	break;
      }
    }else{
      int efsize=0;
      for(int i=0;i<saf.size();i++){
	myMap::iterator it2 = saf[i]->mm.find(it->first);
	if(it2!=saf[i]->mm.end())
	  efsize++;
      }
      if(efsize!=saf.size())
	continue;
    }
    for(int i=0;i<choose(saf.size(),2);i++){
      ares[i].clear();
      bres[i].clear();
    }
    posi.clear();
    while(1) {
      int ret=readdata(saf,gls,nSites,it->first,arg->start,arg->stop,posiToPrint,NULL,arg->fl,1);//read nsites from data
      //  fprintf(stderr,"ret:%d glsx:%lu\n",ret,gls[0]->x);
      //if(gls[0]->x!=nSites&&arg->chooseChr==NULL&&ret!=-3){
	//fprintf(stderr,"continue continue\n");
      //	continue;
      //}
      
      fprintf(stderr,"\t-> Will now do fst temp dump using a chunk of %lu\n",gls[0]->x);
      int inc=0;
      for(int i=0;i<saf.size();i++)
	for(int j=i+1;j<saf.size();j++){
	  //	  fprintf(stderr,"i:%d j:%d inc:%d gls[i]:%p gls[j]:%p sfs:%p a1:%p b1:%p\n",i,j,inc,gls[i],gls[j],sfs[i],&a1[inc][0],&a1[inc][0]);
	  block_coef(gls[i],gls[j],sfs[inc],a1[inc],b1[inc],ares[inc],bres[inc],remaps[inc],remaps_scaling[inc]);
	  inc++;
	}
      for(int i=0;i<gls[0]->x;i++)
	posi.push_back(posiToPrint[i]);

      for(int i=0;i<gls.size();i++)
	gls[i]->x =0;
      if(ret==-2)//no more data in files or in chr, eith way we break;
	break;
    }
    size_t clen = strlen(it->first);
    fwrite(&clen,sizeof(size_t),1,fstfp);
    fwrite(it->first,1,clen,fstfp);
    size_t nit=posi.size();

    assert(1==fwrite(&nit,sizeof(size_t),1,fstfp));
    int64_t tell = bgzf_tell(fstbg);
    fwrite(&tell,sizeof(int64_t),1,fstfp);
    my_bgzf_write(fstbg,&posi[0],posi.size()*sizeof(int));
    int inc =0;
    for(int i=0;i<saf.size();i++)
      for(int j=i+1;j<saf.size();j++){
	my_bgzf_write(fstbg,&(ares[inc][0]),ares[inc].size()*sizeof(double));
	my_bgzf_write(fstbg,&(bres[inc][0]),bres[inc].size()*sizeof(double));
	inc++;
      }
    if(arg->chooseChr!=NULL)
      break;
  }
  delGloc(saf,nSites);
  destroy(gls,nSites);
  destroy_args(arg);
  for(int i=0;i<sfs.size();i++)
    delete [] sfs[i];
#if 0
  fprintf(stderr,"\n\t-> NB NB output is no longer log probs of the frequency spectrum!\n");
  fprintf(stderr,"\t-> Output is now simply the expected values! \n");
  fprintf(stderr,"\t-> You can convert to the old format simply with log(norm(x))\n");
#endif
  bgzf_close(fstbg);
  fclose(fstfp);
  fprintf(stderr,"\t-> fst index finished with no errors!\n");
  return 0;
}

int fst(int argc,char**argv){
  if(argc==0){
    fprintf(stderr,"\t-> Possible options: index print\n");
    return 0;
  }
  if(!strcasecmp(*argv,"index"))  
    fst_index(--argc,++argv);
  else  if(!strcasecmp(*argv,"print"))  
    fst_print(--argc,++argv);
  else if(!strcasecmp(*argv,"stats"))  
    fst_stat(--argc,++argv);
  else if(!strcasecmp(*argv,"stats2"))  
    fst_stat2(--argc,++argv);
  else{
    fprintf(stderr,"unknown option: \'%s\'\n",*argv);
  }
  return 0;
}


void writeAllThetas(BGZF *dat,FILE *idx,char *tmpChr,int64_t &offs,std::vector<int> &p,std::vector<float> *res,int nChr){
  assert(dat!=NULL);
  assert(idx!=NULL);
  assert(p.size()==res[0].size());
  fprintf(stderr,"\t-> Writing %lu sites for chr:%s\n",p.size(),tmpChr);
  for(int i=1;i<5;i++)
    assert(p.size()==res[i].size());//DRAGON, might be discarded during compilation
      
  if(p.size()!=0&&tmpChr!=NULL){
    //write clen and chromoname for both index and bgzf
    size_t clen = strlen(tmpChr);
    fwrite(&clen,sizeof(size_t),1,idx);
    fwrite(tmpChr,1,clen,idx);
    if(sizeof(size_t)!=bgzf_write(dat,&clen,sizeof(size_t))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    if(clen!=bgzf_write(dat,tmpChr,clen)){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }

    //write number of sites for both index and bgzf
    size_t tt = p.size();
    fwrite(&tt,sizeof(size_t),1,idx);
    if(sizeof(size_t)!=bgzf_write(dat,&tt,sizeof(size_t))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    //write nChr for both index and bgzf
    fwrite(&nChr,sizeof(int),1,idx);
    if(sizeof(int)!=bgzf_write(dat,&nChr,sizeof(int))){
      fprintf(stderr,"\t-> Problems writing theta files\n");
      exit(0);
    }
    //write bgzf offset into idx
    fwrite(&offs,sizeof(int64_t),1,idx);
    for(int i=0;i<p.size();i++)
      aio::bgzf_write(dat,&p[i],sizeof(int));

    for(int i=0;i<5;i++)
      for(int j=0;j<p.size();j++)
	aio::bgzf_write(dat,&res[i][j],sizeof(float));
  }

  //reset
  offs = bgzf_tell(dat);
  p.clear();
  for(int i=0;i<5;i++)
    res[i].clear();
}


int saf2theta(int argc,char**argv){
  const char *THETAS =".thetas.gz";
  const char *THETASIDX =".thetas.idx";
  BGZF *theta_dat;
  FILE *theta_idx;
  std::vector<float> theta_res[5];// = new std::vector<float>[5];
  std::vector<int> theta_pos;

  if(argc==0){
    fprintf(stderr,"\t-> Possible options: addoptions\n");
    return 0;
  }
  args *arg = getArgs(argc,argv);
  if(arg->outname==NULL){
    fprintf(stderr,"\t-> Must supply -outname for generating outputfiles\n");
    return 0;
  }
  if(arg->saf.size()!=1){
    fprintf(stderr,"\t-> Must supply one, and only one saf.idx file\n");
    return 0;
  }
  if(arg->sfsfname.size()!=1){
    fprintf(stderr,"\t-> Must supply one, and only one -sfs argument which should contain the prior\n");
    return 0;
  }
  if(arg->nSites==0){
    int block = 4096;
    fprintf(stderr,"\t-> Will read chunks of size: %d\n",block);
    arg->nSites = block;
  }
  arg->saf[0]->kind=2; //<-important orhterwise we dont read positions from saffles\n
  double *prior = new double[arg->saf[0]->nChr+1];
  readSFS(arg->sfsfname[0],arg->saf[0]->nChr+1,prior);
  for(int i=0;i<arg->saf[0]->nChr+1;i++)
    prior[i] = log(prior[i]);
  char buf[8] = "thetav2";
  theta_dat = aio::openFileBG(arg->outname,THETAS);
  theta_idx = aio::openFile(arg->outname,THETASIDX);
  aio::bgzf_write(theta_dat,buf,8);
  fwrite(buf,1,8,theta_idx);
  //  theta_res = new std::vector<float>[5];
  int64_t offs_thetas = bgzf_tell(theta_dat);
  
  
  double aConst=0;
  int nChr = arg->saf[0]->nChr;
  fprintf(stderr,"\t-> nChr:%d\n",nChr);
  for(int i=1;i<nChr;i++)
    aConst += 1.0/i;
  aConst = log(aConst);//this is a1
  
  
  double aConst2 = log((nChr*(nChr-1))/2.0);//choose(nChr,2)
  double aConst3 = log((1.0*nChr-1.0));
  
  double *scalings = new double [nChr+1];
  for(int i=0;i<nChr+1;i++)
    scalings[i] = log(i)+log(nChr-i);

  std::vector<Matrix<float> *> gls;
  for(int i=0;i<arg->saf.size();i++)
    gls.push_back(alloc<float>(arg->nSites,arg->saf[i]->nChr+1));
  
  setGloc(arg->saf,arg->nSites);
  int *posiToPrint = new int[arg->nSites];

  char *tmpChr =NULL;
  static char *curChr=NULL;//why static?

  while(1) {
    int ret=readdata(arg->saf,gls,arg->nSites,arg->chooseChr,arg->start,arg->stop,posiToPrint,&curChr,arg->fl,0);//read nsites from data
    if(arg->chooseChr!=NULL){
      if(curChr==NULL)
	curChr=strdup(arg->chooseChr);
    }
    if(tmpChr==NULL)
      tmpChr = strdup(curChr);
    if(strcmp(tmpChr,curChr)!=0){
      writeAllThetas(theta_dat,theta_idx,tmpChr,offs_thetas,theta_pos,theta_res,nChr);
      free(tmpChr);
      tmpChr=strdup(curChr);
    }
    //calc thetas

    for(int s=0;s<gls[0]->x;s++){
      double workarray[nChr+1];
      for(int i=0;i<nChr+1;i++)//gls->mat is float lets pluginto double
	workarray[i] = gls[0]->mat[s][i];

      //calculate post probs
      double tsum =exp(workarray[0] + prior[0]);
      for(int i=1;i<nChr+1;i++)
	tsum += exp(workarray[i]+prior[i]);
      tsum = log(tsum);
      
      for(int i=0;i<nChr+1;i++){
	workarray[i] = workarray[i]+prior[i]-tsum;
	//      fprintf(stderr,"[%d]:%f\n",i,(workarray[i]));
      }
      //exit(0);
      //First find thetaW: nSeg/a1
      double pv,seq;
      if(arg->fold)
	pv = 1-exp(workarray[0]);
      else
	pv = 1-exp(workarray[0])-exp(workarray[nChr]);
      //      fprintf(stderr,"pv:%f work[0]:%f 2k:%f\n",pv,workarray[0],workarray[nChr]);
      //      exit(0);
      if(pv<0)//catch underflow
	seq=log(0.0);
      else
	seq = log(pv)-aConst;//watterson
      theta_res[0].push_back(seq);
      //     ksprintf(&kb,"%s\t%d\t%f\t",header->target_name[pars->refId],pars->posi[i]+1,seq);
      theta_pos.push_back(posiToPrint[s]);
      // fprintf(stderr,"posiToPrint[s]:%d\n",posiToPrint[s]);
      if(arg->fold==0) {
	double pairwise=0;    //Find theta_pi the pairwise
	double thL=0;    //Find thetaL sfs[i]*i;
	double thH=0;//thetaH sfs[i]*i^2
	for(size_t ii=1;ii<nChr;ii++){
	  
	  pairwise += exp(workarray[ii]+scalings[ii] );
	  double li=log(ii);
	  
	  thL += exp(workarray[ii])*ii;
	  thH += exp(2*li+workarray[ii]);
	}
	theta_res[1].push_back(log(pairwise)-aConst2);
	theta_res[2].push_back(workarray[1]);
	theta_res[3].push_back(log(thH)-aConst2);
	theta_res[4].push_back(log(thL)-aConst3);
      }else{
	double pairwise=0;    //Find theta_pi the pairwise
	for(size_t ii=1;ii<nChr+1;ii++)
	  pairwise += exp(workarray[ii]+scalings[ii] );
	theta_res[1].push_back(log(pairwise)-aConst2);
	for(int i=2;i<=4;i++)
	  theta_res[i].push_back(log(0));
      }
    }
    
#if 0 //just for printout
    for(int s=0;s<gls[0]->x;s++){
      if(arg->chooseChr==NULL)
	fprintf(stdout,"%s\t%d",curChr,posiToPrint[s]+1);
      else
	fprintf(stdout,"%s\t%d",arg->chooseChr,posiToPrint[s]+1);
      for(int i=0;i<arg->saf.size();i++)
	for(int ii=0;ii<gls[i]->y;ii++)
	  fprintf(stdout,"\t%f",log(gls[i]->mat[s][ii]));
      fprintf(stdout,"\n");
    }
#endif
    if(ret==-3&&gls[0]->x==0){//no more data in files or in chr, eith way we break;g
      //      fprintf(stderr,"breaking change of chr\n");
      break;
    }
    for(int i=0;i<gls.size();i++)
      gls[i]->x =0;
    
    if(ret==-2&&arg->chooseChr!=NULL)
      break;
    if(arg->onlyOnce)
      break;
  }
  if(theta_pos.size()>0)
    writeAllThetas(theta_dat,theta_idx,tmpChr,offs_thetas,theta_pos,theta_res,nChr);
  
  fclose(theta_idx);
  bgzf_close(theta_dat);
  if(curChr){
    //  free(curChr);//<- DRAGON leak?
    curChr=NULL;
  }
  if(tmpChr)
    free(tmpChr);
 
  delGloc(arg->saf,arg->nSites);
  destroy(gls,arg->nSites);
  delete [] prior;
  //  delete [] theta_res;
  delete [] scalings;
  delete [] posiToPrint;
  destroy_args(arg);
  return 0;
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
    //    fprintf(stderr, "\t->------------------\n\t-> ./realSFS\n\t->------------------\n");
    // fprintf(stderr,"\t-> This is the new realSFS program which works on the newer binary files from ANGSD!!\n");
    fprintf(stderr, "\t-> ---./realSFS------\n\t-> EXAMPLES FOR ESTIMATING THE (MULTI) SFS:\n\n\t-> Estimate the SFS for entire genome??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx \n");
    fprintf(stderr, "\n\t-> 1) Estimate the SFS for entire chromosome 22 ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -r chr22 \n");
    fprintf(stderr, "\n\t-> 2) Estimate the 2d-SFS for entire chromosome 22 ??\n");
    fprintf(stderr,"\t-> ./realSFS afile1.saf.idx  afile2.saf.idx -r chr22 \n");

    fprintf(stderr, "\n\t-> 3) Estimate the SFS for the first 500megabases (this will span multiple chromosomes) ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -nSites 500000000 \n");

    fprintf(stderr, "\n\t-> 4) Estimate the SFS around a gene ??\n");
    fprintf(stderr,"\t-> ./realSFS afile.saf.idx -r chr2:135000000-140000000 \n");
    fprintf(stderr, "\n\t-> Other options [-P nthreads -tole tolerence_for_breaking_EM -maxIter max_nr_iterations -bootstrap number_of_replications -resample_chr 0/1]\n");

    fprintf(stderr,"\n\t-> See realSFS print for possible print options\n");
    fprintf(stderr,"\t-> Use realSFS print_header for printing the header\n");
    fprintf(stderr,"\t-> Use realSFS cat for concatenating saf files\n");

    fprintf(stderr,"\n\t->------------------\n\t-> NB: Output is now counts of sites instead of log probs!!\n");
    fprintf(stderr,"\t-> NB: You can print data with ./realSFS print afile.saf.idx !!\n");
    fprintf(stderr,"\t-> NB: Higher order SFS's can be estimated by simply supplying multiple .saf.idx files!!\n");
    fprintf(stderr,"\t-> NB: Program uses accelerated EM, to use standard EM supply -m 0 \n");
    fprintf(stderr,"\t-> Other subfunctions saf2theta, cat, check, dadi\n");
    return 0;
  }
  ++argv;
  --argc;
  if(!strcasecmp(*argv,"print"))
    print<float>(--argc,++argv);
  else if(!strcasecmp(*argv,"cat"))
    saf_cat(--argc,++argv);
  else if(!strcasecmp(*argv,"fst"))
    fst(--argc,++argv);
  else if(!strcasecmp(*argv,"dadi"))
    main_dadi<float>(--argc,++argv);
  else if(!strcasecmp(*argv,"saf2theta"))
    saf2theta(--argc,++argv);
  else if(!strcasecmp(*argv,"check"))
    saf_check(--argc,++argv);
  else if(!strcasecmp(*argv,"print_header"))
    print_header(--argc,++argv);
  else {
    args *arg = getArgs(argc,argv);
    if(arg->saf.size()>1)
      fprintf(stderr,"\t-> Multi SFS is 'still' under development. Please report strange behaviour\n");
    if(!arg)
      return 0;

    if(isatty(fileno(stdout))){
      fprintf(stderr,"\t-> You are printing the optimized SFS to the terminal consider dumping into a file\n");
      fprintf(stderr,"\t-> E.g.: \'./realSFS");
      for(int i=0;i<argc;i++)
	fprintf(stderr," %s",argv[i]);
      fprintf(stderr," >sfs.ml.txt\'\n");   
    }
  

    main_opt<float>(arg);
    
  }
  extern size_t *bootstrap;
  if(bootstrap!=NULL)
    delete [] bootstrap;

  if(dumpedFiles.size()){
    fprintf(stderr,"\t-> Output filenames:\n");
    for(int i=0;i<(int)dumpedFiles.size();i++){
      fprintf(stderr,"\t\t->\"%s\"\n",dumpedFiles[i]);
      free(dumpedFiles[i]);
    }
  }
  return 0;
}

