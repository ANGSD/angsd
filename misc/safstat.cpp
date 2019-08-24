#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <htslib/bgzf.h>

#include <cstring>
#include "safstat.h"
#include "realSFS_args.h"
//#include "realSFS.h"
//lazy binomial
int choose(int n,int m){
  if(n==2&&m==2)
    return 1;
  else if(n==3&&m==2)
    return 3;
  else{
    fprintf(stderr,"\t-> Never here\n");
    exit(0);
  }
  return -1;
}
int choose(size_t n,int m){
  return choose((int)n,m);
}
void calcpbs(double fstW[3]){
  double p01 = -log(1.0-fstW[0]);
  double p02 = -log(1.0-fstW[1]);
  double p12 = -log(1.0-fstW[2]);
  fstW[0] = (p01+p02-p12)/2.0;
  fstW[1] = (p12+p01-p02)/2.0;
  fstW[2] = (p12+p02-p01)/2.0;


}

void reynoldFst(int sfs1,int sfs2,double **aMat,double **bMat){
  fprintf(stderr,"\t-> [%s] sfs1:%d sfs2:%d dimspace:%d \n",__FUNCTION__,sfs1,sfs2,(sfs1+1)*(sfs2+1));
  *aMat = new double[(sfs1+1)*(sfs2+1)];
  *bMat = new double[(sfs1+1)*(sfs2+1)];
  int at=0;
  for(int a1=0;a1<=sfs1;a1++)
    for(int a2=0;a2<=sfs2;a2++){
      double p1 = 1.0 * a1/(1.0*sfs1);
      double p2 = 1.0 * a2/(1.0*sfs2);
      double q1 = 1 - p1;
      double q2 = 1 - p2;
      double alpha1 = 1 - (p1*p1 + q1*q1);
      double alpha2 = 1 - (p2*p2 + q2*q2);
      
      double al =  0.5 * ( pow(p1-p2,2.0) + pow(q1-q2,2)) - (sfs1+sfs2) *  (sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      double bal = 0.5 * ( pow(p1-p2,2) + pow(q1-q2,2)) + (4*sfs1*sfs2-sfs1-sfs2)*(sfs1*alpha1 + sfs2*alpha2) / (4*sfs1*sfs2*(sfs1+sfs2-1));
      (*aMat)[at] = al;
      (*bMat)[at] = bal;
      //      fprintf(stderr,"p1:%f p2:%f q1:%f q2:%f alhpa1:%f alpha:%f al:%f bal:%f\n",p1,p2,q1,q2,alpha1,alpha2,al,bal);
      at++;
    }
}
void bhatiaFst(int sfs1,int sfs2,double **aMat,double **bMat){
  fprintf(stderr,"\t-> [%s] sfs1:%d sfs2:%d dimspace:%d \n",__FUNCTION__,sfs1,sfs2,(sfs1+1)*(sfs2+1));
  *aMat = new double[(sfs1+1)*(sfs2+1)];
  *bMat = new double[(sfs1+1)*(sfs2+1)];
  int at=0;
  for(int a1=0;a1<=sfs1;a1++)
    for(int a2=0;a2<=sfs2;a2++){
      double p1 = 1.0 * a1/(1.0*sfs1);
      double p2 = 1.0 * a2/(1.0*sfs2);

      //sample size correction
      (*aMat)[at] = (p1-p2)*(p1-p2)-((p1*(1.0-p1))/((double)sfs1-1))-((p2*(1.0-p2))/((double)sfs2-1.0));
      (*bMat)[at] = p1*(1.0-p2)+p2*(1.0-p1);
      //fprintf(stderr,"(%d,%d) p1:%f p2:%f al:%f bal:%f\n",a1,a2,p1,p2,(*aMat)[at],(*bMat)[at]);
      at++;
    }
}
void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat,int whichFst){
  if(whichFst==0)
    reynoldFst(sfs1,sfs2,aMat,bMat);
  else if(whichFst==1)
    bhatiaFst(sfs1,sfs2,aMat,bMat);
  else{
    fprintf(stderr,"\t-> Fst option is not implemented\n");
    exit(0);
  }

}

void block_coef(Matrix<float > *gl1,Matrix<float> *gl2,double *prior,double *a1,double *b1,std::vector<double> &ares,std::vector<double> &bres,int *remap,int *rescal){
  assert(prior!=NULL);
  double tre[3]={0,0,0};//a/b,sum(a),sum(0)
  for(int s=0;s<gl1->x;s++){
    int inc =0 ;
    double tmp[(gl1->y+1)*(gl2->y+1)];
    for(int jj=0;jj<(gl1->y+1)*(gl2->y+1);jj++)
      tmp[jj] = 0;
    for(int i=0;i<gl1->y;i++)
      for(int j=0;j<gl2->y;j++){
	tmp[remap[inc]] += prior[remap[inc]]* gl1->mat[s][i] *gl2->mat[s][j]*rescal[inc];
	inc++;
      }
    //    exit(0);
    double as=0;
    double bs=0;
    void normalize(double *,size_t);
    normalize(tmp,inc);
    for(int i=0;i<inc;i++){
      as += a1[i]*tmp[i];
      bs += b1[i]*tmp[i];
    }
    tre[0] += as/bs;
    tre[1] += as;
    tre[2] += bs;
    ares.push_back(as);
    bres.push_back(bs);
      
    //    fprintf(stdout,"%f %f\n",as,bs);
    
  }
  //  fprintf(stderr,"bres.size:%lu\n",bres.size());
  // fprintf(stderr,"u:%f w:%f\n",tre[0]/(1.0*gl1->x),tre[1]/tre[2]);
}

#include "fstreader.h"
int fst_print(int argc,char **argv){

  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
  perfst *pf = perfst_init(bname);
  writefst_header(stderr,pf);  
  args *pars = getArgs(--argc,++argv);  
  int *ppos = NULL;
  fprintf(stderr,"choose:%d \n",choose((int)pf->names.size(),2));
  double **ares = new double*[choose((int)pf->names.size(),2)];
  double **bres = new double*[choose((int)pf->names.size(),2)];
  for(myFstMap::iterator it=pf->mm.begin();it!=pf->mm.end();++it){
    if(pars->chooseChr!=NULL){
      it = pf->mm.find(pars->chooseChr);
      if(it==pf->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",pars->chooseChr);
	break;
      }
    }
    if(it->second.nSites==0)
      continue;
    my_bgzf_seek(pf->fp,it->second.off,SEEK_SET);
    ppos = new int[it->second.nSites];
    
    my_bgzf_read(pf->fp,ppos,sizeof(int)*it->second.nSites);
    for(int i=0;i<choose((int)pf->names.size(),2);i++){
      ares[i] = new double[it->second.nSites];
      bres[i] = new double[it->second.nSites];
      my_bgzf_read(pf->fp,ares[i],sizeof(double)*it->second.nSites);
      my_bgzf_read(pf->fp,bres[i],sizeof(double)*it->second.nSites);
    }
    


    size_t first=0;
    if(pars->start!=-1)
      while(ppos[first]<pars->start) 
	first++;
    
    size_t last=it->second.nSites;

    if(pars->stop!=-1&&pars->stop<=ppos[last-1]){
      last=first;
      while(ppos[last]<pars->stop) 
	last++;
    }

    fprintf(stderr,"pars->stop:%d ppos:%d first:%lu last:%lu\n",pars->stop,ppos[last-1],first,last);

    for(size_t s=first;s<last;s++){
      fprintf(stdout,"%s\t%d",it->first,ppos[s]+1);
      for(int i=0;i<choose((int)pf->names.size(),2);i++)
	fprintf(stdout,"\t%f\t%f",ares[i][s],bres[i][s]);
      fprintf(stdout,"\n");
    }
    for(int i=0;i<choose((int)pf->names.size(),2);i++){
      delete [] ares[i];
      delete [] bres[i];
    }
    
    delete [] ppos;
    
    if(pars->chooseChr!=NULL)
      break;
  }
  delete [] ares;
  delete [] bres;
  destroy_args(pars);
  perfst_destroy(pf);
  return 0;
}

int fst_stat2(int argc,char **argv){
  int pS,pE;//physical start,physical end
  int begI,endI;//position in array for pS, pE;
  
  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
  perfst *pf = perfst_init(bname);
  args *pars = getArgs(--argc,++argv);
  fprintf(stderr,"win:%d step:%d\n",pars->win,pars->step);
  int *ppos = NULL;
  int chs = choose((int)pf->names.size(),2);
  // fprintf(stderr,"choose:%d \n",chs);
  double **ares = new double*[chs];
  double **bres = new double*[chs];
  double unweight[chs];
  double wa[chs];
  double wb[chs];
  size_t nObs =0;
 

  //print header
  fprintf(stdout,"region\tchr\tmidPos\tNsites");
  for(int c1=0;c1<chs-1;c1++)
    for(int c2=c1+1;c2<chs;c2++)
      fprintf(stdout,"\tFst%d%d",c1,c2);
  if(chs==3)
    fprintf(stdout,"\tPBS0\tPBS1\tPBS2");
  
  fprintf(stdout,"\n");


  for(myFstMap::iterator it=pf->mm.begin();it!=pf->mm.end();++it){
    if(pars->chooseChr!=NULL){
      it = pf->mm.find(pars->chooseChr);
      if(it==pf->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",pars->chooseChr);
	break;
      }
    }
    fprintf(stderr,"nSites:%lu\n",it->second.nSites);
    if(it->second.nSites==0&&pars->chooseChr!=NULL)
      break;
    else if(it->second.nSites==0&&pars->chooseChr==NULL)
      continue;
    my_bgzf_seek(pf->fp,it->second.off,SEEK_SET);
    ppos = new int[it->second.nSites];
    
    my_bgzf_read(pf->fp,ppos,sizeof(int)*it->second.nSites);
    for(int i=0;i<it->second.nSites;i++)//what? why? dragon!
      ppos[i]++;

    for(int i=0;i<choose((int)pf->names.size(),2);i++){
      ares[i] = new double[it->second.nSites];
      bres[i] = new double[it->second.nSites];
      my_bgzf_read(pf->fp,ares[i],sizeof(double)*it->second.nSites);
      my_bgzf_read(pf->fp,bres[i],sizeof(double)*it->second.nSites);
    }
    

    if(pars->type==0)
      pS = ((pars->start!=-1?pars->start:ppos[0])/pars->step)*pars->step +pars->step;
    else if(pars->type==1)
      pS = (pars->start!=-1?pars->start:ppos[0]);
    else if(pars->type==2)
      pS = 1;
    pE = pS+pars->win;
    begI=endI=0;

    //fprintf(stderr,"ps:%d\n",pS);exit(0);
    if(pE>(pars->stop!=-1?pars->stop:ppos[it->second.nSites-1])){
    fprintf(stderr,"end of dataset is before end of window: end of window:%d last position in chr:%d\n",pE,ppos[it->second.nSites-1]);
    //    return str;
  }

  while(ppos[begI]<pS) begI++;
  
  endI=begI;
  while(endI<it->second.nSites &&ppos[endI]<pE) endI++;

  //fprintf(stderr,"begI:%d endI:%d\n",begI,endI);

  while(1){
    for(int i=0;i<chs;i++)
      unweight[i] = wa[i] = wb[i] =0.0;
    nObs=0;

    for(int s=begI;s<endI;s++){
#if 0
      fprintf(stdout,"%s\t%d",it->first,ppos[s]+1);
      for(int i=0;i<choose(pf->names.size(),2);i++)
	fprintf(stdout,"\t%f\t%f",ares[i][s],bres[i][s]);
      fprintf(stdout,"\n");
#endif
      for(int i=0;i<choose((int)pf->names.size(),2);i++){
	unweight[i] += ares[i][s]/bres[i][s];
	wa[i] += ares[i][s];
	wb[i] += bres[i][s];
      }
      nObs++;
    }
    if(nObs>0)
      fprintf(stdout,"(%d,%d)(%d,%d)(%d,%d)\t%s\t%d\t%d",begI,endI-1,ppos[begI],ppos[endI-1],pS,pE,it->first,pS+(pE-pS)/2,endI-begI+1);
    double fstW[chs];
    for(int i=0;nObs>0&&i<chs;i++){
      fstW[i] = wa[i]/wb[i];
      //      fprintf(stdout,"\t%f\t%f",unweight[i]/(1.0*nObs),fstW[i]);
      fprintf(stdout,"\t%f",fstW[i]);
    }
    if(nObs>0&&chs==3){
      //if chr==3 then we have 3pops and we will also calculate pbs statistics
      calcpbs(fstW);//<- NOTE: the pbs values will replace the fstW values
      for(int i=0;i<3;i++)
	fprintf(stdout,"\t%f",fstW[i]);
    }
    if(nObs>0)
      fprintf(stdout,"\n");

    pS += pars->step;
    pE =pS+pars->win;
    if(pE>(pars->stop!=-1?pars->stop:ppos[it->second.nSites-1]))
      break;
    
    while(ppos[begI]<pS) begI++;
    while(ppos[endI]<pE) endI++;
  }
  for(int i=0;i<choose((int)pf->names.size(),2);i++){
      delete [] ares[i];
      delete [] bres[i];
    }
    
    delete [] ppos;
    
    if(pars->chooseChr!=NULL)
      break;
  }
 
  delete [] ares;
  delete [] bres;
  destroy_args(pars);
  perfst_destroy(pf);
  return 0;
}



int fst_stat(int argc,char **argv){
  
  char *bname = *argv;
  fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
  perfst *pf = perfst_init(bname);
  args *pars = getArgs(--argc,++argv);  
  int *ppos = NULL;
  int chs = choose((int)pf->names.size(),2);
  // fprintf(stderr,"choose:%d \n",chs);
  double **ares = new double*[chs];
  double **bres = new double*[chs];
  double unweight[chs];
  double wa[chs];
  double wb[chs];
  size_t nObs[chs];
  for(int i=0;i<chs;i++){
    unweight[i] = wa[i] = wb[i] =0.0;
    nObs[i] = 0;
  }
  for(myFstMap::iterator it=pf->mm.begin();it!=pf->mm.end();++it){
    if(pars->chooseChr!=NULL){
      it = pf->mm.find(pars->chooseChr);
      if(it==pf->mm.end()){
	fprintf(stderr,"Problem finding chr: %s\n",pars->chooseChr);
	break;
      }
    }
    if(it->second.nSites==0)
      continue;
    my_bgzf_seek(pf->fp,it->second.off,SEEK_SET);
    ppos = new int[it->second.nSites];
    
    my_bgzf_read(pf->fp,ppos,sizeof(int)*it->second.nSites);
    for(int i=0;i<choose((int)pf->names.size(),2);i++){
      ares[i] = new double[it->second.nSites];
      bres[i] = new double[it->second.nSites];
      my_bgzf_read(pf->fp,ares[i],sizeof(double)*it->second.nSites);
      my_bgzf_read(pf->fp,bres[i],sizeof(double)*it->second.nSites);
    }
    


    size_t first=0;
    if(pars->start!=-1)
      while(ppos[first]<pars->start) 
	first++;
    
    size_t last=it->second.nSites;

    if(pars->stop!=-1&&pars->stop<=ppos[last-1]){
      last=first;
      
      while(ppos[last]<pars->stop) 
	last++;
    }

    //  fprintf(stderr,"pars->stop:%d ppos:%d first:%d last:%d\n",pars->stop,ppos[last-1],first,last);

    for(size_t s=first;s<last;s++){
#if 0
      fprintf(stdout,"%s\t%d",it->first,ppos[s]+1);
      for(int i=0;i<choose(pf->names.size(),2);i++)
	fprintf(stdout,"\t%f\t%f",ares[i][s],bres[i][s]);
      fprintf(stdout,"\n");
#endif
      for(int i=0;i<choose((int)pf->names.size(),2);i++){
	if(bres[i][s]!=0){
	  unweight[i] += ares[i][s]/bres[i][s];
	  nObs[i]++;
	}
	wa[i] += ares[i][s];
	wb[i] += bres[i][s];
      }
    }
    for(int i=0;i<choose((int)pf->names.size(),2);i++){
      delete [] ares[i];
      delete [] bres[i];
    }
    
    delete [] ppos;
    
    if(pars->chooseChr!=NULL)
      break;
  }
  double fstUW[chs];
  double fstW[chs];
  for(int i=0;i<chs;i++){
    fstUW[i] = unweight[i]/(1.0*((double)nObs[i]));
    fstW[i] = wa[i]/wb[i];
    fprintf(stderr,"\t-> FST.Unweight[nObs:%lu]:%f Fst.Weight:%f\n",nObs[i],fstUW[i],fstW[i]);
    fprintf(stdout,"%f\t%f\n",fstUW[i],fstW[i]);
  }
  if(chs==3){
    //if chr==3 then we have 3pops and we will also calculate pbs statistics
    calcpbs(fstW);//<- NOTE: the pbs values will replace the fstW values
    for(int i=0;i<3;i++)
      fprintf(stdout,"pbs.pop%d\t%f\n",i+1,fstW[i]);
  }
  delete [] ares;
  delete [] bres;
  destroy_args(pars);
  perfst_destroy(pf);
  return 0;
}

