/*
  jackknifing could be speeded current implemention is O(nsites^2)
*/

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <cassert>
#include <ctype.h>
#include <pthread.h>
#include "kmath.h"

#include "../fet.c"
#define LENS 4096

typedef struct{
  char allele1;
  char allele2;
  double freq;
}hapSite;

void print(FILE *fp,int p,hapSite &h){
  fprintf(fp,"hapmap\t%d\t%d\t%d\t%f\n",p,h.allele1,h.allele2,h.freq);
}

double sd(double *a,int l){
  double ts =0;
  for(int i=0;i<l;i++)
    ts += a[i];

  double u = ts/(1.0*l);
  assert(u!=0);
  ts =0;
  for(int i=0;i<l;i++)
    ts += (a[i]-u)*(a[i]-u);

  double N = l;
  double scal=(N-1.0)/N;
  return sqrt(scal*ts);
}

double ldbinom(int k, int n,double p){
  return lbinom(n,k)+k*log(p)+(n-k)*log(1-p);
}

//    l<-  dbinom(error,d,(1-x)*eps+x*freq)
double likeOld(double x,int len,int *seqDepth,int *nonMajor,double *freq,double eps,int skip){
  double t = 0;
  for(int i=0;i<len;i++)
    if(i!=skip)
      t += ldbinom(nonMajor[i],seqDepth[i],(1-x)*eps+x*freq[i]);
  return t;
}

//    l<-  dbinom(error,d,x*freq*(1-4*eps/3)+eps)
double likeNew(double x,int len,int *seqDepth,int *nonMajor,double *freq,double eps,int skip){
  double t = 0;
  for(int i=0;i<len;i++)
    if(i!=skip)
      t += ldbinom(nonMajor[i],seqDepth[i],x*freq[i]*(1-4*eps/3.0)+eps);
  return t;
  
}


double likeOldMom(int len,int *seqDepth,int *nonMajor,double *freq,double eps,int jump){
  //  return(mean(error/d-eps)/(mean(freq)-eps))
  double top = 0;
  double bot = 0;
  
  for(int i=0;i<len;i++){
    if(i!=jump){
      top += (1.0*nonMajor[i])/(1.0*seqDepth[i]);
      bot += freq[i];
    }
  }
  top = top/(1.0*len)-eps;
  bot = bot/(1.0*len)-eps;
  return top/bot; 
}


double likeNewMom(int len,int *seqDepth,int *nonMajor,double *freq,double eps,int jump){
  //  return(mean(error/d-eps)/(mean(freq)-eps))
  double top = 0;
  double bot = 0;
  
  for(int i=0;i<len;i++){
    if(i!=jump){
      top += (1.0*nonMajor[i])/(1.0*seqDepth[i]);
      bot += freq[i];
    }
  }
  top = top/(1.0*len)-eps;
  bot = bot/(1.0*len)*(1-4.0*eps/3.0);
  return top/bot; 
}


typedef struct{
  int newllh;//<- bool, should we use new likelihood function.
  double eps;//<- error rate
  int len;
  //  int skip;
  int* seqDepth;
  int *nonMajor;
  double *freq;
  int *e1;
  int *d1;
}allPars;


typedef struct{
  allPars *ap;
  double *thetas;
  double *val;
  int from;
  int to;
  int skip;
}tpars;

void print(allPars *ap,char *fname){
  if(fname){
    FILE *fp = fopen(fname,"w");
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"%d %d %d %d %f\n",ap->nonMajor[i],ap->e1[i],ap->seqDepth[i],ap->d1[i],ap->freq[i]);
    fclose(fp);
  }
}
void marshall(allPars *ap,char *fname){
  if(fname){
    FILE *fp = fopen(fname,"w");
    fprintf(fp,"%d\n",ap->newllh);
    fprintf(fp,"%f\n",ap->eps);
    fprintf(fp,"%d\n",ap->len);
    // fprintf(fp,"%d\n",ap->skip);
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"sd[%d]\t%d\n",i,ap->seqDepth[i]);
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"sd[%d]\t%d\n",i,ap->nonMajor[i]);
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"sd[%d]\t%f\n",i,ap->freq[i]);
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"sd[%d]\t%d\n",i,ap->d1[i]);
    for(int i=0;i<ap->len;i++)
      fprintf(fp,"sd[%d]\t%d\n",i,ap->e1[i]);
    fclose(fp);
  }

}


//calculate error rate
double calcEps(int *e1,int *d1,int len,int skip){
  //  fprintf(stderr,"eps len:%d\n",len);
  //fprintf(stderr,"eps skip:%d\n",skip);
  double top=0;double bot=0;
  for(int i=0;i<len;i++){
    if(i!=skip){
      top += e1[i];
      bot += d1[i];
    }
  }
  double eps=top/bot;
  //  fprintf(stderr,"eps:%f\n",eps);
  return eps;
  
}



//wrapper function used for threading the ML opt
double myfun(double x,void *d){
  tpars *tp = (tpars *)d;
  //  fprintf(stderr,"tp->ap->newllh:%d\n",tp->ap->newllh);
  if(tp->ap->newllh)
    return -likeNew(x,tp->ap->len,tp->ap->seqDepth,tp->ap->nonMajor,tp->ap->freq,calcEps(tp->ap->e1,tp->ap->d1,tp->ap->len,tp->skip),tp->skip);
  else
    return -likeOld(x,tp->ap->len,tp->ap->seqDepth,tp->ap->nonMajor,tp->ap->freq,calcEps(tp->ap->e1,tp->ap->d1,tp->ap->len,tp->skip),tp->skip);
}

double jackMom(allPars *ap,int len){
  len=ap->len;
  if(len>0)
    len=std::min(ap->len,len);
  assert(len>-1);
  int *seqDepth=ap->seqDepth;
  int *nonMajor=ap->nonMajor;
  double *freq =ap->freq;
  int isnew = ap->newllh;
  int *e1 = ap->e1;
  int *d1 = ap->d1;
  double *thetas =new double[len];

  for(int i=0;i<len;i++){
    if(isnew==0)
      thetas[i] = likeNewMom(len,seqDepth,nonMajor,freq,calcEps(e1,d1,len,i),i) ;
    else
      thetas[i] = likeOldMom(len,seqDepth,nonMajor,freq,calcEps(e1,d1,len,i),i) ;
    
  }
  double esd = sd(thetas,len);
  delete [] thetas;
  return esd;
}



//each thread will run this.
void *slave(void *ptr){
  tpars *tp = (tpars *)ptr; 
  //  fprintf(stderr,"from:%d to:%d \n",tp->from,tp->to);
  for(int i=tp->from;i<tp->to;i++){
    tp->skip=i;
    tp->val[i] = kmin_brent(myfun,1e-6,0.5,tp,1e-6,&tp->thetas[i]);
    assert(tp->thetas[i]!=1e-6);
  }
  pthread_exit(0);
}



double jackML(allPars *ap,int nthreads,char *fname,int nJack) {
  
  if(nJack==-1)
    nJack = ap->len;
  else
    nJack = std::min(ap->len,nJack);
  assert(nJack>0);
  double *thetas =new double[nJack];
  double *val = new double[nJack];
  if(nthreads>1){
    pthread_t *thd = new pthread_t[nthreads];
    tpars *tp = new tpars[nthreads];
    int block = nJack/nthreads;
    
    for(int i=0;i<nthreads;i++){
      tp[i].thetas = thetas;
      tp[i].val = val;
      tp[i].ap = ap;
      tp[i].from = i==0?0:tp[i-1].to;
      tp[i].to = tp[i].from+block;
    }
    tp[nthreads-1].to = nJack;
    for(int i=0;i<nthreads;i++)
      pthread_create(&thd[i],NULL,slave,&tp[i]);

    for(int i=0;i<nthreads;i++)
      pthread_join(thd[i],NULL);
    
  }else{    //if we do not threads  
    tpars tp;
    tp.ap=ap;
    for(int i=0;i<nJack;i++){
      tp.skip=i;
      val[i]=kmin_brent(myfun,1e-6,0.5-1e-6,&ap,0.0001,thetas+i);
    }
  }
  double esd = sd(thetas,nJack);
  if(fname){
    FILE *fp =fopen(fname,"w");
    for(int i=0;i<nJack;i++){
      fprintf(fp,"%e\t%f\t%e\n",thetas[i],val[i],thetas[i]-1e6);
    }
    fclose(fp);
  }
  delete [] thetas;
  delete [] val;
  return esd;
}


typedef std::map<int,hapSite> aMap;

const char *hapfile=NULL,*mapfile=NULL,*icounts=NULL;
typedef unsigned char uchar;

typedef std::vector<int> iv;
typedef std::vector<int*> iv2;

typedef struct{
  iv pos;//contains position
  iv dist;//0=snpsite;-1=one pos left of snp,+1=one pos right of snp
  std::vector<int *> cn;//four long , counts of C,C,G,T
  aMap myMap;//allele1, allele2, freq
}dat;



int simrbinom(double p){
  if(drand48()<(1-p))
    return 0;
  else
    return 1;
}

void analysis(dat &d,int nThreads,int nJack,int skipML) {
  int *rowSum = new int[d.cn.size()];
  int *rowMax = new int[d.cn.size()];
  int *rowMaxW = new int[d.cn.size()];
  int *error1 = new int[d.cn.size()];//number of non most frequent observed bases 
  int *error2 = new int[d.cn.size()];//sampled 
  size_t mat1[4]={0,0,0,0};//matrix used for fisher for method 1
  size_t mat2[4]={0,0,0,0};//matrix used for fisher for method 2
  size_t tab[2] = {0,0};//used for debug

  for(int i=0;i<d.cn.size();i++) {
    int s =d.cn[i][0];
    int max=s;
    int which=0;
    for(int j=1;j<4;j++){
      s += d.cn[i][j];
      if(d.cn[i][j]>max){
	max=d.cn[i][j];
	which=j;
      }
    }
    rowSum[i] = s;
    rowMax[i]=max;
    rowMaxW[i]=which;
    aMap::iterator it= d.myMap.find(d.pos[i]);
    if(it!=d.myMap.end()){//if site is hapmap site
      // fprintf(stderr,"posi:%d wmax:%d all1:%d freq:%f\n",it->first,rowMaxW[i],it->second.allele1,it->second.freq);
      //if maximum occuring bases is the same as allele1 from hapmap, then set freq to 1-freq
      if(rowMaxW[i]==it->second.allele1)
	//it->first C++ syntax for getting key of iterator
	//it->second C++ syntax for getting value of key of iterator, key->value: key=pos,value=hapSite
 	it->second.freq=1-it->second.freq;
      else
	it->second.freq=it->second.freq;
      // fprintf(stderr,"posi:%d wmax:%d all1:%d freq:%f\n",it->first,rowMaxW[i],it->second.allele1,it->second.freq);
      // exit(0);
    }

    error1[i] = rowSum[i]-rowMax[i];
    error2[i] = simrbinom((1.0*error1[i])/(1.0*rowSum[i]));
    //  fprintf(stdout,"simrbinom\t%d %d %d %d %d %d\n",rowSum[i],rowMax[i],rowMaxW[i],error1[i],error2[i],d.dist[i]);
    if(error1[i]>0)
      tab[1]++;
    else
      tab[0]++;
    if(d.dist[i]==0){//this is a snpsite
      mat1[0] +=error1[i];
      mat1[1] +=rowSum[i]-error1[i];
      mat2[0] +=error2[i];
      mat2[1] +=1-error2[i];
      // fprintf(stdout,"rs %d %d %d %d %d %d %d %f %d %d\n",d.pos[i],rowSum[i],rowMax[i],rowMaxW[i],error1[i],error2[i],d.dist[i],it->second.freq,it->second.allele1,it->second.allele2);
    }else{
      
      mat1[2] +=error1[i];
      mat1[3] +=rowSum[i]-error1[i];
      mat2[2] += error2[i];
      mat2[3] += 1-error2[i];
    }
  }
#if 0
  fprintf(stderr,"tab:%lu %lu\n",tab[0],tab[1]);
  fprintf(stderr,"mat: %lu %lu %lu %lu\n",mat1[0],mat1[1],mat1[2],mat1[3]);
  fprintf(stderr,"mat2: %lu %lu %lu %lu\n",mat2[0],mat2[1],mat2[2],mat2[3]);
#endif
  int n11, n12, n21, n22;
  double left, right, twotail, prob;
  //  fprintf(stderr,"--------\nMAIN RESULTS: Fisher exact test:\n");
  n11=mat1[0];n12=mat1[2];n21=mat1[1];n22=mat1[3];
  prob = kt_fisher_exact(n11, n12, n21, n22, &left, &right, &twotail);
  //  fprintf(stdout,"Method\t n11 n12 n21 n22 prob left right twotail\n");
  //fprintf(stdout,"%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\n", "method1", n11, n12, n21, n22,
  //		prob, left, right, twotail);

  n11=mat2[0];n12=mat2[2];n21=mat2[1];n22=mat2[3];
  prob = kt_fisher_exact(n11, n12, n21, n22, &left, &right, &twotail);
//fprintf(stdout,"%s\t%d\t%d\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\n", "method2", n11, n12, n21, n22,
  //			prob, left, right, twotail);

  //estimate how much contamination
  double c= mat1[2]/(1.0*(mat1[2]+mat1[3]));//this is error for flanking site
  double err= mat1[0]/(1.0*(mat1[0]+mat1[1]));//this is error for snpsite
  fprintf(stderr,"Mismatch_rate_for_flanking:%f MisMatch_rate_for_snpsite:%f \n",c,err);

  int *err0 =new int[d.cn.size()/9];//<-nbases of non frequent most occuring at snpsite
  int *err1 =new int[d.cn.size()/9];//<-nbases of non frequent most occuring at flanking
  int *d0 =new int[d.cn.size()/9];//<-seqdepth for snpsite
  int *d1 =new int[d.cn.size()/9];//<-seqdepth for flanking
  double *freq =new double[d.cn.size()/9];//<- freq for snpsite

  for(int i=0;i<d.cn.size()/9;i++){
    int adj=0;
    int dep=0;
    for(int j=0;j<9;j++){
      
      if(d.dist[i*9+j]!=0){//<- flanking
	adj += error1[i*9+j];
	dep += rowSum[i*9+j];
      }else{ //snpsite
	err0[i] = error1[i*9+j];
	d0[i] = rowSum[i*9+j];
	
	freq[i] =d.myMap.find(d.pos[i*9+j])->second.freq;
	
      }
    }
    err1[i] =adj;
    d1[i] = dep;
#if 0
    if(it==d.myMap.end()){
      fprintf(stderr,"Problem finding:%d\n",d.pos[i]);
      exit(0);
    }
#endif
    //    fprintf(stdout,"cont\t%d\t%d\t%d\t%d\t%f\n",err0[i],err1[i],d0[i],d1[i],freq[i]);
    
  }


  allPars ap;
  ap.len=d.cn.size()/9;
  ap.seqDepth = d0;
  ap.nonMajor = err0;
  ap.freq = freq;
  ap.eps = c;
  ap.newllh =0;
  ap.e1 = err1;
  ap.d1=d1;

  double mom,momJack,ML,mlJack,val;

  ap.newllh =0;
  mom= likeOldMom(d.cn.size()/9,d0,err0,freq,c,-1);
  momJack = jackMom(&ap,nJack);
  tpars tp;tp.ap=&ap;tp.skip=-1;
  //  print(tp.ap,"asdff1");
  if(skipML==0){
    kmin_brent(myfun,1e-6,0.5-1e-6,&tp,0.0001,&ML);
    mlJack= jackML(&ap,nThreads,NULL,nJack);
  }else{
    ML=-1;
    mlJack=-1;
  }
  fprintf(stderr,"\nMethod1: old_llh Version: MoM:%f SE(MoM):%e ML:%f SE(ML):%e",mom,momJack,ML,mlJack);
  
  ap.newllh =1;
  mom=likeNewMom(d.cn.size()/9,d0,err0,freq,c,-1);
  momJack= jackMom(&ap,nJack);
  //marshall(&ap,"prem1");
  if(skipML==0){
    val=kmin_brent(myfun,1e-6,0.5-1e-6,&tp,0.0001,&ML);
    //  fprintf(stderr,"\nM1: ML:%f VAL:%f\n",ML,val);
    mlJack= jackML(&ap,nThreads,NULL,nJack);
  }else{
    val=-1;
    mlJack=-1;
  }
  fprintf(stderr,"\nMethod1: new_llh Version: MoM:%f SE(MoM):%e ML:%f SE(ML):%e",mom,momJack,ML,mlJack);
  // fread(error2,sizeof(int),d.cn.size(),fopen("error2.bin","rb"));
  //for(int i=0;0&&i<d.cn.size();i++)
  //  fprintf(stdout,"pik\t%d\n",error2[i]);
  //exit(0);
  for(int i=0;i<d.cn.size()/9;i++){
    int adj=0;
    for(int j=0;j<9;j++){
      if(d.dist[i*9+j]!=0){
	adj += error2[i*9+j];
      }else{
	err0[i] = error2[i*9+j];
	freq[i] =d.myMap.find(d.pos[i*9+j])->second.freq;
	//	fprintf(stderr,"freq:%f\n",freq[i]);
      }
    }
    err1[i] =adj;
    d0[i] = 1; 
    d1[i] = 8;
#if 0
    if(it==d.myMap.end()){
      fprintf(stderr,"Problem finding:%d\n",d.pos[i]);
      exit(0);
    }
#endif
    //    fprintf(stdout,"cont\t%d\t%d\t%d\t%d\t%f\n",err0[i],err1[i],d0[i],d1[i],freq[i]);
    
  }

  ap.seqDepth = d0;
  ap.nonMajor = err0;
  ap.e1=err1;
  ap.d1=d1;
  ap.freq = freq;

  ap.newllh =0;
  mom= likeOldMom(d.cn.size()/9,d0,err0,freq,c,-1);
  momJack = jackMom(&ap,nJack);
  //print(tp.ap,"asdff2");
  //  exit(0);
  val = kmin_brent(myfun,1e-6,0.5-1e-6,&tp,0.0001,&ML);
  //fprintf(stderr,"\nML2:%f VAL:%f\n",ML,val);
  // exit(0);
  //FILE *fp = fopen("heyaa","w");  print(&ap,fp);fclose(fp);
  //return;
  mlJack= jackML(&ap,nThreads,NULL,nJack);
  fprintf(stderr,"\nMethod2: old_llh Version: MoM:%f SE(MoM):%e ML:%f SE(ML):%e",mom,momJack,ML,mlJack);

  ap.newllh =1;
  mom=likeNewMom(d.cn.size()/9,d0,err0,freq,c,-1);
  momJack= jackMom(&ap,nJack);
  kmin_brent(myfun,1e-6,0.5-1e-6,&tp,0.0001,&ML);
  mlJack= jackML(&ap,nThreads,NULL,nJack);
  fprintf(stderr,"\nMethod2: new_llh Version: MoM:%f SE(MoM):%e ML:%f SE(ML):%e\n",mom,momJack,ML,mlJack);
 

  delete [] rowSum;
  delete [] rowMax;
  delete [] rowMaxW;
  delete [] error1;
  delete [] error2;
  

  delete [] err0;
  delete [] err1;
  delete [] d0;
  delete [] d1;
  delete [] freq;



}

void print(int *ary,FILE *fp,size_t l,char *pre){
  fprintf(fp,"%s\t",pre);
  for(size_t i=0;i<l;i++)
    fprintf(fp,"%d\t",ary[i]);
  fprintf(fp,"\n");
}

#define NVAL -66
dat count(aMap &myMap,std::vector<int> &ipos,std::vector<int*> &cnt){
  int lastP = std::max((--myMap.end())->first,ipos[ipos.size()-1])+5;//<-add five so we dont step out

  char *hit = new char[lastP];
  memset(hit,NVAL,lastP);//-10 just indicate no value..
  
  for(aMap::iterator it=myMap.begin();it!=myMap.end();++it){
      for(int p=0;p<5;p++){
        hit[it->first+p] =hit[it->first-p] = p;
	hit[it->first-p] = -  hit[it->first-p] ;
      }
  }

  //now hit contains long stretches of NVAL and -4,...4.
  //now set NVAL to all the places where we don't have data
  char *aa= new char[lastP];
  memset(aa,NVAL,lastP);
  for(int i=0;i<ipos.size();i++)
    aa[ipos[i]] = 1; 
  //now loop over hit array and set to NVAL if no data
  for(int i=0;i<lastP;i++)
    if(aa[i]==NVAL)
      hit[i] = NVAL;
  delete [] aa;
  //now loop over hitarray and make remove non -4,..4: segments
  int i=0;
  while(i<lastP){
    if(hit[i]!=-4)
      hit[i] = NVAL;
    else{
      int isOk=1;
      for(int j=0;j<9;j++)
	if(hit[i+j]!=j-4){
	  isOk=0;
	  break;
	}
      if(isOk==1)
	i+=9;
      else{
	hit[i]=NVAL;
      }
    }
    i++;
  }

  dat d;
  for(int i=0;i<ipos.size();i++)
    if(hit[ipos[i]]!=NVAL) { 
      d.pos.push_back(ipos[i]);
      d.dist.push_back(hit[ipos[i]]);
      
      d.cn.push_back(cnt[i]);//plugs in pointer to 4ints.
      
      if(hit[ipos[i]]==0){//this is a snpsite
	aMap::iterator it=myMap.find(ipos[i]);
	if(it==myMap.end()){
	  fprintf(stderr,"Problem finding snpsite:%d\n",ipos[i]);
	  exit(0);
	}
	aMap::iterator it2=d.myMap.find(ipos[i]);
	assert(it2==d.myMap.end());
	d.myMap[it->first] = it->second;

      }
    }
  delete [] hit;
  fprintf(stderr,"  After removing SNP sites with no data in 5bp surrounding region`\n  We have nSNP sites: %lu, with flanking: %lu\n",d.myMap.size(),d.cn.size());
  return d;
}


char flip(char c){
  c = toupper(c);
  if(c=='A')
    return 'T';
  if(c=='T')
    return 'A';
  if(c=='G')
    return 'C';
  if(c=='C')
    return 'G';
  if(c==0)
    return 3;
  if(c==1)
    return 2;
  if(c==2)
    return 1;
  if(c==3)
    return 0;
  fprintf(stderr,"Problem interpreting char:%c\n",c);
  return 0;
}

gzFile getgz(const char *fname,const char *mode){
  //  fprintf(stderr,"Trying to open file:%s\n",fname);
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }
  return gz;
}

void printhapsite(hapSite &hs,FILE *fp,int &p){
  fprintf(fp,"p:%d al1:%c al2:%c freq:%f\n",p,hs.allele1,hs.allele2,hs.freq);
}

int charToNum(char c){
  if(c=='A'||c=='a')
    return 0;
  if(c=='C'||c=='c')
    return 1;
  if(c=='G'||c=='g')
    return 2;
  if(c=='T'||c=='t')
    return 3;
  // fprintf(stderr,"Problem with observed char: '%c'\n",c);
  return -1;
}


aMap readhap( char *fname,int minDist,double minMaf,int startPos,int stopPos,int skiptrans){
  //  fprintf(stderr,"[%s] fname:%s\tminDist:%d LENS:%d\n",__FUNCTION__,fname,minDist,LENS);
  gzFile gz=getgz(fname,"rb");
  
  char *buf = new char[LENS];
  int viggo=3;
  aMap myMap;//<- this will be our return object

  gzgets(gz,buf,LENS);
  while(gzgets(gz,buf,LENS)){//loop over sites
    hapSite hs;
    int p = atoi(strtok(buf,"\t\n "));// pos
    if(p<startPos)
      continue;
    if(p>stopPos)
      continue;

    hs.allele1 = charToNum(strtok(NULL,"\t\n ")[0]);
    hs.freq = atof(strtok(NULL,"\t\n "));
    if(hs.freq<minMaf||(1-hs.freq)<minMaf)
      continue;
    char strand= strtok(NULL,"\t\n ")[0];
    hs.allele2 = charToNum(strtok(NULL,"\t\n ")[0]);
    
    if(hs.allele1==-1||hs.allele2==-1)
      continue;
    if(strand=='-'){
      hs.allele1 = flip(hs.allele1);
      hs.allele2 = flip(hs.allele2);
    }
    if(hs.allele1==hs.allele2)
      continue;

    if(skiptrans){
      if(hs.allele1 == 0 && hs.allele2 == 2)
	continue;
      if(hs.allele1 == 2 && hs.allele2 == 0)
	continue;
      if(hs.allele1 == 1 && hs.allele2 == 3)
	continue;
      if(hs.allele1 == 3 && hs.allele2 == 1)
	continue;

    }
    if(myMap.count(p)>0){
      if(viggo>0){
	//	fprintf(stderr,"[%s] Duplicate positions found in file: %s, pos:%d\n",__FUNCTION__,fname,p);
	//fprintf(stderr,"[%s] Will only use first entry\n",__FUNCTION__);
	//fprintf(stderr,"[%s] This message is only printed 3 times\n",__FUNCTION__);
	viggo--;
      }
    }else{
      myMap[p]=hs;
    }
  }
  //  fprintf(stderr,"[%s] We have read: %zu sites from hapfile (after filtering for start/stop pos):%s\n",__FUNCTION__,myMap.size(),fname);
  //fprintf(stderr,"[%s] will remove snp sites to close:\n",__FUNCTION__);

  assert(myMap.size()>0);
  int *vec = new int[myMap.size() -1];
  aMap::iterator it = myMap.begin();
  for(int i=0;i<myMap.size()-1;i++){
    aMap::iterator it2 = it;
    it2++;
    vec[i]=it2->first - it->first;
    it=it2;
  }
  it = myMap.begin();
  aMap newMap;
  for(int i=0;i<myMap.size()-1;i++){
    if(std::abs(vec[i])>=minDist){
      newMap[it->first] = it->second;
      //     fprintf(stdout,"test\t%d\n",it->first);
    }
    it++;
  }
  newMap[it->first] = it->second;
  //..  exit(0);
  delete [] vec;
  fprintf(stderr,"[%s] We now have: %lu snpSites after filtering based on hapMapfile\n",__FUNCTION__,newMap.size());

#if 0
  for(aMap::iterator it=newMap.begin();it!=newMap.end();++it)
    print(stdout,it->first,it->second);
  exit(0);
#endif
  delete [] buf;
  gzclose(gz);
  return newMap;
}




void readicnts(const char *fname,std::vector<int> &ipos,std::vector<int*> &cnt,int minDepth,int maxDepth){
    fprintf(stderr,"[%s] fname:%s minDepth:%d maxDepth:%d\n",__FUNCTION__,fname,minDepth,maxDepth);
  gzFile gz=getgz(fname,"rb");

  int tmp[5];
  int totSite=0;
  while(gzread(gz,tmp,sizeof(int)*5)){
    totSite++;
    int *tmp1=new int[4];
    int d=0;
    for(int i=0;i<4;i++) {
      tmp1[i]=tmp[i+1];
      d += tmp1[i];
    }
   
    if(d>=minDepth&&d<=maxDepth){//if we should use site?
      cnt.push_back(tmp1); //push back [#A,#C,#G,#T]
      ipos.push_back(tmp[0]-1);//push back position
    }else
      delete [] tmp1;
  }
  //  print(cnt[0],stderr,4,"dung");
  fprintf(stderr,"[%s] Has read:%d sites,  %zu sites (after depfilter) from ANGSD icnts file\n",__FUNCTION__,totSite,ipos.size());
  gzclose(gz);
}

int main(int argc,char**argv){
  char *hapfile=NULL;
  char *icounts=NULL;
  double minMaf = 0.05;
  int startPos = 5e6;
  int stopPos = 154900000;
  int minDepth=2;
  int maxDepth=200;
  int skipTrans = 0;
  int nThreads = 1;
  int nJack=-1;
  long int seed = 0;
  int n;
  int skipML = 0;
  int printcounts =0;
  while ((n = getopt(argc, argv, "h:a:m:b:c:d:e:f:p:s:j:l:k:")) >= 0) {
    switch (n) {
    case 'h': hapfile = strdup(optarg); break;
    case 'a': icounts = strdup(optarg); break;
    case 'm': minMaf = atof(optarg); break;
    case 'b': startPos = atoi(optarg); break;
    case 'c': stopPos = atoi(optarg); break;
    case 'd': minDepth = atoi(optarg); break;
    case 'e': maxDepth = atoi(optarg); break;
    case 'f': skipTrans = atoi(optarg); break;
    case 'j': nJack = atoi(optarg); break;
    case 'k': printcounts = atoi(optarg); break;
    case 'p': nThreads = atoi(optarg); break;
    case 's': seed = atol(optarg); break;
    case 'l': skipML = atoi(optarg); break;
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
    }
  }
  if(!hapfile||!icounts){
    fprintf(stderr,"\t-> Must supply -h hapmapfile -a angsd.icnts.gz file\n");
    fprintf(stderr,"\t-> Other options: -m minaf -b startpos -c stoppos -d mindepth -e maxdepth -f skiptrans -p nthreads -s seed -j maxjackknife -k printcounts\n");
    return 0;
  }
  
  int minDist = 10;
  fprintf(stderr,"-----------------\nhapmap:%s counts:%s minMaf:%f startPos:%d stopPos:%d minDepth:%d maxDepth:%d skiptrans:%d nthreads:%d seed:%lu skipML:%d\n",hapfile,icounts,minMaf,startPos,stopPos,minDepth,maxDepth,skipTrans,nThreads,seed,skipML);
  fprintf(stderr,"Method2 is subject to fluctuations due to random sampling\n");
  fprintf(stderr,"Seed value of 0 (zero) will use time as seed\n-----------------\n");
  if(seed==0)
    srand48(time(NULL));
  else
    srand48(seed);
  
  aMap myMap = readhap(hapfile,minDist,minMaf,startPos,stopPos,skipTrans);
  std::vector<int> ipos;
  std::vector<int*> cnt;
  readicnts(icounts,ipos,cnt,minDepth,maxDepth);
  dat d=count(myMap,ipos,cnt);
  if(printcounts){
    aMap dm=d.myMap;
    for(aMap::iterator it=dm.begin();it!=dm.end();it++)
      fprintf(stdout,"hapmap\t%d\t%d\t%d\t%f\n",it->first,it->second.allele1,it->second.allele2,it->second.freq);
    for(int i=0;i<d.pos.size();i++)
      fprintf(stdout,"counts\t%d\t%d\t%d\t%d\t%d\t%d\n",d.pos[i],d.dist[i],d.cn[i][0],d.cn[i][1],d.cn[i][2],d.cn[i][3]);
    return 0;
  }
  analysis(d,nThreads,nJack,skipML);

  //cleanup
  for(int i=0;i<cnt.size();i++)
    delete [] cnt[i];
  for(int i=0;0&&i<d.cn.size();i++)//<-should cleanup here. We have copied pointers not values
    delete [] d.cn[i];
  free(hapfile);
  free(icounts);
  return 0;

}
