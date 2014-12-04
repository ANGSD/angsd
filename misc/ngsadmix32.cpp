/*
  log:
  ngsAdmix32: program now works on osx
  icpc ngsAdmix32.cpp -lz -lpthread  -O3 -o ngsAdmix32

  ngsAdmix31: reads gz frequencies 
  icpc ngsAdmix31.cpp -lz -lpthread  -O3 -o ngsAdmix31
  ngsAdmix30: prints information of the filtered sites
  icpc ngsAdmix30.cpp -lz -lpthread  -O3 -o ngsAdmix30
 ngsAdmix29: more speed. DO_MIS fucked
  icpc ngsAdmix29.cpp -lz -lpthread  -O3 -o ngsAdmix29
  ngsAdmix26: more speed. No longer estimate likelihoods 
   icpc ngsAdmix26.cpp -lz -lpthread  -O3 -o ngsAdmix26domis  -DDO_MIS
  ngsAdmix25: small speed up (expG1 = expG/fstar, expG2 = (2-expG) / ( 1 - fstar))
  icpc ngsAdmix25.cpp -lz -lpthread  -O3 -o ngsAdmix25
  icpc ngsAdmix25.cpp -lz -lpthread  -O3 -o ngsAdmix25domis -DDO_MIS
  ngsAdmix24: clean up. no more genotype calling. CHECK must now be compiled. DO_MIS fixed, nMiss ->minInd
  icpc ngsAdmix24.cpp -lz -lpthread  -O3 -o ngsAdmix24
  ngsAdmix23: prints gz Q matrix. dynamic bounderies
  icpc ngsAdmix23.cpp -lz -lpthread  -O3 -o ngsAdmix23
  ngsAdmix22: prints gz Q matrix
  icpc ngsAdmix22and2.cpp -lz -lpthread  -O3 -o ngsAdmix22and2
  ngsAdmix21: If missing info -maxLike sets GLikes to 0.33
  same for -hweLike but probably never used
  icpc ngsAdmix21.cpp -lz -lpthread  -O3 -o ngsAdmix21

  ngsAdmix20: minLRT was compared to likelihood  - now fixed 
  icpc ngsAdmix20.cpp -lz -lpthread  -O3 -o ngsAdmix20
  
  version 19
  cutoff for HWE fixed - did the same as HWE without cutoff before.
  icpc ngsAdmix19.cpp -lz -lpthread  -O3 -o ngsAdmix19
  still a problem with -DDO_MIS giving -nan likelihoods - not fixed
  
  version 18 line er soed
  cutoff for genotype calling -cut (only works when compiled with -DDO_MIS)
  call genotypes (maxLike, hweLike) and update freq estimate before filtering - because this is what people will do
  icpc ngsAdmix18.cpp -lz -lpthread  -O3 -o ngsAdmix18
  icpc ngsAdmix18.cpp -lz -lpthread -DDO_MIS -O3 -o ngsAdmix18mis
  
  //version 17
  add option to change the difference tolerance for calling missing
  
  //version 16: fox chix
  simple missing implemented
  1) skip missing genotypes
  icpc ngsAdmix16.cpp -lz -lpthread -DDO_MIS -O3
  2) use missing genotypes
  icpc ngsAdmix16.cpp -lz -lpthread  -O3
  //not major allert. if hwe genotyep calling. alleles are swiched. You cannot trust frequency output

  //version 15: Thorfinn fox chix
  added minmaf filter, added hwe genotype calling added lrt filter
  added char** keeps indicating if a site/sample has data
  added nMis which will remove sites with miss>nMissallowed
  
  //version 14: anders sux cux
  after the alpha guess the values are mapped to the domain
  program for estimating admixture from beagle files.
  beagle files contains posterior in normal space
  #introducing squareem
  icpc -O3 -o ngsAdmix ngsAdmix14.cpp -lz -lpthread

  ngsAdmix13: Use random seed. can set seed in options
*/

//optimazation parameteres for maf estimation
#define MAF_START 0.3
#define MAF_ITER 20


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <limits>
#include <zlib.h>
#include <vector>
#include <pthread.h>
#include <signal.h>
#include <vector>
#include <sys/stat.h>

//This is taken from here:
//http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
#ifdef __APPLE__

#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count)
{
    if(count == 0)
    {
        errno = EINVAL;
        return -1;
    }

if(pthread_mutex_init(&barrier->mutex, 0) < 0)
    {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }

    barrier->tripCount = count;
    barrier->count = 0;

    return 0;
}
int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }
    else
    {
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

#endif // PTHREAD_BARRIER_H_
#endif // __APPLE__



//global stuff below, this is very beautifull
char **keeps=NULL;
#define LENS 100000 //this is the max number of bytes perline, should make bigger
int SIG_COND =1;//if we catch signal then quit program nicely
double tol=1e-5; //stopping criteria


double errTolMin=1e-9;
double errTolStart=0.1;
double errTol=errTolStart;//frequencies and admixture coef cannot be less than this or more than 1-this
double misTol=0.05;

//this struct contains nescerray information needed for threading
typedef struct{
  int ID; 
  int threadNumber; 
  int nThreads; 
  double **F_1;
  double **Q_1;
  double **F;
  double **Q;
  int start; 
  int stop;
  int startI; 
  int stopI;
  double ***prodA;
  double ***prodB;
  double nPop;
  double **genos;
  int nInd;
  int nSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
}pars;

//to make life simple we make stuff relating to the threading global
pthread_t *threads = NULL;
pars * myPars= NULL; //den ovenfor definerede struct type

 
double likeFixedMinor(double p,double *likes,int numInds,char *keepInd){
  // should these actually be normed genotype likelihoods? Or not normed?
  // only used for checking minLrt
  // returns minus the log of the likelihood
  double totalLike=0;
  for(int i=0;i<numInds;i++){
    if(keepInd[i])
      totalLike+=log(likes[i*3+0]*(1-p)*(1-p)+likes[i*3+1]*2.0*p*(1-p)+likes[i*3+2]*p*p);
  }
  return -totalLike;
}



double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;


  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }



  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);

    for(int ii=0;0&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	    fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    // exit(0);
  }
  
  return(p);
}



void getExpGandFstar(double** Q, double** F, int nSites_start,int nSites_stop, int nInd, int K,double **genos,int startI,int stopI,char **keeps,double ***prodA,double ***prodB,double **F_1,int thread){

    for(int j=nSites_start;j<nSites_stop;j++){
      double sumAG[K];
      double sumBG[K];
      for(int k=0;k<K;k++){ //time killar
	sumAG[k]=0;
	sumBG[k]=0;
      }
	
      for(int i=0;i<nInd;i++){
#ifdef DO_MIS
	//	if(keeps[j][i]==0){
	//	  expG1[j][i]= std::numeric_limits<double>::quiet_NaN();
	//	  expG2[j][i]= std::numeric_limits<double>::quiet_NaN();
	//	}
	//	else{
#endif
	  double fpart=0;
	  for(int k=0;k<K;k++){ //time killar
	    fpart += F[j][k] * Q[i][k];
	  }
	  fpart=1-fpart;
	  double pp0=fpart*fpart*        genos[j][3*i+0];
	  double pp1=2*(1-fpart)*fpart*  genos[j][3*i+1];
	  double pp2=(1-fpart)*(1-fpart)*genos[j][3*i+2];
	  double sum=pp0+pp1+pp2;
	  double expGG =(pp1+2*pp2)/sum;//range 0-2, this is the expected genotype
	
	 
	  for(int k=0;k<K;k++){
	    prodA[thread][i][k] += expGG/(1-fpart) * F[j][k]; //proteckMe
	    prodB[thread][i][k] += (2-expGG)/fpart *(1- F[j][k]); //proteckMe
	    sumAG[k] += expGG/(1-fpart) * Q[i][k]; //time killar
	    sumBG[k] += (2-expGG)/fpart * Q[i][k];//time killar
  
	  }
#ifdef DO_MIS
	}
#endif
      }
      for(int k=0;k<K;k++){
	sumAG[k] *= F[j][k];
	sumBG[k] *= (1-F[j][k]);
	F_1[j][k]=sumAG[k]/(sumAG[k]+sumBG[k]);
      }

    
    }
  
}



void map2domainF(double** F, int nSites_start,int nSites_stop, int K){
 for(int s=nSites_start;s<nSites_stop;s++)
   for(int k=0;k<K;k++){
     if(F[s][k]<errTol)
       F[s][k] = errTol;
     if(F[s][k]>1-errTol)
       F[s][k] =1-errTol;
     
   }


}



void map2domainQ(double** Q, int nInd, int K){
  for(int i=0;i<nInd;i++){  
    double sum=0;
    for(int k=0;k<K;k++){
      if(Q[i][k]<errTol)
	Q[i][k] = errTol;
      if(Q[i][k]>(1-errTol))
	Q[i][k] =1-errTol;
      sum+=Q[i][k];
    }
    for(int k=0;k<K;k++)
      Q[i][k]=Q[i][k]/sum;
  }
 
}


void updateQ(double** Q_0, double** F_0, int nSites_start,int nSites_stop, int nInd, int K,double **Q_1,int totSites,int startI,int stopI,char **keeps,double ***prodA,double ***prodB,int nThreads){
 for(int i=startI;i<stopI;i++){  
    for(int k=0;k<K;k++){
      int nMisTotSites=totSites;
#ifdef DO_MIS
      nMisTotSites=0;
      for(int j=0;j<totSites;j++){

      	if(keeps[j][i]==0)
      	  continue;

      	nMisTotSites++;
      }
#endif
      //      fprintf(stdout,"nThreads %d\n",nThreads);
      double prod1=0;
      double prod2=0;
      for(int t=0;t<nThreads;t++){
	prod1 += prodA[t][i][k];
	prod2 += prodB[t][i][k];

      }
      //      fprintf(stdout,"nThreads %d %d %f %f\n",nThreads,i,prod1,prod2);
      //      exit(0);
      double sumAGBG = prod1 * Q_0[i][k] + prod2 * Q_0[i][k];
      Q_1[i][k]=sumAGBG/(2*totSites);
     }
  }
}


int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}



std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

gzFile openFileGz(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  gzFile fp = gzopen(c,"w");
  delete [] c;
  return fp;
}



//some struct will all the data from the beagle file
typedef struct{
  double **genos;
  char *major;
  char *minor;
  char **ids;
  int nSites;
  int nInd;
  char **keeps; //matrix containing 0/1 indicating if data or missing
  int *keepInd; //keepInd[nSites] this is the number if informative samples
  float *mafs;
}bgl;

//utility function for cleaning up out datastruct
void dalloc(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];
    free(b.ids[i]);
  }
  delete [] b.minor;
  delete [] b.major;
  delete [] b.genos;
  delete [] b.ids;
}

/*
  plug in random variables

// Q - an nInd x K matrix of population proportions (random)
// F - a nSites x K matrix of frequencies (random)

Q sums to one for every sample
F doesn't do crap
 */

/*
  Returns the bgl struct containing all data from a beagle file.

  It find the nsamples from counting the header
  It finds the number of sites by queing every line in a std::vector
  After the file has been read intotal it reloops over the lines in the vector and parses data
 */

bgl readBeagle(const char* fname) {
  const char *delims = "\t \n";
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  
  bgl ret;
  char buf[LENS];

  //find number of columns
  gzgets(fp,buf,LENS);
  strtok(buf,delims);
  int ncols=1;
  while(strtok(NULL,delims))
    ncols++;
  if(0!=(ncols-3) %3 ){
    fprintf(stderr,"ncols=%d\n",ncols);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;//this is the number of samples
  
  //read every line into a vector
  std::vector<char*> tmp;
  while(gzgets(fp,buf,LENS))
    tmp.push_back(strdup(buf));
  
  //now we now the number of sites
  ret.nSites=tmp.size();
  ret.major= new char[ret.nSites];
  ret.minor= new char[ret.nSites];
  ret.ids = new char*[ret.nSites];
  ret.genos= new double*[ret.nSites];

  //then loop over the vector and parsing every line
  for(int s=0;SIG_COND&& (s<ret.nSites);s++){
    ret.ids[s] = strdup(strtok(tmp[s],delims));
    ret.major[s] =strtok(NULL,delims)[0];
    ret.minor[s] =strtok(NULL,delims)[0];
    ret.genos[s]= new double[3*ret.nInd];
    for(int i=0;i<ret.nInd*3;i++){
      ret.genos[s][i] = atof(strtok(NULL,delims));
      if(ret.genos[s][i]<0){
	fprintf(stderr,"Likelihoods must be positive\n");
	fprintf(stderr,"site %d ind %d geno %d has value %f\n",s,int(i*1.0/3),i%3,ret.genos[s][i]);
	exit(0);
      }
    }
    for(int i=0;i<ret.nInd;i++){
      double tmpS = 0.0;
      for(int g=0;g<3;g++)
	tmpS += ret.genos[s][i*3+g];
      if(!(tmpS>0)){
	fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	fprintf(stderr,"individual %d site %d has sum %f\n",i,s,tmpS);
	exit(0);
      } 
    }
    free(tmp[s]);
  }
  
  // Here additional stuff is calculated from the likelihoods
  // this must be done again while filtering later in main
  ret.keeps=new char*[ret.nSites]; // array nSites x nInd 0 if missing info
  ret.keepInd = new int[ret.nSites];
  ret.mafs = new float[ret.nSites];

  for(int s=0;s<ret.nSites;s++){
    ret.keeps[s] = new char[ret.nInd];
    int nKeep =0;
    for(int i=0;i<ret.nInd;i++){
      double mmin=std::min(ret.genos[s][i*3],std::min(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      double mmax=std::max(ret.genos[s][i*3],std::max(ret.genos[s][i*3+1],ret.genos[s][i*3+2]));
      if(fabs(mmax-mmin)<misTol)
	ret.keeps[s][i] =0;
      else{
	ret.keeps[s][i] =1;
	nKeep++;
      }
    }
    ret.keepInd[s] = nKeep;
    ret.mafs[s] = emFrequency(ret.genos[s],ret.nInd,MAF_ITER,MAF_START,ret.keeps[s],ret.keepInd[s]);
  }
  //  keeps=ret.keeps;
  gzclose(fp); //clean up filepointer
  return ret;
}

void readDouble(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
    if(NULL==fgets(buf,lens,fp)){
      fprintf(stderr,"Increase buffer\n");
      exit(0);
    }
    if(neg)
      d[i][0] = -atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = -atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  fclose(fp);
}

void readDoubleGZ(double **d,int x,int y,const char*fname,int neg){
  fprintf(stderr,"opening : %s with x=%d y=%d\n",fname,x,y);
  const char*delims=" \n";
  gzFile fp = NULL;
  if((fp=gzopen(fname,"r"))==NULL){
    fprintf(stderr,"cont open:%s\n",fname);
    exit(0);
  }
  int lens=1000000 ;
  char buf[lens];
  for(int i=0;i<x;i++){
      if(NULL==gzgets(fp,buf,lens)){
      fprintf(stderr,"Increase buffer\n");
      exit(0);
    }
    if(neg)
      d[i][0] = -atof(strtok(buf,delims));
    else
      d[i][0] = atof(strtok(buf,delims));
    for(int j=1;j<y;j++){
      //      fprintf(stderr,"i=%d j=%d\n",i,j);
      if(neg)
	d[i][j] = -atof(strtok(NULL,delims));
      else
	d[i][j] = atof(strtok(NULL,delims));
    }
  }
  gzclose(fp);
}

double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}

double ***allocDouble3(size_t x,size_t y,size_t z){
  double ***ret= new double**[x];
  for(size_t i=0;i<x;i++){
    ret[i] = new double*[y];
    for(size_t j=0;j<y;j++)
      ret[i][j] = new double[z];
  }
  return ret;
}



void dallocDouble3(double ***ret,size_t x,size_t y){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      delete[] ret[i][j];
  }

  for(size_t i=0;i<x;i++){
    delete ret[i];
  }

  delete ret;
}

void minus(double **fst,double **sec,size_t x,size_t y,double **res){
  //  fprintf(stderr,"x=%lu y=%lu\n",x,y);
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++){
      //  fprintf(stderr,"i=%lu j=%lu\n",i,j);
      res[i][j] = fst[i][j]-sec[i][j];
    }
}

double sumSquare(double **mat,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++)
      tmp += mat[i][j]*mat[i][j];
  return tmp;
}


double sumSquareMinus(double **mat1,double **mat2,size_t x,size_t y){
  double tmp=0;
  for(size_t i=0;i<x;i++)
    for(size_t j=0;j<y;j++){
      double tmp2 = mat1[i][j]-mat2[i][j];
      tmp += tmp2*tmp2;
    }
  return tmp;
}


void dalloc(size_t x,double **ret){
  for(size_t i=0;i<x;i++)
    delete [] ret[i] ;
  delete [] ret;
}
//same as above but swapped.... to lazy to change code
void dalloc(double **ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret[i] ;
  delete [] ret;
}

void printDouble(double **ret,size_t x,size_t y,FILE *fp){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      fprintf(fp,"%.20f ",ret[i][j]);
    fprintf(fp,"\n");
  }
}

void printDoubleGz(double **ret,size_t x,size_t y,gzFile fp){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      gzprintf(fp,"%.20f ",ret[i][j]);
    gzprintf(fp,"\n");
  }
}

void checkFQ(double **F,double **Q,int nSites,int nInd,int K,const char *function){
  
  for(int i = 0; i < nInd; i++) {
    double sum=0;
    for(int k = 0; k < K; k++) {
      //	fprintf(stderr,"%f\t",Q[i][k]);
      sum+=Q[i][k];
      if(Q[i][k]>1 || Q[i][k]<0){
	fprintf(stderr,"[%s] error Q %f\n",function,Q[i][k]);
	exit(0);
      }
    }
    if(sum>1+errTol || sum < 1- errTol){
      fprintf(stderr,"[%s] error sumQ %f\n",function,sum);
      exit(0);
    }
  }
  for(int j = 0; j < nSites; j++) {
    for(int k = 0; k < K; k++) {
      if(F[j][k]<0 || F[j][k] > 1-0){
	fprintf(stderr,"[%s] error freq %f\n",function,F[j][k]);
	exit(0);
      }
    }
  }
}

int printer =0;

//threaded
double likelihood(double** Q, double** F,int nSites, int nInd,int K,double **genos){
  // F is sites times npop
  // Q is nInd times npop
  double prod_sit = 0.0;
  for(int i = 0; i < nInd; i++) {
    double prod_ind = 0.0;
    for(int j = 0; j < nSites; j++) {
      double *gg = genos[j]+i*3;
      double freq = 0.0;
      for(int k = 0; k < K; k++) 
	freq += (F[j][k])*Q[i][k];
      double f =1- freq;
      double sum = gg[0] * f * f;
      sum += gg[1]*2*f*(1-f);
      sum += gg[2]*(1-f)*(1-f);
      prod_ind += log(sum);
    }
    prod_sit += prod_ind;     
  }

  #ifdef CHECK
    checkFQ(F,Q,nSites,nInd,K,__FUNCTION__);
  #endif
  return prod_sit;
}

// log likelihood for single site
// fills it directly into the structure
void likelihood_thd(double** Q, double** F,int nSites,int startI, int stopI,int K,double **genos,double &prod_sit){
  // F is sites times npop
  // Q is nInd times npop
  prod_sit =0;
  for(int i = startI; i < stopI; i++) {
    double prod_ind = 0.0;
    for(int j = 0; j < nSites; j++) {
      double *gg = genos[j]+i*3;
      double freq = 0.0;
      for(int k = 0; k < K; k++) 
	freq += (F[j][k])*Q[i][k];
      double f =1- freq;
      double sum = gg[0] * f * f;
      sum += gg[1]*2*f*(1-f);
      sum += gg[2]*(1-f)*(1-f);
      prod_ind += log(sum); //collect all site for a single individual
    }
    prod_sit += prod_ind;     
  }
  //checkFQ(F,Q,nSites,nInd,K,__FUNCTION__);
}
void *lkWrap(void *a){
  pars *p = (pars *)a; 
  likelihood_thd(p->Q,p->F,p->nSites,p->startI,p->stopI,p->nPop,p->genos,p->lres);
  return NULL;
}

// total log likelihood
double like_tsk(double **Q,double **F,int nThreads){
  for(int i=0;i<nThreads;i++){
    myPars[i].Q=Q;
    myPars[i].F=F;
    if( pthread_create(&threads[i],NULL,lkWrap,&myPars[i]))// here loglike is calculated in threads
      fprintf(stderr,"Problems starting threads\n");
  }
  for(int i=0;i<nThreads;i++)
    if(pthread_join(threads[i], NULL))
      fprintf(stderr,"problems joining\n");
  double res=0;
  for(int i=0;i<nThreads;i++){
    //    fprintf(stderr,"lres[%d]=%f\n",i,myPars[i].lres);
    res += myPars[i].lres;
  }
  return res;
}


void em(double** Q, double** F, int nSites, int nInd, int K,double **genos,double **F_1,double **Q_1) {
  //  fprintf(stderr,"no threads no sqem\n");
  #ifdef CHECK
  checkFQ(F,Q,nSites,nInd,K,"em lorte start ");
  #endif
  //tmp vars
  static double **bNorm=NULL;
  static double **aNorm=NULL;
  static double ***prodA = NULL;
  static double ***prodB = NULL;
  if(bNorm==NULL){
    bNorm= allocDouble(nSites,nInd);
    aNorm = allocDouble(nSites,nInd);
    prodA=allocDouble3(1,nInd,K);
    prodB=allocDouble3(1,nInd,K);
  }
  for(int i=0;i<nInd;i++)
    for(int k=0;k<K;k++){
      prodA[0][i][k]=0;
      prodB[0][i][k]=0;
    }

  if(Q==NULL){//cleanup
    //    fprintf(stderr,"cleanup\n");
    dalloc(bNorm,nSites);
    dalloc(aNorm,nSites);
    dallocDouble3(prodA,1,nInd);
    dallocDouble3(prodB,1,nInd);

    return;
  }

 


  //  getExpGandExpG2(Q,F,0,nSites,nInd,K,genos,expG2,expG1,0,nInd);
  getExpGandFstar(Q,F,0,nSites,nInd,K,genos,0,nInd,keeps,prodA,prodB,F_1,0);
  //  updateF(Q, F,0,nSites,nInd,K,F_1,expG2,expG1,0,nInd);
  //  updateF(Q, F,0,nSites,nInd,K,F_1,expG2,expG1,0,nInd,keeps);
  //  updateQ(Q,F,0,nSites,nInd,K,Q_1,nSites,expG2,expG1,0,nInd);
  updateQ(Q,F,0,nSites,nInd,K,Q_1,nSites,0,nInd,keeps,prodA,prodB,1);

 
  map2domainQ(Q_1,nInd,K);
  #ifdef CHECK
  checkFQ(F_1,Q_1,nSites,nInd,K,__FUNCTION__);
  #endif
  //    fprintf(stderr,"\tq1=%f q2=%f q3=%f\n",Q_1[0][0],Q_1[0][1],Q_1[0][2]);
  // Set step n = n+1
  ///fixup underoverflow
}

pthread_barrier_t barr;
int dumpOld =0;

//Q,F are the old F_1,Q_1 are the next startpoint

void emSQ_threads(double** Q, double** F, int nSites_start,int nSites_stop, int nInd, int K,double **genos,double **F_1,double **Q_1,int totSites,double ***prodA, double*** prodB,int startI,int stopI,int threadNumber,int nThreads) {


  getExpGandFstar(Q,F,nSites_start,nSites_stop,nInd,K,genos,startI,stopI,keeps,prodA,prodB,F_1,threadNumber);
   //we wait untill all threads reach this point.
  int rc = pthread_barrier_wait(&barr);
#ifdef __APPLE__ 
  if(rc != 0 && rc !=1){
#else
  if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD){
#endif
        fprintf(stderr,"If here then we have problems with the barriere, shouldn't happen\n");
      exit(-1);
    }
  //update F and Q 
  updateQ(Q,F,nSites_start,nSites_stop,nInd,K,Q_1,totSites,startI,stopI,keeps,prodA,prodB,nThreads);
  if(1){
    ///fixup underoverflow
    for(int s=nSites_start;s<nSites_stop;s++)
      for(int k=0;k<K;k++){
	if(F_1[s][k]<errTol)
	  F_1[s][k] = errTol;
	if(F_1[s][k]>1-errTol)
	  F_1[s][k] =1-errTol;
      }
    for(int j = startI; j < stopI; j++) 
      for(int k = 0; k < K; k++) {
	if(Q_1[j][k]<errTol)
	  Q_1[j][k] = errTol;
	if(Q_1[j][k]>(1-errTol))
	  Q_1[j][k] =1-errTol;
      }
  }
}


void *emWrap(void *a){
  pars *p = (pars *)a; 
  emSQ_threads(p->Q,p->F,p->start,p->stop,p->nInd,p->nPop,p->genos,p->F_1,p->Q_1,p->nSites,p->prodA,p->prodB,p->startI,p->stopI,p->threadNumber,p->nThreads);
  return NULL;
}

void em_threadStart(double **Q,double **Q_new,double **F,double **F_new,int nThreads){

  for(int t=0;t<nThreads;t++)
    for(int i=0;i<myPars[0].nInd;i++)
      for(int k=0;k<myPars[0].nPop;k++){
	myPars[0].prodA[t][i][k]=0;
	myPars[0].prodB[t][i][k]=0;
      }


  #ifdef CHECK
  checkFQ(F,Q,myPars[0].nSites,myPars[0].nInd,myPars[0].nPop,"lorte start\n");
  #endif
  for(int i=0;i<nThreads;i++){
    myPars[i].Q=Q;myPars[i].F=F;
    myPars[i].Q_1=Q_new;myPars[i].F_1=F_new;
    if( pthread_create(&threads[i],NULL,emWrap,&myPars[i]))
      fprintf(stderr,"Problems starting threads\n");
  }
  for(int i=0;i<nThreads;i++)
    if(pthread_join(threads[i], NULL))
      fprintf(stderr,"problems joining\n");
  
  for(int i=0;i<myPars[0].nInd;i++){  
    double sum=0;
    for(int k=0;k<myPars[0].nPop;k++){
      if(Q_new[i][k]<errTol)
	Q_new[i][k] = errTol;
      if(Q_new[i][k]>(1-errTol))
	Q_new[i][k] =1-errTol;
      sum+=Q_new[i][k];
    }
    for(int k=0;k<myPars[0].nPop;k++)
      Q_new[i][k]=Q_new[i][k]/sum;
  }
  #ifdef CHECK
  checkFQ(F_new,Q_new,myPars[0].nSites,myPars[0].nInd,myPars[0].nPop,__FUNCTION__);
  #endif
}

//returnval =1 continue, rturnval =0, convergence has been achieved
int emAccel(const bgl &d,int nPop,double **F,double **Q,double ***F_new,double ***Q_new,double &lold){
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;

  //we make these huge structures static such that we just allocate them the first time
  static double **F_em1 =NULL;
  static double **Q_em1 =NULL;
  static double **F_diff1 =NULL;
  static double **Q_diff1 =NULL;
  static double **F_em2 =NULL;
  static double **Q_em2 =NULL;
  static double **F_diff2 =NULL;
  static double **Q_diff2 =NULL;
  static double **F_diff3 =NULL;
  static double **Q_diff3 =NULL;
  static double **F_tmp =NULL;
  static double **Q_tmp =NULL;

  if(F_em1==NULL){
    F_em1 =allocDouble(d.nSites,nPop);
    Q_em1 =allocDouble(d.nInd,nPop);
    F_diff1 =allocDouble(d.nSites,nPop);
    Q_diff1 =allocDouble(d.nInd,nPop);
    F_em2 =allocDouble(d.nSites,nPop);
    Q_em2 =allocDouble(d.nInd,nPop);
    F_diff2 =allocDouble(d.nSites,nPop);
    Q_diff2 =allocDouble(d.nInd,nPop);
    F_diff3 =allocDouble(d.nSites,nPop);
    Q_diff3 =allocDouble(d.nInd,nPop);
    F_tmp =allocDouble(d.nSites,nPop);
    Q_tmp =allocDouble(d.nInd,nPop);
  }
  //should cleanup and exit
  if(F==NULL){
    dalloc(F_em1,d.nSites);
    dalloc(F_em2,d.nSites);
    dalloc(F_diff1,d.nSites);
    dalloc(F_diff2,d.nSites);
    dalloc(F_diff3,d.nSites);
    dalloc(F_tmp,d.nSites);
    dalloc(Q_em1,d.nInd);
    dalloc(Q_em2,d.nInd);
    dalloc(Q_diff1,d.nInd);
    dalloc(Q_diff2,d.nInd);
    dalloc(Q_diff3,d.nInd);
    dalloc(Q_tmp,d.nInd);
    return 0;
  }
  #ifdef CHECK
  //  fprintf(stderr,"em Q F \n");
  #endif
  em(Q, F, d.nSites, d.nInd, nPop,d.genos,F_em1,Q_em1);
  minus(F_em1,F,d.nSites,nPop , F_diff1);
  minus(Q_em1,Q,d.nInd,nPop,Q_diff1 );
  
  double sr2 = sumSquare(F_diff1,d.nSites,nPop)+sumSquare(Q_diff1,d.nInd,nPop);
  if(sqrt(sr2)<tol){
    return 0;
    //break;
  }
   #ifdef CHECK
  //fprintf(stderr,"em Q1 F1 \n");
  #endif
  em(Q_em1, F_em1, d.nSites, d.nInd, nPop,d.genos,F_em2,Q_em2);
  
  minus(F_em2,F_em1,d.nSites,nPop , F_diff2);
  minus(Q_em2,Q_em1,d.nInd,nPop,Q_diff2 );
  double sq2 = sumSquare(F_diff2,d.nSites,nPop)+sumSquare(Q_diff2,d.nInd,nPop);
  if(sqrt(sq2)<tol){
    return 0;
    //    break;
  }
  minus(F_diff2,F_diff1,d.nSites,nPop , F_diff3);
  minus(Q_diff2,Q_diff1,d.nInd,nPop,Q_diff3 );
  double sv2 = sumSquare(F_diff3,d.nSites,nPop)+sumSquare(Q_diff3,d.nInd,nPop);
  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));

  //  fprintf(stderr,"alpha=%f %f %f\n",alpha,sq2,sv2);  
  //update with the linear combination
  for(size_t i=0;i<d.nSites;i++)
    for(size_t j=0;j<nPop;j++){
      (*F_new)[i][j] = F[i][j]+2*alpha*F_diff1[i][j]+alpha*alpha*F_diff3[i][j];

     
    }
  map2domainF(*F_new,0,d.nSites,nPop);

  for(size_t i=0;i<d.nInd;i++){
    for(size_t j=0;j<nPop;j++){
      (*Q_new)[i][j] = Q[i][j]+2*alpha*Q_diff1[i][j]+alpha*alpha*Q_diff3[i][j];
    }
  }
  map2domainQ(*Q_new,d.nInd,nPop);
  //  fprintf(stderr,"lik %f\n",-likelihood(*Q_new, *F_new, d.nSites, d.nInd, nPop,d.genos));
  if (fabs(alpha - 1) > 0.01){
    #ifdef CHECK
    //       fprintf(stderr,"em Q2 F2 \n");
    #endif
    em(*Q_new, *F_new, d.nSites, d.nInd, nPop,d.genos,F_tmp,Q_tmp);
    #ifdef CHECK
    checkFQ(F_tmp,Q_tmp,d.nSites,d.nInd,nPop,"bad em\n");
    #endif
    //    double ** tmp;
    //    Q_new = T_tmp
    //    *Q_new = Q_tmp;
    //    *F_new = F_tmp;
    std::swap(*Q_new,Q_tmp);
    std::swap(*F_new,F_tmp);
  }
  double lnew =0;
  #ifdef DO_MIS
  lnew = -likelihood(*Q_new, *F_new, d.nSites, d.nInd, nPop,d.genos);
  if ( (lnew > lold + objfnInc)) {
    fprintf(stderr,"bad guess in %s\n", __FUNCTION__);
    //    *Q_new = Q_em2;
    //    *F_new =F_em2;
    std::swap(*Q_new,Q_em2);
    std::swap(*F_new,F_em2);
    lnew = -likelihood(*Q_new, *F_new, d.nSites, d.nInd, nPop,d.genos);
    if (alpha == stepMax)
      stepMax = std::max(stepMax0, stepMax/mstep);
  } 
  #endif
  if (alpha == stepMax) {
    stepMax = mstep*stepMax;
  }
  #ifdef CHECK
  checkFQ(*F_new,*Q_new,d.nSites,d.nInd,nPop,"bad F new\n");
  #endif
  lold=lnew;
  return 1;
}

//returnval =1 continue, returnval =0, convergence has been achieved
int emAccelThread(const bgl &d,int nPop,double **F,double **Q,double ***F_new,double ***Q_new,double &lold,int nThreads){
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;
  
  //we make these huge structures static such that we just allocate them the first time
  static double **F_em1 =NULL;
  static double **Q_em1 =NULL;
  static double **F_diff1 =NULL;
  static double **Q_diff1 =NULL;
  static double **F_em2 =NULL;
  static double **Q_em2 =NULL;
  static double **F_diff2 =NULL;
  static double **Q_diff2 =NULL;
  static double **F_diff3 =NULL;
  static double **Q_diff3 =NULL;
  static double **F_tmp =NULL;
  static double **Q_tmp =NULL;

  if(F_em1==NULL){
    F_em1 =allocDouble(d.nSites,nPop);
    Q_em1 =allocDouble(d.nInd,nPop);
    F_diff1 =allocDouble(d.nSites,nPop);
    Q_diff1 =allocDouble(d.nInd,nPop);
    F_em2 =allocDouble(d.nSites,nPop);
    Q_em2 =allocDouble(d.nInd,nPop);
    F_diff2 =allocDouble(d.nSites,nPop);
    Q_diff2 =allocDouble(d.nInd,nPop);
    F_diff3 =allocDouble(d.nSites,nPop);
    Q_diff3 =allocDouble(d.nInd,nPop);
    F_tmp =allocDouble(d.nSites,nPop);
    Q_tmp =allocDouble(d.nInd,nPop);
  }

  if(F==NULL){
    dalloc(F_em1,d.nSites);
    dalloc(F_em2,d.nSites);
    dalloc(F_diff1,d.nSites);
    dalloc(F_diff2,d.nSites);
    dalloc(F_diff3,d.nSites);
    dalloc(F_tmp,d.nSites);
    dalloc(Q_em1,d.nInd);
    dalloc(Q_em2,d.nInd);
    dalloc(Q_diff1,d.nInd);
    dalloc(Q_diff2,d.nInd);
    dalloc(Q_diff3,d.nInd);
    dalloc(Q_tmp,d.nInd);
    return 0;
  }

  em_threadStart(Q,Q_em1,F,F_em1,nThreads);
  minus(F_em1,F,d.nSites,nPop,F_diff1);
  minus(Q_em1,Q,d.nInd,nPop,Q_diff1 );
  
  double sr2 = sumSquare(F_diff1,d.nSites,nPop)+sumSquare(Q_diff1,d.nInd,nPop);
  if(sqrt(sr2)<tol){
    return 0;
    //break;
  }
  em_threadStart(Q_em1,Q_em2,F_em1,F_em2,nThreads);
  

  minus(F_em2,F_em1,d.nSites,nPop,F_diff2);
  minus(Q_em2,Q_em1,d.nInd,nPop,Q_diff2 );
  double sq2 = sumSquare(F_diff2,d.nSites,nPop)+sumSquare(Q_diff2,d.nInd,nPop);
  if(sqrt(sq2)<tol){
    return 0;
    //    break;
  }
  minus(F_diff2,F_diff1,d.nSites,nPop,F_diff3);
  minus(Q_diff2,Q_diff1,d.nInd,nPop,Q_diff3);
  double sv2 = sumSquare(F_diff3,d.nSites,nPop)+sumSquare(Q_diff3,d.nInd,nPop);
  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  //  fprintf(stderr,"alpha=%f %f %f\n",alpha,sq2,sv2);  
  //  fprintf(stderr,"some fun\n");
  #ifdef CHECK
  checkFQ(F,Q,d.nSites,d.nInd,nPop,__FUNCTION__);
  //  fprintf(stderr,"some fun2\n");
  checkFQ(F_em1,Q_em1,d.nSites,d.nInd,nPop,__FUNCTION__);
  //  fprintf(stderr,"some fun3\n");
  checkFQ(F_em2,Q_em2,d.nSites,d.nInd,nPop,__FUNCTION__);
  //  fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);
  #endif
 //update with the linear combination
  for(size_t i=0;i<d.nSites;i++)
    for(size_t j=0;j<nPop;j++){
      (*F_new)[i][j] = F[i][j]+2*alpha*F_diff1[i][j]+alpha*alpha*F_diff3[i][j];
      //  if((*F_new)[i][j]<1e-8 || (*F_new)[i][j]>1-1e-8){
      //  fprintf(stderr,"Fnew: %f F: %f Fem1: %f Fem2: %f alpha: %f diff1: %f diff3: %f\n",(*F_new)[i][j],F[i][j],F_em1[i][j],F_em2[i][j],alpha,F_diff1[i][j],F_diff3[i][j]);
      //(*F_new)[i][j]=F_em2[i][j];
      //}
      //	
    }
  map2domainF(*F_new,0,d.nSites,nPop);

  for(size_t i=0;i<d.nInd;i++){
    for(size_t j=0;j<nPop;j++){
      (*Q_new)[i][j] = Q[i][j]+2*alpha*Q_diff1[i][j]+alpha*alpha*Q_diff3[i][j];
      // if((*Q_new)[i][j]<1e-8||(*Q_new)[i][j]>1-1e-8)
      //	fprintf(stderr,"(*Q_new)[i][j]: %f\n",(*Q_new)[i][j]);

      }
  }
  map2domainQ(*Q_new,d.nInd,nPop);
  

  //fprintf(stderr,"lik %f\n",-likelihood(*Q_new, *F_new, d.nSites, d.nInd, nPop,d.genos));
  //  checkFQ(*F_new,*Q_new,d.nSites,d.nInd,nPop,__FUNCTION__);
  if (fabs(alpha - 1) > 0.01){
    em_threadStart(*Q_new,Q_tmp,*F_new,F_tmp,nThreads);
    //*Q_new= Q_tmp;
    //*F_new = F_tmp;
    std::swap(*Q_new,Q_tmp);
    std::swap(*F_new,F_tmp);
  }
  
  //  double lnew = -likelihood(Q_new, F_new, d.nSites, d.nInd, nPop,d.genos);
  double lnew = 1;
  #ifdef DO_MIS
  lnew = -like_tsk(*Q_new, *F_new,nThreads);
  
  if ( (lnew > lold + objfnInc)) {
    //*Q_new= Q_em2;
    // *F_new = F_em2;
    fprintf(stderr,"bad guess in %s\n", __FUNCTION__);
    std::swap(*Q_new,Q_em2);
    std::swap(*F_new,F_em2);
    // lnew = -likelihood(Q_new, F_new, d.nSites, d.nInd, nPop,d.genos);
    lnew = -like_tsk(*Q_new, *F_new, nThreads);
    if ((alpha - stepMax) > -0.001)
      stepMax = std::max(stepMax0, stepMax/mstep);
  } 
  #endif
  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }
  //fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);
 
  lold=lnew;
  return 1;
}


void info(){
  
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\t-likes Beagle likelihood filename\n");
  fprintf(stderr,"\t-K Number of ancestral populations\n"); 
  fprintf(stderr,"Optional:\n");
  fprintf(stderr,"\t-fname Ancestral population frequencies\n"); 
  fprintf(stderr,"\t-qname Admixture proportions\n"); 
  fprintf(stderr,"\t-outfiles Prefix for output files\n"); 
  fprintf(stderr,"\t-printInfo print ID and mean maf for the SNPs that were analysed\n"); 

  fprintf(stderr,"Setup:\n"); 
  fprintf(stderr,"\t-seed Seed for initial guess in EM\n"); 
  fprintf(stderr,"\t-P Number of threads\n"); 
  fprintf(stderr,"\t-method If 0 no acceleration of EM algorithm\n"); 
  fprintf(stderr,"\t-misTol Tolerance for considering site as missing\n");

  fprintf(stderr,"Stop chriteria:\n"); 
  fprintf(stderr,"\t-tolLike50 Loglikelihood difference in 50 iterations\n"); 
  fprintf(stderr,"\t-tol Tolerance for convergence\n"); 
  fprintf(stderr,"\t-dymBound Use dymamic boundaries (1: yes (default) 0: no)\n"); 
  fprintf(stderr,"\t-maxiter Maximum number of EM iterations\n"); 


  fprintf(stderr,"Filtering\n"); 
  fprintf(stderr,"\t-minMaf Minimum minor allele frequency\n"); 
  fprintf(stderr,"\t-minLrt Minimum likelihood ratio value for maf>0\n"); 
  fprintf(stderr,"\t-minInd Minumum number of informative individuals\n");


  exit(0);
}



int VERBOSE =1;
void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"Caught SIGNAL: %d will try to exit nicely (no more threads are created, we will wait for the current threads to finish)\n",s);
  VERBOSE=0;
  SIG_COND=0;
}

float calcThres(double **d1,double **d2, int x,int y){
  // finds the largest difference between 2 arrays
  // arrays has dimention x times y
  float diff=0;
  for(int i=0;i<x;i++)
    for(int j=0;j<y;j++){
      //      fprintf(stderr,"diiffs=%f\n",fabs(d1[i][j]-d2[i][j]));
      if(fabs(d1[i][j]-d2[i][j])>diff)
	diff=fabs(d1[i][j]-d2[i][j]);
    }
  return diff;
}

int whichMax(double *g){
  // equality signs such that when some are equal
  // the lowest genotype is preferred
  if(g[0]>=g[1]&&g[0]>=g[2])
    return 0;
  if(g[1]>g[0]&&g[1]>=g[2])
    return 1;
  else
    return 2;
}

void printLikes(bgl &d){
  // to write likelihoods for debugging
  for(int s=0;s<d.nSites;s++){
    double *g  = d.genos[s];
    for(int i=0;i<d.nInd;i++){
      for(int j=0;j<3;j++)
	fprintf(stdout,"%f\t",g[3*i+j]);
    }
    fprintf(stdout,"\n");
  }
  exit(0);
}

void printKeepSites(bgl &d,FILE *ffilter){
  fprintf(ffilter,"marker\tmajor\tminor\tmaf\tnonMis\n");
 for(int s=0;s<d.nSites;s++){
   fprintf(ffilter,"%s\t%d\t%d\t%f\t%d\n",d.ids[s],d.major[s],d.minor[s],d.mafs[s],d.keepInd[s]);

 }
}

//void modLikesMinMaf(bgl &d,float minMaf){
void filterMinMaf(bgl &d,float minMaf){
  //  fprintf(stderr,"WARNING filtering minMaf=%f \n",minMaf);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    //    fprintf(stderr,"minmaf=%f mafs[%d]=%f lik=%f\n",minMaf,s,d.mafs[s],lik);
    if(d.mafs[s]>minMaf&&d.mafs[s]<1-minMaf){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      //      fprintf(stderr,"skippping\n");
    }
  }
  d.nSites=posi;
}

//void modLikesMiss(bgl &d,int minInd){
void filterMiss(bgl &d,int minInd){
  //  fprintf(stderr,"WARNING filtering mis=%d \n",minInd);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    if(d.keepInd[s]>minInd){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      // fprintf(stderr,"skippping\n");
    }
  }
  d.nSites=posi;
}

//void modLikesMinLrt(bgl &d,float minLrt){
void filterMinLrt(bgl &d,float minLrt){
  //  fprintf(stderr,"WARNING filtering minlrt=%f \n",minLrt);
  int posi =0;
  for(int s=0;s<d.nSites;s++){
    float lik=likeFixedMinor(d.mafs[s],d.genos[s],d.nInd,d.keeps[s]);
    float lik0=likeFixedMinor(0.0,d.genos[s],d.nInd,d.keeps[s]);
    //    fprintf(stderr,"minlrt=%f mafs[%d]=%f lik=%f lik0=%f 2*(lik0-lik)=%f\n",minLrt,s,d.mafs[s],lik,lik0,2.0*(lik0-lik));
    if(2.0*(lik0-lik)>minLrt){
      d.genos[posi] = d.genos[s];
      d.major[posi] = d.major[s];
      d.minor[posi] = d.minor[s];
      d.ids[posi] = d.ids[s];
      d.keeps[posi] = d.keeps[s];
      d.keepInd[posi] = d.keepInd[s];
      d.mafs[posi] = d.mafs[s];
      posi++;   
    }else{
      //  fprintf(stderr,"skippping\n");
    }
  }
  d.nSites=posi;
}


int main(int argc, char **argv){
  if(argc==1){// if no arguments, print info on program
    info();
    return 0;
  }
  //below for catching ctrl+c, and dumping files
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  //initial values
  int dymBound = 0;
  int maxIter = 2000;
  int method = 1;
  int minInd = 0;
  int printInfo = 0;
  float minMaf =0.05;
  float minLrt =0;
  const char* lname = NULL;
  const char* fname = NULL;
  const char* qname = NULL;
  const char* outfiles = NULL;
  int nPop = 3;
  int seed =time(NULL);
  int nThreads = 1;
  float tolLike50=0.1;

  // reading arguments
  argv++;
  while(*argv){
    if(strcmp(*argv,"-likes")==0 || strcmp(*argv,"-l")==0) lname=*++argv; 
    else if(strcmp(*argv,"-K")==0) nPop=atoi(*++argv); 
    // to read start values from output from previous run 
    else if(strcmp(*argv,"-fname")==0 || strcmp(*argv,"-f")==0) fname=*++argv; 
    else if(strcmp(*argv,"-qname")==0 || strcmp(*argv,"-q")==0) qname=*++argv;
    // prefix for output files
    else if(strcmp(*argv,"-outfiles")==0 || strcmp(*argv,"-o")==0) outfiles=*++argv; 
    // settings: seed, threads and if method==0 not accelerated
    else if(strcmp(*argv,"-seed")==0||strcmp(*argv,"-s")==0) seed=atoi(*++argv);
    else if(strcmp(*argv,"-P")==0) nThreads=atoi(*++argv); 
    else if(strcmp(*argv,"-printInfo")==0) printInfo=atoi(*++argv); 
    else if(strcmp(*argv,"-method")==0 || strcmp(*argv,"-m")==0) method=atoi(*++argv); 
    // different stop chriteria
    else if(strcmp(*argv,"-tolLike50")==0||strcmp(*argv,"-lt50")==0) tolLike50=atof(*++argv);
    else if(strcmp(*argv,"-tol")==0||strcmp(*argv,"-t")==0) tol=atof(*++argv);
    else if(strcmp(*argv,"-maxiter")==0 || strcmp(*argv,"-i")==0) maxIter=atoi(*++argv); 
    // different filterings
    else if(strcmp(*argv,"-misTol")==0 || strcmp(*argv,"-mt")==0) misTol=atof(*++argv);
    else if(strcmp(*argv,"-minMaf")==0||strcmp(*argv,"-maf")==0) minMaf=atof(*++argv);
    else if(strcmp(*argv,"-minLrt")==0||strcmp(*argv,"-lrt")==0) minLrt=atof(*++argv);
    else if(strcmp(*argv,"-minInd")==0||strcmp(*argv,"-mis")==0) minInd=atoi(*++argv);
    // different genotype callers
    else if(strcmp(*argv,"-dymBound")==0) dymBound=atoi(*++argv);
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }
  if(lname==NULL){
    fprintf(stderr,"Please supply beagle file: -likes");
    info();
  }
  if(outfiles==NULL){
    fprintf(stderr,"Will use beagle fname as prefix for output\n");
    outfiles=lname;
  }
  FILE *flog=openFile(outfiles,".log");
  FILE *ffilter=openFile(outfiles,".filter");

  fprintf(stderr,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(stderr,"Setup: seed=%d nThreads=%d method=%d\n",seed,nThreads,method);
  fprintf(stderr,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(stderr,"Filters: misTol=%f minMaf=%f minLrt=%f minInd=%d\n",misTol,minMaf,minLrt,minInd);


  fprintf(flog,"Input: lname=%s nPop=%d, fname=%s qname=%s outfiles=%s\n",lname,nPop,fname,qname,outfiles);
  fprintf(flog,"Setup: seed=%d nThreads=%d method=%d\n",seed,nThreads,method);
  fprintf(flog,"Convergence: maxIter=%d tol=%f tolLike50=%f dymBound=%d\n",maxIter,tol,tolLike50,dymBound);
  fprintf(flog,"Filters: misTol=%f minMaf=%f minLrt=%f minInd=%d\n",misTol,minMaf,minLrt,minInd);

  if(dymBound==0){
    errTolStart = errTolMin;
    errTol = errTolMin;
  }
    
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);
  
  bgl d=readBeagle(lname);
  fprintf(stderr,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);

  // filter sites
  if(minMaf!=0.0)
    filterMinMaf(d,minMaf);
  if(minLrt!=0.0)
    filterMinLrt(d,minLrt);
  if(minInd!=0)
    filterMiss(d,minInd);
   if(printInfo)
    printKeepSites(d,ffilter);

  #ifdef DO_MIS
  keeps = d.keeps;
  #endif

  //  printLikes(d);
  fprintf(stderr,"Input file has dim (AFTER filtering): nsites=%d nind=%d\n",d.nSites,d.nInd);
  fprintf(flog,"Input file has dim (AFTER filtering): nsites=%d nind=%d\n",d.nSites,d.nInd);
  fflush(stderr);
  
  //set seed
  srand(seed);
  //unknown parameters
  double **F =allocDouble(d.nSites,nPop);
  double **Q =allocDouble(d.nInd,nPop);
  double **F_new =allocDouble(d.nSites,nPop);
  double **Q_new =allocDouble(d.nInd,nPop);
  

  //get start values
  if(fname==NULL){
    for(int j=0;j<d.nSites;j++)
      for(int k=0;k<nPop;k++){
	F[j][k]=rand()*1.0/RAND_MAX;
	//F[j][k]/1.01+0.005;
	F_new[j][k]=F[j][k];
      }
  }else
    readDoubleGZ(F,d.nSites,nPop,fname,0);
  if(qname==NULL){
    for(int i=0;i<d.nInd;i++) {
      double sum=0;
      for(int k=0;k<nPop;k++){
	Q[i][k]=rand()*1.0/RAND_MAX;
	//Q[i][k]/1.01+0.005;
	sum+=Q[i][k];
      }
      for(int k=0;k<nPop;k++) {
	Q[i][k]= Q[i][k]/sum; 
	Q_new[i][k]=Q[i][k];
      }
    }
  }else
    readDouble(Q,d.nInd,nPop,qname,0);
  
  //  double res =likelihood(Q, F, d.nSites, d.nInd, nPop,d.genos);    
  //  fprintf(stderr,"startres=%f\n",res);
  //emsquare stuff beloq
  
  //  if(1) {
  
  //update the global stuff NOW
  threads = new pthread_t[nThreads];//only used if ntreads!=1
  myPars  = new pars[nThreads];
  //update the internal stuff in the pars for the threading
  if(nThreads!=1) {

    double ***prodA;
    double ***prodB;
    prodA=allocDouble3(nThreads,d.nInd,nPop);
    prodB=allocDouble3(nThreads,d.nInd,nPop);

    int offsets[nThreads+1];
    int offsets2[nThreads+1];
    for(int i=0;i<nThreads;i++){
      prodA[i] = allocDouble(d.nInd,nPop);
      prodB[i] = allocDouble(d.nInd,nPop);
    }
    offsets[0] =0;
    offsets2[0] =0;
    for(int i=1;i<nThreads;i++){

      offsets[i] = d.nSites/nThreads + offsets[i-1];
      offsets2[i] = d.nInd/nThreads + offsets2[i-1];
    }
    offsets[nThreads] = d.nSites;
    offsets2[nThreads] = d.nInd;
    
    for(int i=0;i<nThreads;i++){
      myPars[i].nThreads=nThreads;
      myPars[i].threadNumber=i;
      myPars[i].start=offsets[i];
      myPars[i].stop=offsets[i+1];
      myPars[i].startI=offsets2[i];
      myPars[i].stopI=offsets2[i+1];
      //	fprintf(stderr,"i=%d sites=(%d, %d): nind=(%d, %d)\n",i,myPars[i].start,myPars[i].stop,myPars[i].startI,myPars[i].stopI);
      myPars[i].prodA=prodA;
      myPars[i].prodB=prodB;
      myPars[i].nPop=nPop;
      myPars[i].genos=d.genos;
      myPars[i].nInd=d.nInd;
      myPars[i].nSites=d.nSites;
    }
    //return 0;
    if(pthread_barrier_init(&barr, NULL, nThreads)){
      fprintf(stderr,"Could not create a barrier\n");
      return -1;
    }	
 
  }
  
  double lold = -likelihood(Q, F, d.nSites, d.nInd, nPop,d.genos);
  fprintf(stderr,"iter[start] like is=%f\n",lold);
  
  //below is the main looping trhought the iterations.
  // we have 4 possible ways, threading/nothreading line/noline
  int nit;
  double likeLast= lold;
  #ifdef CHECK
  checkFQ(F,Q,d.nSites,d.nInd,nPop,"bad start guess");
  #endif
  for(nit=1;SIG_COND&& nit<maxIter;nit++) {
    if(nThreads==1){
      if(method==0)//no acceleration
	em(Q, F, d.nSites, d.nInd, nPop,d.genos,F_new,Q_new);
      else{
	if(emAccel(d,nPop,F,Q,&F_new,&Q_new,lold)==0){

	  if(errTol>errTolMin){
	    errTol=errTol/5;
	    if(errTol<errTolMin)
	      errTol=errTolMin;
	    //	    fprintf(stderr,"changing errTol to %f\n",errTol);
	  }
	  else{
	    fprintf(stderr,"EM accelerated has reached convergence with tol %f\n",tol);
	    break; //if we have achieved convergence
	  }
	}
      }
    }else {
      if(method ==0) //no acceleration
	em_threadStart(Q,Q_new,F,F_new,nThreads);
      else{
	if(emAccelThread(d,nPop,F,Q,&F_new,&Q_new,lold,nThreads)==0){
	  if(errTol>errTolMin){
	    errTol=errTol/5;
	    if(errTol<errTolMin)
	      errTol=errTolMin;
	    //	    fprintf(stderr,"changing errTol to %f\n",errTol);
	  }
	  else{
	    fprintf(stderr,"EM accelerated Thread has reached convergence with tol %f\n",tol);
	    break; //if we have achieved convergence
	  }
	}
      }
      
    }
    
    std::swap(Q,Q_new);
    std::swap(F,F_new);
    
    if((nit%50)==0 ){ //stopping criteria
      double lik = likelihood(Q, F, d.nSites, d.nInd, nPop,d.genos);
      // thres is largest differense in admixture fractions
      fprintf(stderr,"iter[%d] like is=%f thres=%f\n",nit,lik,calcThres(Q,Q_new,d.nInd,nPop));
      
      //	fprintf(stderr,"iter[%d] like is=%f like old is=%f %f %f \n",nit,lik,likeLast, lik+likeLast, tolLike50);
      if(errTol>errTolMin){
	errTol=errTol/10;
	if(errTol<errTolMin)
	  errTol=errTolMin;
	//	fprintf(stderr,"50 changing errTol to %f\n",errTol);
      }
      else if(lik+likeLast < tolLike50){
	fprintf(stderr,"Convergence achived because log likelihooditer difference for 50 iteraction is less than %f\n",tolLike50);
	if(lik+likeLast<-1){
	  fprintf(stderr,"Convergence achived because log likelihooditer difference was NEGATIVE\n");
	}
	break;
      }
      
      likeLast=-lik;
      
    }
    
  }
  if(dumpOld){
    FILE *fp=openFile(outfiles,".qopt.ERROR");
    printDouble(Q,d.nInd,nPop,fp);
    fclose(fp);
    fp=openFile(outfiles,".fopt.ERROR");
    printDouble(F,d.nSites,nPop,fp);
    fclose(fp);
    return 0;
    
  }
  //  }
  if(nThreads==1)
    lold = likelihood(Q, F, d.nSites, d.nInd, nPop,d.genos);
  else
    lold = like_tsk(Q, F, nThreads);
  fprintf(stderr,"best like=%f after %d iterations\n",lold,nit);
  
  
  // Print F and Q in files
  FILE *fp=openFile(outfiles,".qopt");
  printDouble(Q,d.nInd,nPop,fp);
  fclose(fp);
  
  gzFile fpGz=openFileGz(outfiles,".fopt.gz");
  printDoubleGz(F,d.nSites,nPop,fpGz);
  gzclose(fpGz);
  
  
  
  //deallocate memory
  
  //if we are using the line approach then cleanup
  double nop;
  if(method &&nThreads==1)
    emAccel(d,0,NULL,NULL,NULL,NULL,nop);
  if(method &&nThreads>1)
    emAccelThread(d,0,NULL,NULL,NULL,NULL,nop,0);
  
  dalloc(F_new,d.nSites);
  dalloc(Q_new,d.nInd);
  
  dalloc(d);
  
  delete [] threads;
  delete [] myPars;
  
  if(nThreads==1){
    em(NULL, NULL, 0, 0, 0,NULL,NULL,NULL);
  }


  for(int j = 0; j < d.nSites; j++) 
    delete[] F[j];
  delete[] F;
  
  for(int i = 0; i < d.nInd; i++)
    delete [] Q[i];
  delete[] Q;
  for(int i=0;1&&i<dumpedFiles.size();i++){
    //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  


  // print to log file

  fprintf(flog, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(flog,"best like=%f after %d iterations\n",lold,nit);
  fclose(flog); 

  fclose(ffilter);
  return 0;

 

}
