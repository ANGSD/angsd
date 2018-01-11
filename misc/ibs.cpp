// g++ -o ibs ibs.cpp -lz -O3
 
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

int model;
int seed;
FILE *flog;
typedef struct
{
  int nInd;
  int theInd;
  size_t nSites;
  size_t totalSites;
  int *keepSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
  double pi[10];
} argu;


typedef struct
{
  int theInd1;
  int theInd2;
  int nInd;
  size_t nSites;
  size_t totalSites;
  int *keepSites;
  double lres;//this is the likelihood for a block of data. total likelihood is sum of lres.
  double pi[100];
} argu2;


void info(){
 fprintf(stderr,"Needed arguments:\n");
 fprintf(stderr,"\t-glf/-f\t\tinput GLF filename:\n");
 fprintf(stderr,"Optional arguments:\n");
 fprintf(stderr,"\t-outFileName/-o\toutput filename(prefix):\n");
 fprintf(stderr,"\t-nInd/-n\tnubmer of individuals in GLF file:\n");
 fprintf(stderr,"\t-ind1/i1\tindividuals 1:\n");
 fprintf(stderr,"\t-ind2/i2\tindividuals 2:\n");
 fprintf(stderr,"\t-allpairs/-a\tanalyse all pairs:\n");
 fprintf(stderr,"\t-maxSites/-m\tmaximum sites to analyze:\n");
 fprintf(stderr,"\t-model\t\tibs model\n\t\t single:\t0 all 10 genotypes, 1 HO/HE\n\t\t pair:\t\t0 all 10 genotypes, 1 HO/HE, 2 ABCD\n");
 fprintf(stderr,"\t-seed\t\tibs -seed 0 use seed for random start\n");
 // fprintf(stderr,"\t-\t\t:\n");


}


double addProtectN(double a[],int len){
  //function does: log(sum(exp(a))) while protecting for underflow
  double maxVal = a[0];

  for(int i=1;i<10;i++)
    if(maxVal<a[i])
      maxVal=a[i];

  double sumVal = 0;
  for(int i=1;i<10;i++)
    sumVal += exp(a[i]-maxVal);

  return log(sumVal) + maxVal;
}

std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  fprintf(stderr,"\t-> Dumping file: %s\n",c);

  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

size_t readGLF(const char* fname,double * &gls,int nInd,size_t maxSites){
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  size_t sizeMax = 100000*nInd;
  size_t nSites = 0;

  //  double *buf = new double[sizeMax];
  double *tmp_gls = new double[sizeMax];

  int temp=0;
  while( (temp = gzread(fp,tmp_gls,sizeof(double)*sizeMax)) ){

    //    fprintf(stdout,"%d\n",temp);
    nSites += temp/sizeof(double)/10/nInd;
    if( nSites > maxSites ){
      nSites = maxSites;
      fprintf(stdout,"max sites reached with %lu sites\n",nSites);
      fflush(stdout);
      break;
    }
  }

  gzclose(fp);
  delete[] tmp_gls;

  gls = new double[nSites*10*nInd];
  fp=gzopen(fname,"r");

  size_t nSites2=0;
  size_t offset=0;
  while(( temp = gzread(fp,gls+offset,sizeof(double)*sizeMax) )){
    nSites2 += temp/sizeof(double)/10/nInd; 
    offset+=sizeMax;
    if(nSites2>=maxSites)
      break;
  }
  
  //  size_t nSites2 = gzread(fp,gls,sizeof(double)*nSites*10*nInd);
  gzclose(fp);
  if(nSites!=nSites2){
    fprintf(stdout,"reading failed nSites%lu != nSites2 %lu\n",nSites,nSites2);   
  }
  return nSites2;

}



void runEM(double *gl,argu *pars){

  float sum;
  
  int maxIter=30;
  int it;
  double W[10];
  double Wtmp[10];
  double p[10];
  double pnew[10];
  for(int w=0;w<10;w++)
    p[w] = 0.1;
  if(seed){
    srand48(seed);
    double ssum=0;
    for(int w=0;w<10;w++){
	p[w] = drand48();
	ssum += p[w];
    }
    for(int w=0;w<10;w++){
      p[w] = p[w]/ssum;
      pnew[w] = p[w];
    }
    
  }


  double tol = 0.00000001;


  for(it=0;it<maxIter;it++){
    sum=0;
    for(int w=0;w<10;w++)
      W[w]=0;
 
    for(size_t i=0;i<pars->totalSites;i++){

      if(pars->keepSites[i]==0)
	continue;
      double sum=0;
      for(int w=0;w<10;w++){

	//Wtmp[w]=exp(gl[pars->theInd*10+w+pars->nInd*10*i])*p[w];
	Wtmp[w]=gl[pars->theInd*10+w+pars->nInd*10*i]*p[w];
	sum += Wtmp[w];
      }
      for(int w=0;w<10;w++)
	W[w] += Wtmp[w]/sum;
    }

    if(model==0){ //full model. all 10 genotypes
      for(int w=0;w<10;w++)
	pnew[w] = W[w]/pars->nSites;

    }
    else if(model==1){// 2 genotype model (HO/HE)
      double HO=W[0]+W[4]+W[7]+W[9];
      int w=0;
      for(int a1=0;a1<4;a1++)
	for(int a2=a1;a2<4;a2++){
	  if(a1==a2)
	    pnew[w] = HO/pars->nSites/4;
	  else
	    pnew[w] =(1-HO/pars->nSites)/6;
	  w++;
	}
    }


    double diff=0;
    for(int w=0;w<10;w++){
      double d = fabs(p[w]-pnew[w]);
      if(d>diff)
	diff=d;
      p[w] = pnew[w];
    }
    if(diff<tol){
      fprintf(stdout,"\ttol reached in %d interactions\tdiff=%f\tmodel=%d\n",it,diff,model);
      fprintf(flog,"\ttol reached in %d interactions\tdiff=%f\tmodel=%d\n",it,diff,model);
      break;
    }

  }
  for(int w=0;w<10;w++)
    pars->pi[w]=p[w];

  double like=0;
  for(size_t i=0;i<pars->totalSites;i++){
    double sum=0;
    for(int w=0;w<10;w++){
      sum += gl[pars->theInd*10+w+pars->nInd*10*i]*p[w];
    }
     like+=log(sum);

  }
    
  pars->lres = like;

}


void model1_2D(double *W,double  pnew[],size_t Nsites){
      double HOHO=0;
      double HOHE=0;
      double HEHO=0;
      double HEHE=0;
      double HOaHO=0;
      int w=0;
      for(int a11=0;a11<4;a11++)
	for(int a12=a11;a12<4;a12++)
	  for(int a21=0;a21<4;a21++)
	    for(int a22=a21;a22<4;a22++){
	      if(a11!=a12){
		if(a21!=a22)
		  HEHE += W[w];
		else
		  HEHO += W[w];
	      }
	      else{
		if(a21!=a22)
		  HOHE += W[w];
		else{
		  if(a11==a22)
		    HOHO += W[w];
		  else
		    HOaHO +=W[w];
		}
	      }
	      w++;
	    }
      //    fprintf(stdout,"%f\t%f\t%f\t%f\t%f\t\n",HOHO,HOHE,HEHO,HEHE,HOaHO);
     w=0;
      for(int a11=0;a11<4;a11++)
	for(int a12=a11;a12<4;a12++)
	  for(int a21=0;a21<4;a21++)
	    for(int a22=a21;a22<4;a22++){
	      if(a11!=a12){
		if(a21!=a22)
		  pnew[w] = HEHE/Nsites/36;
		else
		  pnew[w] = HEHO/Nsites/24;
	      }
	      else{
		if(a21!=a22)
		  pnew[w] = HOHE/Nsites/24;
		else{
		  if(a11==a22)
		    pnew[w] = HOHO/Nsites/4;
		  else
		    pnew[w] = HOaHO/Nsites/12;
		}
	      }
	      w++;
	      //  fprintf(stdout,"w pnew %d\t%f\tmodel %d\n",w,pnew[w],model);
	    }



}









/*
  multi-allelic model


*/


void model2_2D(double *W,double  pnew[],size_t Nsites){
      double AAAA=0;
      double AABA=0;
      double AABB=0;
      double AABC=0;
      double ABAB=0;
      double ABAA=0;
      double ABCC=0;
      double ABAC=0;
      double ABCD=0;
      int w=0;
      for(int a11=0;a11<4;a11++)
	for(int a12=a11;a12<4;a12++)
	  for(int a21=0;a21<4;a21++)
	    for(int a22=a21;a22<4;a22++){
	      if(a11!=a12){ //AB__
		if( (a21==a11 && a22==a12) ||  (a21==a12 && a22==a11)) //ABAB
		  ABAB += W[w];
		else if( (a21==a11 && a22==a11 ) ||  (a21==a12 && a22==a12 ) ) // ABAA
		  ABAA += W[w];
		else if( a21==a22 && a22!=a11 && a22!=a12 ) // ABCC
		  ABCC += W[w];
		else if( a21!=a22 && (a21==a12 || a21==a11 || a22==a12 || a22==a11 ) ) // ABAC  ( (a21==a11 && a22!=a12 ) ||  (a21==a12 && a22!=a11 ) ) 
		  ABAC += W[w];
		else  //ABCD
		  ABCD +=W[w];
	      }
	      else{ //AA__
		if(a21==a11 && a22==a12) //AAAA
		  AAAA += W[w];
		else if( (a21==a11 && a22!=a11 ) ||  (a22==a11 && a21!=a11 ) ) // AABA
		  AABA += W[w];
		else if( a21==a22 && a22!=a11 ) // AABB
		  AABB += W[w];
		else if( a21!=a22 && a22!=a12 ) // AABC
		  AABC += W[w];
		else  
		  exit(0);
	      }
	      w++;
	    }
 
      w=0;
      for(int a11=0;a11<4;a11++)
	for(int a12=a11;a12<4;a12++)
	  for(int a21=0;a21<4;a21++)
	    for(int a22=a21;a22<4;a22++){
	      if(a11!=a12){ //AB__
		if( (a21==a11 && a22==a12) ||  (a21==a12 && a22==a11)) //ABAB
		  pnew[w] = ABAB/Nsites/6;
		else if( (a21==a11 && a22==a11 ) ||  (a21==a12 && a22==a12 ) ) // ABAA
		  pnew[w] = ABAA/Nsites/6/2;
		else if( a21==a22 && a22!=a11 && a22!=a12 ) // ABCC
		  pnew[w] = ABCC/Nsites/6/2;
		else if( (a21==a11 && a22!=a12 ) ||  (a21==a12 && a22!=a11 ) ) // ABAC
		  pnew[w] = ABAC/Nsites/6/4;
		else  //ABCD
		  pnew[w] = ABCD/Nsites/6/1;
	      }
	      else{ //AA__
		if(a21==a11 && a22==a12) //AAAA
		  pnew[w] = AAAA /Nsites/4;
		else if( (a21==a11 && a22!=a11 ) ||  (a22==a11 && a21!=a11 ) ) // AABA
		  pnew[w] = AABA /Nsites/4/3;
		else if( a21==a22 && a22!=a11 ) // AABB
		  pnew[w] = AABB /Nsites/4/3;
		else if( a21!=a22 && a22!=a12 ) // AABC
		  pnew[w] = AABC/Nsites/4/3;
		else  //ABCD
		  exit(0);
	      }
	      w++;
	    }
}






void runEM2D(double *gl,argu2 *pars){


  float sum;
  
  int maxIter=1000;
  int it;
  double W[100];
  double Wtmp[100];
  double p[100];
  double pnew[100];
  for(int w=0;w<100;w++){
    p[w] = 0.01;
    pnew[w] = 0.01;
  }
  if(seed){
    srand48(seed);
    double ssum=0;
    for(int w=0;w<100;w++){
	p[w] = drand48();
	ssum += p[w];
    }
    for(int w=0;w<100;w++){
      p[w] = p[w]/ssum;
      pnew[w] = p[w];
    }
    //    fprintf(stderr,"using seed p[0]=%f \n",p[0]);
  }
  //  fprintf(stderr,"using seed %d, p[0]=%f \n",seed,p[0]);
  double tol = 0.00000001;


  for(it=0;it<maxIter;it++){
    sum=0;
    for(int w=0;w<100;w++)
      W[w]=0;
    for(size_t i=0;i<pars->totalSites;i++){

      if(pars->keepSites[i]==0)
	continue;
      double sum=0;
      double gl1;
      for(int w1=0;w1<10;w1++){
	for(int w2=0;w2<10;w2++){
	  int w=w1 + 10 * w2;

	  // 	Wtmp[w]=exp(gl[pars->theInd1*10+w1+pars->nInd*10*i]+gl[pars->theInd2*10+w2+pars->nInd*10*i])*p[w];
	Wtmp[w]=gl[pars->theInd1*10+w1+pars->nInd*10*i]*gl[pars->theInd2*10+w2+pars->nInd*10*i]*p[w];
	sum += Wtmp[w];
	//	fprintf(stdout,"%f\t%f\t%d\t%d\t%d\t\n",gl[pars->theInd1*10+w1+pars->nInd*10*i],gl[pars->theInd2*10+w2+pars->nInd*10*i],w,w1,w2);
	}
      }

      for(int w=0;w<100;w++)
	W[w] += Wtmp[w]/sum;
    }
    if(model==0){ //full model. all 10 genotypes
      for(int w=0;w<100;w++)
	pnew[w] = W[w]/pars->nSites;

    }
    else if(model==1){// 2 genotype model (HO/HE), diallelic model
      model1_2D(W,pnew,pars->nSites);
    }
    else if(model==2){// 2 genotype model (HO/HE), multiAllelic model
      model2_2D(W,pnew,pars->nSites);
    }
    
    double diff=0;
    for(int w=0;w<100;w++){
      double d = fabs(p[w]-pnew[w]);
      if(d>diff)
	diff=d;
      p[w] = pnew[w];
    }
    if(diff<tol){
      fprintf(stdout,"\ttol reached in %d interations\tdiff=%f\tmodel=%d\n",it,diff,model);
      fprintf(flog,"\ttol reached in %d interations\tdiff=%f\tmodel=%d\n",it,diff,model);
      break;
    }  
    //        fprintf(stdout,"\n");
  }
  
    
  
  for(int w=0;w<100;w++)
    pars->pi[w]=p[w];
  
  double like=0;
  for(size_t i=0;i<pars->totalSites;i++){
    double sum=0;
    if(pars->keepSites[i]==0)
      continue;

      for(int w1=0;w1<10;w1++){
	for(int w2=0;w2<10;w2++){
	  int w=w1 + 10 * w2;
	  sum += gl[pars->theInd1*10+w1+pars->nInd*10*i]*gl[pars->theInd2*10+w2+pars->nInd*10*i]*p[w];
	}
      }
      //    for(int w=0;w<100;w++){
      //      sum += gl[pars->theInd1*10+w+pars->nInd*10*i]*gl[pars->theInd2*10+w+pars->nInd*10*i]*p[w];
      //    }
    like+=log(sum);
    
    }
  
  pars->lres = like;
    
}


int main(int argc, char **argv){
  if(argc==1){// if no arguments, print info on program
    info();
    return 0;
  }

  const char* likeFileName = NULL;
  const char* outFileName = NULL;
  seed =0;
  int nInd=1;
  int p1=-1;
  int p2=-1;
  size_t maxSites=10000000;
  int all=0;
  model=0;
  seed=0;
  // reading arguments
  argv++;
  while(*argv){
    if(strcmp(*argv,"-glf")==0 || strcmp(*argv,"-f")==0) likeFileName=*++argv; 
    else if(strcmp(*argv,"-outFileName")==0 || strcmp(*argv,"-o")==0) outFileName=*++argv; 
    else if(strcmp(*argv,"-nInd")==0 || strcmp(*argv,"-n")==0) nInd=atoi(*++argv);
    else if(strcmp(*argv,"-ind1")==0||strcmp(*argv,"-i1")==0) p1=atoi(*++argv);
    else if(strcmp(*argv,"-ind2")==0||strcmp(*argv,"-i2")==0) p2=atoi(*++argv);
    else if(strcmp(*argv,"-allpairs")==0||strcmp(*argv,"-a")==0) all=atoi(*++argv);
    else if(strcmp(*argv,"-maxSites")==0||strcmp(*argv,"-m")==0) maxSites=atoi(*++argv);
    else if(strcmp(*argv,"-model")==0) model=atoi(*++argv);
    else if(strcmp(*argv,"-seed")==0) seed=atoi(*++argv);
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }

  if(maxSites<10000){
    fprintf(stderr,"maxSites cannot be smaller than 10000\n");
    exit(0);
  }

  fprintf(stderr,"using seed %d\n",seed);
  //
  if(likeFileName==NULL){
      fprintf(stderr,"Please supply glf file: -f");
      info();
  }
  if(outFileName==NULL){
    fprintf(stderr,"Will use glf file name as prefix for output\n");
    outFileName=likeFileName;
  }

  FILE *fibs;
  FILE *fibspair;
  flog=openFile(outFileName,".log");
  if(p2==-1 && all!=1)
    fibs=openFile(outFileName,".ibs");
  if((p1>=0 && p2>=0) || all)
    fibspair=openFile(outFileName,".ibspair");

  double *genoLike=NULL;
  fprintf(stderr,"reading data\n");
  size_t totalSites=readGLF(likeFileName,genoLike,nInd,maxSites);
  fprintf(stdout,"read %lu sites\n",totalSites);
  fprintf(flog,"read %lu sites\n",totalSites);
  

  for(size_t s=0;s<totalSites;s++)
    for(int i=0;i<nInd*10;i++){
      genoLike[s*nInd*10+i]=exp(genoLike[s*nInd*10+i]);
    }

 
  argu * myPars= new argu;
  myPars->totalSites = totalSites;
  myPars->keepSites =  new int[totalSites];
  myPars->nInd = nInd;

  fprintf(stdout,"\tusing model %d \n",model);
  fprintf(flog,"\tusing model %d \n",model);


  if(p2==-1 && all!=1){ // single ind
    //header of ibs files
    if(model==0)
      fprintf(fibs,"ind\tnSites\tLlike\tnAA\tnAC\tnAG\tnAT\tnCC\tnCG\tnCT\tnGG\tnGT\tnTT\n");
    else if(model==1)
      fprintf(fibs,"ind\tnSites\tLlike\tnHO\tnHE\n");
     
    
    for(int theInd=0;theInd<nInd;theInd++){
      fflush(fibs);
      if(p1!=-1 && p1!=theInd)
	continue;
      fprintf(stdout,"analysing individual %d\n",theInd);
       fprintf(flog,"analysing individual %d\n",theInd);
      myPars->theInd = theInd;
      
      size_t nSites=0;
      for(size_t i=0;i<totalSites;i++){
	myPars->keepSites[i] = 0;
	for(int w=0;w<10;w++)
	  if(genoLike[myPars->theInd*10+w+myPars->nInd*10*i]<0.9999){
	    myPars->keepSites[i] = 1;
	    nSites++;
	    break;
	}
      }  
      myPars->nSites = nSites;
      fprintf(stdout,"\tusing nSites= %lu \n",nSites);
    
      runEM(genoLike,myPars);
      
      
      fprintf(fibs,"%d\t%lu\t%f",theInd,myPars->nSites,myPars->lres);

      if(model==0){
	for(int w=0;w<10;w++){
	  fprintf(fibs,"\t%f",myPars->pi[w]*myPars->nSites);
	}
      }
      else if(model==1){
	double HO=(myPars->pi[0]+myPars->pi[4]+myPars->pi[7]+myPars->pi[9])*myPars->nSites;
	fprintf(fibs,"\t%f\t%f",HO,myPars->nSites-HO);
      }
      fprintf(fibs,"\n");
    }
  }

 if(p1>=0 && p2>=0 && all==0){


   fprintf(fibspair,"ind1\tind2\tnSites\tLlike\t");
   //paste(paste0("p",paste(rep(GENO,10),rep(GENO,each=10),sep="_")),collapse="\t")
   if(model==0)
     fprintf(fibspair,"nAA_AA\tnAC_AA\tnAG_AA\tnAT_AA\tnCC_AA\tnCG_AA\tnCT_AA\tnGG_AA\tnGT_AA\tnTT_AA\tnAA_AC\tnAC_AC\tnAG_AC\tnAT_AC\tnCC_AC\tnCG_AC\tnCT_AC\tnGG_AC\tnGT_AC\tnTT_AC\tnAA_AG\tnAC_AG\tnAG_AG\tnAT_AG\tnCC_AG\tnCG_AG\tnCT_AG\tnGG_AG\tnGT_AG\tnTT_AG\tnAA_AT\tnAC_AT\tnAG_AT\tnAT_AT\tnCC_AT\tnCG_AT\tnCT_AT\tnGG_AT\tnGT_AT\tnTT_AT\tnAA_CC\tnAC_CC\tnAG_CC\tnAT_CC\tnCC_CC\tnCG_CC\tnCT_CC\tnGG_CC\tnGT_CC\tnTT_CC\tnAA_CG\tnAC_CG\tnAG_CG\tnAT_CG\tnCC_CG\tnCG_CG\tnCT_CG\tnGG_CG\tnGT_CG\tnTT_CG\tnAA_CT\tnAC_CT\tnAG_CT\tnAT_CT\tnCC_CT\tnCG_CT\tnCT_CT\tnGG_CT\tnGT_CT\tnTT_CT\tnAA_GG\tnAC_GG\tnAG_GG\tnAT_GG\tnCC_GG\tnCG_GG\tnCT_GG\tnGG_GG\tnGT_GG\tnTT_GG\tnAA_GT\tnAC_GT\tnAG_GT\tnAT_GT\tnCC_GT\tnCG_GT\tnCT_GT\tnGG_GT\tnGT_GT\tnTT_GT\tnAA_TT\tnAC_TT\tnAG_TT\tnAT_TT\tnCC_TT\tnCG_TT\tnCT_TT\tnGG_TT\tnGT_TT\tnTT_TT\n");
   else if(model==1) //    HOHO=1, HEHO=2, aHOHO=5, HOHE=11, HEHE 12 
     fprintf(fibspair,"nAA_AA\tnAB_AA\tnAA_BB\tnAA_AB\tnAB_AB\n");
   else if(model==2) //    AAAA=1,, ABAA=2, AABB=41 ,AAAB=11, ABAB=12, ABCC=72, ABAC=22, ABCD=82, AABC=51,  =, =, =, =
     fprintf(fibspair,"nAA_AA\tnAB_AA\tnAA_BB\tnAA_AB\tnAB_AB\tnAB_CC\tnAB_AC\tnAB_CD\tnAA_BC\n");


   argu2 * myPars2D= new argu2;
   myPars2D->totalSites = totalSites;
   myPars2D->keepSites =  new int[totalSites];
   myPars2D->nInd = nInd;
   
   myPars2D->theInd1 = p1;
   myPars2D->theInd2 = p2;
   fprintf(stdout,"analysing pair %d %d\n",p1,p2);
   fprintf(flog,"analysing pair %d %d\n",p1,p2);

   size_t nSites=0;
   for(size_t i=0;i<totalSites;i++){
     myPars2D->keepSites[i] = 0;
     for(int w=0;w<10;w++)
       // if(-genoLike[myPars2D->theInd1*10+w+myPars2D->nInd*10*i] > myPars2D->keepSites[i] && -genoLike[myPars2D->theInd2*10+w+myPars2D->nInd*10*i]>myPars2D->keepSites[i]){
       if(genoLike[myPars2D->theInd1*10+w+myPars2D->nInd*10*i] < 0.9999 && genoLike[myPars2D->theInd2*10+w+myPars2D->nInd*10*i] < 0.9999){
	 myPars2D->keepSites[i] = 1;
	 nSites++;
	 break;
       }
   }
   fprintf(stdout,"\tUsing nSites= %lu \n",nSites);
   myPars2D->nSites = nSites;
   
   runEM2D(genoLike,myPars2D);
   
   if(model==0){
     fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
     for(int w=0;w<100;w++){
       fprintf(fibspair,"\t%f",myPars2D->pi[w]*myPars2D->nSites);
     }
     fprintf(fibspair,"\n");
   }
   else if(model==1){
     
     //   HOHO=1, HEHO=2, aHOHO=5, HOHE=11, HEHE 12 
     fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
     fprintf(fibspair,"\t%f\t%f\t%f\t%f\t%f\n",myPars2D->pi[0]*myPars2D->nSites*4,myPars2D->pi[1]*myPars2D->nSites*24,myPars2D->pi[4]*myPars2D->nSites*12,myPars2D->pi[10]*myPars2D->nSites*24,myPars2D->pi[11]*myPars2D->nSites*36);
     

   }
   else if(model==2){ //    AAAA=1,, ABAA=2, AABB=41 ,AAAB=11, ABAB=12, ABCC=72, ABAC=22, ABCD=82, AABC=51, 
     fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
     fprintf(fibspair,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",myPars2D->pi[0]*myPars2D->nSites*4,myPars2D->pi[1]*myPars2D->nSites*12,myPars2D->pi[40]*myPars2D->nSites*12,myPars2D->pi[10]*myPars2D->nSites*12,myPars2D->pi[11]*myPars2D->nSites*6,myPars2D->pi[71]*myPars2D->nSites*12,myPars2D->pi[21]*myPars2D->nSites*24,myPars2D->pi[81]*myPars2D->nSites*6,myPars2D->pi[50]*myPars2D->nSites*12);
     

   }
   //    runEM2D(genoLike,myPars2D);
   delete[] myPars2D->keepSites;
   delete myPars2D;

 }
 else if(all){
   fprintf(fibspair,"ind1\tind2\tnSites\tLlike\t");
   //paste(paste0("p",paste(rep(GENO,10),rep(GENO,each=10),sep="_")),collapse="\t")
   if(model==0)
     fprintf(fibspair,"pAA_AA\tpAC_AA\tpAG_AA\tpAT_AA\tpCC_AA\tpCG_AA\tpCT_AA\tpGG_AA\tpGT_AA\tpTT_AA\tpAA_AC\tpAC_AC\tpAG_AC\tpAT_AC\tpCC_AC\tpCG_AC\tpCT_AC\tpGG_AC\tpGT_AC\tpTT_AC\tpAA_AG\tpAC_AG\tpAG_AG\tpAT_AG\tpCC_AG\tpCG_AG\tpCT_AG\tpGG_AG\tpGT_AG\tpTT_AG\tpAA_AT\tpAC_AT\tpAG_AT\tpAT_AT\tpCC_AT\tpCG_AT\tpCT_AT\tpGG_AT\tpGT_AT\tpTT_AT\tpAA_CC\tpAC_CC\tpAG_CC\tpAT_CC\tpCC_CC\tpCG_CC\tpCT_CC\tpGG_CC\tpGT_CC\tpTT_CC\tpAA_CG\tpAC_CG\tpAG_CG\tpAT_CG\tpCC_CG\tpCG_CG\tpCT_CG\tpGG_CG\tpGT_CG\tpTT_CG\tpAA_CT\tpAC_CT\tpAG_CT\tpAT_CT\tpCC_CT\tpCG_CT\tpCT_CT\tpGG_CT\tpGT_CT\tpTT_CT\tpAA_GG\tpAC_GG\tpAG_GG\tpAT_GG\tpCC_GG\tpCG_GG\tpCT_GG\tpGG_GG\tpGT_GG\tpTT_GG\tpAA_GT\tpAC_GT\tpAG_GT\tpAT_GT\tpCC_GT\tpCG_GT\tpCT_GT\tpGG_GT\tpGT_GT\tpTT_GT\tpAA_TT\tpAC_TT\tpAG_TT\tpAT_TT\tpCC_TT\tpCG_TT\tpCT_TT\tpGG_TT\tpGT_TT\tpTT_TT\n");
   else if(model==1) //    HOHO=1, HEHO=2, aHOHO=5, HOHE=11, HEHE 12 
     fprintf(fibspair,"pAA_AA\tpAB_AA\tpAA_BB\tpAA_AB\tpAB_AB\n");
    else if(model==2) //    AAAA=1,, ABAA=2, AABB=41 ,AAAB=11, ABAB=12, ABCC=72, ABAC=22, ABCD=82, AABC=51,  =, =, =, =
     fprintf(fibspair,"pAA_AA\tpAB_AA\tpAA_BB\tpAA_AB\tpAB_AB\tpAB_CC\tpAB_AC\tpAB_CD\tpAA_BC\n");

   
   argu2 * myPars2D= new argu2;
   myPars2D->totalSites = totalSites;
   myPars2D->keepSites =  new int[totalSites];
   myPars2D->nInd = nInd;
   
   for(int i1=0;i1<nInd-1;i1++){
     for(int i2=i1+1;i2<nInd;i2++){
       fflush(fibspair);
       p1=i1;
       p2=i2;
       myPars2D->theInd1 = p1;
       myPars2D->theInd2 = p2;
       fprintf(stdout,"analysing pair %d %d\n",p1,p2);
       fprintf(flog,"analysing pair %d %d\n",p1,p2);
       size_t nSites=0;
       for(size_t i=0;i<totalSites;i++){
	 myPars2D->keepSites[i] = 0;
	 for(int w=0;w<10;w++)
	   if(genoLike[myPars2D->theInd1*10+w+myPars2D->nInd*10*i] < 0.9999 && genoLike[myPars2D->theInd2*10+w+myPars2D->nInd*10*i] < 0.9999){
	     myPars2D->keepSites[i] = 1;
	     nSites++;
	     break;
	   }
       }  
       myPars2D->nSites = nSites;
       fprintf(stdout,"\tusing nSites= %lu \n",nSites);
       runEM2D(genoLike,myPars2D);
       if(model==0){
	 fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
	 for(int w=0;w<100;w++){
	   fprintf(fibspair,"\t%f",myPars2D->pi[w]*100);
	 }
	 fprintf(fibspair,"\n");
       }
       else if(model==1){
	 fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
	 fprintf(fibspair,"\t%f\t%f\t%f\t%f\t%f\n",myPars2D->pi[0]*myPars2D->nSites*4,myPars2D->pi[1]*myPars2D->nSites*24,myPars2D->pi[4]*myPars2D->nSites*12,myPars2D->pi[10]*myPars2D->nSites*24,myPars2D->pi[11]*myPars2D->nSites*36);
	 
       }
       else if(model==2){ //    AAAA=1,, ABAA=2, AABB=41 ,AAAB=11, ABAB=12, ABCC=72, ABAC=22, ABCD=82, AABC=51, 
	 fprintf(fibspair,"%d\t%d\t%lu\t%f",p1,p2,myPars2D->nSites,myPars2D->lres);
	 fprintf(fibspair,"\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",myPars2D->pi[0]*myPars2D->nSites*4,myPars2D->pi[1]*myPars2D->nSites*12,myPars2D->pi[40]*myPars2D->nSites*12,myPars2D->pi[10]*myPars2D->nSites*12,myPars2D->pi[11]*myPars2D->nSites*6,myPars2D->pi[71]*myPars2D->nSites*12,myPars2D->pi[21]*myPars2D->nSites*24,myPars2D->pi[81]*myPars2D->nSites*6,myPars2D->pi[50]*myPars2D->nSites*12);
       }
 
     }
   }
   //    runEM2D(genoLike,myPars2D);
   delete[] myPars2D->keepSites;
   delete myPars2D;

 }
 
 

 ////////// clean  
 fclose(flog); 
 if(p2 == -1 && all != 1)
   fclose(fibs); 
 if((p1 >= 0 && p2 >= 0) || all)
   fclose(fibspair); 
 for(int i=0;i<dumpedFiles.size();i++){
   //    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
   free(dumpedFiles[i]);
 }
 
 delete[] genoLike;

 delete[] myPars->keepSites;
 delete myPars;
 
 return 0;
 
}
