#include <math.h>
#define PI 3.141592654
#include <stdio.h>
#include <stdlib.h>
#include <htslib/bgzf.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include <assert.h>
#include <htslib/kstring.h>
int static z_rndu=137;

int simpleRand = 2;
int pileup;
int psmcOut = 0;

void mbgzf_write(BGZF* fp,const void *data,size_t length){
  if(length!=bgzf_write(fp,data,length)){
    fprintf(stderr,"\t-> Problem writing data\n");
    exit(0);
  }

}


#define LENS 100000
#define NBASE_PER_LINE 57
#define PSMC_FOR_WHICH 0
/*U(0,1): AS 183: Appl. Stat. 31:188-190
Wichmann BA & Hill ID.  1982.  An efficient and portable
pseudo-random number generator.  Appl. Stat. 31:188-190
x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.
Suggested to me by Z. Yang who also provided me with the source code used here. */
double uniform_rasmus()
{
  static int x_rndu=11, y_rndu=23;
  double r;

  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  if (x_rndu<0) x_rndu+=30269;
  if (y_rndu<0) y_rndu+=30307;
  if (z_rndu<0) z_rndu+=30323;
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}
 
//I think this one is fine, tsk 5/7 july
double uniform_thorfinn(){
  return((double)rand() / (double)RAND_MAX);
}

// translate from int code to letters, for like computation, maybe it's useless now
char int_to_base(int b) {
  if (b==0) return 'A';
  else if (b==1) return 'C';
  else if (b==2) return 'G';
  else return 'T';
}

double uniform(){
  double sampled;
  if(simpleRand==1)
    sampled = uniform_thorfinn();    
  else if(simpleRand==0)
    sampled = uniform_rasmus();
  else if(simpleRand==2)
    sampled = drand48();
  else{
    fprintf(stderr,"\t-> unknown rand method:%d\n",simpleRand);
    exit(0);
  }
  //  fprintf(stdout,"%f\n",sampled);
  return sampled;
}

double Poisson(double xm)
{
  double lgamma(double xx);
  static double sq,alxm,g,oldm=(-1.0);
  double em, t, y;
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em=-1;
    t=1.0;
    do {
      ++em;
      t *=uniform();
    } while (t>g);
  } 
  else {
    if (xm!=oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-lgamma(xm+1.0);
    }
    do {
      do {
	y=tan(PI*uniform());
	em=sq*y+xm;
      } while (em< 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-lgamma(em+1.0)-g);
    } while (uniform()>t);
  }
  return em;
}

int gv(char c){
  if(c== '0')
    return 0;
  else if(c=='1')
    return 1;
  assert(1==0);
  return -1;
}



int basepick_with_errors(double errate, int inbase){
  int outbase;
  
  if (uniform()<errate){ 
    while ((outbase=(floor(4*uniform()))) == inbase);
    return outbase;
  }
  else return inbase;
}

int pick_a_base(double errate, int genotype[2])	{
  if (uniform()<0.5) 
    return basepick_with_errors(errate, genotype[0]);
  else 
    return basepick_with_errors(errate, genotype[1]);
}


int offsets[4][10]={
  {0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
  {4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
  {7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
  {9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};




void calclike(int base, double errate, double *like)	{
  for(int i=0;0&&i<10;i++)
    fprintf(stderr," %f ",like[i]);
  
  //fprintf(stderr,"base=%d\n",base);
  /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
  static double homTrue,het,homFalse;
  static int preCalc=0;
  if(preCalc==0){
    // fprintf(stderr,"precalc\n");
    homTrue = log(1.0-errate);
    het = log((1.0-errate)/2.0+errate/6.0);
    homFalse = log(errate/3.0);
    preCalc = 1;
  }
  //fprintf(stderr,"offs=%d homeTrue=%f\n",offsets[base][0],homTrue);
  like[offsets[base][0]] += homTrue; //homozygotic hit
  
  for(int o=1;o<4;o++)//heterozygotic hit{
    like[offsets[base][o]] += het;
  
  for(int o=4;o<10;o++)//non hit
    like[offsets[base][o]] += homFalse;
  
}
int sim_invar_site(double errate,int nsam,double *depths,double meandepth,BGZF* glffile,kstring_t *resultfile,int chr, int pos,int onlyPoly){
  //int print_ind_site(double errate, double meandepth, int genotype[2],gzFile glffile,kstring_t *resultfile){
  int genotype[2];
  int i,b, numreads;
  int nInd = (int)(nsam/2);
  int numreadsInd[nInd];
  int totalReads=0;
  int totalBases[4];
  for(int j=0;j<4;j++)
    totalBases[j]=0;

  
  for(int ii=0;ii<nInd;ii++){
    if(depths!=NULL)
      meandepth=depths[ii];
    
    // int tmp =  Poisson(4.0);
      numreadsInd[ii] = Poisson(meandepth);
      totalReads += numreadsInd[ii];
  }

  int bases[totalReads];
  int Qscores[totalReads];
  double newError = errate;
  if(newError<0.00001)
      newError=0.00001;
  for (i=0; i<totalReads; i++){
    b = basepick_with_errors(errate, 0);
    bases[i] = b; 
   
    Qscores[i] = -10 * log10(newError) + 33; // use 0.00001 instead of 0
    
    totalBases[b]++;

  }

  //do not print  sites with less than onlyPoly alternative bases accross samples
  if(onlyPoly){
    if(totalReads - totalBases[0] < onlyPoly)
      return numreads;
    
  }


 


  if(pileup==1){
    ksprintf(resultfile,"%d\t%d\tN",chr,pos);
    int counter=0;
    for(int ii=0;ii<nInd;ii++){
      kputc('\t',resultfile);
      char ch = (char) Qscores[counter];
      ksprintf(resultfile,"%d\t",numreadsInd[ii]);
      for(int iii=0;iii<numreadsInd[ii];iii++)
	kputc(int_to_base(bases[counter + iii]),resultfile);
      kputc('\t',resultfile);
      for (int iii=0; iii<numreadsInd[ii]; iii++){
	kputc(ch,resultfile);
	//	kputc('H',resultfile);
      }
      counter += numreadsInd[ii];
    }
    kputc('\n',resultfile);
    return numreads;
  }
  
  double like[10*nInd];
  for(int ii=0;ii<nInd;ii++){
    for (i=0; i<10; i++)
      like[i] = -0.0;
    
    for (i=0; i<numreads; i++){
      b = pick_a_base(errate,genotype);
    
    
      calclike(b, errate, like);
    }

    //rescale to likeratios
    double mx = like[0];
    for(int i=1;i<10;i++)
      if(like[i]>mx)
	mx=like[i];
    
    for(int i=0;i<10;i++)
      like[i] -= mx;
    if(pileup==0)  
      mbgzf_write(glffile,like,sizeof(double)*10);
    
  }
  return numreads;
}

/*
  /opt/samtools-0.1.19/bcftools/vcfutils.pl

my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
			 GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');


*/



int print_ind_site(double errate, double meandepth, int genotype[2],BGZF* glffile,kstring_t *resultfile){

  int i, b, numreads;
  numreads = Poisson(meandepth);
  char res[numreads]; // char alleles for all the reads


  double newError = errate;
  if(newError<0.00001)
    newError=0.00001;

  if(pileup==1){
    for (i=0; i<numreads; i++){
      b = pick_a_base(errate,genotype);
      res[i] = int_to_base(b); 
    }
   
    int Q = -10 * log10(newError) + 33; //use 0.00001 instead of 0
    ///  char ch=static_cast<char>(Q);
    char ch = (char) Q;
    //fprintf(stdout,"Q=%c,%d,%f,%f\n",Q,Q,log10(errate),errate);
    ksprintf(resultfile,"%d\t",numreads);
    for(int iii=0;iii<numreads;iii++)
      kputc(res[iii],resultfile);
    kputc('\t',resultfile);
    for (i=0; i<numreads; i++){
      kputc(ch,resultfile);
      //kputc('H',resultfile);
    }
  // gzprintf(resultfile,"\t");
    return numreads;

  }
  
  double like[10];

   //fprintf(stderr,"%d (%d,%d)\n",numreads,genotype[0],genotype[1]);

  for (i=0; i<10; i++)
    like[i] = -0.0;
  
  for (i=0; i<numreads; i++){
    b = pick_a_base(errate,genotype);
    res[i] = int_to_base(b); 
    
    calclike(b, errate, like);
  }

    
  if(1) {
    //rescale to likeratios
    double mx = like[0];
    for(int i=1;i<10;i++)
      if(like[i]>mx)
	mx=like[i];
    
    for(int i=0;i<10;i++)
      like[i] -= mx;

  }
  //exit(0);
 
  if(pileup==0)  
    mbgzf_write(glffile,like,sizeof(double)*10);
 
  return numreads;
}

//nsam=2*nind
void test ( int nsam, int segsites, char **list,int *positInt,BGZF* gz,double errate,double meandepth, int regLen,FILE *vSitesFP,double *depths,int dLen,BGZF* gzSeq,int count,int onlyPoly,BGZF* gzPsmc,int do_seq_glf, int printSNP,int simHap){
  kstring_t kpl;kpl.s=NULL;kpl.l=kpl.m=0;
  kstring_t kpsmc;kpsmc.s=NULL;kpsmc.l=kpsmc.m=0;
  int p =0;
  char ** res=malloc(nsam/(2-simHap)*sizeof(char*));
  if(positInt[0]==0)//why
    positInt[0]=1;
  for(int i=1;i<segsites;i++){
    if(positInt[i]<=positInt[i-1])
      positInt[i] = positInt[i-1]+1;
  }
  for(int i=0;0&&i<segsites;i++){
    fprintf(stdout,"pos[%d]\t%d\n",i,positInt[i]);

  }

  //write the positions that are variable to file. First element in length of posit.
  fwrite(&segsites,sizeof(int),1,vSitesFP);
  fwrite(positInt,sizeof(int),segsites,vSitesFP);
  //make genotypes
  for(int h=0;h<nsam;h+=(2-simHap)){
    res[p] = malloc(segsites);
    for(int s =0;s<segsites;s++){
      res[p][s] =gv(list[h][s]) +gv(list[h+1-simHap][s]); 
    }
    p++;
  }
  if(gzPsmc)
    ksprintf(&kpsmc,"@%d\n",count);


  if(printSNP) {
    //only generate likes for truely variable sites
    for(int s=0;s<segsites;s++){
      if(pileup)
	ksprintf(&kpl,"%d\t%d\tN\t",count,s);

      for(int i=0;i<nsam/(2-simHap);i++){
	int genotypes[2] = {0,0};
	if(res[i][s]>=1)
	  genotypes[1] = 1;
	if(res[i][s]==2)
	  genotypes[0] = 1;

	if(gzPsmc!=NULL&&i==PSMC_FOR_WHICH){
	  if(s>0 &&((s % NBASE_PER_LINE) ==0)){
	    ksprintf(&kpsmc,"\n");
	  }
	  int nder=genotypes[0]+genotypes[1];
	  if(nder==0)
	    ksprintf(&kpsmc,"A");
	  if(nder==1)
	    ksprintf(&kpsmc,"M");
	  if(nder==2)
	    ksprintf(&kpsmc,"C");
	  if(kpsmc.l>4096){
	    mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
	    kpsmc.l=0;
	  }

	}
	if(do_seq_glf){
	  if(depths==NULL)
	    print_ind_site(errate,meandepth,genotypes,gz,&kpl);
	  else
	    print_ind_site(errate,depths[i],genotypes,gz,&kpl);
	}
	if(pileup && i<nsam/(2-simHap)-1)
	  kputc('\t',&kpl);
      }

      //all samples for a site
      if(pileup){
	kputc('\n',&kpl);
	if(kpl.l>4096){
	  mbgzf_write(gzSeq,kpl.s,kpl.l);
	  kpl.l=0;
	}
      }
    }
  }else{
    //print for invariable and variable sites
    int s =0;
   
    for(int i=0;i<regLen;i++) {
      //  fprintf(stderr,"i:%d s:%d posit:%d",i,s,positInt[s]);
      if(s<segsites && positInt[s]-1==i) {
	if(pileup)
	  ksprintf(&kpl,"%d\t%d\tN\t",count,i);
	//loop over samples
	for(int n=0;n<nsam/2;n++) {
	  int genotypes[2] = {0,0};
	  if(res[n][s]>=1)
	    genotypes[1] = 1;
	  if(res[n][s]==2)
	    genotypes[0] = 1;

	  if(gzPsmc!=NULL&&n==PSMC_FOR_WHICH){
	    if(i>0&&(i % NBASE_PER_LINE) ==0)
	      ksprintf(&kpsmc,"\n");
	    
	    int nder=genotypes[0]+genotypes[1];
	    if(nder==0)
	      ksprintf(&kpsmc,"A");
	    if(nder==1)
	      ksprintf(&kpsmc,"M");
	    if(nder==2)
	      ksprintf(&kpsmc,"C");
	    if(kpsmc.l>4096){
	      mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
	      kpsmc.l=0;
	    }
	  }
	  if(do_seq_glf){
	    if(depths==NULL)
	      print_ind_site(errate,meandepth,genotypes,gz,&kpl);
	    else
	      print_ind_site(errate,depths[n],genotypes,gz,&kpl);
	  }

	  if(pileup && n<nsam/2-1)
	    kputc('\t',&kpl);
	}
	//after all samples for a site
	if(pileup){
	  kputc('\n',&kpl);
	  if(kpl.l>4096){
	    mbgzf_write(gzSeq,kpl.s,kpl.l);
	    kpl.l=0;
	  }
	}
	s++;
      }else{
	//this 'for' loop is  important, should not be be commented out.
	for(int n=0;1&&n<nsam/2;n++){
	  int genotypes[2] = {0,0};
	  if(gzPsmc!=NULL&&n==PSMC_FOR_WHICH){
	    if(i>0&& (i % NBASE_PER_LINE) ==0){
	      ksprintf(&kpsmc,"\n");
	    }
	    int nder=genotypes[0]+genotypes[1];
	    if(nder==0)
	      ksprintf(&kpsmc,"A");
	    if(nder==1)
	      ksprintf(&kpsmc,"M");
	    if(nder==2)
	      ksprintf(&kpsmc,"C");
	    if(kpsmc.l>4096){
	      mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
	      kpsmc.l=0;
	    }
	    
	  }
	  if(do_seq_glf){
	    if(depths==NULL)
	      print_ind_site(errate,meandepth,genotypes,gz,&kpl);
	    else
	      print_ind_site(errate,depths[n],genotypes,gz,&kpl);
	  }

	  if(pileup && n<nsam/2-1)
	    kputc('\t',&kpl);
	}

	if(pileup){
	  sim_invar_site(errate,nsam,depths,meandepth,gz,&kpl,count,i,onlyPoly);
	  if(kpl.l>4096){
	    mbgzf_write(gzSeq,kpl.s,kpl.l);
	    kpl.l=0;
	  }
	}
      }
    }
  }
  if(pileup){
    mbgzf_write(gzSeq,kpl.s,kpl.l);
    kpl.l=0;free(kpl.s);
  }
  if(psmcOut){
    if(kpsmc.l>0&&kpsmc.s[kpsmc.l-1]!='\n')
      ksprintf(&kpsmc,"\n");
    mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
    kpsmc.l=0;

    ksprintf(&kpsmc,"+\n");
    int up = segsites;
    if(regLen>0)
      up=regLen;
    for(int s=0;s<up;s++){
      if(s>0 &&(s % NBASE_PER_LINE) ==0){
	//fprintf(stderr," s:%d NBA:%d\n",s,NBASE_PER_LINE);
	ksprintf(&kpsmc,"\n");
      }
      ksprintf(&kpsmc,"S");
      if(kpsmc.l>4096){
	mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
	kpsmc.l=0;
      }
    }
    if(kpsmc.l>0&&kpsmc.s[kpsmc.l-1]!='\n')
      ksprintf(&kpsmc,"\n");
    mbgzf_write(gzPsmc,kpsmc.s,kpsmc.l);
    kpsmc.l=0;free(kpsmc.s);
  }
  for(int h=0;h<nsam/2;h++)
    free(res[h]);
  free(res);
}


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

double nucdiv(int, int, char **);
double tajd(int, int, double) ;
double hfay(int, int, char **);
double thetah(int, int, char **);

int maxsites = 1000 ;

void biggerlist(int nsam, unsigned nmax,char ** list ) {
        
  int i;
  maxsites = nmax  ;
  for( i=0; i<nsam; i++){
    list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
    if( list[i] == NULL ) perror( "realloc error. bigger");
  }
}                        


FILE *getFILE(const char*fname,const char* mode){
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;
  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"[%s:%s()]\t->Error opening FILE handle for file:%s exiting\n",__FILE__,__FUNCTION__,fname);
    exit(0);
  }
  return fp;
}

FILE *openFile(const char* a,const char* b){
  if(1)
    fprintf(stderr,"[%s] %s %s\n",__FUNCTION__,a,b);
  //  char *c = new char[strlen(a)+strlen(b)+1];
  char *c = malloc(strlen(a)+strlen(b)+1);
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  FILE *fp = getFILE(c,"w");
  free(c);
  return fp;
}

BGZF* getGz(const char*fname,const char* mode){
  if(1)
    fprintf(stderr,"doing %s with %s\n",fname,mode);
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  BGZF *fp=NULL;
  if(NULL==(fp=bgzf_open(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

BGZF* openFileGz(const char* a,const char* b,const char *mode){
  if(1)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  //  char *c = new char[strlen(a)+strlen(b)+1];
  char *c = malloc(strlen(a)+strlen(b)+1);
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  BGZF *fp = getGz(c,mode);
  //delete [] c;
  free(c);
  return fp;
}

BGZF* openFileGzI(const char* a,const char* b,int i,const char*mode){
  char ary[5000];
  snprintf(ary,5000,"%s%d.gz",b,i);
  return openFileGz(a,ary,mode);  
}


FILE *openFileI(const char* a,const char* b,int i){
  char ary[5000];
  snprintf(ary,5000,"%s%d",b,i);
  return openFile(a,ary);  
}


//return filesize of file
size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


double *getDepths(const char*fname,int len){

  FILE *fp = NULL;
  if((fp=fopen(fname,"r"))==NULL){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }

  char *buf = malloc(sizeof(char)*fsize(fname));
  double *ret = malloc(len*sizeof(double));
  if(fsize(fname)!=fread(buf,sizeof(char),fsize(fname),fp)){
    fprintf(stderr,"Problem reading file=%s\n",fname);
    exit(0);
  }
  int posi=0;
  
  ret[posi++] = atof(strtok(buf,"\n\t "));
  char *tok;
  while((tok=strtok(NULL,"\n\t "))) 
    ret[posi++] = atof(tok);

  return ret;
}



int main(int argc,char **argv){
  double errate = 0.01;
  double meanDepth = 4;
  const char *inS = NULL;
  const char *prefix=NULL;
  int regLen = 0;
  int printSNP = 0; // Added by FGV on 31/08/2018
  int singleOut = 0;
  FILE *in = NULL;
  argv++;
  char *depthFile = NULL;
  double *depths = NULL;
  int onlyPoly = 0;
  int nind = 0;
  char **orig = argv;
  int seed = -1;
  int Nsites=0;
  int do_seq_glf = 1;
  int simHap =0;
  pileup=0;


  while(*argv){
    if(strcasecmp(*argv,"-in")==0) inS = *++argv;
    else if(strcasecmp(*argv,"-out")==0) prefix=*++argv; 
    else if(strcasecmp(*argv,"-err")==0) errate=atof(*++argv); 
    else if(strcasecmp(*argv,"-depth")==0) meanDepth=atof(*++argv); 
    else if(strcasecmp(*argv,"-Nsites")==0) Nsites=atoi(*++argv); 
    else if(strcasecmp(*argv,"-depthFile")==0) depthFile=*++argv; 
    else if(strcasecmp(*argv,"-singleOut")==0) singleOut=atoi(*++argv);
    else if(strcasecmp(*argv,"-simHap")==0) simHap=atoi(*++argv); 
    else if(strcasecmp(*argv,"-regLen")==0) regLen=atoi(*++argv);
    else if(strcasecmp(*argv,"-printSNP")==0) printSNP=atoi(*++argv);
    else if(strcasecmp(*argv,"-onlyPoly")==0) onlyPoly=atoi(*++argv);
    else if(strcasecmp(*argv,"-pileup")==0) pileup=atoi(*++argv);
    else if(strcasecmp(*argv,"-psmc")==0) psmcOut=atoi(*++argv);
    else if(strcasecmp(*argv,"-seed")==0) seed=atoi(*++argv);
    else if(strcasecmp(*argv,"-simpleRand")==0) simpleRand=atoi(*++argv); 
    else if(strcasecmp(*argv,"-nind")==0) nind=atoi(*++argv); 
    else if(strcasecmp(*argv,"-do_seq_glf")==0) do_seq_glf=atoi(*++argv); 
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      return 0;
    }
    ++argv; 
  } 
  if(seed==-1)
    srand48(time(NULL));
  else
    srand48(seed);
  if(inS==NULL||prefix==NULL){
    fprintf(stderr,"Probs with args, supply -in -out\n");
    fprintf(stderr,"also -err -depth -depthFile -singleOut -regLen -printSNP -nind -onlyPoly -seed -pileup -Nsites -psmc -simpleRand -do_seq_glf -simHap\n");
    return 0;
  }
  // If no region length provided, then only print SNPs
  if(regLen==0)
    printSNP=1;

  fprintf(stderr,"-in %s -out %s -err %f -depth %f -Nsites %d -singleOut %d -regLen %d -printSNP %d -onlyPoly %d -pileup %d -simpleRand %d -depthFile %s -seed %d -nind %d -psmc %d -do_seq_glf %d -simHap %d\n",inS,prefix,errate,meanDepth,Nsites,singleOut,regLen,printSNP,onlyPoly,pileup,simpleRand,depthFile,seed,nind,psmcOut,do_seq_glf,simHap);
  //print args file
  FILE *argFP=openFile(prefix,".argg");
  fprintf(argFP,"-in %s -out %s -err %f -depth %f -Nsites %d -singleOut %d -regLen %d -printSNP %d -onlyPoly %d -pileup %d -simpleRand %d -depthFile %s -seed %d -nind %d -psmc %d -do_seq_glf %d -simHap %d\n",inS,prefix,errate,meanDepth,Nsites,singleOut,regLen,printSNP,onlyPoly,pileup,simpleRand,depthFile,seed,nind,psmcOut,do_seq_glf,simHap);
  for(int i=0;i<argc;i++)
    fprintf(argFP,"%s ",orig[i]);
  fclose(argFP);
  
  if(depthFile!=NULL){
    if(nind==0){
      fprintf(stderr,"\nmust supply -nind when supplying -depthFile\n");
      return 0;
    }
    depths=getDepths(depthFile,nind);
    for(int i=0;0&&i<nind;i++)
      fprintf(stderr,"%d %f\n",i,depths[i]);
  }

  in=getFILE(inS,"r");
  
  int nsam, j ,nsites, i,  howmany  ;
  char **list, **cmatrix(), allele,na, line[1001], slashline[1001]  ;
  FILE *pf, *fopen(), *pfin ;
  double *posit   ;
  int *positInt ;
  int   segsites, count  , nadv, probflag  ;
  double pi , h, th  ,prob ;
  char dum[20], astr[100] ;
  char dum2[100], dum3[100] ;
  int  nsegsub, segsub( int nsam, int segsites, char **list ) ;
	

  

  /* read in first two lines of output  (parameters and seed) */
  //  pfin = stdin ;
  pfin = in;

  int ss;
  if(NULL==fgets( line, 1000, pfin))
    fprintf(stderr,"Problem reading from file:\n");
  if(Nsites)
    ss=sscanf(line," %s %s %s  %d %d", dum,dum2,dum3,  &nsam, &howmany);
  else
    ss=sscanf(line," %s  %d %d", dum,  &nsam, &howmany);

  fprintf(stderr,"Number of samples:%d\n",nsam);
  fprintf(stderr,"Number of replications:%d\n",howmany);
  if(ss < 3+2*Nsites){
     fprintf(stderr,"Error -> wrong header for input \n ss=%d\n",ss);
    exit(0);

  }
  if(nsam>1000){
    fprintf(stderr,"Error -> Number of samples:%d , is too high. exit\n Maybe you need -Nsites 1?\n ",nsam);
    exit(0);

  }
  

  double ttt=0;
  if(simHap==1){
    if(nsam % 2 ){
      fprintf(stderr,"\nGL calculation is based on diploid samples, you need to supply an even number of haplotypes\n");
      return 0;
    }
  }
  for(int i=1;i<=nsam-1;i++){
    ttt += 1.0*1/i;
    //    fprintf(stderr,"nsam=%d\tttt=%f,i=%d\n",nsam,ttt,i); 
  }
  if(NULL==fgets( line, 1000, pfin))
    fprintf(stderr,"Problem reading from file:\n");
  if( argc > 1 ) { 
    nadv = atoi( argv[1] ) ; 
  }

  list = cmatrix(nsam,maxsites+1);
  posit = (double *)malloc( maxsites*sizeof( double ) ) ;
  positInt = (int *)malloc( maxsites*sizeof( int ) ) ;
  count=0;
  probflag = 0 ;

  BGZF* gz=NULL;
  BGZF* gzSeq = NULL;
  BGZF* gzPsmc = NULL;
  FILE *vPosFP = NULL;
  //  infoFp = openFile(prefix,".info");
  if(pileup)
    gzSeq = openFileGz(prefix,".pileup.gz","w");
  if(psmcOut)
    gzPsmc = openFileGz(prefix,".fq.gz","w");
  
  if(singleOut==1){
    if(do_seq_glf)
      gz = openFileGz(prefix,".glf.gz","w");
    vPosFP=openFile(prefix,".vPos");
  }
  FILE *pgEst = openFile(prefix,".pgEstH");

  double res[3]={0,0,0};//phi_t,phi_w,D',Dt
  fprintf(stderr,"\n");
 
  while( howmany-count++ > 0 ) {
    // fprintf(stderr,"2nSam %d count %d %d\n",howmany,count,howmany-count);

     fprintf(stderr,"count %d\r",count);
      /* read in a sample */
     do {
       if( fgets( line, 1000, pfin) == NULL ){
	 exit(0);
       }
       if( line[0] == '/' )  strcpy(slashline,line+2);
     }while ( (line[0] != 's') && (line[0] != 'p' ) ) ;
     
     if( line[0] == 'p'){
       sscanf( line, "  prob: %lf", &prob );
       probflag = 1 ;
       if( fgets( line, 1000, pfin) == NULL ){
	 exit(0);
       }
     }
     sscanf( line, "  segsites: %d", &segsites );
     if(regLen!=0&&segsites>regLen){
       fprintf(stderr,"\t-> Must specify a regionsize larger than the number of variable sites\n");
       exit(0);
     }
     if( segsites >= maxsites){
       maxsites = segsites + 10 ;
       posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
       positInt = (int *)realloc( positInt, maxsites*sizeof( int) ) ;
       biggerlist(nsam,maxsites, list) ;
     }
     if( segsites > 0) {
       if(0==fscanf(pfin," %s", astr))
	 fprintf(stderr,"Problem reading stuff:\n");
       for( i=0; i<segsites ; i++){ 
	 if(0==fscanf(pfin," %lf",posit+i))
	   fprintf(stderr,"Problem reading stuff:\n");
	 positInt[i] = posit[i]*regLen;
       }
       for( i=0; i<nsam;i++) 
	 if(0==fscanf(pfin," %s", list[i] ))
	   fprintf(stderr,"Problem reading stuff:\n");
     }
     /* analyse sample ( do stuff with segsites and list) */
     if( argc > 1 ) 
       nsegsub = segsub( nadv, segsites, list) ;
     
     if(singleOut==0){
       if(do_seq_glf)
	 gz = openFileGzI(prefix,".glf",count,"w");
       vPosFP = openFileI(prefix,".vPos",count);
     }
     
     test(nsam, segsites, list,positInt,gz,errate,meanDepth,regLen,vPosFP,depths,nind,gzSeq,count,onlyPoly,gzPsmc,do_seq_glf,printSNP,simHap) ;
     if(singleOut==0){
       if(gz!=NULL){
	 bgzf_close(gz);
	 gz=NULL;
       }
       fclose(vPosFP);
       vPosFP=NULL;
     }
     
     pi = nucdiv(nsam, segsites, list) ;
     h = hfay(nsam, segsites, list) ;
     th = thetah(nsam, segsites, list) ;
     
     //  fprintf(pgEst,"pi:\t%lf\tss:\t%d\tD:\t%lf\tthetaH:\t%lf\tH:\t%lf%s", pi,segsites, tajd(nsam,segsites,pi) , th , h,slashline  ) ;
     fprintf(pgEst,"%d\t%f\t%f\t%f\t%f\n",segsites, pi, tajd(nsam,segsites,pi) , th , h) ;
     
     res[0] += pi;
     //  fprintf(stderr,"%d %f \n",segsites,ttt);
     res[1] += segsites/ttt;
     res[2] +=  tajd(nsam,segsites,pi);
     
   }
   count--;
   if(pileup)
     bgzf_close(gzSeq);
   if(psmcOut)
     bgzf_close(gzPsmc);
   if(gz!=NULL)
     bgzf_close(gz);
   if(vPosFP!=NULL)
     fclose(vPosFP);
   fclose(pgEst);
   fclose(in);
   for(int i=0;i<nsam;i++)
     free(list[i]);
   free(list);
   free(positInt);
   free(posit);
   return 0;

}

	

/* allocates space for gametes (character strings) */
char ** cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}



double nucdiv( int nsam, int segsites, char **list) {
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	return( pi ) ;
}

/*   thetah - pi   */
double hfay( int nsam, int segsites, char **list) {
  int s, frequency( char, int, int, char**);
  double pi, p1, nd, nnm1  ;
  
  pi = 0.0 ;
  
  nd = nsam;
  nnm1 = nd/(nd-1.0) ;
  //  fprintf(stderr,"nnm1=%f\n",nnm1);
  for( s = 0; s <segsites; s++){
    p1 = frequency('1', s,nsam,list)/nd ;
    pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
  }
  return( -pi ) ;
}

/* Fay's theta_H  */
double thetah( int nsam, int segsites, char **list) {
  int s, frequency( char, int, int, char**);
  double pi, p1, nd, nnm1  ;
  
  pi = 0.0 ;
  
  nd = nsam;
  nnm1 = nd/(nd-1.0) ;
  for( s = 0; s <segsites; s++){
    p1 = frequency('1', s,nsam,list) ;
    //    fprintf(stderr,"p1=%f\n",p1);
    pi += p1*p1 ; 
    if(s==1&&0)
      exit(0);
  }
  //  fprintf(stderr,"pi=%f\n",pi);
  return( pi*2.0/( nd*(nd-1.0) )  ) ;
}


int frequency( char allele,int site,int nsam,  char **list) {
  int i, count=0;
  for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
  return( count);
}

int segsub( int nsub, int segsites, char **list ) {
  int i, count = 0 , c1 ;
  int frequency( char, int, int, char**) ;
  
  for(i=0; i < segsites ; i++){
    c1 = frequency('1',i,nsub, list);
    if( ( c1 > 0 ) && ( c1 <nsub )  ) count++;
  }
  return( count ) ;
}
	

/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/


#include <stdio.h>
#include <math.h>



double a1f(int);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);


double tajd(int nsam, int segsites, double sumk){

double  a1, a2, b1, b2, c1, c2, e1, e2; 
 
if( segsites == 0 ) return( 0.0) ;

a1 = a1f(nsam);
a2 = a2f(nsam);
b1 = b1f(nsam);
b2 = b2f(nsam);
c1 = c1f(a1, b1);
c2 = c2f(nsam, a1, a2, b2);
e1 = e1f(a1, c1);
e2 = e2f(a1, a2, c2);
 if(0){
   fprintf(stderr,"e1=%f\te2=%f\n",e1,e2);
   fprintf(stderr,"sqrt=%f top=%f\n",sqrt((e1*segsites) + ((e2*segsites)*(segsites-1))),sumk - (segsites/a1));
   fprintf(stderr,"sumk=%f\t segs/a1=%f\n",sumk, (segsites/a1));
   exit(0);
 }
return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites-1))) ) ;

}

double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam){
  double b1;
  b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
  return (b1);
}


double b2f(int nsam){
  double b2;
  b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
  return (b2);
}


double e1f(double a1, double c1) {
  double e1;
  e1 = c1/a1;
  return (e1);
}

double e2f(double a1, double a2, double c2){ 
  double e2;
  e2 = c2/((a1*a1)+a2);
  return (e2);
}


double c1f(double a1, double b1) {
  double c1;
  c1 = b1 - (1/a1);
  return (c1);
}


double c2f(int nsam, double a1, double a2, double b2) {
  double c2;
  c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
  return (c2);
}
