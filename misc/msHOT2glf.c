#include <math.h>
#define PI 3.141592654
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include <assert.h>
#include <htslib/kstring.h>


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>


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

gzFile getGz(const char*fname,const char* mode){
  if(1)
    fprintf(stderr,"doing %s with %s\n",fname,mode);
  int writeFile = 0;
  for(size_t i=0;i<strlen(mode);i++)
    if(mode[i]=='w')
      writeFile = 1;

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  gzFile fp=Z_NULL;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

gzFile openFileGz(const char* a,const char* b,const char *mode){
  if(1)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  //  char *c = new char[strlen(a)+strlen(b)+1];
  char *c = malloc(strlen(a)+strlen(b)+1);
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  gzFile fp = getGz(c,mode);
  //delete [] c;
  free(c);
  return fp;
}

gzFile openFileGzI(const char* a,const char* b,int i,const char*mode){
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

int static z_rndu=137;

int simpleRand = 2;
int psmcOut = 0;



#define LENS 100000
#define NBASE_PER_LINE 60

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

int print_ind_site(double errate, double meandepth, int genotype[2],gzFile glffile,kstring_t *resultfile){

  int i, b, numreads;
  numreads = Poisson(meandepth);
  char res[numreads]; // char alleles for all the reads


  double newError = errate;
  if(newError<0.00001)
    newError=0.00001;

  
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
  gzwrite(glffile,like,sizeof(double)*10); 
  return numreads;
}

//nsam=2*nind
//only one individual input
void test (int *positInt,gzFile gz,double errate,double meandepth, int regLen,gzFile gzSeq,int count,gzFile gzPsmc,int do_seq_glf){
  fprintf(stderr,"posit:%p errate:%f meandepth:%f regLen:%d count:%d do_seq_glf:%d\n",positInt,errate,meandepth,regLen,count,do_seq_glf);
  kstring_t kpl;kpl.s=NULL;kpl.l=kpl.m=0;
  kstring_t kpsmc;kpsmc.s=NULL;kpsmc.l=kpsmc.m=0;
  int p =0;

  if(gzPsmc)
    ksprintf(&kpsmc,"@%d\n",count);

  for(int i=0;i<regLen;i++) {
    //  fprintf(stderr,"i:%d s:%d posit:%d",i,s,positInt[s]);
    int genotypes[2]={0,0};
    if(positInt[i])
      genotypes[1] = 1;
	//loop over samples
    if(gzPsmc!=Z_NULL){
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
	gzwrite(gzPsmc,kpsmc.s,kpsmc.l);
	kpsmc.l=0;
      }
    }
    if(do_seq_glf){
      print_ind_site(errate,meandepth,genotypes,gz,&kpl);

    }
  }
  
  if(psmcOut){
    if(kpsmc.l>0&&kpsmc.s[kpsmc.l-1]!='\n')
      ksprintf(&kpsmc,"\n");
    gzwrite(gzPsmc,kpsmc.s,kpsmc.l);
    kpsmc.l=0;
    
    ksprintf(&kpsmc,"+\n");
    for(int s=0;s<regLen;s++){
      if(s>0 &&(s % NBASE_PER_LINE) ==0){
	//fprintf(stderr," s:%d NBA:%d\n",s,NBASE_PER_LINE);
	ksprintf(&kpsmc,"\n");
      }
      ksprintf(&kpsmc,"S");
      if(kpsmc.l>4096){
	gzwrite(gzPsmc,kpsmc.s,kpsmc.l);
	kpsmc.l=0;
      }
    }
    if(kpsmc.l>0&&kpsmc.s[kpsmc.l-1]!='\n')
      ksprintf(&kpsmc,"\n");
    gzwrite(gzPsmc,kpsmc.s,kpsmc.l);
    kpsmc.l=0;free(kpsmc.s);
  }
}

int main(int argc,char **argv){
  double errate = 0.01;
  double meanDepth = 4;
  const char *inS = NULL;
  const char *prefix=NULL;
  int singleOut = 0;
  FILE *in = NULL;
  argv++;
  int nind = 0;
  char **orig = argv;
  int seed = -1;
  int do_seq_glf = 1;
  int nsam,howmany,theta,rho, regLen;
  int j ,nsites, i;
  char line[1001], slashline[1001]  ;
  FILE *pf, *fopen(), *pfin ;
  int *positInt ;
  int   segsites;
  
  while(*argv){
    if(strcasecmp(*argv,"-in")==0) inS = *++argv;
    else if(strcasecmp(*argv,"-out")==0) prefix=*++argv; 
    else if(strcasecmp(*argv,"-err")==0) errate=atof(*++argv); 
    else if(strcasecmp(*argv,"-depth")==0) meanDepth=atof(*++argv); 
    else if(strcasecmp(*argv,"-singleOut")==0) singleOut=atoi(*++argv); 
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
    fprintf(stderr,"also -err -depth -depthFile -singleOut -regLen --nind -seed -pileup -Nsites -psmc -do_seq_glf\n");
    return 0;
  }

  fprintf(stderr,"-in=%s -out=%s -err=%f -depth=%f -singleOut=%d -regLen=%d -seed %d -nind %d -psmc %d -do_seq_glf: %d\n",inS,prefix,errate,meanDepth,singleOut,regLen,seed,nind,psmcOut,do_seq_glf);
  //print args file
  FILE *argFP=openFile(prefix,".argg");
  fprintf(argFP,"-in=%s -out=%s -err=%f -depth=%f -singleOut=%d -regLen=%d -seed %d -nind %d -psmc %d -do_seq_glf: %d\n",inS,prefix,errate,meanDepth,singleOut,regLen,seed,nind,psmcOut,do_seq_glf);
  for(int i=0;i<argc;i++)
    fprintf(argFP,"%s ",orig[i]);
  fclose(argFP);
  
  in=getFILE(inS,"r");
  
  
  /* read in first two lines of output  (parameters and seed) */
  //  pfin = stdin ;
  pfin = in;

  int ss;
  if(NULL==fgets( line, 1000, pfin))
    fprintf(stderr,"Problem reading from file:\n");

  ss=sscanf(line,"msHOT-lite %d %d -t %d -r %d %d\n", &nsam, &howmany,&theta,&rho,&regLen);
  assert(ss==5);
  fprintf(stderr,"\t-> Number of samples:%d\n",nsam);
  fprintf(stderr,"\t-> Number of replications:%d\n",howmany);
  fprintf(stderr,"\t-> theta: %d rho: %d regLen:%d\n",theta,rho,regLen);
  if(nsam!=2){
    fprintf(stderr,"Only defined for a single diploid\n");
    return 0;
  }

  positInt = (int *)malloc( regLen*sizeof( int ) ) ;


  gzFile gz=Z_NULL;
  gzFile gzSeq = Z_NULL;
  gzFile gzPsmc = Z_NULL;

  if(psmcOut)
    gzPsmc = openFileGz(prefix,".fq.gz","w");
  
  if(singleOut==1){
    if(do_seq_glf)
      gz = openFileGz(prefix,".glf.gz","w");
  }
  FILE *pgEst = openFile(prefix,".pgEstH");

  fprintf(stderr,"\n");
  for(int i=0;i<howmany;i++){
      /* read in a sample */

    while(1){
      if( fgets( line, 1000, pfin) == NULL ){
	fprintf(stderr,"\t-> Problem reading file\n");
	break;
      }
      ss=sscanf(line,"@begin %d\n",&segsites);
      if(ss==1)
	break;
    }
    fprintf(stderr,"\t-> Parsing replicate:%d which should have nseg: %d ",i,segsites);
    memset(positInt,0,sizeof(int)*regLen);
    if(singleOut==0){
      if(do_seq_glf)
	gz = openFileGzI(prefix,".glf",i,"w");
    }
    if( fgets( line, 1000, pfin) == NULL ){
      fprintf(stderr,"\t-> Problem reading file\n");
      break;
    }
    int newlen;
    int rv = sscanf(line,"%d\n",&newlen);
    assert(rv==1&&newlen==regLen);
    for(int n=0;n<segsites;n++){
      if( fgets( line, 1000, pfin) == NULL ){
	fprintf(stderr,"\t-> Problem reading file\n");
	break;
      }
      if(strcmp(line,"@end\n")==0){
	fprintf(stderr,"breaking at:%d (lost: %d snps, due to binning of positions...)\n",n,segsites-n);
	break;
      }      
      int at;
      ss=sscanf(line,"%d ",&at);
      assert(ss=1);
      positInt[at]=1;
    }

     test(positInt,gz,errate,meanDepth,regLen,gzSeq,i,gzPsmc,do_seq_glf) ;
     if(singleOut==0){
       if(gz!=Z_NULL){
	 gzclose(gz);
	 gz=Z_NULL;
       }
     }
   }
   
   if(psmcOut)
     gzclose(gzPsmc);
   if(gz!=Z_NULL)
     gzclose(gz);
   fclose(pgEst);
   fclose(in);
   free(positInt);
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
