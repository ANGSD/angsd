#include <math.h>
#define PI 3.141592654
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include <assert.h>
int static z_rndu=137;

int simpleRand = 0;


#define LENS 100000

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

double uniform(){
  double sampled;
  if(simpleRand)
    sampled = uniform_thorfinn();    
  else
    sampled = uniform_rasmus();
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

int print_ind_site(double errate, double meandepth, int genotype[2],gzFile glffile){

  int i, b, numreads;
  double like[10];

  numreads = Poisson(meandepth);
  //fprintf(stderr,"%d (%d,%d)\n",numreads,genotype[0],genotype[1]);

  for (i=0; i<10; i++)
    like[i] = -0.0;
  
  for (i=0; i<numreads; i++){
    b = pick_a_base(errate,genotype);
    for(int j=0;0&&j<10;j++)
      fprintf(stderr," %f ",like[j]);
    //    fprintf(stderr," pre\n ");
    calclike(b, errate, like);
    for(int j=0;0&&j<10;j++)
      fprintf(stderr," %f ",like[j]);
    //fprintf(stderr," post\n\n");
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
  gzwrite(glffile,like,sizeof(double)*10);
  return numreads;
}

//nsam=2*nind
void test ( int nsam, int segsites, char **list,int *positInt,gzFile gz,double errate,double meandepth, int regLen,FILE *vSitesFP,double *depths,int dLen){
  fprintf(stderr,"segsites=%d\n",segsites);
  
  int p =0;
  char ** res=malloc(nsam/2*sizeof(char*));
  if(positInt[0]==0)
    positInt[0]=1;
  for(int i=1;i<segsites;i++){
    if(positInt[i]<=positInt[i-1])
      positInt[i] = positInt[i-1]+1;
  }
  //write the positions that are variable to file. First element in length of posit.
  fwrite(&segsites,sizeof(int),1,vSitesFP);

  fwrite(positInt,sizeof(int),segsites,vSitesFP);
  fflush(vSitesFP);
  //make genotypes
  for(int h=0;h<nsam;h+=2){
    res[p] = malloc(segsites);
    for(int s =0;s<segsites;s++){
      res[p][s] =gv(list[h][s]) +gv(list[h+1][s]); 
    }
    p++;
  }

  if(regLen==0){
    //only generate likes for truely variable sites
    for(int s=0;s<segsites;s++){
      for(int i=0;i<nsam/2;i++){
	int genotypes[2] = {0,0};
	if(res[i][s]>=1)
	  genotypes[1] = 1;
	if(res[i][s]==2)
	  genotypes[0] = 1;
	if(depths==NULL)
	  print_ind_site(errate,meandepth,genotypes,gz);
	else
	  print_ind_site(errate,depths[i],genotypes,gz);
      }

    }
  }else{

    int s =0;
    //shifted with one, to match the positions. This shouldbn't matte for correctness
    for(int i=1;i<=regLen;i++) {
      if(s<segsites&&positInt[s]==i) {
	for(int i=0;i<nsam/2;i++){
	  int genotypes[2] = {0,0};
	  if(res[i][s]>=1)
	    genotypes[1] = 1;
	  if(res[i][s]==2)
	    genotypes[0] = 1;
	  if(depths==NULL)
	    print_ind_site(errate,meandepth,genotypes,gz);
	  else
	    print_ind_site(errate,depths[i],genotypes,gz);
	}
	s++;
      }else{
	for(int i=0;i<nsam/2;i++){
	  int genotypes[2] = {0,0};
	  if(depths==NULL)
	    print_ind_site(errate,meandepth,genotypes,gz);
	  else
	    print_ind_site(errate,depths[i],genotypes,gz);
	}
      }
      
    }
  }
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
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
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
  int singleOut = 0;
  FILE *in = NULL;
  argv++;
  char *depthFile = NULL;
  double *depths = NULL;
  int nind = 0;
  char **orig = argv;
  while(*argv){
    if(strcmp(*argv,"-in")==0) inS = *++argv;
    else if(strcmp(*argv,"-out")==0) prefix=*++argv; 
    else if(strcmp(*argv,"-err")==0) errate=atof(*++argv); 
    else if(strcmp(*argv,"-depth")==0) meanDepth=atof(*++argv); 
    else if(strcmp(*argv,"-depthFile")==0) depthFile=*++argv; 
    else if(strcmp(*argv,"-singleOut")==0) singleOut=atoi(*++argv); 
    else if(strcmp(*argv,"-regLen")==0) regLen=atoi(*++argv); 
    else if(strcmp(*argv,"-nind")==0) nind=atoi(*++argv); 
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      return 0;
    }
    ++argv;
  }

  if(inS==NULL||prefix==NULL){
    fprintf(stderr,"Probs with args, supply -in -out\n");
    fprintf(stderr,"also -err -depth -depthFile -singleOut -regLen -nind\n");
    return 0;
  }

  fprintf(stderr,"-in=%s -out=%s -err=%f -depth=%f -singleOut=%d -regLen=%d \n",inS,prefix,errate,meanDepth,singleOut,regLen);
  //print args file
  FILE *argFP=openFile(prefix,".argg");
  fprintf(argFP,"-in=%s -out=%s -err=%f -depth=%f -singleOut=%d -regLen=%d \n",inS,prefix,errate,meanDepth,singleOut,regLen);
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
  int  nsegsub, segsub( int nsam, int segsites, char **list ) ;
	

  

  /* read in first two lines of output  (parameters and seed) */
  //  pfin = stdin ;
  pfin = in;
  if(NULL==fgets( line, 1000, pfin))
    fprintf(stderr,"Problem reading from file:\n");
  sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
  double ttt=0;
  
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

  gzFile gz=Z_NULL;
  FILE *vPosFP = NULL;
  //  infoFp = openFile(prefix,".info");

  if(singleOut==1){
    gz=openFileGz(prefix,".glf.gz","w");
    vPosFP=openFile(prefix,".vPos");
  }
  FILE *pgEst = openFile(prefix,".pgEstH");

  double res[3]={0,0,0};//phi_t,phi_w,D',Dt
  while( howmany-count++ ) {

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
    gz = openFileGzI(prefix,".glf",count,"w");
    vPosFP = openFileI(prefix,".vPos",count);
  }
  //  if(1||count==58)
  test(nsam, segsites, list,positInt,gz,errate,meanDepth,regLen,vPosFP,depths,nind) ;
  if(singleOut==0){
    gzclose(gz);
    fclose(vPosFP);
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
  if(singleOut==1)
    gzclose(gz);
  fclose(pgEst);
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
