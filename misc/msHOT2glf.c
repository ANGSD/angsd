#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <pthread.h>
#include <sys/stat.h>
#include <assert.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int block = 100;

#define LENS 100000
#define NBASE_PER_LINE 60

size_t adepth[1024];

double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
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
  if(0)
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

BGZF *getGz(const char*fname,const char* mode){
  if(0)
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

BGZF *openFileGz(const char* a,const char* b,const char *mode){
  if(0)
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

BGZF *openFileGzI(const char* a,const char* b,int i,const char*mode){
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
	y=tan(M_PI*uniform());
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

int print_ind_site(double errate, double meandepth, int genotype[2],BGZF *glffile,kstring_t *resultfile,BGZF *outfileSAF){

  int i, b, numreads;
  numreads = Poisson(meandepth);
  char res[numreads]; // char alleles for all the reads
  if(numreads<1024)
    adepth[numreads]++;
  else
    adepth[1023]++;

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
  if(glffile)
    assert(sizeof(double)*10*bgzf_write(glffile,like,sizeof(double)*10));
  int noTrans =0;
  if(outfileSAF){
    // AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    //  0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    //  =. -. *. -. =. -. *. =. -. =
    double homo =log(0);
    double het= log(0);

    if(numreads!=0){
      const char homoz[] = {0,4,7,9};//=
      const char transverz[] = {1,3,5,8};//-
      const char transitz[] ={2,6};//* 
      
      for(int i=0;i<4;i++){
	homo = addProtect2(homo,like[homoz[i]]);
	het =  addProtect2(het,like[transverz[i]]);
      }
      for(int i=0;(noTrans==0&&i<2);i++)
	het = addProtect2(het,like[transitz[i]]);
    }
    if(sizeof(double)!=bgzf_write(outfileSAF,&homo,sizeof(double))){
      fprintf(stderr,"\t-> Problem writing hom\n");
      exit(0);
    }
    if(sizeof(double)!=bgzf_write(outfileSAF,&het,sizeof(double))){
      fprintf(stderr,"\t-> Problem writing het\n");
      exit(0);
    }
  }
  
  return numreads;
}

//nsam=2*nind
//only one individual input
void test (int *positInt,BGZF *gz,double errate,double meandepth, int regLen,BGZF *gzSeq,int count,BGZF *gzPsmc,int do_seq_glf,BGZF*outfileSAF,BGZF*outfileSAFPOS){
  //  fprintf(stderr,"posit:%p errate:%f meandepth:%f regLen:%d count:%d do_seq_glf:%d\n",positInt,errate,meandepth,regLen,count,do_seq_glf);
  kstring_t kpl;kpl.s=NULL;kpl.l=kpl.m=0;
  kstring_t kpsmc;kpsmc.s=NULL;kpsmc.l=kpsmc.m=0;
  int p =0;


  for(int i=0;i<regLen;i++) {
    int genotypes[2]={0,0};
    if(positInt[i])
      genotypes[1] = 1;
    if(do_seq_glf||outfileSAF!=NULL)
      print_ind_site(errate,meandepth,genotypes,gz,&kpl,outfileSAF);
    if(outfileSAFPOS)
      if(sizeof(int)!=bgzf_write(outfileSAFPOS,&i,sizeof(int))){
	fprintf(stderr,"\t-> Problem writing i\n");
	exit(0);
      }
  }
  

  if(gzPsmc!=Z_NULL){

    if(gzPsmc)
      ksprintf(&kpsmc,">%d\n",count+1);
    int at=0;
    for(int i=0;i<regLen;i+=block) {
      int isHet =0;
      for(int b=0;b<block;b++)
	if(positInt[i+b])
	  isHet++;
      if(at==NBASE_PER_LINE){
	ksprintf(&kpsmc,"\n");
	at=0;
      }
      if(isHet)
	ksprintf(&kpsmc,"K");
      else
	ksprintf(&kpsmc,"T");
      at++;
    }
    if(kpsmc.l>0&&kpsmc.s[kpsmc.l-1]!='\n')
      ksprintf(&kpsmc,"\n");
    assert(kpsmc.l==bgzf_write(gzPsmc,kpsmc.s,kpsmc.l));
    kpsmc.l=0;
  }
  free(kpsmc.s);
}

int main(int argc,char **argv){
  for(int i=0;i<1024;i++)
    adepth[i] = 0;
  double errate = 0.01;
  double meanDepth = 4;
  const char *inS = NULL;
  const char *prefix=NULL;
  FILE *in = NULL;
  argv++;
  int nind = 0;
  char **orig = argv;
  int seed = -1;\
  int do_seq_glf = 1;
  int nsam,howmany,theta,rho;
  int regLen=0;
  int j ,nsites, i;
  char line[1001] ;
  FILE *pf, *fopen(), *pfin ;
  int *positInt ;
  int   segsites;
  int psmcfa = 0;
  int psmc2 = 0;
  int nChr =-1;
  kstring_t tree_kstring;tree_kstring.s=NULL;tree_kstring.l=tree_kstring.m=0;
  while(*argv){
    if(strcasecmp(*argv,"-in")==0) inS = *++argv;
    else if(strcasecmp(*argv,"-out")==0) prefix=*++argv; 
    else if(strcasecmp(*argv,"-err")==0) errate=atof(*++argv); 
    else if(strcasecmp(*argv,"-depth")==0) meanDepth=atof(*++argv); 
    else if(strcasecmp(*argv,"-nChr")==0) nChr=atoi(*++argv); 
    else if(strcasecmp(*argv,"-win")==0) block=atoi(*++argv); 
    else if(strcasecmp(*argv,"-psmcfa")==0) psmcfa=atoi(*++argv);
    else if(strcasecmp(*argv,"-psmc2")==0) psmc2=atoi(*++argv);
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
    seed = time(NULL); 
  
  srand48(seed);
  if(inS==NULL||prefix==NULL){
    fprintf(stderr,"Probs with args, supply -in -out\n");
    fprintf(stderr,"also -err -depth -depthFile -regLen -nind -seed -pileup -psmc -do_seq_glf -win -psmcfa 1 -psmc2 -nChr \n");
    return 0;
  }
  fprintf(stderr,"-in=%s -out=%s -err=%f -depth=%f -regLen=%d -seed %d -nind %d -psmcfa %d -do_seq_glf: %d -win:%d -nChr: %d\n",inS,prefix,errate,meanDepth,regLen,seed,nind,psmcfa,do_seq_glf,block,nChr);
  //print args file
  FILE *argFP=openFile(prefix,".argg");
  fprintf(argFP,"-in=%s -out=%s -err=%f -depth=%f -regLen=%d -seed %d -nind %d -psmc %d -do_seq_glf: %d -nChr: %d\n",inS,prefix,errate,meanDepth,regLen,seed,nind,psmcfa,do_seq_glf,nChr);

  for(int i=0;i<argc;i++)
    fprintf(argFP,"%s ",orig[i]);
  fclose(argFP);
  
  const char *SAF = ".psmc.gz";
  const char *SAFPOS =".psmc.pos.gz";
  const char *SAFIDX =".psmc.idx";
  int64_t offs[2];
  BGZF *outfileSAF = NULL;
  BGZF *outfileSAFPOS = NULL;
  FILE *outfileSAFIDX = NULL;
  if(psmc2){
    outfileSAF = openFileGz(prefix,SAF,"wb");
    outfileSAFPOS =  openFileGz(prefix,SAFPOS,"wb");
    outfileSAFIDX = openFile(prefix,SAFIDX);
    char buf[8]="psmcv1";
    if(8!=bgzf_write(outfileSAF,buf,8)){
      fprintf(stderr,"\t-> Problem writing buf\n");
      return 0;
    }
    if(8!=bgzf_write(outfileSAFPOS,buf,8)){
      fprintf(stderr,"\t-> Problem writing buf\n");
      return 0;
    }
    fwrite(buf,1,8,outfileSAFIDX);
    offs[0] = bgzf_tell(outfileSAFPOS);
    offs[1] = bgzf_tell(outfileSAF);
  }
  
  in=getFILE(inS,"r");
  
  
  /* read in first two lines of output  (parameters and seed) */
  //  pfin = stdin ;
  pfin = in;

  int ss;
  if(NULL==fgets( line, 1000, pfin))
    fprintf(stderr,"Problem reading from file:\n");
  char *pch = strstr(line,"msHOT-lite");
  while(1){
    char *tmppch=strstr(pch+1,"msHOT-lite");
    if(tmppch=='\0')
      break;
    pch=tmppch;
  }
  ss=sscanf(pch,"msHOT-lite %d %d -t %d -r %d %d\n", &nsam, &howmany,&theta,&rho,&regLen);
  assert(ss==5);
  fprintf(stderr,"\t-> Number of samples:%d (haploids)\n",nsam);
  fprintf(stderr,"\t-> Number of replications:%d\n",howmany);
  fprintf(stderr,"\t-> theta: %d rho: %d regLen:%d\n",theta,rho,regLen);
  if(nsam!=2){
    fprintf(stderr,"Only defined for a single diploid\n");
    return 0;
  }

  positInt = (int *)malloc( regLen*sizeof( int ) ) ;


  BGZF *gz=NULL;
  BGZF *gzSeq = NULL;
  BGZF *gzPsmc = NULL;
  BGZF *tree = NULL;
  
  if(psmcfa)
    gzPsmc = openFileGz(prefix,".fa.gz","w");
  
  if(do_seq_glf)
    gz = openFileGz(prefix,".glf.gz","w");
  if(1)
    tree = openFileGz(prefix,".tree.gz","w");
    
  FILE *pgEst = openFile(prefix,".pgEstH");
  FILE *aDepthFP = openFile(prefix,".acounts.txt");

  fprintf(stderr,"\n");
  if(nChr!=-1)
    howmany=nChr;
  for(int i=0;i<howmany;i++){
      /* read in a sample */

    while(1){
      if( fgets( line, 1000, pfin) == NULL ){
	fprintf(stderr,"\t-> Problem reading file\n");
	break;
      }
      if(strcmp(line,"//\n")==0){
	//fprintf(stderr,"skipping line\n");
	continue;
      }
      //fprintf(stderr,"pre[%d]: %s",i,line);
      int span,en,to;
      float fl1,fl2;
      int treentok=sscanf(line,"[%d](%d:%f,%d:%f);\n",&span,&en,&fl1,&to,&fl2);
      if(treentok==5){
	//fprintf(stderr,"treentok:%d span:%d en:%d to:%d fl1:%f\n",treentok,span,en,to,fl1);
	//      fprintf(stderr,"treentok:%d span:%d en:%d to:%d fl1:%f fl2:%f\n",treentok,span,en,to,fl1,fl2);
	for(int j=0;j<span;j++)
	  ksprintf(&tree_kstring,"%d\t%d\t%f\n",i,j,fl1);
      }
      ss=sscanf(line,"@begin %d\n",&segsites);
      if(ss==1)
	break;
    }
    if(tree!=NULL & tree_kstring.l>4096){
      assert(tree_kstring.l==bgzf_write(tree,tree_kstring.s,tree_kstring.l));
      tree_kstring.l=0;
    }
    fprintf(stderr,"\t-> Parsing replicate:%d which should have nseg: %d ",i,segsites);
    memset(positInt,0,sizeof(int)*regLen);

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
    
    test(positInt,gz,errate,meanDepth,regLen,gzSeq,i,gzPsmc,do_seq_glf,outfileSAF,outfileSAFPOS);
    if(outfileSAFIDX){
      char tmpChr[1024];
      sprintf(tmpChr,"%d",i+1 );
      size_t clen = strlen(tmpChr);
      fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
      fwrite(tmpChr,1,clen,outfileSAFIDX);
      size_t tt = regLen;
      fwrite(&tt,sizeof(size_t),1,outfileSAFIDX);
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
      fflush(outfileSAFIDX);
      offs[0] = bgzf_tell(outfileSAFPOS);
      offs[1] = bgzf_tell(outfileSAF);
    }

   }


  for(int i=0;i<1024-1;i++)
    fprintf(aDepthFP,"%lu\t",adepth[i]);
  fprintf(aDepthFP,"%lu\n",adepth[i]);
  fclose(aDepthFP);
  if(tree){
    assert(tree_kstring.l==bgzf_write(tree,tree_kstring.s,tree_kstring.l));
    bgzf_close(tree);
  }
   if(psmcfa)
     bgzf_close(gzPsmc);
   if(gz!=Z_NULL)
     bgzf_close(gz);
   fclose(pgEst);
   fclose(in);
   free(positInt);
   if(outfileSAF)
     bgzf_close(outfileSAF);
   if(outfileSAFPOS)
     bgzf_close(outfileSAFPOS);
   if(outfileSAFIDX)
     fclose(outfileSAFIDX);
   return 0;
}
