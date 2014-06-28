
// program to simulate NGS data

// author Rasmus Nielsen

// Thorfinn:

// copied from simnextgen2.c
// changed code such that ancestral is always A
// 
// moved minfreq after argparsing
// 10 changed .geno output such that it reflects the "real" genotypes, an doesnt depend on what I've sampled. Added a much faster uniform sampling function. The old one can still be used by -simpleRand 0
// 9 added inbreeding
// 8 added a variable errate, mean value is as supplied by user, but range=[0,2*errate]
// 7 added inbreeding
// 6 added true genotype output
// 
// Matteo:
  // 2P version: 1 or 2 populations allowed:
  // main changes:
  //	function for generating random values from Beta distribution
  //	independent drawns from Balding-Nichols distribution for estimating pop alle freq for 2 subpopulations given an ancestral pop alle freq and a FST
  // use the unfolded spectrum (0 is ancestral)
  // output for each population and for the whole sample of pooled populations
  // output joint-SFS
  // use the unfolded spectrum


// to compile: cc -lm -lz -O simnextgenpop.c -o simnextgenpop

/*Note 0=A, 1=C, 2=G, 3=T*/
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <cstring>
#include <zlib.h>
#include <sys/stat.h>
#define PI 3.141592654 // for Poisson distribution

#include "rbeta.cpp" // random generator of values beta-distributed

// define some constants or global variables

int static z_rndu=137;
double minfreq, myConst;

int dumpBinary = 1; // output as a binary file?
int simpleRand = 0; // use Rasmus's or Thorfinn's uniform function?

// declare functions

void SetSeed(int seed) {
  z_rndu = 170*(seed%178) + 137;
}

/// QUESTION: SetSeed is never used through the program. Should we include this parameter?

/*U(0,1): AS 183: Appl. Stat. 31:188-190
Wichmann BA & Hill ID.  1982.  An efficient and portable pseudo-random number generator.  Appl. Stat. 31:188-190 x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.

Suggested to me by Z. Yang who also provided me with the source code used here. */
double uniform_rasmus() {
  static int x_rndu=11, y_rndu=23;
                                    // set seed here?
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

double uniform_thorfinn() {
  return((double)rand() / (double)RAND_MAX);
}

// which uniform function to use
double uniform() {
  double sampled;
  if(simpleRand)
    sampled = uniform_thorfinn();    
  else
    sampled = uniform_rasmus();
  //  fprintf(stdout,"%f\n",sampled);
  return sampled;
}

// function to check for file existence, using stat.h.
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

// functions to get gz files
gzFile getGz(const char*fname,const char* mode){
  if(fexists(fname)){
    fprintf(stderr,"\t->File exists: %s exiting...\n",fname);
    exit(0);
  }
  gzFile fp;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

FILE *getFile(const char*fname,const char *mode){
  if(fexists(fname)){
    fprintf(stderr,"\t-> File exists: %s exiting... Terminate\n",fname);
    exit(0);
  }
  FILE *fp=NULL;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening File handle for file:%s\n",fname);
    exit(0);
  }
  return fp;
}

// gamma distribution
double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947,-86.50532032942, 24.01409824083,-1.231739572450, 0.1208650973866e-2,-0.5395239385e-5};

  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -=(x+0.5)*log(tmp);
  ser=1.00000000019015;
  for (j=0; j<5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

//??????????? ERROR HERE !!!!!!!!!!!!!!!!!!!!!!!!

// // beta distribution using gamma distribution (Matteo)
// double beta(double a, double b) {
//   //  fprintf(stderr,"%f\t%f\n",a,b);
//   double x = uniform();
// //  double matteo = (tgamma(a+b)/(tgamma(a)*tgamma(b)))*pow(x, a-1)*pow(1-x, b-1);
//   double thorfinn = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(x)+(b-1)*log(1-x);
// //  fprintf(stderr,"%f\t%f\n",matteo,exp(thorfinn));
//   return(exp(thorfinn));
// //  return(matteo);
// }

// // beta distribution using gamma distribution (Matteo)
// double beta(double a, double b) {
//   double X = gammln(a);
//   double Y = gammln(b);
//   return (X/(X+Y));
// }




// Poisson distribution
double Poisson(double xm) {
  double gammln(double xx);
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
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(PI*uniform());
	em=sq*y+xm;
      } while (em< 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (uniform()>t);
  }
  return em;
}

// generate sequencing errors by picking a different base
int basepick_with_errors(double errate, int inbase) {
  // from error rate and a starting base
  int outbase;
  if (uniform()<errate){ // if error occurs then change base
    while ((outbase=(floor(4*uniform()))) == inbase); // then take a different base
    return outbase;
  }
  else return inbase;
}

// pick a genotype, but with errors for each allele of the genotype[2]
int pick_a_base(double errate, int genotype[2])	{
  if (uniform()<0.5) 
    return basepick_with_errors(errate, genotype[0]);
  else 
    return basepick_with_errors(errate, genotype[1]);
}

// calculate the likelihood for each genotype configuration
void calclike(int base, double errate, double *like)	{
  /*0=AA, 1=AC, 2=AG, 3=AT, 4=CC, 5=CG, 6=CT, 7=GG, 8= GT, 9= TT*/
  double  prob;
  int i, j, k=0;
  for (i=0; i<4; i++){  // all 4 possible alleles for 1st base of genotype
    for (j=i; j<4; j++){ // all 4 possible alleles for 2nd base of genotype 
      if (base==i) 
	prob = (1.0-errate)/2.0; // if right base
      else 
	prob = errate/6.0; // if wrong base
      if (base==j) 
	prob = prob +(1.0-errate)/2.0;
      else 
	prob = prob + errate/6.0;
       if (prob <= 0.0) 
	like[k] = -1000000.0;
      else 
	like[k] = like[k] + log(prob); // add all log(prob) as likelihood
      k++;
    }
  }
}

// translate from int code to letters, for like computation, maybe it's useless now
char int_to_base(int b) {
  if (b==0) return 'A';
  else if (b==1) return 'C';
  else if (b==2) return 'G';
  else return 'T';
}

// compute and print results into files
int print_ind_site(double errate, double meandepth, int genotype[2], gzFile resultfile, gzFile glffile) {

  int i, b, numreads;
  double like[10];

  numreads = Poisson(meandepth);  // mumber of reads, poisson distributed with lambda=meandepth
  char res[numreads]; // char alleles for all the reads

  for (i=0; i<10; i++) like[i] = 0.0; // initialize GL values for all 10 possible genotypes
  
  // compute like
  for (i=0; i<numreads; i++) {
    b = pick_a_base(errate, genotype); // pick a base including errors
    res[i] = int_to_base(b); 
    calclike(b, errate, like); // compute likelihood
  }

  // write into files
  if(dumpBinary) { // if binary 
    gzwrite(resultfile, res, numreads);
    gzwrite(glffile, like, sizeof(double)*10); // print likelihoods values
  } else { // text
    fprintf(stderr,"textout not implemented\n"); // not yet...
    exit(0);
    //    for (i=0; i<10; i++) //fprintf(glffile," %f",like[i]);
    //  fprintf(glffile," ,");
  }

  return numreads; // output is the number fo reads

}

// randomly pick a base from prior basefreq
int pick_base_from_prior(double basefreq[4])	{
  int i = 0;
  double U, p;
  
  U = uniform();
  p = basefreq[0];
  while (U>p)
    {
      i++;
      p = p + basefreq[i];
    }
  return i;
}

// simulate sfs for 1 population (or the ancestral one)
double simfreq() {
  return exp(uniform()*myConst-myConst);
}

// simulate sfs for 2 subpopulations using Balding-Nichols distribution (Matteo)
double simfreqBN(double F, double p_anc) {
  // FST values for the 2 subpops
  // p_anc: ancestral population allele frequency, drawn from a truncated exponential distribution (some other authors use a uniform distribution, I like using an exponential better).
  return rbeta( ((1-F)/F)*p_anc, ((1-F)/F)*(1-p_anc)) ;
}

// ???
char *append(const char* a,const char *b){
  char *c =(char *) malloc((strlen(a)+strlen(b)+1)*sizeof(char));
  strcpy(c,a);
  strncat(c,b,strlen(b));
  return c;
}

// help printout
void info() {
  fprintf(stderr,"\t -> Required arg:\n\t\t-outfiles PREFIX\t PREFIX.seq PREFIX.glf PREFIX.frq PREFIX.arg\n");
  fprintf(stderr,"\t -> Optional arg:\n\t\t-npop\tNumber of populations. This MUST be set before -nind [1]\n");
  fprintf(stderr,"\t\t-nind\tNumber of diploid individuals for each population [10]\n");
  fprintf(stderr,"\t\t-nsites\tNumber of sites [500000]\n");
  fprintf(stderr,"\t\t-errate\tThe sequencing error rate [0.0075]\n");
  fprintf(stderr,"\t\t-depth\tMean sequencing depth [5]\n");
  fprintf(stderr,"\t\t-pvar\tProbability that a site is variable in the population [0.015]\n");
  fprintf(stderr,"\t\t-mfreq\tMinimum population frequency [0.0001]\n");
  fprintf(stderr,"\t\t-F\tinbreeding coefficient for each population [0]\n");
  fprintf(stderr,"\t\t-model\t0=fixed errate 1=variable errate [1]\n");
  fprintf(stderr,"\t\t-simpleRand\tboolean [1]\n");
  fprintf(stderr,"\t\t-base_freq\tBackground allele frequencies for A,C,G,T [0.25 0.25 0.25 0.25]\n");
}

///		///
/// MAIN	///
///		///

// to compile: g++ -lm -lz -Wall -O simnextgenpop.c -o simnextgenpop

int main(int argc, char *argv[]) { // read input parameters

  if(argc==1) { // if no argument (call of the program is considered as 1)
    info(); // retunr info
    return 0; // terminate
  }

  /// define and initialize the variables (with default values)
  
  int i=0, j=0, k=0, b1=0, b2=0, var=0, nsites = 500000, nind = 10, npop=1, model=1, nind1=0, nind2=0, increment=0;
  static int genotype[2], genotype1[2], genotype2[2]; // array /matrix of genotypes for all pops
  double pfreq=0.0, pfreq1=0.0, pfreq2=0.0, pvar= 0.015, meandepth = 5, errate = 0.0075, F=0.0, F1=0.0, F2=0.0, minfreq=0.0001;
  double basefreq[4] = {0.25, 0.25, 0.25, 0.25}; // background frequencies
  // as pointers
/*  int *freqspec = NULL; // whole sample
  int *freqspec1 = NULL; // 1st pop
  int *freqspec2 = NULL; // 2nd pop*/
 
// for debug
  static int basecheck[4], basecheck1[4], basecheck2[4];

  //filehandles and their associated names
  // resultfile: reads
  // glffile: genotype likelihoods per site
  // freqfile: whole SFS
  // argfile: input parameters
  gzFile glffile, resultfile, glffile1, glffile2, resultfile1, resultfile2;
  FILE  *freqfile, *argfile, *genofile, *freqfile1, *freqfile2, *genofile1, *genofile2, *freqfile12;
  char *fGlf=NULL,*fFreq=NULL,*fSeq=NULL,*fArg=NULL,*fGeno=NULL;
  char *fGlf1=NULL,*fFreq1=NULL,*fSeq1=NULL,*fGeno1=NULL;
  char *fGlf2=NULL,*fFreq2=NULL,*fSeq2=NULL,*fGeno2=NULL;
  char *fFreq12=NULL;
  char *outfiles=NULL;
  
  /// read command line parameters and assign
  
  // read input parameters and assign to values variables (if defined)
  int argPos=1; // we checked before if there were no inputs
  increment=0; // increment over the input parameters (one input may have more values than 1)
  
  while (argPos<argc) { // for all inputs

    // increment of 2 is for only 1 value input, otherwise is larger 
    int increment=0; // increment over the input parameters (one input may have more values than 1)

    if(strcmp(argv[argPos],"-npop")==0) npop = atoi(argv[argPos+1]); 

    else if(strcmp(argv[argPos],"-nind")==0) {
      increment = increment + npop -1; // argPos will be then 1+1=2 if 2 pops
      if (npop==1) {
	nind = atoi(argv[argPos+1]);
	nind1 = nind2 =0;
      } else if (npop==2) {
	nind1 = atoi(argv[argPos+1]);
	nind2 = atoi(argv[argPos+2]);
	nind = nind1+nind2;
      } else {
	printf("\n-npop must be 1 or 2. Terminate.\n");
	return 0;
      } // end if npop
    }

    else if(strcmp(argv[argPos],"-F")==0) {
      increment = increment + npop -1; // argPos will be then 1+1=2 if 2 pops
      if (npop==1) {
	F = atof(argv[argPos+1]);
      } else if (npop==2) {
	F1 = atof(argv[argPos+1]);
	F2 = atof(argv[argPos+2]);
	F = 0.0; // no "extra inbreeding"
      } // end if npop
    }

    // argPos+1 because you count the -flag too (argPos is 1!)
    // argPos+1+npop-1 == argPos+npop; es. 2 pops, increment is 3 for the next inpute read 

    else if(strcmp(argv[argPos],"-errate")==0)   errate  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-depth")==0) meandepth  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-pvar")==0)  pvar  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-mfreq")==0)  minfreq  = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-outfiles")==0) outfiles  = (argv[argPos+1]);
    else if(strcmp(argv[argPos],"-nsites")==0)  nsites  = atoi(argv[argPos+1]);
//    else if(strcmp(argv[argPos],"-FST")==0)  FST = atof(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-model")==0)  model  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-simpleRand")==0) simpleRand  = atoi(argv[argPos+1]);
    else if(strcmp(argv[argPos],"-base_freq")==0) {
     increment = increment + 3;
     basefreq[0] = atof(argv[argPos+1]); basefreq[1] = atof(argv[argPos+2]); basefreq[2] = atof(argv[argPos+3]); basefreq[3] = atof(argv[argPos+4]); 
    }
        
    else { // input is not a valid one 
     
      printf("\tUnknown arguments: %s\n",argv[argPos]);
      info();
      return 0; // terminate
      
    }
    
    argPos = argPos + 2 + increment; // increment next value (+2 because un count the input flags too)
  
  } // end while all inputs
  
  // check if there is the output file
   if(outfiles == NULL) {
    fprintf(stderr,"\nMust supply -outfiles. Terminate\n");
    return 0;
  }

  // to adjust the exponential function according to the lowest allele frequency detectable
  myConst = -log(minfreq);
  
  // print input arguments
  fprintf(stderr,"\t->Using args: -npop %d -nind %d -nind1 %d -nind2 %d -errate %f -depth %f -pvar %f -mfreq %f -nsites %d -F %f -F1 %f -F2 %f -model %d -simpleRand %d -base_freq %f %f %f %f\n", npop, nind, nind1, nind2, errate, meandepth, pvar, minfreq, nsites, F, F1, F2, model, simpleRand, basefreq[0], basefreq[1], basefreq[2], basefreq[3]); 
  
  /// output files
  //   and  get gz and file
  // files
  fArg = append(outfiles,".args");
  // whole
  fGlf = append(outfiles,".glf.gz");
  fFreq = append(outfiles,".frq");
  fSeq = append(outfiles,".seq.gz");
  fGeno = append(outfiles,".geno");
  // print message
  fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile: %s\ttruefreq: %s args:%s geno:%s\n",fSeq,fGlf,fFreq,fArg,fGeno);

  argfile=getFile(fArg,"w");
  // whole
  resultfile=getGz(fSeq,"w");
  glffile = getGz(fGlf,"w");
  freqfile= getFile(fFreq,"w"); 
  genofile =getFile(fGeno,"w");
  

  // eventually subpops and messages prints
  if (npop==2) {
    
    fGlf1 = append(outfiles,"1.glf.gz");
    fFreq1 = append(outfiles,"1.frq");
    fSeq1 = append(outfiles,"1.seq.gz");
    fGeno1 = append(outfiles,"1.geno");
    fGlf2 = append(outfiles,"2.glf.gz");
    fFreq2 = append(outfiles,"2.frq");
    fSeq2 = append(outfiles,"2.seq.gz");
    fGeno2 = append(outfiles,"2.geno");
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile1: %s\ttruefreq1: %s geno1:%s\n", fSeq1, fGlf1, fFreq1, fGeno1);
    fprintf(stderr,"\t->Dumping files: sequencefile: %s\tglffile2: %s\ttruefreq2: %s geno2:%s\n", fSeq2, fGlf2, fFreq2, fGeno2);

    resultfile1=getGz(fSeq1, "w");
    glffile1 = getGz(fGlf1, "w");
    freqfile1= getFile(fFreq1, "w"); 
    genofile1 =getFile(fGeno1, "w");
    resultfile2=getGz(fSeq2, "w");
    glffile2 = getGz(fGlf2, "w");
    freqfile2= getFile(fFreq2, "w"); 
    genofile2 =getFile(fGeno2, "w");

        // joint-SFS
    fFreq12 = append(outfiles,"12.frq");    
    freqfile12 = getFile(fFreq12, "w"); 
  
  }

  // write args file
  fprintf(argfile,"\t->Using args: -npop %d -nind %d -nind1 %d -nind2 %d -errate %f -depth %f -pvar %f -mfreq %f -nsites %d -F %f -F1 %f -F2 %f -model %d -simpleRand %d -base_freq %f %f %f %f\n", npop, nind, nind1, nind2, errate, meandepth, pvar, minfreq, nsites, F, F1, F2, model, simpleRand, basefreq[0], basefreq[1], basefreq[2], basefreq[3]); 

  /// COMPUTE
  
  // initialize SFS (whole locus), as the unfolded spectrum (ancestral and derived)
   // as arrays
  int freqspec[nind*2+1], freqspec1[nind1*2+1], freqspec2[nind2*2+1];
  int freqspec12[nind1*2+1][nind2*2+1]; // joint-SFS, matrix

//    , sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));


  for (i=0; i<nind*2+1; i++) freqspec[i]=0; // initialize sfs (from 0 to 2*nind)
  for (i=0; i<nind1*2+1; i++) freqspec1[i]=0;
  for (i=0; i<nind2*2+1; i++) freqspec2[i]=0;
  for (i=0; i<nind1*2+1; i++) {
    for (j=0; j<nind2*2+1; j++) {
      freqspec12[i][j]=0;
    }
  } // end initialzie joint-SFS

//     printf("\n\n 0   size of %d %d %d :", sizeof(freqspec)/sizeof(freqspec[0]), sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));

  /// FOR EACH SITE
  
  for (i=0; i<nsites; i++) {
    
    /*debug code*/
    for (k=0; k<4; k++) basecheck[k]=basecheck1[k]=basecheck2[k]=0; // initialize base checks      
    
    // basechecks are arrays of 4 elements counting the nr of alleles 0,1,2,3
    
    /// test if site is variable or not
    
    if (uniform() >= pvar) { // if it is NOT variable... 

      var=0; // not variable (boolean)
      // assign alleles to genotype
      genotype[0]=genotype[1]=0; // all ancestral A

      /*debug code*/
      basecheck[genotype[0]] = 2*nind; // all individuals are monorphic for ancestral A, basechecks are arrays of 4 elements counting the nr of alleles 0,1,2,3
      if (npop==2) {
	genotype1[0]=genotype1[1]=genotype2[0]=genotype2[1]=0;
	basecheck1[genotype1[0]] = 2*nind1; // all individuals are monorphic
	basecheck2[genotype2[0]] = 2*nind2; // all individuals are monorphic
      }
      
      
    } else { // if it IS variable: finally a SNP!
      
      var = 1; // TRUE 
      
      // choose the 2 alleles
      b2 = 0; //changed such that reference is always A
      while ((b1=pick_base_from_prior(basefreq))==b2); // // take second allele, different from the first A
      
      // we will assign genotypes on the basis of the 2 alleles later, not here
      // here it is just a check if the site is variable and which alleles to take.
	
      // simulate population allele frequency (or ancestral if 2 subpops)
      pfreq=simfreq(); // if site is not variable, you don't need to compute pfreq
      
    } // end test if it is variable or not (pvar)
    
    /// WHOLE POPULATIONS (or the only one if -npop=1)
 
    // now for each individual, assigned genotypes at each individual as couples of alleles/bases
   if (npop==1) {
    for (j=0; j<nind; j++) {

      if (var==1) {

	if (uniform()>=F) { //no inbreeding case
	  for (k=0; k<2; k++) {
	    if (uniform()<=pfreq) 
	      genotype[k] = b1;
	    else genotype[k] = b2; 
	  }
	} else { //inbreeding case
	  if (uniform()<=pfreq) {
	    genotype[0] = b1;
	    genotype[1] = b1;
	  } else {
	    genotype[0] = b2; 
	    genotype[1] = b2;
	  }
	} // end if uniform<=> F
	
	/*debug code*/ 
	// increment the count for each allele type on the basis of just assigned genotype
	basecheck[genotype[0]]++; basecheck[genotype[1]]++;
	
      } // end assignment of genotypes of site is variable
      
      // write genotypes in the output file (append it)
      fprintf(genofile,"%d %d\t",genotype[0],genotype[1]);      
      int has_reads =0;
      // compute and print likelihoods
      if (model==0)
	has_reads=print_ind_site(errate, meandepth, genotype, resultfile, glffile);
      else
	has_reads=print_ind_site(2*errate*uniform(), meandepth, genotype, resultfile, glffile);

      // now write binary files
      if (j<nind-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile,sep,1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	  //fprintf(glffile,",");
	  // fprintf(resultfile,", ");
	}
      } // for nind

    } // end for j in nind, all individuals
   } // if npop==1
    
    /// MULTI POPULATIONS
    
    if (npop==2) {

      // then assign subpop freq based on time-split? theory of pop genet
  
      // use now balding-nichols and write one set of files for each population
      // simulate distinct pops allele freq from Balding-Nickols distribution, given an ancestral population allele frequency and an FST value
      // 2 independent draws
      // pfreq already simulated!
      pfreq1=simfreqBN(F1, pfreq);
      pfreq2=simfreqBN(F2, pfreq);
      
      
/*           printf("\nPANC: %f",pfreq);
	    printf("\tP1: %f",pfreq1);
          printf("\tP2: %f",pfreq2);*/
	  
// 	  //.poi togli
// 	  pfreq1=0.1;
// 	  	  pfreq2=0.1;
      
      /// 1st pop (same as the whole)
      for (j=0; j<nind1; j++) {
	if (var==1) {
	  if (uniform()>F) { // this will alwasy be TRUE now 
	    for (k=0; k<2; k++) {
	      if (uniform()<=pfreq1) 
		genotype1[k] = b1;
	    else genotype1[k] = b2; 
	  }
	} else {
	  if (uniform()<=pfreq1	) {
	    genotype1[0] = b1;
	    genotype1[1] = b1;
	  } else {
	    genotype1[0] = b2; 
	    genotype1[1] = b2;
	  }
	}	
	basecheck1[genotype1[0]]++; basecheck1[genotype1[1]]++;
	basecheck[genotype1[0]]++; basecheck[genotype1[1]]++;
	}
      fprintf(genofile1,"%d %d\t",genotype1[0],genotype1[1]);     
      // write also into whole genotype file
      fprintf(genofile,"%d %d\t",genotype1[0],genotype1[1]);
      int has_reads1 =0, has_reads=0;
      if (model==0) {
	has_reads1=print_ind_site(errate, meandepth, genotype1, resultfile1, glffile1);
	has_reads=print_ind_site(errate, meandepth, genotype1, resultfile, glffile);
      } else {
	has_reads1=print_ind_site(2*errate*uniform(), meandepth, genotype1, resultfile1, glffile1);
	has_reads=print_ind_site(2*errate*uniform(), meandepth, genotype1, resultfile, glffile);
      }
      if (j<nind1-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile1, sep, 1);
	  gzwrite(resultfile, sep, 1); // write also into whole genotype file
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);
	}
      }
      }

      /// 2nd pop (same as above)
      for (j=0; j<nind2; j++) {
	if (var==1) {
	  if (uniform()>F) { // this will always be TRUE now
	    for (k=0; k<2; k++) {
	      if (uniform()<=pfreq2) 
		genotype2[k] = b1;
	    else genotype2[k] = b2; 
	  }
	} else {
	  if (uniform()<=pfreq2	) {
	    genotype2[0] = b1;
	    genotype2[1] = b1;
	  } else {
	    genotype2[0] = b2; 
	    genotype2[1] = b2;
	  }
	}	
	basecheck2[genotype2[0]]++; basecheck2[genotype2[1]]++;
	basecheck[genotype2[0]]++; basecheck[genotype2[1]]++;
      }

      fprintf(genofile2,"%d %d\t",genotype2[0],genotype2[1]);      
      fprintf(genofile,"%d %d\t",genotype2[0],genotype2[1]);      
      int has_reads2 =0, has_reads=0;
      if (model==0) {
	has_reads2=print_ind_site(errate, meandepth, genotype2, resultfile2, glffile2);
	has_reads=print_ind_site(errate, meandepth, genotype2, resultfile, glffile);
      } else {
	has_reads2=print_ind_site(2*errate*uniform(), meandepth, genotype2, resultfile2, glffile2);
	has_reads=print_ind_site(2*errate*uniform(), meandepth, genotype2, resultfile, glffile);
      }
      if (j<nind2-1) {
	if(dumpBinary) {
	  char sep[1]={'\t'};
	  gzwrite(resultfile2, sep, 1);
	  gzwrite(resultfile, sep, 1);
	} else {
	  fprintf(stderr,"non binary output disabled\n");
	  exit(0);

	}
      }
    }
    
    } // end if multi populations npop==2


    // EOL
    fprintf(genofile,"\n");
    if (npop==2) { 
      fprintf(genofile1,"\n"); 
      fprintf(genofile2,"\n"); 
    }

    /// DEBUG CODE AND COMPUTE SFS

    if (1) { /// QUESTION: what does this sentence (1) mean? What's TRUE? And why?
            
/*    printf("\n\n 1   size of %d %d %d :", sizeof(freqspec)/sizeof(freqspec[0]), sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));*/
      
      
      int k=0, j=0, j1=0, j2=0;
      for (k=1; k<4; k++) { // count only derived alleles
	if (basecheck[k]>j) j=basecheck[k];
      } // j is the count of the derived allele   
      if ((j+basecheck[0])!=nind*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	printf("error in freqspec calculation (%i %i %i %i)\n", basecheck[0], basecheck[1], basecheck[2], basecheck[3]); 
	exit(-1);
      }
      // then now j is the count of non-ancestral alleles (not 0)
      freqspec[j]++; // this if unfolded

//     printf("\n\n 2   size of %d %d %d :", sizeof(freqspec)/sizeof(freqspec[0]), sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));


      if (npop==2) {

	k=0, j1=0;
	for (k=1; k<4; k++) { // count only derived alleles
	  if (basecheck1[k]>j1) j1=basecheck1[k];
	} // j is the count of the derived allele     
	if ((j1+basecheck1[0])!=nind1*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	  printf("error in freqspec1 calculation (%i %i %i %i)\n", basecheck1[0], basecheck1[1], basecheck1[2], basecheck1[3]); 
	  exit(-1);
	}
	// then now j is the count of non-ancestral alleles (not 0)
	freqspec1[j1]++; // this if unfolded

//     printf("\n\n 3   size of %d %d %d :", sizeof(freqspec)/sizeof(freqspec[0]), sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));

	k=0, j2=0;
	for (k=1; k<4; k++) { // count only derived alleles
	  if (basecheck2[k]>j2) j2=basecheck2[k];
	} // j is the count of the derived allele     
	if ((j2+basecheck2[0])!=nind2*2) { // if sum of allele-counts is not equal to nind*2 (nr chromosomes)
	  printf("error in freqspec2 calculation (%i %i %i %i)\n", basecheck2[0], basecheck2[1], basecheck2[2], basecheck2[3]); 
	  exit(-1);
	}
	// then now j is the count of non-ancestral alleles (not 0)
	freqspec2[j2]++; // this if unfolded
	
/*    printf("\n\n 4   size of %d %d %d :", sizeof(freqspec)/sizeof(freqspec[0]), sizeof(freqspec1)/sizeof(freqspec1[0]), sizeof(freqspec2)/sizeof(freqspec2[0]));*/
	
	// joint-SFS
	freqspec12[j1][j2]++;
      }
    } // end if (1) 
     
    /// write reads and GLF
    if (dumpBinary) {
      char sep[1] = {'\n'};
      gzwrite(resultfile, sep, 1);
    } else {
      fprintf(stderr, "non binary output disabled\n");
      exit(0);
    }
    if (npop==2) {
      // 1st
      if (dumpBinary) {
      char sep[1] = {'\n'};
      gzwrite(resultfile1, sep, 1);
    } else {
      fprintf(stderr, "non binary output disabled\n");
      exit(0);

    }
      // 2nd
      if (dumpBinary) {
      char sep[1] = {'\n'};
      gzwrite(resultfile2, sep, 1);
    } else {
      fprintf(stderr, "non binary output disabled\n");
      exit(0);

    }
    } // end dump binary if 2 pops

    // computed and written for a site
    
    /// QUESTION: should we provide the true SFS-per-site? Kind of useless for further analyses? I think so.

  } // end for in i each sites

  /// COMPUTE WHOLE SFS
  k=0;
  for (i=0; i<nind*2+1; i++)
    k=k+freqspec[i]; // how many alleles in the whole sample
  for (i=0; i<nind*2+1; i++) { // compute relative frequency for each allele freq bin
    fprintf(freqfile,"%f\t",(double)freqspec[i]/(double)k);
  }
  fprintf(freqfile,"\n"); // EOL
       
  // npop==2
  if (npop==2) {
    k=0;
    for (i=0; i<nind1*2+1; i++)
      k=k+freqspec1[i]; // how many alleles in the whole sample
    for (i=0; i<nind1*2+1; i++) { // compute relative frequency for each allele freq bin
      fprintf(freqfile1,"%f\t",(double)freqspec1[i]/(double)k);
    }
    fprintf(freqfile1,"\n"); // EOL

    k=0;
    for (i=0; i<nind2*2+1; i++)
    k=k+freqspec2[i]; // how many alleles in the whole sample
    for (i=0; i<nind2*2+1; i++) { // compute relative frequency for each allele freq bin
      fprintf(freqfile2,"%f\t",(double)freqspec2[i]/(double)k);
    }
    fprintf(freqfile2,"\n"); // EOL 
    
    k=0; // joint
    for (i=0; i<nind1*2+1; i++) {
      for (j=0; j<nind2*2+1; j++) {
	k=k+freqspec12[i][j];
      }
    }

    for (i=0; i<nind1*2+1; i++) {
      for (j=0; j<nind2*2+1; j++) {
	fprintf(freqfile12,"%f\t",(double)freqspec12[i][j]/(double)k);
      }
      fprintf(freqfile12,"\n"); // EOL
    }
    
  } // end if 2 pops

  /// FREE MEMORY, FLUSH AND CLOSE
  
//  free(freqspec);
  free(fGlf);
  free(fFreq);
  free(fSeq);
  free(fArg);

  gzclose(resultfile);//fclose flushed automaticly
  gzclose(glffile);
  fclose(argfile);
  fclose(freqfile);
  fclose(genofile);

  if (npop==2) {

//    free(freqspec1);
    free(fGlf1);
    free(fFreq1);
    free(fSeq1);

//    free(freqspec2);
    free(fGlf2);
    free(fFreq2);
    free(fSeq2);
   
    free(fFreq12);
    
    gzclose(resultfile1);//fclose flushed automaticly
    gzclose(glffile1);
    fclose(freqfile1);
    fclose(genofile1);
    gzclose(resultfile2);//fclose flushed automaticly
    gzclose(glffile2);
    fclose(freqfile2);
    fclose(genofile2);

    fclose(freqfile12);
    
  }

  return 0; // return value

} // end main












