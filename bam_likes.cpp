/*
  This code is very much copied from bam2bcf.c and errmod.c from samtools 1.18

  credit to auther of SAMtools

 */

#include <cstdio>
#include <stdint.h>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <ctype.h>
#include "ksort.h"
#include "mUpPile.h"


KSORT_INIT_GENERIC(uint16_t)

#define ERR_DEP 0.83f

extern int refToInt[256];

#if 0
typedef struct __errmod_coef_t {
  


  double *fk;//fk function of n: fk(n)=0.83^n*0.97+0.03
  double *beta; 
  double *lhet;//lhet function of n,k : lhet(n,k)=(1/2^n)*choose(n,k)
} errmod_coef_t;



typedef struct {
	double depcorr;
	errmod_coef_t *coef;
} errmod_t;


typedef struct {
  double bsum[16];
  uint32_t c[16];
} call_aux_t;


 errmod_coef_t *cal_coef(double depcorr, double eta)
{
  
  int k, n, q;
  long double sum, sum1;
  double *lC;
  errmod_coef_t *ec;
  
  ec =(errmod_coef_t *) calloc(1, sizeof(errmod_coef_t));
  
  ec->fk = (double*)calloc(256, sizeof(double));
  ec->fk[0] = 1.0;
  for (n = 1; n != 256; ++n){
    
    ec->fk[n] = pow(1. - depcorr, n) * (1.0 - eta) + eta;
    
  }

  ec->beta = (double*)calloc(256 * 256 * 64, sizeof(double));
  lC = (double*)calloc(256 * 256, sizeof(double));
  for (n = 1; n != 256; ++n) {
    double lgn = lgamma(n+1);
    for (k = 1; k <= n; ++k)
      lC[n<<8|k] = lgn - lgamma(k+1) - lgamma(n-k+1);
  }
  for (q = 1; q != 64; ++q) {
    double e = pow(10.0, -q/10.0);//phredprob
    double le = log(e);//log(phredprob)
    double le1 = log(1.0 - e); //log(1-phredborg)
    for (n = 1; n <= 255; ++n) {
      double *beta = ec->beta + (q<<16|n<<8);
      sum1 = sum = 0.0;
      for (k = n; k >= 0; --k, sum1 = sum) {
	sum = sum1 + expl(lC[n<<8|k] + k*le + (n-k)*le1);
	beta[k] = -10. / M_LN10 * logl(sum1 / sum);
      }
    }
  }

  ec->lhet = (double*)calloc(256 * 256, sizeof(double));
  for (n = 0; n < 256; ++n)
    for (k = 0; k < 256; ++k)
      ec->lhet[n<<8|k] = lC[n<<8|k] - M_LN2 * n;
  free(lC);
  return ec;
}

//depcorr = 1-0.83=0.17
errmod_t *errmod_init(float depcorr) {
  //  fprintf(stderr,"[%s]",__FUNCTION__);
  errmod_t *em;
  em = (errmod_t*)calloc(1, sizeof(errmod_t));
  em->depcorr = depcorr;
  em->coef = cal_coef(depcorr, 0.03);
  return em;
}

void errmod_destroy(errmod_t *em)
{
  if (em == 0) return;
  free(em->coef->lhet); free(em->coef->fk); free(em->coef->beta);
  free(em->coef); free(em);
}



// qual:6, strand:1, base:4
int errmod_cal(const errmod_t *em, int n, int m, uint16_t *bases, float *q)
{
  //  fprintf(stderr,"n=%d\n",n);
  if (n == 0) 
    return 0; //no data

  //reset results array;
  memset(q, 0, m * m * sizeof(float));

  //set "counts" array
  int w[32];
  memset(w, 0, 32 * sizeof(int));
  for(int i=0;0&&i<32;i++)
    fprintf(stderr,"pre w[%d]=%d\n",i,w[i]);
      
  call_aux_t aux;
  memset(&aux, 0, sizeof(call_aux_t));

  //sample 255 if depth >255
  if ((n > 255)) { // then sample 255 bases //THIS MAKES
    ks_shuffle(uint16_t, n, bases);
    n = 255;
  }
  ks_introsort(uint16_t, n, bases);
  
  for (int j = n - 1; j >= 0; --j) {
    uint16_t b = bases[j];
    //    fprintf(stderr,"j=%d q=%d strand=%d base=%d\n",j,b>>5,0x1&(b>>4),b&0xf);
    //    fprintf(stderr,"j=%d q=%d\tc=%c\t",j,(b>>5) + 33,(b>>5) + 33);
    //cap quality at [4,63]
    int q = b>>5 < 4? 4 : b>>5;
    if (q > 63) 
      q = 63;

    int k = b&0x1f;
    
    
    aux.bsum[k&0xf] += em->coef->fk[w[k]] * em->coef->beta[q<<16|n<<8|aux.c[k&0xf]];
    ++aux.c[k&0xf];
    ++w[k];
    for(int i=0;0&&i<32;i++)
      fprintf(stderr,"w[%d]=%d\n",i,w[i]);
  }
  if(0){
    //floating point inprecision compared with samtools binary output. But is correct
    
    //the genotype like p(data|A1=g1,A2=g2) = p(data|A1=g2,A2=g1)
    for (int g1 = 0; g1 <5; ++g1) {//allele1=0,1,2,3,4
      for (int g2 = g1; g2<5; ++g2) {//allele2=0,1,2,3,4
	if(g1!=g2){ // A1!=A2 - heterozygoues
	  int cjk = aux.c[g1] + aux.c[g2];//total depth for allele g1+allele g2
	  //binomial when ignoring non A1/A2 alleles: Bin (n,k,p);n=cjk=#g1+#g2 , k= aux.c[g2] = #g2 ; p=0.5. returns log
	  q[g1*5+g2] = -4.343 * em->coef->lhet[cjk<<8|aux.c[g2]]; 
	}
	for (int k = 0; k <5; ++k){
	  if(k!=g1 && k!=g2) //if a read has a non A1/A2 alleles it is an error. add the log of the prob of these reads
	    q[g1*5+g2] += aux.bsum[k];
	  
	}

	//mirror
	if(g1!=g2)
	  q[g2*5+g1] =  q[g1*5+g2];
	
	if (q[g1*5+g2] < 0.0) 
	  q[g1*5+g2] = 0.0;
	
      }
    }
    return 0;
  }
  // generate likelihood THIS WORKS PERFECTLY june 4 ande
  for (int g1 = 0; g1 <5; ++g1) {//j=0,1,2,3,4
    for (int g2 = g1; g2<5; ++g2) {//j=0,1,2,3,4
      if(g1==g2){ 
	for (int k = 0; k <5; ++k){
	  if(k!=g1)
	    q[g1*5+g2] += aux.bsum[k];
	}
      }
      else{
	int other=0;
	float tmp1=0;
	for (int k = 0; k < 5; ++k) 
	  if (k != g1 && k != g2) {
	    tmp1 += aux.bsum[k]; 
	    other = 1; 
	}
	int cjk = aux.c[g1] + aux.c[g2];
	if (other) 
	  q[g1*5+g2] = q[g2*5+g1] = -4.343 * em->coef->lhet[cjk<<8|aux.c[g2]] + tmp1;
	else 
	  q[g1*5+g2] = q[g2*5+g1] = -4.343 * em->coef->lhet[cjk<<8|aux.c[g2]]; // all the bases are either j or k

      }
      if (q[g1*5+g2] < 0.0) 
	q[g1*5+g2] = 0.0;

    }
  }
 
  return 0;
  

  //old original almost
  for (int j = 0; j != m; ++j) {//j=0,1,2,3,4
    // homozygous
    for (int k = 0; k != m; ++k){//only updates if k!=j and aux.c[k]!=0
      fprintf(stderr,"\t-> j=%d aux.c[%d]=%d aux.bsum[%d]=%f res=%d ",j,k,aux.c[k],k,aux.bsum[k],j*m+j);
      if (k != j && aux.c[k]) {
	fprintf(stderr,"USING \n");
	q[j*m+j] += aux.bsum[k];
      }else
	fprintf(stderr,"skipping\n");
    }
    
    // heterozygous
    for (int k = j + 1; k < m; ++k) {//k=1,...,4
      float tmp1=0.0;
      int isHe=0;
      for (int i = 0; i < m; ++i) 
	if (i != j && i != k) {
	  tmp1 += aux.bsum[i]; 
	  isHe += aux.c[i]; 
	}
      

      int cjk = aux.c[j] + aux.c[k];
      fprintf(stderr,"j=%d k=%d RES=%d\n",j,k,j*m+k);
      if(1) {

	if (isHe) 
	  q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]] + tmp1;
	else 
	  q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]]; // all the bases are either j or k
      }else{

	q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk<<8|aux.c[k]];
	if(isHe)
	  q[j*m+k] = q[k*m+j] =  q[k*m+j]+tmp1;

      }

    }

    //set to zero if negative, shoulnd't happen
    for (int k = 0; k != m; ++k) 
      if (q[j*m+k] < 0.0) 
	q[j*m+k] = 0.0;
  }
  //  exit(0);
  return 0;
}
#endif
//this is global this is nasty
errmod_t *mod = NULL;

//err_dep = 0.83
void bam_likes_init(){
  mod = errmod_init(1-ERR_DEP);
}

void bam_likes_destroy(){
  //  fprintf(stderr,"bam_likes_destroy\n");
  if(mod!=NULL)
    errmod_destroy(mod);    
}


void call_bam(chunkyT *chk,double **lk,int trim){
  //  fprintf(stderr,"trim=%d\n",trim);
  if(chk==NULL){
    if(mod!=NULL)
      bam_likes_destroy();
    return;
  }
  assert(mod!=NULL);
  for(int s=0;s<chk->nSites;s++){
    //fprintf(stderr,"posi=%d\n",chk->nd[s][0].refPos);
    lk[s] = new double[10*chk->nSamples];
    for(int ii=0;ii<10*chk->nSamples;ii++)
      lk[s][ii] = -0.0;
    //fprintf(stderr,"chr=%s minbaseq=%d\n",chk->refName,minQ);
    for(int i=0;i<chk->nSamples;i++){

      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;
      uint16_t bases[nd->l];
      int numItems =0;
      for(int j=0;j<nd->l;j++) {
	//	fprintf(stderr,"posi=%d isop=%d +=%d\n",nd->posi[j],nd->isop[j],nd->posi[j]+nd->isop[j]);
	int q=nd->qs[j];
	//	fprintf(stderr,"q=%c\n",q+33);
	//fprintf(stderr,"skipping q=%d mapQ=%d posi=%d isop=%d\n",q,nd.mapQ[j],nd.posi[j],nd.isop[j]);


	extern int minQ;
	if(q<minQ || nd->posi[j]<trim||nd->isop[j]<trim){
	  //fprintf(stderr,"skipping q=%d mapQ=%d posi=%d isop=%d ind=%d\n",q,nd->mapQ[j],nd->posi[j],nd->isop[j],i);
	  continue;
	}


	/*
	  extern int minQ;
	  if(q<minQ || nd.posi[j]<trim||nd.isop[j]<trim){
	  // fprintf(stderr,"skipping q=%d mapQ=%d posi=%d isop=%d ind=%d\n",q,nd.mapQ[j],nd.posi[j],nd.isop[j],i);
	  continue;
	  }
	 */
	if (q>nd->mapQ[j])//dont allow qscores in read to exceed mapping quality
	  q=nd->mapQ[j];
	if (q > 63) q = 63;
	if (q < 4) q = 4;
	//	fprintf(stderr,"q=%d qshift=%d\n",q,q<<5);
	//plug in values
	//	bca->bases[n++] = q<<5 | (int)bam1_strand(p->b)<<4 | b;
	int isLower = islower(nd->seq[j])>0;
	//fprintf(stderr,"strandstuff=%d\tshiftstrandstuff=%d\n",isLower, isLower<<4 );
	bases[numItems++] = q<<5 | isLower<<4 | refToInt[nd->seq[j]];
	//	fprintf(stderr,"bases[%d]=%d\n",j,bases[numItems-1]);
      }
      float likes[25];
      for(int j=0;0&&j<25;j++) 
	likes[j]=0;
      if(numItems!=0) {
	errmod_cal(mod,numItems,5,bases,likes);
	int ord[] = {0,1,2,3,6,7,8,12,13,18};
	for(int gg=0;gg<10;gg++){
	  //fprintf(stderr,"%f\n",likes[ord[gg]]);
	  lk[s][i*10+gg] = -log(10.0) *likes[ord[gg]]/10.0;
	}

      }
    }
    //    break;
  }
}
