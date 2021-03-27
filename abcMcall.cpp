#include "abc.h"
#include "analysisFunction.h"
#include "shared.h"
#include "abcMcall.h"
#include <cfloat>

/** log(exp(a)+exp(b)) */
static inline double logsumexp2(double a, double b)
{
    if ( a>b )
        return log(1 + exp(b-a)) + a;
    else
        return log(1 + exp(a-b)) + b;
}


void abcMcall::printArg(FILE *fp){
  fprintf(fp,"\t-> doMcall=%d\n",domcall);
  
}

double theta = 1.1e-3;
void abcMcall::run(funkyPars *pars){
  int trim=0;
  if(!domcall)
    return;

  // Watterson factor, here aM_1 = aM_2 = 1
  double aM = 1;
  for (int i=2; i<pars->nInd*2; i++) aM += 1./i;
  theta *= aM;
  if ( theta >= 1 )
    {
      fprintf(stderr,"The prior is too big (theta*aM=%.2f), going with 0.99\n", theta);
      theta = 0.99;
    }
  theta = log(theta);
  //fprintf(stderr,"theta: %f\n",theta);
  
  
  double QS_ind[4][pars->nInd];
  
  chunkyT *chk = pars->chk;
  for(int s=0;s<chk->nSites;s++) {
    float QS_glob[5]={0,0,0,0,0};
    // fprintf(stderr,"\t-> s:%d REF: %d\n",s,pars->ref[s]);
    for(int i=0;i<chk->nSamples;i++) {
      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;

      for(int j=0;j<4;j++)
	QS_ind[j][i] = 0;
      for(int j=0;j<nd->l;j++){
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
	if(qs>nd->mapQ[j])
	  qs = nd->mapQ[j];
	//filter qscore, mapQ,trimming, and always skip n/N
	if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
	  continue;
	}
	QS_ind[allele][i] += qs;
      }
      for(int j=0;0&&j<4;j++)
	fprintf(stderr,"QS[%d][%d] %f\n",i,j,QS_ind[j][i]);
    }
    for(int i=0;i<chk->nSamples;i++){
      double partsum = 0;
      for(int j=0;j<4;j++){
	partsum += QS_ind[j][i];
	//	fprintf(stderr,"partsum: %f QS_ind[%d][%d]: %f\n",partsum,j,i,QS_ind[j][i]);
      }
      for(int j=0;(partsum>0)&&j<4;j++){
	QS_glob[j] += QS_ind[j][i]/partsum;
	//	fprintf(stderr,"qs_glob[%d]: %f partsum: %f\n",j,QS_glob[j],partsum);
      }
    }

    double partsum = QS_glob[0]+QS_glob[1]+QS_glob[2]+QS_glob[3]+QS_glob[4];
    for(int i=0;1&&i<5;i++){
      QS_glob[i] = QS_glob[i]/partsum;
      //    fprintf(stderr,"qsum global %d) %f\n",i,QS_glob[i]);
    }
    if(partsum==0){
      fprintf(stderr,"Q1: nan Q2: nan\n");
      continue;
    }
      
    //   exit(0);
  
#if 1
    //this macro will discard precision to emulate PL, it can be removed whenever things work
    for(int i=0;i<pars->nInd;i++){
      float min = FLT_MAX;
      for(int j=0;j<10;j++)
	if (min > -10*log10(exp(pars->likes[s][i*10+j])))
	  min = -10*log10(exp(pars->likes[s][i*10+j]));
      for(int j=0;j<10;j++){
	pars->likes[s][i*10+j]  = (int)(-10*log10(exp(pars->likes[s][i*10+j])) - min + .499);
      }
      //now pars->likes is in phred PL scale and integerized
      for(int j=0;j<10;j++)
	pars->likes[s][i*10+j] = log(pow(10,-pars->likes[s][10*i+j]/10.0));

      //now pars->likes is in logscale but we have had loss of praecision
    }
#endif

    // sort qsum in ascending order (insertion sort), copied from bam2bcf.c bcftools 21march 2021
    float *ptr[5], *tmp;
    for (int i=0; i<5; i++)
      ptr[i] = &QS_glob[i];
    for (int i=1; i<4; i++)
      for (int j=i; j>0 && *ptr[j] < *ptr[j-1]; j--)
	tmp = ptr[j], ptr[j] = ptr[j-1], ptr[j-1] = tmp;


 // Set the reference allele and alternative allele(s)
    int calla[5];
    float callqsum[5];
    int callunseen;
    int calln_alleles;
    for (int i=0; i<5; i++)
      calla[i] = -1;
    for (int i=0; i<5; i++)
      callqsum[i] = 0;
    callunseen = -1;
    calla[0] = pars->ref[s];
    int i,j;
    for (i=3, j=1; i>=0; i--)   // i: alleles sorted by QS; j, a[j]: output allele ordering
    {
        int ipos = ptr[i] - QS_glob;   // position in sorted qsum array
        if ( ipos==pars->ref[s] )
            callqsum[0] = QS_glob[ipos];    // REF's qsum
        else
        {
            if ( !QS_glob[ipos] ) break;       // qsum is 0, this and consequent alleles are not seen in the pileup
            callqsum[j] = QS_glob[ipos];
            calla[j++]  = ipos;
        }
    }
    if (pars->ref[s] >= 0)
    {
        // for SNPs, find the "unseen" base
        if (((pars->ref[s] < 4 && j < 4) || (pars->ref[s] == 4 && j < 5)) && i >= 0)
            callunseen = j, calla[j++] = ptr[i] - QS_glob;
        calln_alleles = j;
    }
    else
    {
      fprintf(stderr,"NEVER HERE\n");exit(0);
        calln_alleles = j;
        if (calln_alleles == 1) exit(0);//return -1; // no reliable supporting read. stop doing anything
    }
    for(int i=0;0&&i<5;i++)
      fprintf(stderr,"%d: = %d %f unseen: %d nallele: %d\n",i,calla[i],callqsum[i],callunseen,calln_alleles);
   
    double liks[10*pars->nInd];//<- this will be the work array
    //this block of code will plug in the relevant gls and rescale to normal
    //double newlik[10*pars->nInd];
    for(int i=0;i<pars->nInd;i++)
      for(int j=0;j<10;j++)
	liks[i*10+j] = 0;
      for(int i=0;i<pars->nInd;i++)
      for(int ii=0;ii<4;ii++)
	for(int iii=ii;iii<4;iii++){
 	  int b1=calla[ii];
	  int b2=calla[iii];
	  if(b1==-1||b2==-1)
	    continue;
	  // fprintf(stderr,"b1: %d b2: %d\n",b1,b2);
	  if(b1!=-1&&b2!=-1)
	    liks[i*10+angsd::majorminor[b1][b2] ] = exp(pars->likes[s][i*10+angsd::majorminor[b1][b2]]);
	}
      for(int i=0;i<pars->nInd;i++){
	double tsum = 0.0;
	for(int j=0;j<10;j++)
	  tsum += liks[i*10+j];
	//	fprintf(stderr,"tsum: %f\n",tsum);
	for(int j=0;j<10;j++)
	  liks[i*10+j] = liks[i*10+j]/tsum;

	for(int j=0;0&&j<10;j++)
	  fprintf(stderr,"lk[%d][%d]: %f\n",i,j,liks[10*i+j]);
      }

    
  
    //monomophic
    double monollh[4]={HUGE_VAL,HUGE_VAL,HUGE_VAL,HUGE_VAL};
    for(int b1=0;b1<4;b1++){
      int b = calla[b1];
      if(b==-1)
	continue;
      double llh =0;
      for(int i=0;i<pars->nInd;i++)
	llh += log(liks[10*i+angsd::majorminor[b][b]]);
      //      fprintf(stderr,"CALLING MOHO llh[%d]: llh: %f\n",b,llh);
      monollh[b] = llh;
    }
    //diallelic
    double dillh[4][4];
    for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	dillh[i][ii] = HUGE_VAL;
    //    fprintf(stderr,"calling diallelic\n");
    for(int a=0;a<4;a++){
      for(int b=a+1;b<4;b++){
	int b1=calla[a];
	int b2=calla[b];
	if(b1==-1||b2==-1)
	  continue;
	if(QS_glob[b1]==0||QS_glob[b2]==0)
	  continue;
	double tot_lik = 0;
	//	fprintf(stderr,"CALL %d %d QS_glob: %e %e %e %e\n",b1,b2,QS_glob[0],QS_glob[1],QS_glob[2],QS_glob[3]);
	double fa  = QS_glob[b1]/(QS_glob[b1] + QS_glob[b2]);
	double fb  = QS_glob[b2]/(QS_glob[b1] + QS_glob[b2]);
	double fa2 = fa*fa;
	double fb2 = fb*fb;
	double fab = 2*fa*fb;
	double val = 0;
	//	fprintf(stderr,"fa: %f fb: %f fa2: %f fb2: %f fab: %f\n",fa,fb,fa2,fb2,fab);
	for(int i=0;i<pars->nInd;i++){
	  val= fa2*liks[i*10+angsd::majorminor[b1][b1]] + fab*liks[i*10+angsd::majorminor[b1][b2]] + fb2*liks[i*10+angsd::majorminor[b2][b2]];
	  //	  fprintf(stderr,"ABC val: %f\n",val);
	  tot_lik += log(val);
	}

	dillh[b1][b2] = tot_lik;
	//	fprintf(stderr,"DIALLELIC b1: %d b2: %d tot_lik: %f\n",b1,b2,tot_lik);
      }
    }
  
    //triallelic
    
    double trillh[4][4][4];
    for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	for(int iii=0;iii<4;iii++)
	trillh[i][ii][iii] = HUGE_VAL;
    //    fprintf(stderr,"calling triallelic\n");
    for(int a=0;a<4;a++){
      for(int b=a+1;b<4;b++){
	for(int c=b+1;c<4;c++){
	  int b1=calla[a];
	  int b2=calla[b];
	  int b3=calla[c];
	if(b1==-1||b2==-1||b3==-1)
	  continue;
	if(QS_glob[b1]==0||QS_glob[b2]==0||QS_glob[b3]==0)
	  continue;
	double tot_lik = 0;
	//	fprintf(stderr,"CALL %d %d %d QS_glob: %e %e %e %e\n",b1,b2,b3,QS_glob[0],QS_glob[1],QS_glob[2],QS_glob[3]);
	double fa  = QS_glob[b1]/(QS_glob[b1] + QS_glob[b2]+QS_glob[b3]);
	double fb  = QS_glob[b2]/(QS_glob[b1] + QS_glob[b2]+QS_glob[b3]);
	double fc  = QS_glob[b3]/(QS_glob[b1] + QS_glob[b2]+QS_glob[b3]);
	double fa2 = fa*fa;
	double fb2 = fb*fb;
	double fc2 = fc*fc;
	double fab = 2*fa*fb,fac=2*fa*fc,fbc=2*fb*fc;
	double val = 0;
	//	fprintf(stderr,"fa: %f fb: %f fa2: %f fb2: %f fab: %f fc %f fc2 %f fbc %f\n",fa,fb,fa2,fb2,fab,fc,fc2,fbc);
	for(int i=0;i<pars->nInd;i++){
	  val= fa2*liks[i*10+angsd::majorminor[b1][b1]] + fab*liks[i*10+angsd::majorminor[b1][b2]]+ fac*liks[i*10+angsd::majorminor[b1][b3]]+ fbc*liks[i*10+angsd::majorminor[b2][b3]] + fb2*liks[i*10+angsd::majorminor[b2][b2]]+ fc2*liks[i*10+angsd::majorminor[b3][b3]];
	  //  fprintf(stderr,"val: %f\n",val);
	  tot_lik += log(val);
	}
	//	fprintf(stderr,"trIALLELIC b1: %d b2: %d  b3:%dtot_lik: %f\n",b1,b2,b3,tot_lik);
	trillh[b1][b2][b3] = tot_lik;
	}
      }
    }
    std::map<double,char*> llh_als;
    
    for(int i=0;i<4;i++){
      if(monollh[i]!=HUGE_VAL){
	if(i!=pars->ref[s])
	  monollh[i] += theta;
	char *val = new char[4];
	memset(val,'N',4);
	val[0]=i+'0';
	llh_als[monollh[i]] = val;
      }
    }
    
    for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	if(dillh[i][ii]!=HUGE_VAL){
	  if(i!=pars->ref[s])
	    dillh[i][ii] += theta;
	  if(ii!=pars->ref[s])
	    dillh[i][ii] += theta;

	  char *val = new char[4];
	  memset(val,'N',4);
	  val[0]=i+'0';
	  val[1]=ii+'0';
	  llh_als[dillh[i][ii]] = val;
	}
    
    
    for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	for(int iii=0;iii<4;iii++)
	  if(dillh[i][ii]!=HUGE_VAL){
	    if(i!=pars->ref[s])
	      trillh[i][ii][iii] += theta;
	    if(ii!=pars->ref[s])
	      trillh[i][ii][iii] += theta;
	    if(iii!=pars->ref[s])
	      trillh[i][ii][iii] += theta;
	      char *val = new char[4];
	      memset(val,'N',4);
	      val[0]=i+'0';
	      val[1]=ii+'0';
	      val[2]=iii+'0';
	      if(trillh[i][ii][iii]!=HUGE_VAL)
		llh_als[trillh[i][ii][iii]] = val;
	  }

    double totlik1=log(0);//this value DOES NOT CONTAIN llh of the ref (to avoid underlow)
    int isvar = 0;
    for(std::map<double,char*>::reverse_iterator it=llh_als.rbegin();it!=llh_als.rend();it++){
      //fprintf(stderr,"%s %f\n",it->second,it->first);
      if(it->second[1]!='N')
	isvar++;
      if(it->second[0]==pars->ref[s]+'0'&&it->second[1]=='N'){
	//	fprintf(stderr,"skipping: %f\n",it->first);
	continue;
      }
      totlik1 = logsumexp2(totlik1,it->first);
      //fprintf(stderr,"totlik1: %f\n",totlik1);
    }
     
    double ref_llh =  monollh[pars->ref[s]];
    // fprintf(stderr,"ref_llh: %f totlik1: %f isvar: %d \n",ref_llh,totlik1,isvar);

    //this is abit strange we shouldnt put the ref_llh in the denominator the ratio
    double Q1 =  -4.343*(ref_llh - logsumexp2(totlik1,ref_llh));
    double Q2 =  -4.343*(totlik1 - logsumexp2(totlik1,ref_llh));
    fprintf(stderr,"Q1: %f Q2: %f\n",Q1,Q2);
#if 0
    for(int i=0;i<4;i++)
      if(monollh[i]!=HUGE_VAL)
	fprintf(stderr,"mono[%d] llh: %f\n",i,monollh[i]);
    for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	if(dillh[i][ii]!=HUGE_VAL)
	  fprintf(stderr,"di[%d][%d] llh: %f\n",i,ii,dillh[i][ii]);
     for(int i=0;i<4;i++)
      for(int ii=0;ii<4;ii++)
	for(int iii=0;iii<4;iii++)
	  if(trillh[i][ii][iii]!=HUGE_VAL)
	    fprintf(stderr,"tri[%d][%d][%d] llh: %f\n",i,ii,iii,trillh[i][ii][iii]);
#endif 
     //  exit(0);
  }
}


void abcMcall::clean(funkyPars *fp){
  if(!domcall)
    return;

  
}

void abcMcall::print(funkyPars *fp){
  if(!domcall)
    return;
      
}


void abcMcall::getOptions(argStruct *arguments){
  //  fprintf(stderr,"asdfadsfadsf\n");
  //default
  domcall=0;

  //from command line
  domcall=angsd::getArg("-domcall",domcall,arguments);
  if(domcall==-999){
    domcall=0;
    printArg(stderr);
    exit(0);
  }
  if(domcall==0)
    return;

  printArg(arguments->argumentFile);

}


abcMcall::abcMcall(const char *outfiles,argStruct *arguments,int inputtype){
  getOptions(arguments);
  printArg(arguments->argumentFile);
}

abcMcall::~abcMcall(){


}
