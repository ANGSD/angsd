#include <cassert>
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

void offset2allele(int a,int &b,int &c){
  assert(a>=0&&a<10);
  if(a==0){//AA
    b=0;
    c=0;
  }
  if(a==1){//AC
    b=0;
    c=1;
  }
  if(a==2){//AG
    b=0;
    c=2;
  }
  if(a==3){//AT
    b=0;
    c=3;
  }
  if(a==4){//CC
    b=1;
    c=1;
  }
  if(a==5){//CG
    b=1;
    c=2;
  }
  if(a==6){//CT
    b=1;
    c=3;
  }
  if(a==7){//GG
    b=2;
    c=2;
  }
  if(a==8){//GT
    b=2;
    c=3;
  }
  if(a==9){//TT
    b=3;
    c=3;
  }
}


void abcMcall::run(funkyPars *pars) {
  double theta = 1.1e-3;
  if(!domcall)
    return;

  angsd_mcall *dat = new angsd_mcall;
  dat->quals = new float[2*pars->numSites];
  dat->QS = new float[5*pars->numSites];
  dat->gcdat = new int*[pars->numSites];
  dat->als = new char[4*pars->numSites];
  memset(dat->als,4,4*pars->numSites);
  for(int s=0;s<pars->numSites;s++)
    dat->gcdat[s] = new int[2*pars->nInd];

  //genoCall struct to pars
  pars->extras[index] = dat;
  
  int trim=0;
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
    if(pars->keepSites[s]==0)
      continue;
    //    fprintf(stderr,"SSSSSS: %d\n",s);
    float *QS_glob = dat->QS + 5*s;
    QS_glob[0] = QS_glob[1] = QS_glob[2] = QS_glob[3] = QS_glob[4] = 0.0;
    
    // fprintf(stderr,"\t-> s:%d REF: %d\n",s,pars->ref[s]);
    for(int i=0;i<chk->nSamples;i++) {
      for(int j=0;j<4;j++)
	QS_ind[j][i] = 0;
      tNode *nd = chk->nd[s][i];
      if(nd==NULL){
	//fprintf(stderr,"nodata[%d]\n",i);
	continue;
      }
      //  fprintf(stderr,"nodata[%d].l:%d\n",i,nd->l);
      for(int j=0;j<nd->l;j++){
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
	if(qs>nd->mapQ[j])
	  qs = nd->mapQ[j];
	if(qs<4)
	  qs=4;
	//	fprintf(stderr,"i: %d j: %d mapq: %d allele: %d qs: %d\n",i,j,nd->mapQ[j],allele,qs);
	//filter qscore, mapQ,trimming, and always skip n/N
	if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
	  continue;
	}
	QS_ind[allele][i] += qs;
      }
      for(int j=0;0&&j<4;j++)
	fprintf(stderr,"QSind[allele=%d][ind=%d] %f\n",j,i,QS_ind[j][i]);
    }
    for(int i=0;i<chk->nSamples;i++){
      double partsum = 0;
      for(int j=0;j<4;j++){
	partsum += QS_ind[j][i];
	//	fprintf(stderr,"partsum: %f QS_ind[%d][%d]: %f\n",partsum,j,i,QS_ind[j][i]);
      }
      for(int j=0;(partsum>0)&&j<4;j++){
	QS_glob[j] += QS_ind[j][i]/partsum;
	///	fprintf(stderr,"qs_glob[%d]: %f partsum: %f\n",j,QS_glob[j],partsum);
      }
    }

    double partsum = QS_glob[0]+QS_glob[1]+QS_glob[2]+QS_glob[3]+QS_glob[4];
    for(int i=0;0&&i<5;i++){
      QS_glob[i] = QS_glob[i]/partsum;
      //      fprintf(stderr,"qsum global %d) %f\n",i,QS_glob[i]);
    }
    if(partsum==0){
      dat->quals[s*2] = log(0);
      dat->quals[s*2+1] = log(0);
      //      fprintf(stderr,"Q1: nan Q2: nan\n");
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
	//	fprintf(stderr,"ADSF: %f\n",pars->likes[s][i*10+j]);
	if(pars->likes[s][i*10+j]>255)
	  pars->likes[s][i*10+j]=255;
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
    char hithit[pars->nInd];//<- if set there is missing data
    //this block of code will plug in the relevant gls and rescale to normal
    //double newlik[10*pars->nInd];

    for(int i=0;i<pars->nInd;i++){
      hithit[i] =0;
      for(int j=0;j<10;j++)
	liks[i*10+j] = 0;
      for(int j=0;0&&j<10;j++)
	fprintf(stderr,"RAW: lk[%d][%d]: %f llk: %f\n",i,j,exp(pars->likes[s][i*10+j]),pars->likes[s][i*10+j]);

      int hasdata = 0;
      for(int j=1;j<10;j++){
	//	fprintf(stderr," %f %f diff=%e\n",liks[i*10+j],liks[i*10],liks[i*10+j]-liks[i*10]);
	if(pars->likes[s][i*10+j]!=pars->likes[s][i*10]){
	  //fprintf(stderr,"NOT same\n");
	  hasdata = 1;
	  break;
	}
      }
      if(hasdata==1)
	hithit[i] = 1;
    }
    	 
   
   
   
    for(int i=0;i<pars->nInd;i++){
      for(int ii=0;ii<4;ii++)
	for(int iii=ii;iii<4;iii++){
 	  int b1=calla[ii];
	  int b2=calla[iii];
	  if(b1==-1||b2==-1)
	    continue;
	  
	  liks[i*10+angsd::majorminor[b1][b2] ] = exp(pars->likes[s][i*10+angsd::majorminor[b1][b2]]);
	}
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
      for(int i=0;i<pars->nInd;i++){
	//	fprintf(stderr,"b: %d i: %d lik:%f llh: %f\n",b,i,liks[10*i+angsd::majorminor[b][b]],log(liks[10*i+angsd::majorminor[b][b]]));
	llh += log(liks[10*i+angsd::majorminor[b][b]]);
      }
      //   fprintf(stderr,"CALLING MOHO llh[%d]: llh: %f\n",b,llh);
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
      //      fprintf(stderr,"ALLALE: %d llh: %f ref: %d\n",i,monollh[i],pars->ref[s]);
      if(monollh[i]!=HUGE_VAL){

	if(i!=pars->ref[s])
	  monollh[i] += theta;
	char *val = new char[4];
	memset(val,'N',4);
	val[0]=i+'0';
	llh_als[monollh[i]] = val;
	//fprintf(stderr,"ALLALE: %d llh: %f ref: %d\n",i,monollh[i],pars->ref[s]);
      }
    }
    //    fprintf(stderr,"SIZE: %lu\n",llh_als.size());//exit(0);
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
    
    double totlik1=log(0);//this value DOES NOT CONTAIN llh of the ref (to avoid underflow)
    // fprintf(stderr,"llh_als.size(): %lu\n",llh_als.size());
    char *DAS_BEST=NULL;
    for(std::map<double,char*>::reverse_iterator it=llh_als.rbegin();it!=llh_als.rend();it++){
      if(DAS_BEST==NULL)
	DAS_BEST = it->second;
      //   fprintf(stderr,"%s %f\n",it->second,it->first);
      if(it->second[0]==pars->ref[s]+'0'&&it->second[1]=='N'){
	//	fprintf(stderr,"skipping: %f\n",it->first);
	continue;
      }
      totlik1 = logsumexp2(totlik1,it->first);
      //fprintf(stderr,"totlik1: %f\n",totlik1);
    }
    if(DAS_BEST==NULL){
      fprintf(stderr,"Couldnt find best alleleic configurartion");
      exit(0);
    }
 
    double ref_llh =  monollh[pars->ref[s]];

    //this is abit strange we shouldnt put the ref_llh in the denominator the ratio
    double Q1 =  -4.343*(ref_llh - logsumexp2(totlik1,ref_llh));
    double Q2 =  -4.343*(totlik1 - logsumexp2(totlik1,ref_llh));
    dat->quals[2*s] = Q1;
    dat->quals[2*s+1] = Q2;
    //    fprintf(stderr,"Q1: %f Q2: %f\n",Q1,Q2);
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
     if(DAS_BEST[2]!='N'){
       fprintf(stderr,"Multi allelic not tested for >2alleles will skip site\n");
       pars->keepSites[s] = 0;
     }
     
     
     //     fprintf(stderr,"CALLING genotypes with BEST alleles: %s\n",DAS_BEST);
     for(int i=0;i<4;i++)
       DAS_BEST[i] = refToInt[DAS_BEST[i]];
#if 1
     int refpos = -1;
     int lastmissing=-1;
     for(int i=0;i<4;i++){
       DAS_BEST[i] = refToInt[DAS_BEST[i]];
       if(DAS_BEST[i]==pars->ref[s])
	 refpos = i;
       if(DAS_BEST[i]==4){
	 lastmissing = i;
	 break;
       }
     }
     if(refpos==-1){
       fprintf(stderr,"\t-> Big problem reference not among alelleset: last missing is: %d\n",lastmissing);
       for(int i=lastmissing;i>=1;i--){
	 //	 fprintf(stderr,"i: %d\n",i);
	 DAS_BEST[i] = DAS_BEST[i-1];
       }
       DAS_BEST[0] = pars->ref[s];
     }else if(refpos>0){
       fprintf(stderr,"\t-> Big problem reference not the most frequently occuring allele last missing is: %d refpos: %d\n",lastmissing,refpos);
       int tmptmp = DAS_BEST[0];
       DAS_BEST[0] = DAS_BEST[refpos];
       DAS_BEST[refpos] = tmptmp;
     }
#endif
     for(int jj=0;jj<4;jj++){
       dat->als[s*4+jj] = refToInt[DAS_BEST[jj]];
       //       fprintf(stderr,"YOYO dat->als[%d]: %d\n",jj,dat->als[s*4+jj]);
       //fprintf(stderr,"ASDFASDF\n");
     }
     //continue;
     partsum = QS_glob[0]+QS_glob[1]+QS_glob[2]+QS_glob[3]+QS_glob[4];
     float FREQ[5];
     for(int i=0;1&&i<5;i++){
       FREQ[i] = QS_glob[i]/partsum;
       // fprintf(stderr,"FREQ qsum global %d) %f\n",i,FREQ[i]);
     }
     double gc_gls[10*pars->nInd];//this will contain a copy and pp of the gls perind
     double gc_llh[10*pars->nInd];//this will contain the llh of the different genotypes
   
     for(int i=0;i<pars->nInd;i++) {
       //   fprintf(stderr,"hit[%d]: %d\n",i,hithit[i]);
       if(hithit[i]==0){
	 //	 fprintf(stderr,"\t-> Skipping hit[%d]: %d\n",i,hithit[i]);
	 dat->gcdat[s][2*i] =dat->gcdat[s][2*i+1]= -1;
	 continue;
       }


       for(int j=0;j<10;j++){
	 gc_gls[i*10+j]  = gc_llh[i*10+j] = 0;
       }
      
       //these two for loops just copies over the relevant gls and normalizes
       double tsum = 0;
       for(int a=0;a<4;a++){
	 int b1 = refToInt[DAS_BEST[a]];
	 if(b1==4)
	   continue;
	 for(int b=a;b<4;b++){
	   int b2 = refToInt[DAS_BEST[b]];
	   if(b2==4)
	     continue;
	   //   fprintf(stderr,"YOYO b1: %d b2: %d\n",b1,b2);
	   tsum += gc_gls[i*10+angsd::majorminor[b1][b2]] = liks[i*10+angsd::majorminor[b1][b2]];
	 }
       }
       //    fprintf(stderr,"tsum: %f\n",tsum);
       for(int j=0;0&&j<10;j++){
	 gc_gls[i*10+j]  /= tsum;
	 fprintf(stderr,"gc_gls[%d][%d]: %f\n",i,j,gc_gls[i*10+j]);
       }
       /*
	 these for loops just calculates the llh for the different diploid genotype configurations.
	 The prior is the allele frequency estimated from the sum of qscores
	 this should be updated to take into account deviations from hwe
       */
       char called[2] = {0,0};
       double maxgcllh = -1;
       for(int a=0;a<4;a++){
	 int b1 = refToInt[DAS_BEST[a]];
	 if(b1==4)
	   continue;
	 for(int b=a;b<4;b++){
	   int b2 = refToInt[DAS_BEST[b]];

	   if(b2==4)
	     continue;
	   //  fprintf(stderr,"b1: %d b2: %d\n",b1,b2);
	   int offs = angsd::majorminor[b1][b2];

	   
	   if(b1==b2)
	     gc_llh[10*i+offs] = gc_gls[10*i+offs] * FREQ[b1] * FREQ[b1];
	   else
	     gc_llh[10*i+offs] = gc_gls[10*i+offs] * 2 * FREQ[b1] * FREQ[b2];

	   double tmpgcllh = gc_llh[10*i+offs];
	   if(tmpgcllh>maxgcllh){
	     maxgcllh = tmpgcllh;
	     called[0] = b1;
	     called[1] = b2;
	   }
	   //  fprintf(stderr,"QS_glob[%d]: %f QS_glob[%d]: %f offs: %d llh: %f gc_gls: %f\n",b1,FREQ[b1],b2,FREQ[b2],offs,gc_llh[10*i+offs],gc_gls[10*i+offs]);
	 }
       }
       dat->gcdat[s][2*i] = called[0];
       dat->gcdat[s][2*i+1] = called[1];
#if 0 
       for(int j=0;1&&j<10;j++)
	 fprintf(stderr,"gc_llh[%d][%d]: %f\n",i,j,gc_llh[i*10+j]);
      
       int whichmax = 0;
       partsum = gc_llh[i*10+whichmax];
       for(int j=1;1&&j<10;j++){
	 //fprintf(stderr,"gc_llh[%d]: %f comparing with: %f \n",j,gc_llh[i*10+j],gc_llh[i*10+whichmax]);
	 partsum += gc_llh[i*10+j];
	 if(gc_llh[i*10+j]>gc_llh[i*10+whichmax])
	   whichmax = j;
       }
       int a1,a2;
       offset2allele(whichmax,a1,a2);
       int b1=refToInt[DAS_BEST[0]];
       int b2=refToInt[DAS_BEST[1]];
       //    fprintf(stderr,"offset2allele: which=%d a1=%d a2=%d DAS: b1=%d b2=%d lk: %f nd->l: %d\n",whichmax,a1,a2,b1,b2,gc_llh[i*10+whichmax],pars->chk->nd[s][i]?pars->chk->nd[s][i]->l:0);
       if(b2!=4){
	 assert(a1==b1||a1==b2);
	 assert(a2==b1||a2==b2);
       }
       dat->gcdat[s][i] = -1;
       if(pars->chk->nd[s][i]&&pars->chk->nd[s][i]->l>0){
	 if(pars->ref[s]==b1||pars->ref[s]==b2){
	   if(a1==a2&&a1==b1)
	     dat->gcdat[s][i] = 0;
	   if(a1==a2&&a1==b2)
	     dat->gcdat[s][i] = 2;
	   
	   if(a1==b1&&a2==b2||a1==b2&&a2==b1)
	     dat->gcdat[s][i] = 1;
	 }else{
	   if(a1==a2&&pars->ref[s]!=a1)
	     dat->gcdat[s][i] = 2;
	 }
	 
       //   fprintf(stderr,"whichmax: %d a1: %d a2: %d pp: %f das_best: %s gcdat: %d\n",whichmax,a1,a2,gc_llh[i*10+whichmax]/partsum, DAS_BEST,dat->gcdat[s][i]);
       }
#endif
     }

  }
}


void abcMcall::clean(funkyPars *fp){
  if(!domcall)
    return;

  angsd_mcall *dat =(angsd_mcall *) fp->extras[index];
  delete [] dat->quals;
  delete [] dat->QS;
  for(int s=0;s<fp->numSites;s++)
    delete [] dat->gcdat[s];
  delete [] dat->gcdat;
  delete [] dat->als;
  delete dat;
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
  int gl = 0;
  gl=angsd::getArg("-gl",gl,arguments);
  if(gl==0){
    fprintf(stderr,"\t-> Error, you must supply a genotype likelihood model for -doMcall\n");
    exit(0);
  }
  char *ref = NULL;
  ref=angsd::getArg("-ref",ref,arguments);
  if(ref==NULL){
    fprintf(stderr,"\t-> Error, mcall requires reference (-ref)\n");
    exit(0);
  }
  printArg(arguments->argumentFile);

}


abcMcall::abcMcall(const char *outfiles,argStruct *arguments,int inputtype){
  getOptions(arguments);
  printArg(arguments->argumentFile);
}

abcMcall::~abcMcall(){


}
