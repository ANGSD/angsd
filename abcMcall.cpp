#include "abc.h"
#include "analysisFunction.h"
#include "shared.h"
#include "abcMcall.h"
#include <cfloat>
void abcMcall::printArg(FILE *fp){
  fprintf(fp,"\t-> doMcall=%d\n",domcall);
  
}

double theta = 1.1e-3;
void abcMcall::run(funkyPars *pars){
  int trim=0;
  if(!domcall)
    return;
  float QS_glob[5]={0,0,0,0,0};
  double QS_ind[4][pars->nInd];
  
  chunkyT *chk = pars->chk;
  for(int s=0;s<chk->nSites;s++) {
    fprintf(stderr,"REF: %d\n",pars->ref[s]);
    for(int i=0;i<chk->nSamples;i++) {
      tNode *nd = chk->nd[s][i];
      if(nd==NULL)
	continue;

      for(int j=0;j<4;j++)
	QS_ind[j][i] = 0;
      for(int j=0;j<nd->l;j++){
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
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
      for(int j=0;j<4;j++)
	partsum += QS_ind[j][i];
      for(int j=0;j<4;j++)
	QS_glob[j] += QS_ind[j][i]/partsum;
    }
    for(int i=0;1&&i<5;i++)
      fprintf(stderr,"%d) %f\n",i,QS_glob[i]);

    double liks[10*pars->nInd];//<- this will be the work array


    //switch to PL just to compare numbers
    /*
      p = 10^(-q/10)
      log10(p) = -q/10
      -10*log10(p) = q;
    */
#if 1
    //this macro will discard precision to emulate PL, it can be removed whenever things work
    for(int i=0;i<pars->nInd;i++){
      float min = FLT_MAX;
      for(int j=0;j<10;j++)
	if (min > -10*log10(exp(pars->likes[s][i*10+j])))
	  min = -10*log10(exp(pars->likes[s][i*10+j]));
      for(int j=0;j<10;j++){
	pars->likes[s][i*10+j]  = (int)(-10*log10(exp(pars->likes[s][i*10+j])) - min + .499);
	//fprintf(stderr,"PL[%d][%d]: %f\n",i,j,pars->likes[s][i*10+j]);
      }
      //now pars->likes is in phred PL scale and integerized
      for(int j=0;j<10;j++){
	//	fprintf(stderr,"PRINTER %f\n",pars->likes[s][i*10+j]);
	pars->likes[s][i*10+j] = log(pow(10,-pars->likes[s][10*i+j]/10.0));

      }
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
    for(int i=0;i<5;i++)
      fprintf(stderr,"%d: = %d %f unseen: %d nallele: %d\n",i,calla[i],callqsum[i],callunseen,calln_alleles);
   

    double newlik[10*pars->nInd];
    for(int i=0;i<pars->nInd;i++)
      for(int j=0;j<10;j++)
	newlik[i*10+j] = 0;
  
    for(int i=0;i<pars->nInd;i++)
      for(int ii=0;ii<4;ii++)
	for(int iii=ii;iii<4;iii++){
	  int b1=calla[ii];
	  int b2=calla[iii];
	  // fprintf(stderr,"b1: %d b2: %d\n",b1,b2);
	  if(b1!=-1&&b2!=-1)
	    newlik[i*10+angsd::majorminor[b1][b2] ] = exp(pars->likes[s][i*10+angsd::majorminor[b1][b2]]);
	}
    for(int i=0;0&&i<10*pars->nInd;i++)
      fprintf(stderr,"newlik[%d]: %f\n",i,newlik[i]);
    //exit(0);
    

    //   fprintf(stderr,"pars->nInd: %d\n",pars->nInd);
    for(int i=0;i<pars->nInd;i++){
      double tsum = 0.0;
      for(int j=0;j<10;j++)
	tsum += newlik[i*10+j];
      //    fprintf(stderr,"tsum: %f\n",tsum);
      for(int j=0;j<10;j++){
	//	fprintf(stderr,"PRE[%d][%d]: %f\n",i,j,newlik[i*10+j]);
	liks[i*10+j] = newlik[i*10+j]/tsum;
      }
      for(int j=0;j<10;j++)
	fprintf(stderr,"lk[%d][%d]: %f\n",i,j,liks[10*i+j]);
     
    }
    exit(0);
     exit(0);
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
    fprintf(stderr,"theta: %f\n",theta);
    exit(0);
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
  fprintf(stderr,"asdfadsfadsf\n");
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
