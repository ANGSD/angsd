/*
  little class that does 
1) hwe using genotype likelihoods
2) a) fisher exact
   b) GATK approach
   c) SB
   These are based on guo 2012 BMC genomics 2013 13:666
3) 2 sample wilcox rank sum for baseQ bias
   This file should be completely rewritten, much to slow and stupid. But it should give correct results

 */


#include <cmath>
#include <ctype.h>
#include "analysisFunction.h"
#include "shared.h"
#include "fet.c"
#include "chisquare.h" //<- stuff related to calculating pvalues from LRT tests
#include "abc.h"
#include "abcFilterSNP.h"
#include "abcHWE.h"
#include "abcCallGenotypes.h"
#include <htslib/kstring.h>

/*
  wilcox manwhitney u test, whatever ,implementation is from wiki
  validated with wilcox.test (qscores of major,qscores of minor,exact=T,correct=F)
  returns Zscore, for some reason...
 */

double mann(int majD[255],int minD[255]){
  double U[255];
  double cumLess1 =0;
  for(int i=0;i<255;i++){
    U[i]=(double) minD[i]*(majD[i]/2.0+cumLess1);
    cumLess1 += majD[i];
  }
  double U1 =0;
  for(int i=0;i<255;i++) 
    U1+= U[i];
  //below is a check
  double cumLess2 =0;
  for(int i=0;i<255;i++){
    U[i]=(double) majD[i]*(minD[i]/2.0+cumLess2);
    cumLess2 += minD[i];
  }
  double U2 =0;
  for(int i=0;i<255;i++) 
    U2+= U[i];
  

  double mu=cumLess1*cumLess2/2.0;
  double sigu=sqrt((cumLess1*cumLess2*(cumLess1+cumLess2+1))/12.0);
  double Z=(std::min(U1,U2)-mu)/sigu;

 
  //  fprintf(stderr,"U1:%f U2:%f U1+U2:%f nObs:%f Z:%f p.value:%f\n",U1,U2,U1+U2,cumLess1*cumLess2,Z,2*phi(Z));

  return Z;

}



double edgebias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[std::min(nd->posi[l],nd->isop[l])]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[std::min(nd->isop[l],nd->posi[l])]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}



double mapQbias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[nd->mapQ[l]]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[nd->mapQ[l]]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}



double baseQbias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[255];
  int minD[255];
  memset(majD,0,sizeof(int)*255);
  memset(minD,0,sizeof(int)*255);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      if(obB==maj){
	majD[nd->qs[l]]++;
	//	fprintf(stdout,"maj\t%d\n",nd->qs[l]);
      }else if(obB==min){
	minD[nd->qs[l]]++;
	//fprintf(stdout,"min\t%d\n",nd->qs[l]);
      }
    }
  }

  return mann(majD,minD);
}

Chisqdist chi(1);

//guo 2012 mutat res 2012, 744(2):154-160
double sb1(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double top=b/(a+b)-d/(c+d);
  double bot=(b+d)/(a+b+c+d);
  return top/bot;
}

//the gatk way
double sb2(int cnts[4]){
  double a=cnts[0];double b=cnts[1];double c=cnts[2];double d=cnts[3];
  double en=(b/(a+b))*(c/(c+d));
  double to=(a+c)/(a+b+c+d);
  double tre=(d/(c+d))*(a/(a+b));
  return std::max(en/to,tre/to);
}

//strandbias using fisher

double sb3(int cnts[4]){

  double left,right,twotail,prob;
  kt_fisher_exact(cnts[0], cnts[1], cnts[2], cnts[3], &left, &right, &twotail);
  return twotail;
}


void abcFilterSNP::printArg(FILE *argFile){
   fprintf(argFile,"-----BETA---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSnpStat %d\n",doSnpStat);
  fprintf(argFile,"\t-edge_pval %f\n",edge_pval);
  fprintf(argFile,"\t-mapQ_pval %f\n",mapQ_pval);
  fprintf(argFile,"\t-sb_pval %f\n",sb_pval);
  fprintf(argFile,"\t-hwe_pval %f\n",hwe_pval);
  fprintf(argFile,"\t-qscore_pval %f\n",qscore_pval);
  fprintf(argFile,"\t-hetbias_pval %f\n",hetbias_pval);

}
void abcFilterSNP::run(funkyPars *pars){
  if(!doSnpStat)
    return;
  chunkyT *chk = pars->chk;
  
  if(doSnpStat==1){
    kstring_t *bufstr = new kstring_t;
    bufstr->s=NULL;bufstr->l=bufstr->m=0;
    //loop over sites;
    kstring_t persite;
    persite.s=NULL;persite.l=persite.m=0;
    //pull 
   
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue;
      //loop over samples
      int cnts[4]={0,0,0,0};
      for(int i=0;i<pars->nInd;i++){
	tNode *nd = chk->nd[s][i];
	if(nd==NULL)
	  continue;
	for(int l=0;l<nd->l;l++){
	  int obB = refToInt[nd->seq[l]];
	  //	    fprintf(stderr,"%c ",nd.seq[l]);
	  int strand = (isupper(nd->seq[l])==0)<<1;
	  //  fprintf(stderr,"strand:%d\n",strand);
	  if(obB==4)
	    continue;
	  if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	    continue;
	  if(obB!=pars->major[s])
	    strand +=1;
	  //fprintf(stderr,"strand=%d\n",strand);
	  cnts[strand]++;
	}
      }
      ksprintf(&persite,"%s\t%d\t%d %d %d %d\t",header->target_name[pars->refId],pars->posi[s]+1, cnts[0],cnts[1],cnts[2],cnts[3]);
      ksprintf(&persite,"%f:%f:%f\t",sb1(cnts),sb2(cnts),sb3(cnts));
      //      fprintf(stderr,"sb4:%f\n",sb3(cnts));
      if(sb_pval!=-1 && sb3(cnts)<sb_pval)
	pars->keepSites[s] = 0;

      funkyHWE *hweStruct = (funkyHWE *) pars->extras[8];//THIS IS VERY NASTY! the ordering within general.cpp is now important
      double lrt = 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
      double pval;
      if(std::isnan(lrt))
	pval=lrt;
      else if(lrt<0)
	pval =1;
      else
	pval =1- chi.cdf(lrt);
      ksprintf(&persite,"%f:%e\t",lrt,pval);
      //   ksprintf(&persite,"%f:%e\t",lrt,pval);
      if(hwe_pval!=-1 && pval<hwe_pval)
	pars->keepSites[s] = 0;

      double Z = baseQbias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e\t",Z,2*phi(Z));
      if(qscore_pval!=-1 && 2*phi(Z)<qscore_pval)
	pars->keepSites[s] = 0;
    

      Z = mapQbias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e\t",Z,2*phi(Z));

      if(mapQ_pval!=-1 && 2*phi(Z)<mapQ_pval)
	pars->keepSites[s] = 0;

      Z = edgebias(chk->nd[s],pars->nInd,refToInt[pars->major[s]],refToInt[pars->minor[s]]);
      ksprintf(&persite,"%f:%e",Z,2*phi(Z));
      if(2*phi(Z)<edge_pval)
	pars->keepSites[s] = 0;

      genoCalls *gcw =(genoCalls *) pars->extras[10];
      int **gc=NULL;
      if(gcw)
	gc = gcw->dat;
      if(gc){
	cnts[0]=cnts[1]=cnts[2]=cnts[3]=0;

	int nsampleswithdata =0;
	for(int i=0;i<pars->nInd;i++){
	  if(gc[s][i]!=1)
	    continue;
	  tNode *nd = chk->nd[s][i];
	  if(nd==NULL)
	    continue;
	  nsampleswithdata++;
	  for(int l=0;l<nd->l;l++){
	    int obB = refToInt[nd->seq[l]];
	    //	    fprintf(stderr,"%c ",nd.seq[l]);
	    int strand = (isupper(nd->seq[l])==0)<<1;
	    //  fprintf(stderr,"strand:%d\n",strand);
	    if(obB==4)
	      continue;
	    if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	      continue;
	    if(obB!=pars->major[s])
	      strand +=1;
	    //fprintf(stderr,"strand=%d\n",strand);
	    cnts[strand]++;
	  }
	}

	double tsum= cnts[0]+cnts[1]+cnts[2]+cnts[3];
	double fA=(cnts[0]+cnts[2])/tsum;
	double fa=(cnts[1]+cnts[3])/tsum;
	int n = tsum;
	ksprintf(&persite,"\t%d %d %d %d %d\t",cnts[0],cnts[1],cnts[2],cnts[3],n);

	double lrt = 2*tsum*(fA-0.5)*(fA-0.5) + 2*tsum*(fa-0.5)*(fa-0.5);
	double pval;
	if(std::isnan(lrt))
	  pval=lrt;
	else if(lrt<0)
	  pval =1;
	else
	  pval =1- chi.cdf(lrt);
	ksprintf(&persite,"%f:%e\t",lrt,pval);
	if(hetbias_pval!=-1&&pval<hetbias_pval)
	  pars->keepSites[s]=0;
      }

      
      if(pars->keepSites[s]!=0&&persite.l>0){
	ksprintf(&persite,"\n");
	ksprintf(bufstr,"%s",persite.s);
      }
      persite.l=0;
   }
    free(persite.s);
    pars->extras[index] = bufstr;
  }
}

void abcFilterSNP::clean(funkyPars *fp){
  if(!doSnpStat)
    return;

  
}

void abcFilterSNP::print(funkyPars *pars){
  if(!doSnpStat)
    return;
  kstring_t *bufstr =(kstring_t*) pars->extras[index];
  aio::bgzf_write(outfileZ,bufstr->s,bufstr->l);bufstr->l=0;
  free(bufstr->s);
  delete bufstr;
}


void abcFilterSNP::getOptions(argStruct *arguments){
  //default


  //from command line
  doSnpStat=angsd::getArg("-doSnpStat",doSnpStat,arguments);


  if(doSnpStat==0)
    return;
  int domajorminor=0;
  domajorminor = angsd::getArg("-domajorminor",domajorminor,arguments);
  if(domajorminor==0){
    fprintf(stderr,"\t-> Must supply -doMajorMinor for running dosnpstat (needs to look a distributions of major and minor alleles)\n");
    exit(0);
  }
  //from command line

  edge_pval=angsd::getArg("-edge_pval",edge_pval,arguments);      
  mapQ_pval=angsd::getArg("-mapQ_pval",mapQ_pval,arguments);      
  sb_pval=angsd::getArg("-sb_pval",sb_pval,arguments);    
  hwe_pval=angsd::getArg("-hwe_pval",hwe_pval,arguments);    
  qscore_pval=angsd::getArg("-qscore_pval",qscore_pval,arguments);
  hetbias_pval=angsd::getArg("-hetbias_pval",hetbias_pval,arguments);    
  int doHWE=0;
  doHWE=angsd::getArg("-doHWE",doHWE,arguments);    
  if(doHWE==0){
    fprintf(stderr,"must use -doHWE 1, to test for HWE\n");
    exit(0);
  }  
}


abcFilterSNP::abcFilterSNP(const char *outfiles,argStruct *arguments,int inputtype){
  doSnpStat=0;
  outfileZ = NULL;
  edge_pval=mapQ_pval=sb_pval=hwe_pval=qscore_pval=hetbias_pval-1;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSnpStat")||!strcasecmp(arguments->argv[1],"-doPost")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  
  if(doSnpStat==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);  
  int dogeno=0;
  dogeno=angsd::getArg("-dogeno",dogeno,arguments);
  if(doSnpStat){
    //    fprintf(stderr,"running doSnpStat=%d\n",doSnpStat);
    const char *postfix=".snpStat.gz";
    outfileZ = aio::openFileBG(outfiles,postfix);
    kstring_t bufstr;bufstr.s=NULL;bufstr.l=bufstr.m=0;
    ksprintf(&bufstr,"Chromo\tPosition\t+Major +Minor -Major -Minor\tSB1:SB2:SB3\tHWE_LRT:HWE_pval\tbaseQ_Z:baseQ_pval\tmapQ_Z:mapQ_pval\tedge_z:edge_pval");
    if(dogeno){
      ksprintf(&bufstr,"\t+MajorHet +MinorHet -MajorHet -MinorHet nHet\thetStat:hetStat_pval");

    }
    ksprintf(&bufstr,"\n");
    aio::bgzf_write(outfileZ,bufstr.s,bufstr.l);bufstr.l=0;
    free(bufstr.s);
  }

}

abcFilterSNP::~abcFilterSNP(){
  if(outfileZ!=NULL)
    bgzf_close(outfileZ);

}
