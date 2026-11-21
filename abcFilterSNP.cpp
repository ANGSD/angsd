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
#include "aio.h"

/*
  wilcox manwhitney u test, whatever ,implementation is from wiki
  validated with wilcox.test (qscores of major,qscores of minor,exact=T,correct=F)
  returns Zscore, for some reason...
  The values are counted in NBIN bins, base and mapping qualities are unsigned char (0-255)
 */
#define NBIN 256

double mann(int majD[NBIN],int minD[NBIN]){
  double U[NBIN];
  double cumLess1 =0;
  for(int i=0;i<NBIN;i++){
    U[i]=(double) minD[i]*(majD[i]/2.0+cumLess1);
    cumLess1 += majD[i];
  }
  double U1 =0;
  for(int i=0;i<NBIN;i++) 
    U1+= U[i];
  //below is a check
  double cumLess2 =0;
  for(int i=0;i<NBIN;i++){
    U[i]=(double) majD[i]*(minD[i]/2.0+cumLess2);
    cumLess2 += minD[i];
  }
  double U2 =0;
  for(int i=0;i<NBIN;i++) 
    U2+= U[i];
  

  double mu=cumLess1*cumLess2/2.0;
  double sigu=sqrt((cumLess1*cumLess2*(cumLess1+cumLess2+1))/12.0);
  double Z=(std::min(U1,U2)-mu)/sigu;

 
  //  fprintf(stderr,"U1:%f U2:%f U1+U2:%f nObs:%f Z:%f p.value:%f\n",U1,U2,U1+U2,cumLess1*cumLess2,Z,2*phi(Z));

  return Z;

}



double edgebias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[NBIN];
  int minD[NBIN];
  memset(majD,0,sizeof(int)*NBIN);
  memset(minD,0,sizeof(int)*NBIN);

  for(int i=0;i<nInd;i++){
    tNode *nd = tn[i];
    if(nd==NULL)
      continue;
    for(int l=0;l<nd->l;l++){
      int obB = refToInt[nd->seq[l]];

      //distance to the nearest read end, long reads are put in the last bin
      int dist = std::min((int)nd->posi[l],(int)nd->isop[l]);
      if(dist>NBIN-1)
	dist = NBIN-1;
      if(obB==maj)
	majD[dist]++;
      else if(obB==min)
	minD[dist]++;
    }
  }

  return mann(majD,minD);
}



double mapQbias(tNode **tn,int nInd,int maj,int min){
  //  fprintf(stderr,"phi:%f\n",phi(3));
  int majD[NBIN];
  int minD[NBIN];
  memset(majD,0,sizeof(int)*NBIN);
  memset(minD,0,sizeof(int)*NBIN);

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
  int majD[NBIN];
  int minD[NBIN];
  memset(majD,0,sizeof(int)*NBIN);
  memset(minD,0,sizeof(int)*NBIN);

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
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doSnpStat\t%d\t(1: write per site statistics to .snpStat.gz)\n",doSnpStat);
  fprintf(argFile,"\t-edge_pval\t%g\t(remove sites with an edge bias p-value below)\n",edge_pval);
  fprintf(argFile,"\t-mapQ_pval\t%g\t(remove sites with a mapping quality bias p-value below)\n",mapQ_pval);
  fprintf(argFile,"\t-sb_pval\t%g\t(remove sites with a strand bias p-value below)\n",sb_pval);
  fprintf(argFile,"\t-hwe_pval\t%g\t(remove sites with a HWE p-value below, needs -doHWE 1)\n",hwe_pval);
  fprintf(argFile,"\t-qscore_pval\t%g\t(remove sites with a base quality bias p-value below)\n",qscore_pval);
  fprintf(argFile,"\t-hetbias_pval\t%g\t(remove sites with a heterozygote bias p-value below, needs -doGeno and -doPost)\n",hetbias_pval);
  fprintf(argFile,"\n\tA p-value of -1 means the filter is not used.\n");
  fprintf(argFile,"\tNeeds the reads (bam/cram or -pileup input) and -doMajorMinor. The statistics are\n");
  fprintf(argFile,"\tcalculated from the bases that are the major or the minor allele.\n");
  fprintf(argFile,"\n\tThe columns in .snpStat.gz are:\n");
  fprintf(argFile,"\t  chromosome, position\n");
  fprintf(argFile,"\t  number of bases: +Major +Minor -Major -Minor (+/- is the strand)\n");
  fprintf(argFile,"\t  strand bias SB1:SB2:SB3, SB3 is the Fisher exact p-value used by -sb_pval\n");
  fprintf(argFile,"\t  HWE_LRT:HWE_pval (nan:nan without -doHWE 1)\n");
  fprintf(argFile,"\t  Z:p-value of the Mann-Whitney tests between major and minor bases for\n");
  fprintf(argFile,"\t  base quality, mapping quality and distance to the nearest read end (edge).\n");
  fprintf(argFile,"\t  Distances above 255 are counted as 255.\n");
  fprintf(argFile,"\t  With -doGeno: the base counts as above but only in the heterozygotes, their sum, and\n");
  fprintf(argFile,"\t  hetStat:p-value of a chi-square test that the two alleles are 50/50 in the heterozygotes\n");
  fprintf(argFile,"\tOnly the sites that pass all filters are written.\n");
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

      //without -doHWE the HWE column is nan:nan
      double lrt = NAN;
      double pval = NAN;
      if(doHWE){
	funkyHWE *hweStruct = (funkyHWE *) pars->extras[9];//THIS IS VERY NASTY! the ordering within general.cpp is now important
	lrt = 2*hweStruct->like0[s]-2*hweStruct->likeF[s];
	if(std::isnan(lrt))
	  pval=lrt;
	else if(lrt<0)
	  pval =1;
	else
	  pval =1- chi.cdf(lrt);
      }
      ksprintf(&persite,"%f:%e\t",lrt,pval);
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
      if(edge_pval!=-1 && 2*phi(Z)<edge_pval)
	pars->keepSites[s] = 0;

      genoCalls *gcw =(genoCalls *) pars->extras[11];
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
  if(doSnpStat!=1){
    fprintf(stderr,"\t-> -doSnpStat %d is not an option, use -doSnpStat 1\n",doSnpStat);
    exit(1);
  }
  //the statistics are calculated from the reads
  if(arguments->inputtype!=INPUT_BAM&&arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"\t-> -doSnpStat needs the reads, so it only works with bam/cram (-bam/-i) or -pileup input\n");
    exit(1);
  }
  int domajorminor=0;
  domajorminor = angsd::getArg("-domajorminor",domajorminor,arguments);
  if(domajorminor==0){
    fprintf(stderr,"\t-> Must supply -doMajorMinor for running dosnpstat (needs to look a distributions of major and minor alleles)\n");
    exit(1);
  }
  //from command line

  edge_pval=angsd::getArg("-edge_pval",edge_pval,arguments);
  mapQ_pval=angsd::getArg("-mapQ_pval",mapQ_pval,arguments);
  sb_pval=angsd::getArg("-sb_pval",sb_pval,arguments);
  hwe_pval=angsd::getArg("-hwe_pval",hwe_pval,arguments);
  qscore_pval=angsd::getArg("-qscore_pval",qscore_pval,arguments);
  hetbias_pval=angsd::getArg("-hetbias_pval",hetbias_pval,arguments);
  doHWE=angsd::getArg("-doHWE",doHWE,arguments);
  if(hwe_pval!=-1&&doHWE==0){
    fprintf(stderr,"\t-> -hwe_pval needs -doHWE 1\n");
    exit(1);
  }
  doGeno=angsd::getArg("-doGeno",doGeno,arguments);
  if(hetbias_pval!=-1&&doGeno==0){
    fprintf(stderr,"\t-> -hetbias_pval needs the called genotypes from -doGeno\n");
    exit(1);
  }
}


abcFilterSNP::abcFilterSNP(const char *outfiles,argStruct *arguments,int inputtype){
  doSnpStat=0;
  outfileZ = NULL;
  edge_pval=mapQ_pval=sb_pval=hwe_pval=qscore_pval=hetbias_pval=-1;
  doHWE=doGeno=0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doSnpStat")){
      printArg(stderr);
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
  if(doSnpStat){
    //    fprintf(stderr,"running doSnpStat=%d\n",doSnpStat);
    const char *postfix=".snpStat.gz";
    outfileZ = aio::openFileBG(outfiles,postfix);
    kstring_t bufstr;bufstr.s=NULL;bufstr.l=bufstr.m=0;
    ksprintf(&bufstr,"Chromo\tPosition\t+Major +Minor -Major -Minor\tSB1:SB2:SB3\tHWE_LRT:HWE_pval\tbaseQ_Z:baseQ_pval\tmapQ_Z:mapQ_pval\tedge_z:edge_pval");
    if(doGeno){
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
