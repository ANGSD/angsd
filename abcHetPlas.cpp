#include <cmath>
#include <ctype.h>
#include <assert.h>
#include "abc.h"
#include "analysisFunction.h"
#include "shared.h"
#include "abcHetPlas.h"
#define MAX_QS 256

//struct which contains the results for a single chunk
typedef struct{
  double **freq;//contains the estimate of the freq of the 4 alleles and the unknown minor freq;
  double **llh;//contains the LRT and the llhNull llhAlt;
  int **nItr;
  double **diff;
}resStruct;


double **gen_probs(){
  double **ret =new double*[MAX_QS];
  for(int i=0;i<MAX_QS;i++)
    ret[i] = new double[16];

  for(int i=0;i<MAX_QS;i++){
    double p1 = log(1-pow(10,-i/10.0));
    double p3 = log((1-exp(p1))/3.0);
    double val= std::max(p1,p3);
    val=std::min(p1-val,p3-val);
    for(int ii=0;ii<16;ii++){
      ret[i][ii] = exp(val);
    }
    ret[i][0]=ret[i][5]=ret[i][10]=ret[i][15] =exp(0);
  }
  
  return ret;
}


void abcHetPlas::printArg(FILE *fp){
  fprintf(fp,"------------------------\n%s:\n",__FILE__); 
  fprintf(fp,"\t-doHetPlas=%d (Perform hetplasmid analysis)\n",doHetPlas);
  fprintf(fp,"\t-maxIter=%d\t(Max number of iterations)\n",maxIter);
  fprintf(fp,"\t-minLRT=%f\n",minLRT);
}

int abcHetPlas::makellhs(tNode *tn,double **liks,int *oldsize){
  if(*oldsize<tn->l){
    if(*oldsize>0){
      for(int i=0;i<4;i++)
	delete [] liks[i];
    }
    for(int i=0;i<4;i++)
      liks[i] = new double[tn->l];
    *oldsize = tn->l;
  }
  int depth =0;
  for(int l=0;l<tn->l;l++){
    if(refToInt[tn->seq[l]]!=4){
      for(int ll=0;ll<4;ll++){
     	liks[ll][depth] = probs[tn->qs[l]][refToInt[tn->seq[l]]*4+ll];
      }
      depth++;
    }
  }
  return depth;
}

//program will calculate the likelihood of the x' given the liks, seqdepth is the dimension
double like(double *x,double **liks,int seqdepth){
  double totLik=0;
  for(int s=0;s<seqdepth;s++){
    double tmp=0;
    for(int i=0;i<4;i++)
      tmp+=x[i]*liks[i][s];
    totLik += log(tmp);
  }
  return totLik;
}

// function will perform em estimation of the frequencines
//function will modfy x,itr,and diff

double em3(double *x,double **liks,int seqdepth,int maxIter,double tole,int &itr,double &diff){
  itr=0;
  double likev=-1e12;
  
  while(itr<maxIter){
    double newx[4]={0,0,0,0};
    memcpy(newx,x,sizeof(double)*4);
    for(int s=0;s<seqdepth;s++){
      double tsum=0;
      for(int i=0;i<4;i++)
	tsum += x[i]*liks[i][s];
      for(int i=0;i<4;i++)
	newx[i] += x[i]*liks[i][s]/tsum;
    }
    double tsum =0;
    for(int i=0;i<4;i++)
      tsum +=newx[i];
    for(int i=0;i<4;i++)
      newx[i] /=tsum;
	
    double newlike=like(newx,liks,seqdepth);
    diff=fabs(newlike-likev);
    if(diff<tole){
      likev=newlike;
      memcpy(x,newx,sizeof(double)*4);
      //      fprintf(stderr,"breaking\n");
      break;
    }

    if(1&&(newlike<likev)){
      fprintf(stderr,"Problems: newlike=%f oldlike=%f diff=%e itr=%d\n",newlike,likev,diff,itr);
      break;
      exit(0);
    }
    assert(newlike>=likev);
    //fprintf(stderr,"lik=%f newlike=%f\n",likev,newlike);

    likev=newlike;
    memcpy(x,newx,sizeof(double)*4);
    itr++;
  }
  //fprintf(stderr,"itr=%d like=%f par=(%f,%f,%f,%f)\n",itr,likev,x[0],x[1],x[2],x[3]);
  return 0;
}




void abcHetPlas::doNew(funkyPars *pars){
  //fprintf(stderr,"pars->numSites=%d\n",pars->numSites);
  //  fflush(stderr);
  resStruct *rs = new resStruct;
  rs->freq=  new double*[pars->numSites]; 
  rs->llh= new double*[pars->numSites]; 
  rs->nItr = new int*[pars->numSites];
  rs->diff = new double *[pars->numSites];


  double **liks = new double*[4];
  for(int i=0;i<4;i++) liks[i] =NULL;
  int oldsize =0;
  for(int s=0;s<pars->numSites;s++){
    //    fprintf(stderr,"AAAA %d %d\n",pars->posi[s]+1,pars->keepSites[s]);
    if(pars->keepSites[s]) {
      rs->freq[s] = new double[5*pars->nInd];  
      rs->llh[s] = new double [3*pars->nInd];  
      rs->diff[s] = new double [pars->nInd];  
      rs->nItr[s] = new int [pars->nInd];  

      for(int i=0;i<pars->nInd;i++){
	if(pars->chk->nd[s][i]==NULL)
	  continue;
	int seqdepth=makellhs(pars->chk->nd[s][i],liks,&oldsize);
	if(seqdepth==0){
	  pars->keepSites[s]=0;
	  continue;//DRAGON should only skip ind, is ok with 1 sample
	}
	double par[4]={0.25,0.25,0.25,0.25};
	//¯	fprintf(stderr,"lik=%f\n",like(par,liks,seqdepth));
	em3(par,liks,seqdepth,maxIter,1e-6,rs->nItr[s][i],rs->diff[s][i]);
	int maxBase=angsd::whichMax(par,4);
	//assert(maxBase!=-1);//<- if vals are equal...
	if(maxBase==-1){
	  fprintf(stderr,"CHECK SITE: %d  par=(%f,%f,%f,%f) seqdepth=%d\n",pars->posi[s]+1,par[0],par[1],par[2],par[3],seqdepth);
	  pars->keepSites[s] =0;
	  continue;
	}
	double unknownMinorSum =0;
	for(int um=0;um<4;um++)
	  if(um!=maxBase)
	    unknownMinorSum += par[um];
	memcpy(rs->freq[s]+4*i,par,sizeof(double)*4);
	rs->freq[s][4*i+4]=unknownMinorSum;
	double llhAlt = like(par,liks,seqdepth);
	par[0]=par[1]=par[2]=par[3]=0.0;
	par[maxBase]=1.0;
	double llhNeu=like(par,liks,seqdepth);
	rs->llh[s][i*3]=2*(llhAlt-llhNeu);
	rs->llh[s][i*3+1]=llhNeu;
	rs->llh[s][i*3+2]=llhAlt;

      }
      
    }
  }
  if(liks!=NULL){
    for(int i=0;i<4;i++)
      delete [] liks[i];
    delete [] liks;
  }
  pars->extras[index] = rs;
}

void abcHetPlas::run(funkyPars *pars){
  //  fprintf(stderr,"pars->numSites=%d\n",pars->numSites);
  //exit(0);
  if(doHetPlas==0)
    return ;
  if(!doHetPlas||pars->numSites==0){
    fprintf(stderr,"skipping due to nos sites\n");
    return;
  }
  assert(pars->numSites>0);
  assert(pars->chk!=NULL);
  if(doHetPlas){
    //exit(0);
    doNew(pars);
  }
  
   
}

void abcHetPlas::clean(funkyPars *pars){
  if(doHetPlas==0)
    return;
  if(doHetPlas){
    resStruct *rs = (resStruct *)pars->extras[index];
    for(int i=0;i<pars->numSites;i++)
      if(pars->keepSites[i]){
	delete [] rs->freq[i];
	delete [] rs->llh[i];
	delete [] rs->diff[i];
	delete [] rs->nItr[i];
      }
    delete [] rs->freq;
    delete [] rs->llh;
    delete [] rs->nItr;
    delete [] rs->diff;
    delete rs;
  }
    

}


void abcHetPlas::print(funkyPars *pars){
  if(doHetPlas==0)
    return;
  if(doHetPlas){
    resStruct *rs =(resStruct *) pars->extras[index];
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;

      double maxLRT=0;
      for(int i=0;i<pars->nInd;i++)
	if(rs->llh[s][i*3+0]>maxLRT)
	  maxLRT=rs->llh[s][i*3+0];

      if(maxLRT<minLRT)
	 continue;

	fprintf(outputFile,"%s\t%d",header->target_name[pars->refId],pars->posi[s]+1);
	if(pars->ref!=NULL)
	  fprintf(outputFile,"\t%c",intToRef[pars->ref[s]]);
	else
	  fprintf(outputFile,"\tN");
	for(int i=0;i<pars->nInd;i++) {
	  double *dbl = rs->freq[s]+i*4;
	  int max1 = angsd::whichMax(dbl,4);
	  //	  fprintf(stderr,"max1:%d\n",max1);
	  int max2 = max1==0?1:0;
	  //	  fprintf(stderr,"max2:%d\n",max2);
	  for(int ii=0;ii<4;ii++)
	    //fprintf(stderr,"dbl[%d]:%f dbl[%d]:%f\n",ii,dbl[ii],max2,dbl[max2]);
	    if(ii!=max1 && dbl[ii]>dbl[max2])
	      max2=ii;
	    
	    //fprintf(stderr,"max2:%d\n",max2);
	  
	  fprintf(outputFile,"\t%c",intToRef[max2]);
	  //write freq
	  for(int j=0;j<5;j++)
	    fprintf(outputFile,"\t%f",rs->freq[s][i*4+j]);
	  for(int j=0;j<3;j++)
	    fprintf(outputFile,"\t%f",rs->llh[s][i*3+j]);
	  fprintf(outputFile,"\t%d\t%e",rs->nItr[s][i],rs->diff[s][i]);
	}
	fprintf(outputFile,"\n");
  }
  }
}

void abcHetPlas::getOptions(argStruct *arguments){
  //from command line
  

  doHetPlas=angsd::getArg("-doHetPlas",doHetPlas,arguments);
  
  if(doHetPlas==0)
    return;

  maxIter = angsd::getArg("-maxIter",maxIter,arguments);
  minLRT = angsd::getArg("-minLRT",minLRT,arguments);


}


abcHetPlas::abcHetPlas(const char *outfiles,argStruct *arguments,int inputtype){
  doHetPlas =0;
  maxIter =100;
  minLRT =-1;
  probs=NULL;
  outputFile=NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doHetPlas")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);

  if(doHetPlas>0&&arguments->nInd!=1){
    fprintf(stderr,"Heteroplasmy analysis only implemented for single bamfiles\n");
    exit(0);
  }

  printArg(arguments->argumentFile);

  if(doHetPlas){
    probs=gen_probs();
    outputFile=aio::openFile(outfiles,".hetGL");
    fprintf(outputFile,"Chr\tPos\tRef\tAlt\t");
    fprintf(outputFile,"fA\tfC\tfG\tfT\tfU\t");
    fprintf(outputFile,"LRT\tllhNeu\tllhAlt\t");
    fprintf(outputFile,"nIter\tdiff\n");
  }else
    shouldRun[index] =0;
}


void printStuff(FILE *fp,size_t mat[MAX_QS][MAX_QS][2]){
  for(int st=0;st<2;st++){
    for(int i=0;i<MAX_QS;i++){
      for(int j=0;j<MAX_QS;j++)
	fprintf(fp,"%zu\t",mat[i][j][st]);
      fprintf(fp,"%d\n",st);
    }
  }

}


abcHetPlas::~abcHetPlas(){
  if(outputFile) fclose(outputFile);
  
}


