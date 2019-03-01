 /*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd
  Class that works with counts
*/

#include <cassert>
#include "analysisFunction.h"
#include "abc.h"
#include "abcCounts.h"
#include <htslib/kstring.h>

void abcCounts::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doCounts\t%d\t(Count the number A,C,G,T. All sites, All samples)\n",doCounts);
  fprintf(argFile,"\t-minQfile\t%s\t file with individual quality score thresholds)\n",minQfile);
  fprintf(argFile,"\t-setMaxDepth\t%d\t(If total depth is larger then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMaxDepth);
  fprintf(argFile,"\t-setMinDepth\t%d\t(If total depth is smaller then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMinDepth);
  fprintf(argFile,"\t-setMaxDepthInd\t%d\t(If depth persample is larger then individual is removed from analysis (from site).\n\t\t\t\t -1 indicates no filtering)\n",setMaxDepthInd);
  fprintf(argFile,"\t-setMinDepthInd\t%d\t(If depth persample is smaller then individual is removed from analysis (from site).\n\t\t\t\t -1 indicates no filtering)\n",setMinDepthInd);

  fprintf(argFile,"\t-minInd\t\t%d\t(Discard site if effective sample size below value.\n\t\t\t\t 0 indicates no filtering)\n",minInd);
  fprintf(argFile,"\t-setMaxDiffObs\t%d\t(Discard sites where we observe to many different alleles.\n\t\t\t\t 0 indicates no filtering)\n",setMaxObs);
  fprintf(argFile,"Filedumping:\n");
  fprintf(argFile,"\t-doDepth\t%d\t(dump distribution of seqdepth)\t%s,%s\n",doDepth,postfix4,postfix5);
  fprintf(argFile,"\t  -maxDepth\t%d\t(bin together high depths)\n",maxDepth);
  
  fprintf(argFile,"\t-doQsDist\t%d\t(dump distribution of qscores)\t%s\n",doQsDist,postfix3);
  fprintf(argFile,"\t-minQ\t%d\t(minimumQ)\n",minQ);
  fprintf(argFile,"\t-dumpCounts\t%d\n",dumpCounts);
  fprintf(argFile,"\t  1: total seqdepth for site\t%s\n",postfix1);
  fprintf(argFile,"\t  2: seqdepth persample\t\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  3: A,C,G,T sum over samples\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  4: A,C,G,T sum every sample\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t-iCounts\t%d (Internal format for dumping binary single chrs,1=simple,2=advanced)\n",iCounts);
  fprintf(argFile,"\t-qfile\t%s\t(Only for -iCounts 2)\n",qfileFname);
  fprintf(argFile,"\t-ffile\t%s\t(Only for -iCounts 2)\n",ffileFname);
}

int calcSum(suint *counts,int len){
  int tmp=0;
  for(int i=0;i<len;i++)
    tmp += counts[i];
  return tmp;
}

void printCounts(char *chr,int *posi,suint **counts,int nSites,size_t nInd,kstring_t &bpos,kstring_t &bbin,int dumpType,int *keepSites){
  bpos.l=bbin.l=0;

  for(int s=0;s<nSites;s++){
    if(keepSites[s]==0)
      continue;
    ksprintf(&bpos, "%s\t%d\t%d\n",chr,posi[s]+1,calcSum(counts[s],4*nInd));
    
    //if we need per sample info
    if(dumpType>1) {
      if(dumpType==4)//count A,C,G,T
	for(int i=0;i<4*nInd;i++)
	  ksprintf(&bbin,"%u\t",counts[s][i]);
      else if(dumpType==2){//print A+C+G+T
	for(int n=0;n<nInd;n++)
	  ksprintf(&bbin,"%u\t",counts[s][n*4]+counts[s][n*4+1]+counts[s][n*4+2]+counts[s][n*4+3]);
      }else{//overall sum of A,C,G,T
	size_t tsum[4]={0,0,0,0};
	for(int i=0;i<4*nInd;i++)
	  tsum[i%4] +=counts[s][i];
	ksprintf(&bbin,"%zu\t%zu\t%zu\t%zu",tsum[0],tsum[1],tsum[2],tsum[3]);
      }
      kputc('\n',&bbin);	
    }
  }

}


void abcCounts::getOptions(argStruct *arguments){


  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  if(doCounts==0)
    return;
  //from command line
  minQfile=angsd::getArg("-minQfile",minQfile,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  dumpCounts=angsd::getArg("-dumpCounts",dumpCounts,arguments);
  doQsDist=angsd::getArg("-doQsDist",doQsDist,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);

  setMaxObs = angsd::getArg("-setMaxDiffObs",setMaxObs,arguments);
  setMaxDepth = angsd::getArg("-setMaxDepth",setMaxDepth,arguments);
  setMinDepth = angsd::getArg("-setMinDepth",setMinDepth,arguments);
  setMaxDepthInd = angsd::getArg("-setMaxDepthInd",setMaxDepthInd,arguments);
  setMinDepthInd = angsd::getArg("-setMinDepthInd",setMinDepthInd,arguments);

  doDepth=angsd::getArg("-doDepth",doDepth,arguments);
  maxDepth=angsd::getArg("-maxDepth",maxDepth,arguments);
  iCounts=angsd::getArg("-iCounts",iCounts,arguments);
  

  if(doCounts && arguments->inputtype!=INPUT_BAM && arguments->inputtype!=INPUT_PILEUP){
    fprintf(stderr,"\t-> Cannot calculate -doCounts 1 from GL/GP data\n");
    exit(0);
  }
  if(dumpCounts&&doCounts==0){
    fprintf(stderr,"\t-> You must supply -doCounts if you want to dumpcounts\n");
    exit(0);
  }

  if(doDepth!=0&&doCounts==0){
    fprintf(stderr,"\t-> Must supply -doCounts 1 if you want depth distribution\n");
    exit(0);
  }

 if(doQsDist!=0&&doCounts==0){
    fprintf(stderr,"\t-> Must supply -doCounts 1 if you want qscore distribution\n");
    exit(0);
  }
 if(iCounts!=0&&doCounts==0){
   fprintf(stderr,"\t-> Must supply -doCounts 1 if you want to write iCounts\n");
   exit(0);
 }
 if(iCounts!=0&&arguments->nInd!=1){
   fprintf(stderr,"\t-> iCounts only for single chromosomes and single samples\n");
   exit(0);
 }
 ffileFname = angsd::getArg("-ffile",ffileFname,arguments);
 qfileFname = angsd::getArg("-qfile",qfileFname,arguments);
 
 if(ffileFname&&iCounts!=2){
   fprintf(stderr,"\t-> -ffile is only used for -iCounts 2\n");
   exit(0);
 }
 if(qfileFname&&iCounts!=2){
   fprintf(stderr,"\t-> -qfile is only used for -iCounts 2\n");
   exit(0);
 }
 int tmp=-1;
 tmp = angsd::getArg("-minQ",tmp,arguments);
 if(iCounts==2 && qfileFname && tmp==-1){
   fprintf(stderr,"\t-> Must fix -minQ 0 when using -qfile \n");
   exit(0);
 }
 if(doQsDist)
   if(minQ>0)
     fprintf(stderr,"\t-> NB: you are doing qscore distribution, set -minQ to zero if you want to count qscores below: %d\n",minQ);
}
//constructor
abcCounts::abcCounts(const char *outfiles,argStruct *arguments,int inputtype){
  minQ = MINQ;
  globCount = NULL;
  const char *delim = "\t\n ";
  ffileFname=qfileFname=NULL;
  for(int i=0;i<255;i++)
    lookup[i]=-1;
  lookup['A'] = 0;
  lookup['C'] = 1;
  lookup['G'] = 2;
  lookup['T'] = 3;
  lookup['N'] = 4;
  lookup['a'] = 5;
  lookup['c'] = 6;
  lookup['g'] = 7;
  lookup['t'] = 8;
  lookup['n'] = 4;

  oFileIcounts = NULL;
  iCounts =0;
  nInd=arguments->nInd;
  minInd = 0;
  setMinDepth = -1;
  setMaxDepth = -1;
  setMinDepthInd = -1;
  setMaxDepthInd = -1;

  dumpCounts =0;
  doCounts = 0;
  doQsDist = 0;
  setMaxObs = 0;
  doDepth = 0;
  maxDepth = 100;

  minQfile=NULL;

  //make output files
  postfix1=".pos.gz";
  postfix2=".counts.gz";
  postfix3=".qs";
  postfix4=".depthSample";
  postfix5=".depthGlobal";
  bpos.s=NULL;bpos.l=bpos.m=0;
  bbin.s=NULL;bbin.l=bbin.m=0;
  bufstr.s=NULL;bufstr.l=bufstr.m=0;

  //from command line
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doCounts")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  if(doCounts==0)
    shouldRun[index] = 0;




  oFiles=strdup(outfiles);
  printArg(arguments->argumentFile);
  
  if(qfileFname){
    FILE *infileP=NULL;
    if(!((infileP=fopen(qfileFname,"r")))){
      fprintf(stderr,"\t-> Problem openingfile: %s\n",qfileFname);
      exit(0);
    }
    int line=0;
    char buf[LENS];
    while(fgets(buf,LENS,infileP)){
      if(line==4){ // new: to make sure qCutoff order corresponds with lookup
	qCutoff[line]=0;
	line++;
      }
      qCutoff[line]=atoi(strtok(buf,delim));
      fprintf(stdout,"%d\n",qCutoff[line]);
      fflush(stdout);
      line++;
    }
    fclose(infileP);
  }else{
    for(int i=0;i<9;i++)
      qCutoff[i] = MINQ;//<-from abc.h

  }
  if(ffileFname){
    FILE *infileP=NULL;
    if(!((infileP=fopen(ffileFname,"r")))){
      fprintf(stderr,"\t-> Problem openingfile: %s\n",qfileFname);
      exit(0);
    }

    int line=0;
    char buf[LENS];
    while(fgets(buf,LENS,infileP)){
      if(line==4){ // new: to make sure qCutoff order corresponds with lookup
	fCutoff[line]=1;
	line++;
      }
      fCutoff[line]=atof(strtok(buf,delim));
      fprintf(stdout,"%f\n",fCutoff[line]);
      fflush(stdout);
      line++;
    }
    fclose(infileP);

  }else{
    for(int i=0;i<9;i++)
      fCutoff[i] = 1;
  }

  //  oFileCountsPos = oFileCountsBin = oFileQs = NULL;
  oFileCountsPos = oFileCountsBin =  NULL;

  if(dumpCounts){
    oFileCountsPos = aio::openFileBG(outfiles,postfix1);
    bufstr.l=0;ksprintf(&bufstr,"chr\tpos\ttotDepth\n");
    aio::bgzf_write(oFileCountsPos,bufstr.s,bufstr.l);bufstr.l=0;
    if(dumpCounts>1)
      oFileCountsBin = aio::openFileBG(outfiles,postfix2);
    if(dumpCounts==2){
      bufstr.l=0;
      for(int i=0;i<arguments->nInd;i++)
	ksprintf(&bufstr,"ind%dTotDepth\t",i);
      aio::bgzf_write(oFileCountsBin,bufstr.s,bufstr.l);bufstr.l=0;
    }
    if(dumpCounts==3){
      bufstr.l=0;
      ksprintf(&bufstr,"totA\ttotC\ttotG\ttotT");
      aio::bgzf_write(oFileCountsBin,bufstr.s,bufstr.l);bufstr.l=0;
    }
    if(dumpCounts==4){
      bufstr.l=0;
      for(int i=0;i<arguments->nInd;i++)
	ksprintf(&bufstr,"ind%d_A\tind%d_C\tind%d_G\tind%d_T\t",i,i,i,i);
      aio::bgzf_write(oFileCountsBin,bufstr.s,bufstr.l);bufstr.l=0;
    }
    if(dumpCounts>1){
      bufstr.l=0;ksprintf(&bufstr,"\n");
      aio::bgzf_write(oFileCountsBin,bufstr.s,bufstr.l);bufstr.l=0;
    }
  }
  if(iCounts){
    oFileIcounts = aio::openFileBG(outfiles,".icnts.gz");

  }
  if(doQsDist){
    //datastructures needed
    qsDist = new size_t[256];
    memset(qsDist,0,256*sizeof(size_t));
    //prepare outputfile
    
  }
  if(doDepth){
    depthCount=new size_t *[arguments->nInd];
    for(int i=0;i<nInd;i++)
      depthCount[i]=new size_t[maxDepth+1];
    for(int i=0;i<nInd;i++)
      for(int j=0;j<maxDepth+1;j++)
	depthCount[i][j]=0;
    
    globCount = new size_t[maxDepth+1];
    memset(globCount,0,sizeof(size_t)*(maxDepth+1));

    
  }
  if(minQfile!=NULL){
    minQmat = angsd::getMatrix(minQfile,0,100000);
    if(minQmat.x!=nInd){
      fprintf(stderr,"Number of lines in the minQfile does not match the number of individuals \n");
      exit(0);
    }
    if(!(minQmat.y==1||minQmat.y==4)){
      fprintf(stderr,"Number of colums in the minQfile has to be 1 or 4 \n");
      exit(0);
    }
  }
  
  
}


void printQs(FILE *fp,size_t *ary){
  int firstidx=0;
  for(int i=0;i<256;i++)
     if(ary[i]!=0){
      firstidx=i;
      break;
    }
  int lastidx=255;
  for(int i=255;i>=0;i--)
    if(ary[i]!=0){
      lastidx=i;
      break;
    }
  for(int i=firstidx;i<=lastidx;i++)
    fprintf(fp,"%d\t%lu\n",i,ary[i]);
  
}


abcCounts::~abcCounts(){
  if(oFileCountsBin!=NULL)
    bgzf_close(oFileCountsBin);
  if(oFileCountsPos!=NULL)
    bgzf_close(oFileCountsPos);
  if(doQsDist){
    FILE *oFileQs = NULL;
    oFileQs = aio::openFile(oFiles,postfix3);
    fprintf(oFileQs,"qscore\tcounts\n");
    printQs(oFileQs,qsDist);
    if(oFileQs) fclose(oFileQs);
    delete[] qsDist;

  }

  if(doDepth){
    FILE *oFileSamplDepth = aio::openFile(oFiles,postfix4);
    FILE *oFileGlobDepth = aio::openFile(oFiles,postfix5);
    for(int i=0;i<nInd;i++){
      for(int j=0;j<maxDepth+1;j++){
	fprintf(oFileSamplDepth,"%lu\t",depthCount[i][j]);
      }
      fprintf(oFileSamplDepth,"\n");
    }
    //thorfinn
    for(int j=0;j<maxDepth+1;j++)
      fprintf(oFileGlobDepth,"%lu\t",globCount[j]);
    fprintf(oFileGlobDepth,"\n");
  

    //clean depthCount
    for(int i=0;i<nInd;i++)
      delete[]  depthCount[i];
    delete[] depthCount; 
    
    if(oFileSamplDepth) fclose(oFileSamplDepth);
    if(oFileSamplDepth) fclose(oFileGlobDepth);
  }
  
  if(minQfile!=NULL){
    //  angsd::printMatrix(minQmat,stderr);
    angsd::deleteMatrix(minQmat);
  }
  if(oFileIcounts!=NULL)
    bgzf_close(oFileIcounts);
  
  free(oFiles);
  free(bpos.s);
  free(bbin.s);
  free(bufstr.s);
  if(globCount)
    delete [] globCount;
}

void countQs(const chunkyT *chk,size_t *ret,int *keepSites,int minQ){
  // fprintf(stderr,"chk->nSites:%d\n",chk->nSites);
  for(int s=0;s<chk->nSites;s++){
    if(keepSites[s]==0)
      continue;
    //loop over sites
    for(int n=0;n<chk->nSamples;n++){
      //loop over samples
      for(int l=0;chk->nd[s][n]&&(l<chk->nd[s][n]->l);l++){
	//loop over persample reads for this position/sample
	if(chk->nd[s][n]->qs[l]>=minQ)
	  ret[chk->nd[s][n]->qs[l]]++;
	
      }
    }
  }
}



void abcCounts::print(funkyPars *pars){
  if(pars->numSites==0)
    return;
  if(dumpCounts)
    printCounts(header->target_name[pars->refId],pars->posi,pars->counts,pars->numSites,pars->nInd,bpos,bbin,dumpCounts,pars->keepSites);
  if(bbin.l>0)
    aio::bgzf_write(oFileCountsBin,bbin.s,bbin.l);bbin.l=0;
  if(bpos.l>0)
    aio::bgzf_write(oFileCountsPos,bpos.s,bpos.l);bpos.l=0;

  if(doQsDist)
    countQs(pars->chk,qsDist,pars->keepSites,minQ);
  
  if(doDepth!=0){
    assert(pars->counts!=NULL);
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      for(int i=0;i<pars->nInd;i++){
	int sum=0;
	for(int a=0;a<4;a++)
	  sum+=pars->counts[s][i*4+a];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
	depthCount[i][sum]++;	
      }
    }
    //thorfinn below
    
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      int sum=0;
      for(int i=0;i<4*pars->nInd;i++){
	sum+=pars->counts[s][i];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
      }
      globCount[sum]++;
    } 
  }
  if(iCounts==1){
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      int cnt[4];
      for(int i=0;i<4;i++)
	cnt[i] = pars->counts[s][i];
     
      int p = pars->posi[s]+1;
      aio::bgzf_write(oFileIcounts,&p,sizeof(int)); // new
      aio::bgzf_write(oFileIcounts,cnt,sizeof(int)*4); // new

    }
  }
  if(iCounts==2){

    for(int s=0;s<pars->chk->nSites;s++){
      //      fprintf(stderr,"TT: %d %d\n",pars->chk->nd[s][0].l2,pars->chk->nd[s][0].deletion);
      if(pars->keepSites[s]==0||pars->chk->nd==NULL||pars->chk->nd[s][0]->l2||pars->chk->nd[s][0]->deletion)
	continue;

      int count[4]={0,0,0,0};
      //loop over persample reads
      for(int l=0;pars->chk->nd[s][0]&&l<pars->chk->nd[s][0]->l;l++){
	int allele = refToInt[pars->chk->nd[s][0]->seq[l]];
	int tmp=lookup[pars->chk->nd[s][0]->seq[l]];
#if 0
	fprintf(stderr,"l:%d b:%c qs:%d\n",l,pars->chk->nd[s][0].seq[l],pars->chk->nd[s][0].qs[l]);
#endif
	if(allele==4)//skip of 'n'/'N'
	  continue;

	if(qfileFname&&pars->chk->nd[s][0]->qs[l]<=qCutoff[tmp]) //skip if we are using qscores and we are below
	  continue;
	  
	// if q score bin to sample from: count if sampled 
	if(pars->chk->nd[s][0]->qs[l]==qCutoff[tmp]+1){
	  double x =drand48();
	  //fprintf(stderr,"skipping\n");
	  if(x>=fCutoff[tmp])
	    count[tmp%5]++;
	}else // if in q score bigger: always count 
	  count[tmp%5]++;
	
      }
      if(1||count[0]+ count[1]+count[2]+count[3]){
	int p = pars->posi[s]+1;
	aio::bgzf_write(oFileIcounts,&p,sizeof(int)); // new
	aio::bgzf_write(oFileIcounts,count,sizeof(int)*4); // new
      }
    }
  }
}



void abcCounts::clean(funkyPars *pars){
  if(doCounts||dumpCounts){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->counts[i];
    delete [] pars->counts;
  }

}


//dragon update with keeplist so we only count necessary sites
suint **abcCounts::countNucs(const chunkyT *chk,int *keepSites,int mmin,int mmax){
  //  fprintf(stderr,"nsamples:%d\n",chk->nSamples);
  suint **cnts = new suint*[chk->nSites];
  if(minQfile==NULL) {
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;chk->nd[s][n]&&l<chk->nd[s][n]->l;l++){
	  //	  fprintf(stderr,"s:%d n:%d ret.l:%d l:%d seq:%c\n",s,n,chk->nd[s][n].l,l,chk->nd[s][n].seq[l]);
	  int allele = refToInt[chk->nd[s][n]->seq[l]];
	
	  if(allele==4){
	    continue;
	  }
	  cnts[s][4*n+allele]++;
	}
	if(chk->nd[s][n]&&(mmin!=-1||mmax!=-1)){
	  int ndep=0;
	  for(int jj=0;jj<4;jj++) 
	    ndep += cnts[s][4*n+jj];
	  if(mmin!=-1&&ndep<mmin){
	    for(int jj=0;jj<4;jj++)
	      cnts[s][4*n+jj] = 0;
	    chk->nd[s][n]->l=0;
	  }
	  if(mmax!=-1&&ndep>mmax){
	    for(int jj=0;jj<4;jj++)
	      cnts[s][4*n+jj] = 0;
	    chk->nd[s][n]->l=0;
	  }    
	}
      }
    }
  }
  else if(minQmat.y>1){
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;chk->nd[s][n]&&l<chk->nd[s][n]->l;l++){
	  int allele = refToInt[chk->nd[s][n]->seq[l]];	  
	  if(allele==4)
	    continue;
	  if(chk->nd[s][n]->qs[l] < minQmat.matrix[n][allele]){
	  continue;
	  }
	  cnts[s][4*n+allele]++;
	}
      }
    }

  }
 else{
    for(int s=0;s<chk->nSites;s++){
      cnts[s] = new suint[4*chk->nSamples];
      if(keepSites[s]==0)
	continue;
      memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
      //loop over samples
      for(int n=0;n<chk->nSamples;n++){
	//loop over persample reads
	for(int l=0;chk->nd[s][n]&&l<chk->nd[s][n]->l;l++){
	  int allele = refToInt[chk->nd[s][n]->seq[l]];
	  if(chk->nd[s][n]->qs[l] < minQmat.matrix[n][0]||allele==4){
 
	  continue;
	  }
	  cnts[s][4*n+allele]++;
	}
      }
    }

  }
  return cnts;
}




void abcCounts::run(funkyPars *pars){
  if(doCounts==0)
    return;
  assert(pars->chk!=NULL&&pars->counts==NULL);
  pars->counts = countNucs(pars->chk,pars->keepSites,setMinDepthInd,setMaxDepthInd);
  // fprintf(stderr,"%d\n",pars->keepSites[0]);
  for(int s=0;s<pars->numSites;s++){// Why is this loop needed? is it to remove sites with no data above minQ filters?
    size_t ss=0;
    for(int i=0;i<4*pars->nInd;i++)
      if(pars->keepSites[s]) 
	if(pars->counts[s][i]){
	  ss++;
	  break;
	}
    if(ss==0)
      pars->keepSites[s]=0;
  }
  //modify keepsites;
  if(minInd!=0) {
    for(int i=0;i<pars->numSites;i++){
      if(pars->keepSites[i]==0)
	continue;
      //THIS IS REALLY STUPID but lets count number of samples wiht info
      int nDep =0;
      for(int s=0;s<pars->nInd;s++){
	int dep=0;
	for(int j=0;j<4;j++)
	  dep += pars->counts[i][s*4+j];
	if(dep)
	  nDep++;
      }
      //nDep is now the number of sapmles wiht info
      if(nDep<minInd)
	pars->keepSites[i] = 0;
      else
	pars->keepSites[i] = nDep;
      
    }
  }
  if(setMaxDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      size_t totSum = calcSum(pars->counts[s],4*nInd);
      if(totSum>setMaxDepth)
	pars->keepSites[s]=0;
    }
  }
  //fprintf(stderr,"minDepth=%d nInd=%d pars->numsites=%d\n",minDepth,nInd,pars->numSites);
  if(setMinDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      size_t totSum = calcSum(pars->counts[s],4*pars->nInd);
      if(totSum<setMinDepth)
	pars->keepSites[s]=0;
      else{
	int nDep =0;
	for(int i=0;i<pars->nInd;i++){
	  int dep=0;
	  for(int j=0;j<4;j++)
	    dep += pars->counts[s][i*4+j];
	  if(dep)
	    nDep++;
	}
	pars->keepSites[s]= nDep;

      }
    }
  }
  if(setMaxObs!=0){
    //only for first sample...
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      int ndiff=0;
      for(int i=0;i<4;i++)
	if(pars->counts[s][i])
	  ndiff++;
      if(ndiff>setMaxObs)
	pars->keepSites[s]=0;

    }
  }
}
