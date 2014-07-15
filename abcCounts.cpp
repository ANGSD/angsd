 /*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  part of angsd
  Class that works with counts
*/

#include <zlib.h>
#include <cassert>
#include "analysisFunction.h"
#include "abc.h"
#include "abcCounts.h"
#include "kstring.h"

void abcCounts::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doCounts\t%d\t(Count the number A,C,G,T. All sites, All samples)\n",doCounts);
  
  
  fprintf(argFile,"\t-minQfile\t%s\t file with individual quality score thresholds)\n",minQfile);
  fprintf(argFile,"\t-setMaxDepth\t%d\t(If total depth is larger then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMaxDepth);
  fprintf(argFile,"\t-setMinDepth\t%d\t(If total depth is smaller then site is removed from analysis.\n\t\t\t\t -1 indicates no filtering)\n",setMinDepth);
  fprintf(argFile,"\t-minInd\t\t%d\t(Discard site if effective sample size below value.\n\t\t\t\t 0 indicates no filtering)\n",minInd);
  fprintf(argFile,"Filedumping:\n");
  fprintf(argFile,"\t-doDepth\t%d\t(dump distribution of seqdepth)\t%s,%s\n",doDepth,postfix4,postfix5);
  fprintf(argFile,"\t  -maxDepth\t%d\t(bin together high depths)\n",maxDepth);
  
  fprintf(argFile,"\t-doQsDist\t%d\t(dump distribution of qscores)\t%s\n",doQsDist,postfix3);  
  fprintf(argFile,"\t-dumpCounts\t%d\n",dumpCounts);
  fprintf(argFile,"\t  1: total seqdepth for site\t%s\n",postfix1);
  fprintf(argFile,"\t  2: seqdepth persample\t\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  3: A,C,G,T sum over samples\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  4: A,C,G,T sum every sample\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t-iCounts\t%d (Internal format for dumping binary single chrs,1=simple,2=advanced)\n",iCounts);
  fprintf(argFile,"\t-qfile\t%s\t(Onlyfor -iCounts 2)\n",qfileFname);
  fprintf(argFile,"\t-ffile\t%s\t(Onlyfor -iCounts 2)\n",ffileFname);
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

  //from command line
  minQfile=angsd::getArg("-minQfile",minQfile,arguments);
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  dumpCounts=angsd::getArg("-dumpCounts",dumpCounts,arguments);
  doQsDist=angsd::getArg("-doQsDist",doQsDist,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);
  setMaxDepth = angsd::getArg("-setMaxDepth",setMaxDepth,arguments);
  setMinDepth=angsd::getArg("-setMinDepth",setMinDepth,arguments);
  doDepth=angsd::getArg("-doDepth",doDepth,arguments);
  maxDepth=angsd::getArg("-maxDepth",maxDepth,arguments);
  iCounts=angsd::getArg("-iCounts",iCounts,arguments);
  

  if(dumpCounts&&doCounts==0){
    fprintf(stderr,"You must supply -doCounts if you want to dumpcounts\n");
    exit(0);
  }

  if(doDepth!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want depth distribution");
    exit(0);
  }

 if(doQsDist!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want qscore distribution");
    exit(0);
  }
 if(iCounts!=0&&doCounts==0){
   fprintf(stderr,"Must supply -doCounts 1 if you want to write iCounts\n");
   exit(0);
 }
 if(iCounts!=0&&arguments->nInd!=1){
   fprintf(stderr,"iCounts only for single chromosomes and single samples\n");
   exit(0);
 }
 ffileFname = angsd::getArg("-ffile",ffileFname,arguments);
 qfileFname = angsd::getArg("-qfile",qfileFname,arguments);
 
 if(ffileFname&&iCounts!=2){
   fprintf(stderr,"-ffile is only used for -iCounts 2\n");
   exit(0);
 }
 if(qfileFname&&iCounts!=2){
   fprintf(stderr,"-qfile is only used for -iCounts 2\n");
   exit(0);
 }
 int tmp=-1;
 tmp = angsd::getArg("-minQ",tmp,arguments);
 if(iCounts==2 && qfileFname && tmp==-1){
   fprintf(stderr,"Must fix -minQ 0 when using -qfile \n");
   exit(0);
 }

}
//constructor
abcCounts::abcCounts(const char *outfiles,argStruct *arguments,int inputtype){
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

  oFileIcounts = Z_NULL;
  iCounts =0;
  nInd=arguments->nInd;
  minInd = 0;
  setMinDepth =-1;
  dumpCounts =0;
  doCounts = 0;
  doQsDist = 0;
  
  doDepth = 0;
  maxDepth = 100;
  setMaxDepth = -1;
  minQfile=NULL;

  //make output files
  postfix1=".pos.gz";
  postfix2=".counts.gz";
  postfix3=".qs";
  postfix4=".depthSample";
  postfix5=".depthGlobal";
  bpos.s=NULL;bpos.l=bpos.m=0;
  bbin.s=NULL;bbin.l=bbin.m=0;

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
  oFileCountsPos = oFileCountsBin =  Z_NULL;

  if(dumpCounts){
    oFileCountsPos = aio::openFileGz(outfiles,postfix1,GZOPT);
    gzprintf(oFileCountsPos,"chr\tpos\ttotDepth\n");
    if(dumpCounts>1)
      oFileCountsBin = aio::openFileGz(outfiles,postfix2,GZOPT);
    if(dumpCounts==2)
      for(int i=0;i<arguments->nInd;i++)
	gzprintf(oFileCountsBin,"ind%dTotDepth\t",i);
    if(dumpCounts==3)
      gzprintf(oFileCountsBin,"totA\ttotC\ttotG\ttotT");
    if(dumpCounts==4)
      for(int i=0;i<arguments->nInd;i++)
	gzprintf(oFileCountsBin,"ind%d_A\tind%d_C\tind%d_G\tind%d_T\t",i,i,i,i);
    if(dumpCounts>1)
      gzprintf(oFileCountsBin,"\n");
  }
  if(iCounts){
    oFileIcounts = aio::openFileGz(outfiles,".icnts.gz",GZOPT);

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
  if(oFileCountsBin!=Z_NULL)    gzclose(oFileCountsBin);
  if(oFileCountsPos!=Z_NULL)    gzclose(oFileCountsPos);
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
  if(oFileIcounts!=Z_NULL) gzclose(oFileIcounts);
  
  free(oFiles);
  free(bpos.s);
  free(bbin.s);
}

void countQs(const chunkyT *chk,size_t *ret,int *keepSites){
  // fprintf(stderr,"chk->nSites:%d\n",chk->nSites);
  for(int s=0;s<chk->nSites;s++){
    if(keepSites[s]==0)
      continue;
    //loop over sites
    for(int n=0;n<chk->nSamples;n++){
      //loop over samples
      for(int l=0;l<chk->nd[s][n].l;l++){
	//loop over persample reads for this position/sample
	ret[chk->nd[s][n].qs[l]]++;
	
      }
    }
  }
}



void abcCounts::print(funkyPars *pars){
  if(pars->numSites==0)
    return;
  if(dumpCounts)
    printCounts(header->name[pars->refId],pars->posi,pars->counts,pars->numSites,pars->nInd,bpos,bbin,dumpCounts,pars->keepSites);
  gzwrite(oFileCountsBin,bbin.s,bbin.l);
  gzwrite(oFileCountsPos,bpos.s,bpos.l);

  if(doQsDist)
    countQs(pars->chk,qsDist,pars->keepSites);
  
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
      gzwrite(oFileIcounts,&p,sizeof(int)); // new
      gzwrite(oFileIcounts,cnt,sizeof(int)*4); // new

    }
  }
  if(iCounts==2){

    for(int s=0;s<pars->chk->nSites;s++){
      //      fprintf(stderr,"TT: %d %d\n",pars->chk->nd[s][0].l2,pars->chk->nd[s][0].deletion);
      if(pars->keepSites[s]==0||pars->chk->nd[s][0].l2||pars->chk->nd[s][0].deletion)
	continue;
      ////BELOW IS SOLELY TO KEEP format similar
      if(0&&s<pars->chk->nSites-1&&pars->chk->nd[s+1][0].deletion)
	continue;
      ////THE ABOVE IS SOLELY TO KEEP FORMAT SIMILAR

      int count[4]={0,0,0,0};
      //loop over persample reads
      for(int l=0;l<pars->chk->nd[s][0].l;l++){
	int allele = refToInt[pars->chk->nd[s][0].seq[l]];
	int tmp=lookup[pars->chk->nd[s][0].seq[l]];
#if 0
	fprintf(stderr,"l:%d b:%c qs:%d\n",l,pars->chk->nd[s][0].seq[l],pars->chk->nd[s][0].qs[l]);
#endif
	if(allele==4)//skip of 'n'/'N'
	  continue;

	if(qfileFname&&pars->chk->nd[s][0].qs[l]<=qCutoff[tmp]) //skip if we are using qscores and we are below
	  continue;
	  
	// if q score bin to sample from: count if sampled 
	if(pars->chk->nd[s][0].qs[l]==qCutoff[tmp]+1){
	  double x =drand48();
	  //fprintf(stderr,"skipping\n");
	  if(x>=fCutoff[tmp])
	    count[tmp%5]++;
	}else // if in q score bigger: always count 
	  count[tmp%5]++;
	
      }
      if(1||count[0]+ count[1]+count[2]+count[3]){
	int p = pars->posi[s]+1;
	gzwrite(oFileIcounts,&p,sizeof(int)); // new
	gzwrite(oFileIcounts,count,sizeof(int)*4); // new
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
suint **abcCounts::countNucs(const chunkyT *chk,int *keepSites){
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
	for(int l=0;l<chk->nd[s][n].l;l++){
	  //	  fprintf(stderr,"s:%d n:%d ret.l:%d l:%d seq:%c\n",s,n,chk->nd[s][n].l,l,chk->nd[s][n].seq[l]);
	  int allele = refToInt[chk->nd[s][n].seq[l]];
	
	  if(allele==4){
	    continue;
	  }
	  cnts[s][4*n+allele]++;
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
	for(int l=0;l<chk->nd[s][n].l;l++){
	  int allele = refToInt[chk->nd[s][n].seq[l]];	  
	  if(allele==4)
	    continue;
	  if(chk->nd[s][n].qs[l] < minQmat.matrix[n][allele]){
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
	for(int l=0;l<chk->nd[s][n].l;l++){
	  int allele = refToInt[chk->nd[s][n].seq[l]];
	  if(chk->nd[s][n].qs[l] < minQmat.matrix[n][0]||allele==4){
 
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
  pars->counts = countNucs(pars->chk,pars->keepSites);
  // fprintf(stderr,"%d\n",pars->keepSites[0]);
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
      else{
	suint *ps = pars->counts[s];
	for(int i=0;i<pars->nInd;i++){
	  int iSum = ps[i*4]+ps[i*4+1]+ps[i*4+2]+ps[i*4+3];
	  if(totSum>iSum){
	    pars->keepSites[s]=0;
	    break;
	  }
	}
      }
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
  
}
