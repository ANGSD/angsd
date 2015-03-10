#include "mpileup.h"
#include "mUpPile.h"
#include <ctype.h>
#include <cassert>
mpileup::mpileup(int nInd_a,gzFile gz_a,int bpl,const aMap* revMap_a,int minQ_a){
  nInd = nInd_a;
  gz = gz_a;
   len = bpl;
  revMap = revMap_a;
  minQ = minQ_a;
  buffer=original = new char [bpl];
}
mpileup::~mpileup(){
  fprintf(stderr,"mpileup:\n");
 delete [] original;
}

tNode **parseNd(char *line,int nInd,const char *delims,int minQ){

  tNode **ret = new tNode*[nInd];

  for(int i=0;i<nInd;i++) {

    char *tok = strtok_r(NULL,delims,&line);
    assert(tok);
    int seqDepth = atoi(tok);
    ret[i] = tNode_init1(seqDepth);
    ret[i]->l = seqDepth;
    if(ret[i]->l==0){
      strtok_r(NULL,delims,&line);
      strtok_r(NULL,delims,&line);
    }else{
      char *tok = strtok_r(NULL,delims,&line);
      int at =0;
      unsigned inner =0;
      while(inner<strlen(tok)){
	//	fprintf(stderr,"inner:%u len:%zu rest:%s\n",inner,strlen(tok),tok+inner);
	if(tok[inner]=='$')
	  inner++;
	else if (tok[inner]=='^'){
	  inner+=2;
	}else if(isalpha(tok[inner])){
	  ret[i]->seq[at] =tok[inner];
	  at++;inner++;
	}else if(tok[inner]=='-'||tok[inner]=='+'){
	  static int show =1;
	  if(show-->0){
	    fprintf(stderr,"\t-> Problem with indels in -pileup: please use raw binary for bam files and NOT pileup files for inserts:%s\n",tok);
	    fprintf(stderr,"\t-> Parsing will only use non-indel information. This msg will appear only once\n");
	  }
	  int opLen=0;
	  char *bb = new char[strlen(tok)];
	  sscanf(tok+inner+1,"%d%s",&opLen,bb);
	  //fprintf(stderr,"opLen:%d\n",opLen);
	  if(opLen<10)
	    inner +=2+opLen;
	  else if(opLen<100)
	    inner +=3+opLen;
	  delete [] bb;
	}else if(tok[inner]=='*'){
	  ret[i]->seq[at++]='N';
	  inner++;
	}else{
	  fprintf(stderr,"This will never happen, ever...:%s sub:%s\n",tok,tok+inner);
	}

      }
#if 0
      for(int ii=0;ii<ret[i].l;ii++)
	fprintf(stderr,"ii:%d %c\n",ii,ret[i].seq[ii]);
#endif
      tok =strtok_r(NULL,delims,&line);
      for(unsigned inner=0;inner<strlen(tok);inner++){
	ret[i]->qs[inner] = tok[inner]-33;
	//	fprintf(stderr,"qs[%d]:%d\n",inner,ret[i].qs[inner]);
	if(ret[i]->qs[inner]<minQ)
	  ret[i]->seq[inner] = 'n';
      }

    }

  }
  return ret;
}

funkyPars *mpileup::fetch(int chunksize){
  static const char *delims = "\t \n";
  
  funkyPars * myfunky =allocFunkyPars();
  myfunky->posi = new int[chunksize];
  myfunky->ref = new char[chunksize];
  myfunky->chk = new chunkyT;
  myfunky->chk->nd = new tNode **[chunksize];
  myfunky->chk->refPos = NULL;
  int nSites=0;
  static int lastRefId =-1;
  static int changed =0;
  aMap ::const_iterator it;
  
  if(changed){
    //parse an entire site:
    myfunky->refId = lastRefId;
    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims,&buffer))-1;
    myfunky->ref[nSites] = strtok_r(NULL,delims,&buffer)[0];
    myfunky->chk->nd[nSites++] = parseNd(buffer,nInd,delims,minQ);
    
    changed =0;
  }
  buffer=original;
  while(gzgets(gz,buffer,len)) {
    //    fprintf(stderr,"buf:%s strlen:%zu\n",buffer,strlen(buffer));
    if(strlen(buffer)==len-1){
      fprintf(stderr,"\t-> Increase -bytesPerLine value\n");
      exit(0);
    }
    char *tok = strtok_r(buffer,delims,&buffer);
    it=revMap->find(tok);
    if(it==revMap->end()){
      fprintf(stderr,"Probleme finding chr \'%s\' in fai file\n",tok);
      exit(0);
    }
    if(lastRefId==-1)
      lastRefId = it->second;
    if(lastRefId!=it->second){
      changed =1;
      lastRefId = it->second; 
      break;
    }
    lastRefId = it->second;
    myfunky->refId = lastRefId;

    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims,&buffer))-1;
    myfunky->ref[nSites] = strtok_r(NULL,delims,&buffer)[0];
    myfunky->chk->nd[nSites++] = parseNd(buffer,nInd,delims,minQ);
    
    buffer=original;
    if(nSites>=chunksize)
      break;
  }
  //  fprintf(stderr,"afterloop\n");
  myfunky->nInd=myfunky->chk->nSamples= nInd;
  myfunky->numSites = myfunky->chk->nSites=nSites;

  if(nSites==0){

    deallocFunkyPars(myfunky);
    return(NULL);

  }
 return(myfunky);

}
