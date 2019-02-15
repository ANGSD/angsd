#include "mpileup.h"
#include "mUpPile.h"
#include <ctype.h>
#include <cassert>
#include "pooled_alloc.h"
extern tpool_alloc_t *tnodes;

mpileup::mpileup(int nInd_a,gzFile gz_a,const aMap* revMap_a,int minQ_a){
  fprintf(stderr,"\t-> You are using -pileup, this means:\n");
  fprintf(stderr,"\t-> 1) Internal positions both from front and back is coded to 255\n");
  fprintf(stderr,"\t-> 2) All mapping qualities (mapQ) are set to 30\n");
  fprintf(stderr,"\t-> 3) Program will not represent insertions, use raw BAM/CRAM for that\n");

  nInd = nInd_a;
  gz = gz_a;
  l = 8;
  revMap = revMap_a;
  minQ = minQ_a;
  buffer=original = (char*) calloc(l,sizeof(char));
}
mpileup::~mpileup(){
  free(original);
}

tNode **parseNd(char *line,int nInd,const char *delims,int minQ,char ref){

  tNode **ret =(tNode**) calloc(nInd,sizeof(tNode*));

  for(int i=0;i<nInd;i++) {

    char *tok = strtok_r(NULL,delims,&line);
    assert(tok);
    int seqDepth = atoi(tok);
    ret[i] = initNodeT(seqDepth);
    ret[i]->l = seqDepth;
    if(line[0]=='\t')
      continue;
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
	}
	else if(tok[inner]=='.' ){
	  ret[i]->seq[at++]=ref;
	  inner++;
	}
	else if(tok[inner]==','){
	  //  fprintf(stderr,"This wil...:%s sub:%c\n",tok,ref);
	  ret[i]->seq[at++]=ref;
	  inner++;
	}
	else{
	  fprintf(stderr,"This will never happen, ever...:%s sub:%s\n",tok,tok+inner);
	  exit(0);
	}

      }
#if 0
      for(int ii=0;ii<ret[i].l;ii++)
	fprintf(stderr,"ii:%d %c\n",ii,ret[i].seq[ii]);
#endif
      tok =strtok_r(NULL,delims,&line);
      for(unsigned inner=0;inner<strlen(tok);inner++){
	ret[i]->qs[inner] = tok[inner]-33;
	ret[i]->posi[inner] = 255;
	ret[i]->isop[inner] = 255;
	ret[i]->mapQ[inner] = 30;
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
  
  funkyPars * myfunky =funkyPars_init();
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
    myfunky->ref[nSites] = refToInt[strtok_r(NULL,delims,&buffer)[0]];
    myfunky->chk->nd[nSites] = parseNd(buffer,nInd,delims,minQ,myfunky->ref[nSites]);
    nSites++;
    changed =0;
  }
  buffer=original;
  while(aio::tgets(gz,&buffer,&l)) {
    if(buffer!=original) 
      original=buffer;
    char *tok = strtok_r(buffer,delims,&buffer);
    it=revMap->find(tok);
    if(it==revMap->end()){
      fprintf(stderr,"\t-> Problems finding reference \'%s\' in fai file\n",tok);
      exit(0);
    }
    if(lastRefId==-1)
      lastRefId = it->second;
    if(lastRefId!=it->second){
      changed =1;
      lastRefId = it->second; 
      myfunky->refId = lastRefId;
      break;
    }
    lastRefId = it->second;
    myfunky->refId = lastRefId;

    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims,&buffer))-1;
    char ref = strtok_r(NULL,delims,&buffer)[0];
    myfunky->ref[nSites] = refToInt[ref];

    myfunky->chk->nd[nSites++] = parseNd(buffer,nInd,delims,minQ,ref);
    
    buffer=original;
    if(nSites>=chunksize)
      break;
  }
  //  fprintf(stderr,"afterloop\n");
  myfunky->nInd=myfunky->chk->nSamples= nInd;
  myfunky->numSites = myfunky->chk->nSites=nSites;
  //fprintf(stderr,"\nchange2 %d\tnSites %d\tlastRefId\t%d %d %d\n",lastRefId,nSites,lastRefId,myfunky->refId,it->second);
  //  fflush(stderr);
  if(nSites==0 && changed == 0){
    funkyPars_destroy(myfunky);
    return(NULL);

  }
 return(myfunky);

}
