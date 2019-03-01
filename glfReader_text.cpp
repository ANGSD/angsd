#include <cassert>
#include "analysisFunction.h"
#include "glfReader_text.h"

glfReader_text::~glfReader_text(){
  free(buffer);
}


glfReader_text::glfReader_text(int nInd_a,gzFile gz_a,const aMap *revMap_a){
  fprintf(stderr,"\t-> NB text glf files are still in the developing stage, please report strange behaviour to thorfinn@binf.ku.dk\n\t-> NB So far only tested with -dosaf -domaf -domajorminor\n");
  gz=gz_a;
  revMap=revMap_a;
  l=128;
    
  original=buffer =(char *) calloc(l,sizeof(char));
  nInd = nInd_a;
}

void parselikes10(char *buffer,double *likes,int nInd,const char *delims){
  for(int i=0;i<nInd*10;i++){
    char *tsk = strtok_r(NULL,delims,&buffer);
    assert(tsk!=NULL);
    likes[i] = atof(tsk);
  }
  
}


funkyPars *glfReader_text::fetch(int chunksize){
  //  fprintf(stderr,"[glfreader_text]nind:%d\n",nInd);

  static const char *delims = "\t \n";
  static const char *delims2 = "_\t \n";
  double **likes = new double*[chunksize];
  
  funkyPars * myfunky =funkyPars_init();
  myfunky->posi = new int[chunksize];

  for(int s=0;s<chunksize;s++)
    likes[s] = new double[nInd*10];
  
  int nSites=0;
  static int positions =0;//every site is a new position across different chunks
  static int lastRefId =-1;
  static int changed =0;

  if(changed){
    //parse an entire site:

    myfunky->refId = lastRefId;
    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims2,&buffer))-1;
    parselikes10(buffer,likes[nSites],nInd,delims);    
    nSites++;
    changed =0;
  }

  buffer=original;
  while(aio::tgets(gz,&buffer,&l)) {
    if(buffer!=original)
      original=buffer;

    char *tok = strtok_r(buffer,delims2,&buffer);
    aMap ::const_iterator it = revMap->find(tok);
    if(it==revMap->end()){
      fprintf(stderr,"\t-> Problem finding chr:%s from faifile\n",tok);
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
    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims2,&buffer))-1;
 
    parselikes10(buffer,likes[nSites],nInd,delims);
    buffer=original;
    
    nSites++;
    if(nSites>=chunksize)
      break;
  }

  if(nSites<chunksize){
    for(int s=nSites;s<chunksize;s++)
      delete[] likes[s];
  }
  myfunky->nInd=nInd;
  myfunky->likes=likes;
  myfunky->numSites = nSites;
  
  if(nSites==0){
    funkyPars_destroy(myfunky);
    return(NULL);

  }
 return(myfunky);
}
