
#include <cassert>
#include "analysisFunction.h"
#include "beagleReader.h"

beagle_reader::~beagle_reader(){
  free(buffer);
}


beagle_reader::beagle_reader(gzFile gz_a,const aMap *revMap_a,int intName_a,int &nInd_a){
  gz=gz_a;
  revMap=revMap_a;
  l=128;
  intName = intName_a;
  
  original=buffer =(char *) calloc(l,sizeof(char));
  const char *delims = "\t \n";

  int nCol=1;

  aio::tgets(gz,&buffer,&l);
  if(buffer!=original)
      original=buffer;
  strtok_r(buffer,delims,&buffer);
  while(strtok_r(NULL,delims,&buffer))
    nCol++;
  if(nCol % 3 ){
    fprintf(stderr,"\t-> Number of columns should be a multiple of 3, nCol=%d\n",nCol);
    exit(0);
  } 
  nInd_a=nCol/3-1;
  nInd = nInd_a;
}

void parsepost(char *buffer,double *post,int nInd,const char *delims){
  for(int i=0;i<nInd*3;i++){
    char *tsk = strtok_r(NULL,delims,&buffer);
    assert(tsk!=NULL);
    post[i] = atof(tsk);
  }
  
}


funkyPars *beagle_reader::fetch(int chunksize){
  char refToChar[256] = {
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
  };
  

  static const char *delims = "\t \n";
  static const char *delims2 = "_\t \n";
  double **post = new double*[chunksize];
  
  funkyPars * myfunky =funkyPars_init();
  myfunky->posi = new int[chunksize];
  myfunky->major = new char[chunksize];
  myfunky->minor = new char[chunksize];
  
  
  for(int s=0;s<chunksize;s++)
    post[s] = new double[nInd*3];
  
  int nSites=0;
  static int positions =0;//every site is a new position across different chunks
  static int lastRefId =-1;
  static int changed =0;
 READAGAIN:
  if(changed){
    //parse an entire site:
    //fprintf(stdout,"nSites %d\n",nSites);
    myfunky->refId = lastRefId;
    //  fprintf(stdout,"BUF= %s\n",buffer);
    myfunky->posi[nSites] = atoi(strtok_r(NULL,delims2,&buffer))-1;

    myfunky->major[nSites] = refToChar[strtok_r(NULL,delims,&buffer)[0]];
    myfunky->minor[nSites] = refToChar[strtok_r(NULL,delims,&buffer)[0]];
 
    parsepost(buffer,post[nSites],nInd,delims);    

    nSites++;
    changed =0;
  }
  buffer=original;
  while(aio::tgets(gz,&buffer,&l)) {
    if(buffer!=original)
      original=buffer;
    if(intName){

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
	if(nSites==0) // if chromosome if finish then read next one
	  goto READAGAIN; 
	break;
      }
      lastRefId = it->second;
      myfunky->refId = lastRefId;
      myfunky->posi[nSites] = atoi(strtok_r(NULL,delims2,&buffer))-1;
    }
    else{
      char *tok = strtok_r(buffer,delims,&buffer);
      myfunky->refId = 0;
      myfunky->posi[nSites] = positions++;
    }
    
    myfunky->major[nSites] = refToChar[strtok_r(NULL,delims,&buffer)[0]];
    myfunky->minor[nSites] = refToChar[strtok_r(NULL,delims,&buffer)[0]];
    
    parsepost(buffer,post[nSites],nInd,delims);
    buffer=original;
    
    nSites++;
    if(nSites>=chunksize)
      break;
  }

  if(nSites<chunksize){
    for(int s=nSites;s<chunksize;s++)
      delete[] post[s];
  }
  myfunky->nInd=nInd;
  myfunky->post=post;
  myfunky->numSites = nSites;
  
  if(nSites==0 & changed==0){
    fprintf(stdout,"Done reading beagle\n");
    funkyPars_destroy(myfunky);
    return(NULL);

  }
 return(myfunky);
}
