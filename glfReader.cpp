#include <cassert>
#include "glfReader.h"

glfReader::~glfReader(){

};

glfReader::glfReader(int &nInd_a,gzFile gz_a,int isSim_a){
  //  fprintf(stderr,"nind:%d from:%d to:%d gz:%p isSim:%d\n",nInd_a,from_a,to_a,gz_a,isSim_a);
  nInd = nInd_a;
   gz=gz_a;
  isSim = isSim_a;
}


funkyPars *glfReader::fetch(int chunkSize){
  funkyPars *r = allocFunkyPars();  
  r->likes=new double*[chunkSize];
  r->posi=new int[chunkSize];
  r->anc = new char[chunkSize];
  memset(r->anc,0,chunkSize);
  r->refId = 0;
  static int pos = 0;
  int i;
  for(i=0;i<chunkSize;i++){
    double *lk = new double [10*nInd];
    size_t bytesRead = gzread(gz,lk,sizeof(double)*10*nInd);
    if(bytesRead==0){
      delete [] lk;
      break;
    }else if(bytesRead!=sizeof(double)*10*nInd){
      fprintf(stderr,"\n[%s:%s():%d] Error reading full chunk: bytesRead:%zu expected:%zu will exit\n",__FILE__,__FUNCTION__,__LINE__,bytesRead,sizeof(double)*10*nInd);
      exit(0);
    }else{ 
      r->likes[i] = lk;
      r->posi[i] = i+pos;
    }
  }
  if(i==0){
    delete [] r->likes;
    delete [] r->posi;
    delete [] r->anc;
    delete r;
    return NULL;
  }
  r->nInd=nInd;
  
  r->numSites=i;
  pos += i;
  
  return r;
}

