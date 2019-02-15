#include <cassert>
#include "glfReader.h"
#include <cfloat>

glfReader::~glfReader(){

};

glfReader::glfReader(int &nInd_a,gzFile gz_a,int nGL_a,int isSim_a){
  //  fprintf(stderr,"nind:%d from:%d to:%d gz:%p isSim:%d\n",nInd_a,from_a,to_a,gz_a,isSim_a);
  nInd = nInd_a;
  gz = gz_a;
  nGL = nGL_a;
  isSim = isSim_a;
}


funkyPars *glfReader::fetch(int chunkSize){
  funkyPars *r = funkyPars_init();  
  r->likes=new double*[chunkSize];
  r->posi=new int[chunkSize];
  r->anc = new char[chunkSize];
  memset(r->anc,0,chunkSize);
  r->refId = 0;
  static int pos = 0;
  int l;
  for(l=0;l<chunkSize;l++){
    double *lk = new double [10*nInd];
    size_t bytesRead = gzread(gz,lk,sizeof(double)*nGL*nInd);
    if(bytesRead==0){
      delete [] lk;
      break;
    }else if(bytesRead!=sizeof(double)*nGL*nInd){
      fprintf(stderr,"\n[%s:%s():%d] Error reading full chunk: bytesRead:%zu expected:%zu will exit\n",__FILE__,__FUNCTION__,__LINE__,bytesRead,sizeof(double)*nGL*nInd);
      exit(0);
    }else{
      if(nGL == 10){
	r->likes[l] = lk;
	r->posi[l] = l+pos;
      }else if(nGL == 3){
	r->likes[l] = new double [10*nInd];
	for(int i=0; i<nInd; i++){
	  for(int g=0; g<10; g++)
	    r->likes[l][i*10+g] = -DBL_MAX;

	  r->likes[l][i*10+angsd::majorminor[refToInt['A']][refToInt['A']]] = lk[i*3+0];
	  r->likes[l][i*10+angsd::majorminor[refToInt['A']][refToInt['C']]] = lk[i*3+1];
	  r->likes[l][i*10+angsd::majorminor[refToInt['C']][refToInt['C']]] = lk[i*3+2];
	}
	r->posi[l] = l+pos;
	delete [] lk;
      }else{
	fprintf(stderr,"\n[%s:%s():%d] Invalid number of GLs: %d. Will exit\n",__FILE__,__FUNCTION__,__LINE__,nGL);
	exit(0);
      }
    }
  }
  if(l==0){
    delete [] r->likes;
    delete [] r->posi;
    delete [] r->anc;
    delete r;
    return NULL;
  }
  r->nInd=nInd;
  
  r->numSites=l;
  pos += l;
  
  return r;
}

