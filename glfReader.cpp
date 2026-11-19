#include "glfReader.h"
#include <cfloat>

glfReader::~glfReader(){
}

glfReader::glfReader(int &nInd_a,gzFile gz_a,int nGL_a){
  nInd = nInd_a;
  gz = gz_a;
  nGL = nGL_a;
  pos = 0;
  if(nGL!=10&&nGL!=3){
    fprintf(stderr,"\t-> [%s] Invalid number of GLs: %d. Will exit\n",__FILE__,nGL);
    exit(1);
  }
}


funkyPars *glfReader::fetch(int chunkSize){
  const size_t siteBytes = sizeof(double)*nGL*nInd;
  funkyPars *r = funkyPars_init();
  r->likes=new double*[chunkSize];
  r->posi=new int[chunkSize];
  r->refId = 0;
  int l;
  for(l=0;l<chunkSize;l++){
    double *lk = new double [nGL*nInd];
    int bytesRead = gzread(gz,lk,siteBytes);
    if(bytesRead==0){
      delete [] lk;
      break;
    }
    if(bytesRead<0){
      int err;
      fprintf(stderr,"\t-> Problem reading glf file: %s. Will exit\n",gzerror(gz,&err));
      exit(1);
    }
    if((size_t)bytesRead!=siteBytes){
      fprintf(stderr,"\t-> glf file ends with a partial site: read %d bytes, expected %zu (-nInd %d x %d likelihoods x %zu bytes). Is the file truncated or -nInd wrong? Will exit\n",bytesRead,siteBytes,nInd,nGL,sizeof(double));
      exit(1);
    }
    if(nGL == 10){
      r->likes[l] = lk;
    }else{
      //glf3 holds AA, AC, CC
      r->likes[l] = new double [10*nInd];
      for(int i=0; i<nInd; i++){
	for(int g=0; g<10; g++)
	  r->likes[l][i*10+g] = -DBL_MAX;
	r->likes[l][i*10+angsd::majorminor[refToInt['A']][refToInt['A']]] = lk[i*3+0];
	r->likes[l][i*10+angsd::majorminor[refToInt['A']][refToInt['C']]] = lk[i*3+1];
	r->likes[l][i*10+angsd::majorminor[refToInt['C']][refToInt['C']]] = lk[i*3+2];
      }
      delete [] lk;
    }
    r->posi[l] = pos+l;
  }
  r->nInd=nInd;
  r->numSites=l;
  if(l==0){
    funkyPars_destroy(r);
    return NULL;
  }
  pos += l;
  return r;
}
