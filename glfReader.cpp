#include "glfReader.h"

glfReader::~glfReader(){

};

glfReader::glfReader(int &nInd_a,int nInd2_a,int from_a, int to_a,gzFile gz_a,int isSim_a){
  //  fprintf(stderr,"nind:%d nind2:%d from:%d to:%d gz:%p isSim:%d\n",nInd_a,nInd2_a,from_a,to_a,gz_a,isSim_a);
  nInd = nInd_a;
  nInd2 = nInd2_a;
  from = from_a;
  to = to_a;
  gz=gz_a;
  isSim = isSim_a;
  if(from!=-1&&to!=-1){
    nInd2=to-from+1;

  }

  if(nInd2!=-1)
    nInd_a = nInd2;
  else
    nInd_a = nInd;  
  nInd = nInd_a;
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
      fprintf(stderr,"\t-> Error reading full chunk\n");
      exit(0);
    }else{ 
      if(from==-1)
	r->likes[i] = lk;
      else{ 
	double *lk2=new double[10*nInd2];
	memcpy(lk2,lk+10*from,sizeof(double)*10*nInd2);
	r->likes[i] = lk2;
	delete [] lk;
      }
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
  if(nInd2==-1)
    r->nInd=nInd;
  else
    r->nInd=nInd2;
  r->numSites=i;
  pos += i;
  
  return r;
}

