#include <pthread.h>
#include "printRes.h"


//print all chunks that are next in order, must be called with printMut held
void printRes::loopyLoop(){
  while(!myMap.empty()){
    std::map<int,funkyPars*>::iterator it=myMap.begin();
    if(it->first<=lastChunk){
      fprintf(stderr,"Error, corruption in chunkordering in the printout function:\n");
      fprintf(stderr,"Last chunk=%d\t, smallest chunk now:%d\n",lastChunk,it->first);
      fflush(stderr);
      exit(1);
    }
    if(it->first!=lastChunk+1)
      break;//the next chunk in order is still being computed
    if(it->first!=it->second->chunkNumber){
      fprintf(stderr,"Error, corruption in chunkordering in the printout function:\n");
      fprintf(stderr,"[%s] first=%d second->chunkNumber=%d\n",__FUNCTION__,it->first,it->second->chunkNumber);
      fflush(stderr);
      exit(1);
    }
    funkyPars *fp = it->second;
    myMap.erase(it);
    lastChunk++;
    printFunky(fp);//fp doesn't exists anymore from this point
  }
}


//myMap is modified by worker threads in prop() under printMut,
//so reading its size must also hold printMut
size_t printRes::contains(){
  pthread_mutex_lock(&printMut);
  size_t n = myMap.size();
  pthread_mutex_unlock(&printMut);
  return n;
}

//called by the worker threads when a chunk is done.
//printing is kept under printMut so output is written by one thread at a time and in order
void printRes::prop(funkyPars *fp){
  pthread_mutex_lock(&printMut);
  if(!myMap.insert(std::pair<int,funkyPars *>(fp->chunkNumber,fp)).second){
    fprintf(stderr,"Error, corruption in chunkordering in the printout function:\n");
    fprintf(stderr,"[%s] chunkNumber=%d was received twice\n",__FUNCTION__,fp->chunkNumber);
    fflush(stderr);
    exit(1);
  }
  loopyLoop();
  pthread_mutex_unlock(&printMut);
}
