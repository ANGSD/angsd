#include <pthread.h>
#include <iostream>
#include "printRes.h"


void printRes::popPars(funkyPars *fp){
  /*
  if(fp->killSig)
    fprintf(stderr,"received killsig in [%s]\n",__FUNCTION__);
  else
    fprintf(stderr,"\t->printing chunkNumber:%d current number of chunks waiting to be printed: %lu\n",fp->chunkNumber,myMap.size());
  */

  // fflush(stderr);
  lastChunk++;

  //unlock the mutex while printing and allowing new data to be added
      /*
#ifdef USE_SPINLOCK  
  pthread_spin_unlock(&printMut);
#else
  pthread_mutex_unlock(&printMut);
#endif
    */

  printFunky(fp);//fp doesn't exists anymore from this point
  
  //lock again to avoid strange things happening. (actually this migth not be nescessary)
     /*
#ifdef USE_SPINLOCK  
  pthread_spin_lock(&printMut);
#else
  pthread_mutex_lock(&printMut);
#endif
     */
  //fprintf(stderr,"")

}

void printRes::loopyLoop(){
  if(!myMap.empty()){
    std::map<int,funkyPars*>::iterator it=myMap.begin();
    if(it->first<lastChunk){
      fprintf(stderr,"Error, corruption in chunkordering in the printout function:\n");
      fprintf(stderr,"Last chunk=%d\t, smallest chunk now:%d\n",lastChunk,it->first);      
      fflush(stderr);
      exit(0);
    }
    if(it->first==(lastChunk+1)){
      if(it->first!=it->second->chunkNumber){
	fprintf(stderr,"Error, corruption in chunkordering in the printout function:\n");
	fprintf(stderr,"[%s] first=%d second->chunkNumber=%d\n",__FUNCTION__,it->first,it->second->chunkNumber);
	fflush(stderr);
	exit(0);
      }
      popPars(it->second);
      myMap.erase(it);
      loopyLoop();
    }
    
  }
}


void printRes::prop(funkyPars *fp){
  //  fprintf(stderr,"[%s]. propping chunking number:%d\n",__FUNCTION__,fp->chunkNumber);
#ifdef USE_SPINLOCK  
  pthread_spin_lock(&printMut);
#else
  pthread_mutex_lock(&printMut);
#endif
  //  fprintf(stderr,"after lock\n");
  // fflush(stderr);
  // contains();
  // if(fp->chunkNumber-1==lastChunk){
  // popPars(fp);
    //}else{
    //  fprintf(stderr,"inserting :%d\n",fp->chunkNumber);
    //    fflush(stderr);
    myMap.insert(std::pair<int,funkyPars *>(fp->chunkNumber,fp));
    //}
  loopyLoop();
#ifdef USE_SPINLOCK  
  pthread_spin_unlock(&printMut);
#else
  pthread_mutex_unlock(&printMut);
#endif
}

