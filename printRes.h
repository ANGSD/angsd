//#pragma once
#include <map>
#include "shared.h"

//#define USE_SPINLOCK 1

extern void printFunky(funkyPars *p);

class printRes{

  pthread_mutex_t printMut;

  int lastChunk;
  std::map<int, funkyPars *> myMap;
  void popPars(funkyPars *fp);
  void loopyLoop();

  
public:
  size_t contains() {return myMap.size();};
  int getlast() {return lastChunk;}
  int getfirst() {
    if(myMap.size()==0) 
      return 0; 
    else {
      std::map<int,funkyPars*>::iterator it=myMap.begin(); 
      return it->first;
    }
  }
  //this is anders mindfuck
  printRes() {
    lastChunk=0;  
    pthread_mutex_init(&printMut,NULL);
  };

  void prop(funkyPars *fp);
  int isEmpty() { 
    return (myMap.size()==0);
  }
};
