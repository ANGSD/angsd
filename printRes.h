#pragma once
#include <pthread.h>
#include <map>
#include "shared.h"

extern void printFunky(funkyPars *p);

//collects finished chunks from the worker threads and prints them in chunkNumber order
class printRes{

  pthread_mutex_t printMut;//protects myMap and lastChunk, and makes printing exclusive

  int lastChunk;//chunkNumber of the last printed chunk
  std::map<int, funkyPars *> myMap;//finished chunks waiting for an earlier chunk
  void loopyLoop();


public:
  size_t contains();//number of chunks waiting to be printed, takes printMut
  //this is anders mindfuck
  printRes() {
    lastChunk=0;
    pthread_mutex_init(&printMut,NULL);
  };

  void prop(funkyPars *fp);
};
