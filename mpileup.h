#include <zlib.h>
#include "shared.h"
class mpileup{
  int nInd;
  gzFile gz;
  char *original;
  char *buffer;
  int l;
  const aMap *revMap;
  int minQ;
public:
  funkyPars *fetch(int chunkSize);
  mpileup(int nInd_a,gzFile gz,const aMap *revMap,int minQ_a);
  ~mpileup();
};

