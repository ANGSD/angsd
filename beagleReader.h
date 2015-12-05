
#include <zlib.h>
#include "shared.h"
class beagle_reader{
private:
  const aMap *revMap;
  void getOptions(argStruct *arguments);
  void printArg(FILE *fp);
  int l;
  char *buffer;
  char *original;
  gzFile gz;
  int intName;
  int nInd;
public:
  beagle_reader(gzFile gz_a,const aMap *revMap_a,int,int&);
  ~beagle_reader();

  funkyPars *fetch(int chunkSize);
};

