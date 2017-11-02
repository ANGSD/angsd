
#include <zlib.h>
#include "shared.h"
class glfReader_text{
private:
  const aMap *revMap;
  void getOptions(argStruct *arguments);
  void printArg(FILE *fp);
  int l;
  char *buffer;
  char *original;
  gzFile gz;
  int nInd;
public:
  glfReader_text(int nInd,gzFile gz_a,const aMap *revMap_a);
  ~glfReader_text();

  funkyPars *fetch(int chunkSize);
};

