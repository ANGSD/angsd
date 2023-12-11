#pragma once
#include <map>
#include "glfReader.h"
#include "vcfReader.h"
#include "bcfReader.h"
#include "beagleReader.h"
#include "bgenReader.h"
#include "argStruct.h"
#include "mpileup.h"
#include "parseArgs_bambi.h"
#include "bammer_main.h"
#include "glfReader_text.h"
class multiReader{
private:
  int nLines;
  glfReader *myglf;
  glfReader_text *myglf_text;
  vcfReader *myvcf;
  bcfReader *mybcf;
  beagle_reader *bglObj;
  bgenReader *bgenObj;
  mpileup *mpil;
  char *fname;
  int intName; // intrepret SNP name as chr_pos
  gzFile gz;
  //  int type;//0=soap;1=glf;2=glfclean,3=tglf,4=simfiles
  void getOptions(argStruct *arguments);
  void printArg( FILE *fp,argStruct *);
  std::map<char*,int> mMap;
  int isSim;
  int nInd;
  int minQ;
  argStruct *args;
  aMap *revMap;
  int pl_or_gl; // <-pl: pl_or_gl=0,gl:pl_or_gl=1
  uint32_t bcf_tags_to_use;
public:
  multiReader(argStruct *arg);
  funkyPars *fetch();
  ~multiReader();
};
