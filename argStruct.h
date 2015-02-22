#pragma once
#include <cstring>


//little struct for keeping information of regions to extract
typedef struct{
  int refID;
  int start;
  int stop;
}regs;


#include "bams.h"


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef std::map<const char *,int,ltstr> aMap;


typedef struct{
  int nInd;//number of inds inferred from filelists
  int argc;
  char **argv;
  int inputtype;//
  int *usedArgs; //array of ints telling if args been used
  FILE *argumentFile; //logfile
  char *outfiles; //prefix output
  const aHead *hd;
  const aMap *revMap;
  char *infile;//contains, the -bam fname,-glf fname, -pileup fname
  std::vector<char *> nams;//contains either the above or the contents of -bam;
  std::vector<regs> regions;//regions to use -r/-rf when using bamfiles
  int nLines;//nLines;
  int show;
}argStruct;
