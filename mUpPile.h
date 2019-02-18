#pragma once
#include "argStruct.h"
#include "bambi_interface.h"

#include <htslib/kstring.h>
#include <htslib/hts.h>
#include "makeReadPool.h"
#include "abcGetFasta.h"


int uppile(int show,int nThreads,bufReader *rd,int NLINES,int nFiles,std::vector<regs> &regions,abcGetFasta *gf);
tNode *initNodeT(int l);


//typedef struct{

//}slist;//<simple list, just an array.
