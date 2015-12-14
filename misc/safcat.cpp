/*

  The functionality of this file, has replaced the old emOptim and testfolded.c programs.

  part of ANGSD

  GNU license or whetever its called

  thorfinn@binf.ku.dk

  fixme: minor leaks in structures related to the thread structs, and the append function.
  
  Its july 13 2013, it is hot outside

  april 13, safv3 added, safv2 removed for know. Will be reintroduced later.
  april 20, removed 2dsfs as special scenario
  april 20, split out the safreader into seperate cpp/h
  may 5, seems to work well now
*/

const char *SAF = ".saf.gz";
const char *SAFPOS =".saf.pos.gz";
const char *SAFIDX =".saf.idx";




#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <unistd.h>
#include <zlib.h>
#include <htslib/tbx.h>
#include "Matrix.hpp"
#include "safstat.h"
#include <libgen.h>
#include <algorithm>
#include "safreader.h"

int saf_cat(int argc,char **argv){
  fprintf(stderr,"\t-> This will cat together .saf files from angsd\n");
  fprintf(stderr,"\t-> regions has to be disjoint between saf files. This WONT be checked (alot) !\n");
  fprintf(stderr,"\t-> This has only been tested on safs for different chrs !\n");
  char *outnames = NULL;
  std::vector<persaf *> saf;
  while(*argv){
    if(!strcasecmp(*argv,"-outnames"))
      outnames = strdup(*(++argv));
    else
      saf.push_back(persaf_init<float>(*argv));
    argv++;
  }
  fprintf(stderr,"\t-> outnames: \'%s\' number of safs:%lu\n",outnames,saf.size());
  if(!outnames)
    return 0;
  BGZF *outfileSAF =  openFileBG(outnames,SAF);
  BGZF *outfileSAFPOS =  openFileBG(outnames,SAFPOS);
  FILE *outfileSAFIDX = openFile(outnames,SAFIDX);

  char buf[8]="safv3";
  bgzf_write(outfileSAF,buf,8);
  bgzf_write(outfileSAFPOS,buf,8);
  fwrite(buf,1,8,outfileSAFIDX);
  int64_t offs[2];
  offs[0] = bgzf_tell(outfileSAFPOS);
  offs[1] = bgzf_tell(outfileSAF);
  
  for(uint i=0;i<saf.size();i++){
    fprintf(stderr,"\r\t-> Merging %d/%lu ",i,saf.size());fflush(stderr);
    if(i==0)
      fwrite(&saf[0]->nChr,sizeof(size_t),1,outfileSAFIDX);
    if(i>0 && saf[0]->nChr!=saf[i]->nChr){
      fprintf(stderr,"\t-> Different number of samples in the different saf files, will exit\n");
      break;
    }
    for(myMap::iterator it=saf[i]->mm.begin();it!=saf[i]->mm.end();++it){
      int *ppos = new int[it->second.nSites];
      bgzf_seek(saf[i]->pos,it->second.pos,SEEK_SET);
      bgzf_seek(saf[i]->saf,it->second.saf,SEEK_SET);
      bgzf_read(saf[i]->pos,ppos,sizeof(int)*it->second.nSites);
      bgzf_write(outfileSAFPOS,ppos,sizeof(int)*it->second.nSites);
      float flt[saf[0]->nChr+1];
      for(uint s=0;s<it->second.nSites;s++){
	bgzf_read(saf[i]->saf,flt,sizeof(float)*(saf[0]->nChr+1));
	bgzf_write(outfileSAF,flt,sizeof(float)*(saf[0]->nChr+1));
      }
      delete [] ppos;

      size_t clen = strlen(it->first);
      fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
      fwrite(it->first,1,clen,outfileSAFIDX);
      fwrite(&it->second.nSites,sizeof(size_t),1,outfileSAFIDX);
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);

      offs[0] = bgzf_tell(outfileSAFPOS);
      offs[1] = bgzf_tell(outfileSAF);
    }
    fprintf(stderr,"\n");
  }

  for(int i=0;i<saf.size();i++)
    persaf_destroy(saf[i]);
  
  fclose(outfileSAFIDX);
  bgzf_close(outfileSAF);
  bgzf_close(outfileSAFPOS);
  free(outnames);
  return 0;
}
