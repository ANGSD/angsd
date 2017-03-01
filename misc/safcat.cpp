/*
  small utility functions for saf files
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
#include "realSFS_args.h"

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
  my_bgzf_write(outfileSAF,buf,8);
  my_bgzf_write(outfileSAFPOS,buf,8);
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
      my_bgzf_seek(saf[i]->pos,it->second.pos,SEEK_SET);
      my_bgzf_seek(saf[i]->saf,it->second.saf,SEEK_SET);
      my_bgzf_read(saf[i]->pos,ppos,sizeof(int)*it->second.nSites);
      my_bgzf_write(outfileSAFPOS,ppos,sizeof(int)*it->second.nSites);
      float flt[saf[0]->nChr+1];
      for(uint s=0;s<it->second.nSites;s++){
	my_bgzf_read(saf[i]->saf,flt,sizeof(float)*(saf[0]->nChr+1));
	my_bgzf_write(outfileSAF,flt,sizeof(float)*(saf[0]->nChr+1));
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


int saf_check(int argc,char **argv) {
  fprintf(stderr,"\t-> ./realSFS check your.saf.idx\n\t-> This will check that the positions are ordered correctly\n");
  if(argc==0)
    return 0;
  args *pars = getArgs(argc,argv);
  assert(pars->saf.size()==1);
  for(myMap::iterator it=pars->saf[0]->mm.begin();it!=pars->saf[0]->mm.end();++it){
    int *ppos = new int[it->second.nSites];
    my_bgzf_seek(pars->saf[0]->pos,it->second.pos,SEEK_SET);
    my_bgzf_read(pars->saf[0]->pos,ppos,sizeof(int)*it->second.nSites);
    for(int i=1;i<it->second.nSites;i++){
      if(ppos[i]<=ppos[i-1])
	fprintf(stderr,"\t-> problems with unsorted saf file chromoname: \'%s\' pos[%d]:%d vs posd[%d-1]:%d\n",it->first,i,ppos[i],i,ppos[i-1]);
    }
      
    delete [] ppos;
  }
  destroy_args(pars);
  return 0;
}
