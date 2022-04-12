/*
  small utility functions for saf files
*/



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
#include "safcat.h"
#include "realSFS_args.h"
#include <fstream>


std::vector<char*> getFilenames(const char * name,int nInd){
  fprintf(stderr,"\t-> Reading saf.idx from: %s\n",name);
  if(strchr(name,'\r')){
    fprintf(stderr,"\t\t-> Filelist contains carriage return. Looks like a windows file please remove hidden \'\r\' from filelist\n");
    exit(0);
  }
  
  if(!fexists(name)){
    fprintf(stderr,"[%s]\t-> Problems opening file: %s\n",__FUNCTION__,name);
    exit(0);
  }
  const char* delims = " \t";
  std::vector<char*> ret;
  std::ifstream pFile(name,std::ios::in);

  char buffer[LENS];
  while(!pFile.eof()){
    pFile.getline(buffer,LENS);
    char *tok = strtok(buffer,delims);
    while(tok!=NULL){
      if(tok[0]!='#')
	ret.push_back(strdup(buffer));
      tok = strtok(NULL,delims);
    }
  }
  if(nInd>0) {
     if(ret.size()<nInd)
      fprintf(stderr,"\t-> Number of samples is smaller than subset requested %lu vs %d\n",ret.size(),nInd);
    else{
      //   fprintf(stderr,"\t-> Will remove tail of filename list\n");
      for(int ii=nInd;ii<ret.size();ii++)
	free(ret[ii]);
      ret.erase(ret.begin()+nInd,ret.end());//we don't free this memory, it doesn't really matter
      // fprintf(stderr,"\t->  filename list now contains only: %lu\n",ret.size());
    }
#if 0
     for(size_t ii=0;ii<ret.size();ii++)
       fprintf(stderr,"%zu->%s\n",ii,ret[ii]);
     fprintf(stderr,"\n");
#endif
  }

  return ret;
}

//this function will cat together saf files. It is not super clever since it decompresses and compresses. We should be able to just copy
//ideally this should functionality should be included in the saf_init context such that the analyses parts doesnt care about the actual backend of the data.
int saf_cat(int argc,char **argv){
  fprintf(stderr,"\t-> This will cat together .saf files from angsd\n");
  fprintf(stderr,"\t-> regions has to be disjoint between saf files. This WONT be checked (alot) !\n");
  fprintf(stderr,"\t-> This has only been tested on safs for different chrs !\n");
  fprintf(stderr,"\t Examples:\n\t\t \'realSFS cat chr1.saf.idx chr2.saf.idx -outnames merged -P 10\'\n");
  fprintf(stderr,"\t\t \'realSFS cat -b file.list -outnames merged -P 10\'\n");
  fprintf(stderr,"\t\t Possible options: -b file.list -P nThreads \n");
  char *outnames = NULL;
  std::vector<persaf *> saf;
  std::vector<char *> saffiles; //<- when user is using -b
  int nThreads = 1;
  while(*argv){
    if(!strcasecmp(*argv,"-outnames"))
      outnames = strdup(*(++argv));
    else if(!strcasecmp(*argv,"-P"))
      nThreads = atoi(*(++argv));
    else if(!strcasecmp(*argv,"-b")){
      char *shifted = *(++argv);
      fprintf(stderr,"\t-> Will read in filenames frome file: %s\n",shifted);
      saffiles = getFilenames(shifted,0);
      fprintf(stderr,"\t-> Has added: %lu\n",saffiles.size());
    }else
      saffiles.push_back(strdup(*argv));
    argv++;
  }
  if(outnames==NULL){
    fprintf(stderr,"\t-> You need to supply -outnames PREFIX\n");
    exit(0);
  }
  for(int i=0;i<saffiles.size();i++)
    saf.push_back(persaf_init<float>(saffiles[i],0));

  fprintf(stderr,"\t-> outnames: \'%s\' number of safs:%lu nThreads: %d\n",outnames,saf.size(),nThreads);
  
  if(!outnames)
    return 0;
  BGZF *outfileSAF =  openFileBG(outnames,SAF);
  BGZF *outfileSAFPOS =  openFileBG(outnames,SAFPOS);
  fprintf(stderr,"\t-> Number of threads: %d\n",nThreads);
  if(nThreads>1){
    fprintf(stderr,"\t-> Setting threads to: %d for both .saf.gz and saf.pos.gz\n",nThreads);
    bgzf_mt(outfileSAF,nThreads,256);
    bgzf_mt(outfileSAFPOS,nThreads,256);
  }
  FILE *outfileSAFIDX = openFile(outnames,SAFIDX);
  int safversion = saf[0]->version;
  char buf[8]="safv3";
  if(safversion==3)
    buf[4]='4';
  my_bgzf_write(outfileSAF,buf,8);
  my_bgzf_write(outfileSAFPOS,buf,8);
  fwrite(buf,1,8,outfileSAFIDX);
  int64_t offs[2];
  assert(0==bgzf_flush(outfileSAFPOS));assert(0==bgzf_flush(outfileSAF));
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
    if(i>1 &&safversion!=saf[i]->version){
      fprintf(stderr,"\t-> Different safversion in the different saf files, will exit\n");
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
	if(safversion==2){
	  my_bgzf_read(saf[i]->saf,flt,sizeof(float)*(saf[0]->nChr+1));
	  my_bgzf_write(outfileSAF,flt,sizeof(float)*(saf[0]->nChr+1));
	}else if(safversion==3){
	  int band[2];
	  my_bgzf_read(saf[i]->saf, band, 2*sizeof(int));
	  my_bgzf_read(saf[i]->saf, flt, band[1]*sizeof(float));
	  my_bgzf_write(outfileSAF,band,2*sizeof(int));
	  my_bgzf_write(outfileSAF,flt,band[1]*sizeof(float));
	}
      }
      delete [] ppos;

      size_t clen = strlen(it->first);
      fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
      fwrite(it->first,1,clen,outfileSAFIDX);
      fwrite(&it->second.nSites,sizeof(size_t),1,outfileSAFIDX);
      if(safversion==3)
	fwrite(&it->second.sumBand,sizeof(size_t),1,outfileSAFIDX);
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
      assert(0==bgzf_flush(outfileSAFPOS));assert(0==bgzf_flush(outfileSAF));
      offs[0] = bgzf_tell(outfileSAFPOS);
      offs[1] = bgzf_tell(outfileSAF);
    }
    fprintf(stderr,"\n");
  }

  for(int i=0;i<saf.size();i++)
    persaf_destroy(saf[i]);
  for(int i=0;i<saffiles.size();i++)
    free(saffiles[i]);
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
