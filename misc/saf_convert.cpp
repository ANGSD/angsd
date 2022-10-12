#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "safcat.h"
#include "header.h"
#include "misc.h"

int main_text2safv4(int argc,char **argv){
  fprintf(stderr,"\t-> This will convert a banded text saf file into safv3 file\n");
  fprintf(stderr,"\t Examples:\n\t\t \'realSFS text2safv3 file.txt -outnames converted -P 10\'\n");
  fprintf(stderr,"\t\t \'realSFS cat -b file.list -outnames merged -P 10\'\n");
  fprintf(stderr,"\t\t Possible options: -P nThreads \n");
  char *outnames = NULL;
  char *infile = NULL;
  int nThreads = 1;
  size_t nchr = 0;
  int isbanded = -1;
  while(*argv){
    if(!strcasecmp(*argv,"-outnames"))
      outnames = strdup(*(++argv));
    else if(!strcasecmp(*argv,"-nchr"))
      nchr = atoi(*(++argv));
    else if(!strcasecmp(*argv,"-isbanded"))
      isbanded = atoi(*(++argv));
    else if(!strcasecmp(*argv,"-P"))
      nThreads = atoi(*(++argv));
    else
      infile = strdup(*argv);
    argv++;
  }
  if(outnames==NULL || nchr==0 ||isbanded==-1 ){
    fprintf(stderr,"\t-> You need to supply -outnames PREFIX -nchr INT -isbanded [0/1]\n");
    exit(0);
  }
  fprintf(stderr,"\t-> outnames: \'%s\' infile: %s nThreads: %d nchr: %lu isbanded: %d\n",outnames,infile,nThreads,nchr,isbanded);
  if(!outnames)
    return 0;
  
  BGZF *infilefp = NULL;
  infilefp = bgzf_open(infile,"rb");
  ASSERT(infilefp!=NULL);
   
  BGZF *outfileSAF =  openFileBG(outnames,SAF);
  BGZF *outfileSAFPOS =  openFileBG(outnames,SAFPOS);
  fprintf(stderr,"\t-> Number of threads: %d\n",nThreads);
  if(nThreads>1){
    fprintf(stderr,"\t-> Setting threads to: %d for both .saf.gz and saf.pos.gz\n",nThreads);
    bgzf_mt(outfileSAF,nThreads,256);
    bgzf_mt(outfileSAFPOS,nThreads,256);
  }
  FILE *outfileSAFIDX = openFile(outnames,SAFIDX);
  char buf[8]="safv4";
  my_bgzf_write(outfileSAF,buf,8);
  my_bgzf_write(outfileSAFPOS,buf,8);
  fwrite(buf,1,8,outfileSAFIDX);
  fwrite(&nchr,sizeof(size_t),1,outfileSAFIDX);
  int64_t offs[2];
  ASSERT(0==bgzf_flush(outfileSAFPOS));ASSERT(0==bgzf_flush(outfileSAF));
  offs[0] = bgzf_tell(outfileSAFPOS);
  offs[1] = bgzf_tell(outfileSAF);

  kstring_t *kstr = new kstring_t;
  kstr->s = NULL;kstr->l=kstr->m = 0;

  char *chr = NULL;
  std::vector<int> pos;
  size_t sumband = 0;
  int atline =0;
  while(bgzf_getline(infilefp,'\n',kstr)>0){
    atline++;
    //    fprintf(stderr,"pos.size(): %lu\n",pos.size());
    int bin_len[2];
    float flt[nchr+1];
    //  fprintf(stdout,"kstr.s: %s\n",kstr->s);
    char *tmpChr = strtok(kstr->s,"\t\n ");
    if(chr == NULL)
      chr = strdup(tmpChr);
    if(strcmp(tmpChr,chr)!=0&&sumband>0){
      //switch chr
      size_t clen = strlen(chr);
      fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
      fwrite(chr,1,clen,outfileSAFIDX);
      size_t NSITES = pos.size();
      fwrite(&NSITES,sizeof(size_t),1,outfileSAFIDX);
      fwrite(&sumband,sizeof(size_t),1,outfileSAFIDX);
      fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
      sumband =0;
      free(chr);
      chr= strdup(tmpChr);
      ASSERT(bgzf_write(outfileSAFPOS, pos.data(), sizeof(int)*pos.size())==sizeof(int)*pos.size());
      pos.clear();
      ASSERT(0==bgzf_flush(outfileSAF));
      ASSERT(0==bgzf_flush(outfileSAFPOS));
      offs[0] = bgzf_tell(outfileSAFPOS);
      offs[1] = bgzf_tell(outfileSAF);
    }
    int tmppos = atoi(strtok(NULL,"\t\n "))-1;
    //  fprintf(stderr,"tmppos: %d\n",tmppos);
    if(pos.size()>1&&tmppos<=pos[pos.size()-1]){
      fprintf(stderr,"\t-> Inputfile looks unsorted, please check tmppos: %d pos[pre]: %d atline: %d\n",tmppos,pos[pos.size()-1],atline);
      exit(0);
    }
    pos.push_back(tmppos);
    
    if(isbanded){
      bin_len[0] = atoi(strtok(NULL,"\t\n "));
      bin_len[1] = atoi(strtok(NULL,"\t\n "));
    }else{
      bin_len[0] = 0;
      bin_len[1] = nchr;
    }
    //    fprintf(stderr,"chr:%s pos: %d bin_len: (%d,%d) \n",chr,tmppos,bin_len[0],bin_len[1]);
    sumband += bin_len[1];
    int howmany = bin_len[1];
    if(isbanded==0)
      howmany++;
    for(int i=0;i<howmany;i++){
      flt[i] = atof(strtok(NULL,"\t\n "));
      //      fprintf(stderr,"flt[%d]: %f\n",i,flt[i]);
    }
    ASSERT(bgzf_write(outfileSAF, bin_len, sizeof(int)*2)==sizeof(int)*2);
    ASSERT(bgzf_write(outfileSAF, flt, sizeof(float)*bin_len[1])==sizeof(float)*bin_len[1]);
  }

  if(sumband>0){
    //switch chr
    size_t clen = strlen(chr);
    //  fprintf(stderr,"clen: %lu\n",clen);
    fwrite(&clen,sizeof(size_t),1,outfileSAFIDX);
    fwrite(chr,1,clen,outfileSAFIDX);
    size_t NSITES = pos.size();
    fwrite(&NSITES,sizeof(size_t),1,outfileSAFIDX);
    fwrite(&sumband,sizeof(size_t),1,outfileSAFIDX);
    fwrite(offs,sizeof(int64_t),2,outfileSAFIDX);
    //    fprintf(stderr,"nsites: %lu\n",pos.size());
    ASSERT(bgzf_write(outfileSAFPOS, pos.data(), sizeof(int)*pos.size())==sizeof(int)*pos.size());
    sumband =0;
  }
  
  if(outfileSAF) bgzf_close(outfileSAF);;
  if(outfileSAFPOS) bgzf_close(outfileSAFPOS);
  if(outfileSAFIDX) fclose(outfileSAFIDX);
  return 0;
}
