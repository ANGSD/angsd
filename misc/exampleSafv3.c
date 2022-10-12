/*

norm<-function(x) x/sum(x)
sim1<-function(x,ndip){
    d<-c(0,1/(1:(2*ndip)))
    d <- c(0,1/(1:(2*ndip)))+ rnorm(2*ndip+1,sd=0.01)+1000
    d <- norm(d)
    log(d)-max(log(d))
}

d <- t(sapply(1:1000,sim1,ndip=20))

write.table(d,file='small.txt',row.names=F,col.names=F,quote=F)

gcc convert.c -lhts
cat small.txt|./a.out
 */

#include <htslib/bgzf.h>
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <stdlib.h>

#define NMAX 96 

int main(){
  FILE *fpindex = fopen("test.saf.idx","wb");
  BGZF *pos = bgzf_open("test.saf.pos.gz","wb");
  BGZF *saf = bgzf_open("test.saf.gz","wb");

  const char *magic="safv3";

  //write magic
  fwrite(magic,sizeof(char),8,fpindex);
  bgzf_write(pos,magic,8*sizeof(char));
  bgzf_write(saf,magic,8*sizeof(char));
  bgzf_flush(pos);bgzf_flush(saf);
  int64_t offs[2]={bgzf_tell(pos),bgzf_tell(saf)};
  char buffer[4096];
  float data[NMAX];
  size_t nbin=-1;
  size_t nsites =0;
  while(fgets(buffer,4096,stdin)){
    int at=0;
    data[at++] =atof(strtok(buffer,"\t\n "));
    char *tok=NULL;
    while(((tok=strtok(NULL,"\t\n ")))){
        data[at++] = atof(tok);
    }
    if(nbin==-1)
      nbin=at;
    if(at!=nbin) exit(1);
    bgzf_write(saf,data,sizeof(float)*nbin);
    nsites++;
  }
  for(int i=0;i<nsites;i++)
    bgzf_write(pos,&i,sizeof(int)*1);

  //write index
  nbin--;
  fwrite(&nbin,sizeof(size_t),1,fpindex);
  const char* nam="chr1X";
  size_t clen = strlen(nam);
  fwrite(&clen,sizeof(size_t),1,fpindex);
  fwrite(nam,sizeof(char),clen,fpindex);
  fwrite(&nsites,sizeof(size_t),1,fpindex);
  fwrite(offs,sizeof(int64_t),2,fpindex);
  fclose(fpindex);
  bgzf_close(pos);
  bgzf_close(saf);
  return 0;
}
