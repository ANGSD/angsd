#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>

gzFile getGz(const char*fname,const char* mode){
  gzFile fp=Z_NULL;
  if(Z_NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}


int main(int argc,char**argv){
  
  if(argc!=5){
    fprintf(stderr,"\nProgram extract ranges from output of msToGlf\n./splitgl file.glf.gz nindTotal firstInd lastInd\n");
    fprintf(stderr,"\nTo extract first the GLS for the first 10 samples in an glf.gz that contains 25 samples\n");
    fprintf(stderr,"Examples\n\n\t1) Extract sample 1 to 12 from a glf.gz file containing 20samples:\n\t\t\t./splitgl raw.glf.gz 20 1 12\n");
    fprintf(stderr,"\t2) Extract sample 13 to 20 from a glf.gz file containing 20samples:\n\t\t\t./splitgl raw.glf.gz 20 13 20\n");
    return 0;
  }
  
  const char *fname = argv[1];
  int tot = atoi( argv[2]);
  int first = atoi( argv[3]);
  int last = atoi( argv[4]);
  first--;last--;

  fprintf(stderr,"fname:%s tot:%d first:%d last:%d\n",fname,tot,first,last);

  gzFile gz=getGz(fname,"r");
  double *lk = malloc(sizeof(double)*10*tot);
  int nsites =0;
  while(1){
    int br=gzread(gz,lk,10*sizeof(double)*tot);
    if(br==0)
      break;
    nsites++;
    assert(br==10*sizeof(double)*tot);
    int beg=10*first;
    int nDoubles=(last-first+1)*10;
    //   fprintf(stderr,"beg: %d nDoubles:%d\n",beg,nDoubles);return 0;
    br=fwrite(lk+beg,sizeof(double),nDoubles,stdout);
    if(br!=nDoubles){
      fprintf(stderr,"Problem writing full chunk\n");
      return 0;
    }
  }
  fprintf(stderr,"nsites:%d processed\n",nsites);
  free(lk);
  gzclose(gz);
}
