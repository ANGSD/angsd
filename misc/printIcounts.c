#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>


int main(int argc,char **argv){
  int i;
  if(argc==1){
    fprintf(stderr,"Supply icounts.gz file\n");
    return 0;
  }
  int buf[5];
  gzFile gz = Z_NULL;
  if(!((gz=gzopen(argv[1],"r")))){
      fprintf(stderr,"problems opening file:%s\n",argv[1]);
      return 0;
  }
  while(gzread(gz,buf,sizeof(int)*5)){
    for( i=0;i<4;i++)
      fprintf(stdout,"%d\t",buf[i]);
    fprintf(stdout,"%d\n",buf[4]);
  }
  gzclose(gz);
}
