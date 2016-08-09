#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <cstdlib>

unsigned len = 65;


char *my_gets(char *buf,gzFile gz){
  unsigned clen = 0;
  for(;;){
    gzgets(gz,buf+clen,len-clen);
    clen = strnlen(buf,len);
    if(buf[clen-1]=='\n')
      break;
    len *=2;
    buf = (char*) realloc(buf,len);
    clen = strlen(buf);
  }
  return buf;
}

void make_tped(char *buf,gzFile gz){
  for(;;){
    buf=my_gets(buf,gz);
    fprintf(stderr,"srlen:%lu buf:%s %p\n",strlen(buf),buf,gz);
  }


}

int main(int argc,char **argv){
  fprintf(stderr,"\t-> ./haploToPlink input.haplo.gz outputname\n");
  if(argc!=3)
    return 0;
  gzFile gz = Z_NULL;
  if(((gz=(gzopen(argv[1],"rb"))))==Z_NULL){
    fprintf(stderr,"\t-> file: \'%s\' does not exists will exit\n",argv[1]);
    return 0;
  }
  char *buf = (char*)calloc(len,sizeof(char));
  buf = my_gets(buf,gz);
  unsigned ncol=1; //we assume we dont have a null line
  for(unsigned i=0;i<strlen(buf);i++)
    if(buf[i]=='\t')
      ncol++;
  fprintf(stderr,"\t-> We have %d number of fields, which means we have %d samples\n",ncol,ncol-3);

  make_tped(buf,gz);
  
  gzclose(gz);
  return 0;
}
