#include <cstdio>
#include <zlib.h>
#include <cstring>
#include <cstdlib>

unsigned len = 4096;


char *my_gets(char *buf,gzFile gz){
  buf[0] = '\0';
  unsigned clen = 0;
  for(;;){
    gzgets(gz,buf+clen,len-clen);
    clen = strnlen(buf,len);
    if(clen==0||buf[clen-1]=='\n')
      break;
    len *=2;
    buf = (char*) realloc(buf,len);
    clen = strlen(buf);
  }
  return buf;
}

void make_tped(char *buf,gzFile gz,FILE *of){
  while(((buf=my_gets(buf,gz)))[0]!='\0'){
    char *chr = strtok(buf,"\n\t ");
    char *pos = strtok(NULL,"\n\t");
    strtok(NULL,"\n\t");	
    fprintf(of,"%s\t%s_%s\t0\t%s",chr,chr,pos,pos);
    char *tok=NULL;
    while(((tok=strtok(NULL,"\t\n ")))){
      fprintf(of,"\t%s %s",tok,tok);
    }
    fprintf(of,"\n");
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
  FILE *of = NULL;
  char *onam = (char*)calloc(strlen(argv[2])+6,sizeof(char));
  onam = strcat(onam,argv[2]);
  onam = strcat(onam,".tped");
  fprintf(stderr,"\t-> Will write: \'%s\'\n",onam);
  if((of=fopen(onam,"w"))==NULL){
    fprintf(stderr,"\t-> Problem opening filehandle: %s %p\n",onam,(void*)of);
    return 0;
  }
  make_tped(buf,gz,of);
  fclose(of);
  free(onam);
  onam = (char*)calloc(strlen(argv[2])+6,sizeof(char));
  onam = strcat(onam,argv[2]);
  onam = strcat(onam,".tfam");
  fprintf(stderr,"\t-> Will write: \'%s\'\n",onam);
  if((of=fopen(onam,"w"))==NULL){
    fprintf(stderr,"\t-> Problem opening filehandle: %s %p\n",onam,(void*)of);
    return 0;
  }
  for(unsigned i=0;i<ncol-3;i++)
    fprintf(of,"ind%u\tind%u\t0\t0\t0\t0\n",i,i);

  
  fclose(of);
  of=NULL;
  gzclose(gz);
  return 0;
}
