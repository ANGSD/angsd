
#include "header.h"

void normalize(double *tmp,int len){
  double s=0;
  for(int i=0;i<len;i++)
    s += tmp[i];
  for(int i=0;i<len;i++)
    tmp[i] /=s;
}



size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}



BGZF *openFileBG(const char* a,const char* b){

  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  BGZF *fp = bgzf_open(c,"wb");
  delete [] c;
  return fp;
}
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

int fexists(const char* str){
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.                             
}

