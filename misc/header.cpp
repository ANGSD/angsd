
#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "header.h"

void normalize(double *tmp,size_t len){
  double s=0;
  for(size_t i=0;i<len;i++)
    s += tmp[i];
  for(size_t i=0;i<len;i++)
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


void my_bgzf_write(BGZF *fp, const void *data, size_t length){
  if(bgzf_write(fp,data,length)!=length){
    fprintf(stderr,"\t-> Problem writing bgzf block of size: %lu\n",length);
    exit(0);
  }

}
void my_bgzf_seek(BGZF *fp, int64_t pos, int whence){
  if(bgzf_seek(fp,pos,whence)<0){
    fprintf(stderr,"\t-> Problems seeking in bgzf_seek");
    exit(0);
  }
}
void my_bgzf_read(BGZF *fp, void *data, size_t length){
  if(length!=bgzf_read(fp,data,length)){
    fprintf(stderr,"\t-> Problem reading chunk in bgzf_read\n");
    exit(0);
  }

}


#ifdef __APPLE__
size_t getTotalSystemMemory(){
  uint64_t mem;
  size_t len = sizeof(mem);
  sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
  return mem;
}
#else
size_t getTotalSystemMemory(){
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
#endif


