
/*
  loads first 8bytes and checks magic value
  //returns 1 if it is new format
  //return 0 if it is old format
 */

#include <zlib.h>


int isNewFormat(const char *fname){
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file: \'%s\'",fname);
    exit(0);
  }
  char buf[8];
  gzread(gz,buf,8);
  return (0==strcmp(buf,"safv2"));
}

/*
  maybe not this format, just a suggestion

 */
typedef struct{
  int start;
  int len;
  double *llh;
}aSite;



int main_1dsfs_v2(const char * fname, int nchr,int nsites,int nthreads,char *startsfs,double tole,int maxiter){
  



}
