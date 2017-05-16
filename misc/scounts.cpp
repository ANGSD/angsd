#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <zlib.h>


typedef unsigned char uchar;

typedef struct{
  unsigned int rel_pos;
  uchar A:2;
  uchar C:2;
  uchar G:2;
  uchar T:2;    
}counts;

void set(counts &cnt,int type,int val){
  if(type==0)
    cnt.A=val;
  if(type==1)
    cnt.C=val;
  if(type==2)
    cnt.G=val;
  if(type==3)
    cnt.T=val;
}

void printit(counts &cnt,FILE *fp){
  fprintf(stdout,"%u\t%d\t%d\t%d\t%d\n",cnt.rel_pos,cnt.A,cnt.C,cnt.G,cnt.T);
}
const char *bit_rep[16] = {
    [ 0] = "0000", [ 1] = "0001", [ 2] = "0010", [ 3] = "0011",
    [ 4] = "0100", [ 5] = "0101", [ 6] = "0110", [ 7] = "0111",
    [ 8] = "1000", [ 9] = "1001", [10] = "1010", [11] = "1011",
    [12] = "1100", [13] = "1101", [14] = "1110", [15] = "1111",
};

void print_byte(uint8_t byte)
{
  fprintf(stderr,"%s%s\n", bit_rep[byte >> 4], bit_rep[byte & 0x0F]);
}

FILE *myfopen(const char *file,char *mode){
  FILE *fp = NULL;
  fp = fopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(0);
  }
  return fp;
}
gzFile mygzopen(const char *file,char *mode){
  gzFile fp = Z_NULL;
  fp = gzopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(0);
  }
  return fp;
}


int main(int argc,char **argv){
  
  unsigned nsites = 70e6;
  
  if(argc==1){
    fprintf(stderr,"./a.out -c \n./a.out -d \n");
    return 0;
  }
  char *c=argv[1];
  if(!strcasecmp(c,"-c")){
    char *keep = new char[nsites];
    keep=(char*) memset(keep,0,nsites*sizeof(char));
    FILE *dat = myfopen("data4.txt",(char*)"wb");
    fprintf(stderr,"encrypying\n");
    int nkeep=0;
    for(int i=0;i<nsites;i++){
      counts cnt={0,0,0,0};
      int val =0;
      double r=drand48();
      if(r<=0.06){
	cnt.rel_pos = i;
	//	fprintf(stderr,"hit at:%d\n",i);
	nkeep++;
        int which =lrand48()%4;
	for(int i=0;i<4;i++){
	  int val = 0;
	  if(i==which){
	    
	    val++;
	  }if(drand48()>0.9999)
	    val++;
	  set(cnt,i,val);
	}
	val=1;
	//	print_byte((uint8_t,cnt);
	fwrite(&cnt,sizeof(counts),1,dat);
	keep[i] = 1;
      }
      //printit(cnt,stdout);
      //      fprintf(stderr,"val si:%d\n",val);
    }
    fprintf(stderr,"nkeep:%d\n",nkeep);
    fclose(dat);
  }else if(!strcasecmp(c,"-d")){
    assert(argc ==3);
    gzFile dat = mygzopen(argv[2],(char*)"rb");
    counts cnts;
    while(gzread(dat,&cnts,sizeof(counts))){
      printit(cnts,stdout);

    }
  }
  return 0;
}
