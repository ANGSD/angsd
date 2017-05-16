#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <zlib.h>
#include <vector>
#include <map>

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

typedef std::map<int,int*> smap;


void counter(smap &sm,char *fname){
  gzFile dat = mygzopen(fname,(char*)"rb");
  counts cnts;
  while(gzread(dat,&cnts,sizeof(counts))){
    smap::iterator it=sm.find(cnts.rel_pos);
    
    if(it==sm.end()){
      int *tmp = new int[4];
      tmp[0]=tmp[1]=tmp[2]=tmp[3]=0;
      sm[cnts.rel_pos] = tmp;
      it = sm.find(cnts.rel_pos);
    }
    int *total = it->second;
    int tmp[4] = {cnts.A,cnts.C,cnts.G,cnts.T};
    int tmp2[4] = {0,0,0,0};
    int at=0;
    for(int i=0;i<4;i++)
      if(tmp[i]>0)
	tmp2[at++] =i;
    //    fprintf(stderr,"informative sites:%d\n",at);
    int which = lrand48() % at;
    it->second[tmp2[which]]++;
  }
}

int main(int argc,char **argv){
  long seed =0;
  unsigned nsites = 70e6;
  
  if(argc==1){
    fprintf(stderr,"./a.out filelist population \n");
    return 0;
  }
  if(0==strcmp(argv[1],"print")){
    fprintf(stderr,"printing\n");
    gzFile dat = mygzopen(argv[2],(char*)"rb");
    counts cnts;
    while(gzread(dat,&cnts,sizeof(counts)))
      printit(cnts,stdout);
    return 0;
    fprintf(stderr,"done printing\n");
  }
  
  if(argc==4){
    fprintf(stderr,"\t-> Setting seed: \n");
    seed = atol(argv[3]);
  }
  srand48(seed);
  char *c=argv[1];
  smap sm;
  gzFile gz = mygzopen(argv[1],(char*)"rb");
  char buf[4096];
  
  while(gzgets(gz,buf,4096)){
    char *file = strtok(buf," \t");
    char *pop = strtok(NULL," \t");
    if(strcmp(pop,argv[2])){
      counter(sm,file);
    }
  }
  for(smap::iterator it=sm.begin();it!=sm.end();it++){
    int *ary = it->second;
    fprintf(stdout,"%d\t%d\t%d\t%d\t%d\n",it->first,ary[0],ary[1],ary[2],ary[3]);
  }
  return 0;
}
