#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <cassert>
#include <ctype.h>
#define LENS 4096
int minDist = 10;
typedef struct{
  char allele1;
  char allele2;
  double freq;
}hapSite;

std::vector<int> ipos;
std::vector<int*> cnt;

typedef std::map<int,hapSite> aMap;
aMap myMap;

//counts of[dist][allele]
//counts[5] is snpsites,counts[5-1] is position -1 left of snpsites 
size_t counts[9][2];


const char *hapfile=NULL,*mapfile=NULL,*icounts=NULL;
typedef unsigned char uchar;

void count(){
  int lastP = std::max((--myMap.end())->first,ipos[ipos.size()-1])+5;//<-add five so we dont step out
    char *hit = new char[lastP];
  memset(hit,0,lastP);
  int pp=-1;
  for(aMap::iterator it=myMap.begin();it!=myMap.end();++it){
    if(pp!=-1&&(it->first-pp)<=minDist)
      continue;
    pp=it->first;
    for(int p=0;p<5;p++){
      hit[it->first+p] =hit[it->first-p] = 1+p;
    }
  }
  #if 1
  for(int i=0;i<lastP;i++)
    if(hit[i])
      fprintf(stdout,"%d\t%d\n",i,(int)hit[i]);
  #endif
}


char flip(char c){
  c = toupper(c);
  if(c=='A')
    return 'T';
  if(c=='T')
    return 'A';
  if(c=='G')
    return 'C';
  if(c=='C')
    return 'G';
  fprintf(stderr,"Problem interpreting char:%c\n",c);
  return 0;
}

gzFile getgz(const char *fname,const char *mode){
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }
  return gz;
}

void printhapsite(hapSite &hs,FILE *fp,int &p){
  fprintf(fp,"p:%d al1:%c al2:%c freq:%f\n",p,hs.allele1,hs.allele2,hs.freq);
}

void readhap(const char *fname){
  gzFile gz=getgz(fname,"rb");
  char buf[LENS];
  while(gzgets(gz,buf,LENS)){
    hapSite hs;
    int p = atoi(strtok(buf,"\t\n "));
    hs.allele1 = strtok(NULL,"\t\n ")[0];
    hs.freq = atof(strtok(NULL,"\t\n "));
    char strand= strtok(NULL,"\t\n ")[0];
    hs.allele2 = strtok(NULL,"\t\n ")[0];

    //    printhapsite(hs,stdout,p);
    if(strand=='-'){
      hs.allele1 = flip(hs.allele1);
      hs.allele2 = flip(hs.allele2);
    }
    if(myMap.count(p)>0){
      fprintf(stderr,"Duplicate positions found in file: %s, pos:%d\n",fname,p);
      fprintf(stderr,"Will only use first entry\n");
    }else{
      myMap[p]=hs;

    }
  }
  fprintf(stderr,"We have read: %zu sites from hapfile:%s\n",myMap.size(),fname);
}


void readicnts(const char *fname){
  gzFile gz=getgz(fname,"rb");

  int tmp[5];
  while(gzread(gz,tmp,sizeof(int)*5)){
    ipos.push_back(tmp[0]);
    int tmp1[4];
    memcpy(tmp1,tmp+1,sizeof(int)*4);
    cnt.push_back(tmp1);
  }
  
  fprintf(stderr,"Has read: %zu sites from icnts file\n",ipos.size());
}




int main(int argc,char**argv){
  hapfile="../RES/HapMapChrX.gz";
  icounts="../angsdput.icnts.gz";
  mapfile="../RES/chrX.unique.gz";
  readhap(hapfile);
  readicnts(icounts);
  count();
  return 0;
}
