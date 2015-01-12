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

const char *hapfile=NULL,*mapfile=NULL,*icounts=NULL;
typedef unsigned char uchar;

typedef std::vector<int> iv;

typedef struct{
  iv pos;
  iv dist;//0=snpsite;-1=one pos left of snp,+1=one pos right of snp
  std::vector<int *> cn;//four long , counts of C,C,G,T
  aMap myMap;
}dat;

#define NVAL -66
dat count(){
  int lastP = std::max((--myMap.end())->first,ipos[ipos.size()-1])+5;//<-add five so we dont step out
  char *hit = new char[lastP];
  memset(hit,NVAL,lastP);//-10 just indicate no value..
  int pp=-1;
  for(aMap::iterator it=myMap.begin();it!=myMap.end();++it){
    if(pp!=-1&&(it->first-pp)<=minDist)
      continue;
    pp=it->first;
    for(int p=0;p<5;p++){
    
      hit[it->first+p] =hit[it->first-p] = p;
      hit[it->first-p] = -  hit[it->first-p] ;
      //      fprintf(stderr,"p:%d\n",p);
      //fprintf(stderr,"+p %d\n",hit[it->first+p]);
      //fprintf(stderr,"-p %d\n",hit[it->first-p]);
    }
    //    exit(0);
  }
#if 0
  for(int i=0;i<lastP;i++)
    if(hit[i]!=NVAL)
      fprintf(stdout,"%d\t%d\n",i,(int)hit[i]);
  exit(0);
#endif

  dat d;
  for(int i=0;i<ipos.size();i++)
    if(hit[ipos[i]]!=NVAL) { 
      d.pos.push_back(ipos[i]);
      d.dist.push_back(hit[i]);
      d.cn.push_back(cnt[i]);//plugs in pointer to 4ints.
     
      if(hit[ipos[i]]==0){//this is a snpsite
	aMap::iterator it=myMap.find(ipos[i]);
	if(it==myMap.end()){
	  fprintf(stderr,"Problem finding snpsite:%d\n",ipos[i]);
	  exit(0);
	}
	aMap::iterator it2=d.myMap.find(ipos[i]);
	assert(it2==d.myMap.end());
	d.myMap[it->first] = it->second;
	//	fprintf(stdout,"%d\n",it->first);
      }
    }
  fprintf(stderr,"nSNP sites: %lu, with flanking: %lu\n",d.myMap.size(),d.cn.size());
  return d;
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

void readhap(const char *fname,int minDist=10){
  fprintf(stderr,"[%s] fname:%s minDist:%d\n",__FUNCTION__,fname,minDist);
  gzFile gz=getgz(fname,"rb");
  char buf[LENS];
  int viggo=3;
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
      if(viggo>0){
	fprintf(stderr,"Duplicate positions found in file: %s, pos:%d\n",fname,p);
	fprintf(stderr,"Will only use first entry\n");
	fprintf(stderr,"This message is only printed 3 times\n");
	viggo--;
      }
    }else{
      myMap[p]=hs;

    }
  }
  fprintf(stderr,"[%s] We have read: %zu sites from hapfile:%s\n",__FUNCTION__,myMap.size(),fname);
  if(1){
    fprintf(stderr,"[%s] will remove snp sites to close:\n",__FUNCTION__);
    for(aMap::iterator it = myMap.begin();it!=myMap.end();++it){
      aMap::iterator it2=it++;
      if(it2==myMap.end())
	break;
      if(it->first-it2->first<minDist){
	myMap.erase(it);
	it=it2;
      }
    }
    fprintf(stderr,"[%s] We know have: %lu snpSites\n",__FUNCTION__,myMap.size());
  }
  
}


void readicnts(const char *fname,int minDepth=2,int maxDepth=20){
  fprintf(stderr,"[%s] fname:%s minDepth:%d maxDepth:%d\n",__FUNCTION__,fname,minDepth,maxDepth);
  gzFile gz=getgz(fname,"rb");

  int tmp[5];
  int totSite=0;
  while(gzread(gz,tmp,sizeof(int)*5)){
    totSite++;
    int tmp1[4];
    memcpy(tmp1,tmp+1,sizeof(int)*4);
    
    int d=tmp1[0]+tmp1[1]+tmp1[2]+tmp1[3];
    if(d>=minDepth&&d<=maxDepth){
      cnt.push_back(tmp1);
      ipos.push_back(tmp[0]);
    }
  }
  
  fprintf(stderr,"Has read:%d sites,  %zu sites (after depfilter) from ANGSD icnts file\n",totSite,ipos.size());
}




int main(int argc,char**argv){
  hapfile="../RES/hapMapCeuXlift.map.gz";
  icounts="../angsdput.icnts.gz";
  mapfile="../RES/chrX.unique.gz";
  readhap(hapfile);
  return 0;
  readicnts(icounts);
  count();
  return 0;
}
