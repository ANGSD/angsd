#include <cstdlib>
#include <vector>
#include <map>
#include <zlib.h>

struct cmp_str
{
  bool operator()(char const *a, char const *b) const
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<char*,std::vector<int>,cmp_str> mm;

void print(std::vector<char*> &v){
  for(int i=0;i<v.size();i++)
    fprintf(stderr,"%d %s\n",i,v[i]);

}
void print(std::vector<int> &v){
  for(int i=0;i<v.size();i++)
    fprintf(stderr,"%d %d\n",i,v[i]);

}


mm* read_pops(char *fname){
  mm *p =new mm;
  gzFile fp = Z_NULL;
  if(((fp=gzopen(fname,"rb")))==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(0);
  }
  char buf[1024];
  std::vector<char*> popnames;
  while(gzgets(fp,buf,1024)){
    popnames.push_back(strdup(strtok(buf,"\n\t\r ")));
  }
  gzclose(fp);

  for(int i=0;i<popnames.size();i++){
    mm::iterator it=p->find(popnames[i]);
    if(it==p->end()){
      std::vector<int> tmp;
      tmp.push_back(i);
      p->insert(std::pair<char *,std::vector<int> >(popnames[i],tmp));
    }else
      it->second.push_back(i);
  }
  
  for(mm::iterator it = p->begin();it!=p->end();++it){
    fprintf(stderr,"%s/%lu\n",it->first,it->second.size());
    print(it->second);
  }
  
  return p;
}

#ifdef __WITH_MAIN__
int main(int a,char **b){
  fprintf(stderr,"a:%d b[1]:%s\n",a,b[1]);
  if(a>0)
    read_pops(b[1]);
  return 0;
}
#endif

