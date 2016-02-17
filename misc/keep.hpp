#pragma once

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

template <typename T>
struct keep{
  size_t m;//length of d array
  size_t l;//number of elements with a set value
  size_t first;//lowest index
  size_t last;//highest index
  T *d;
};

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

template<typename T>
void keep_info(keep<T> *k,FILE *fp,int full,int hit){
  fprintf(fp,"[%s] \t k->m:%lu hit:%d\n",__FUNCTION__,k->m,hit);
  for(size_t i=0;full&&i<k->m;i++)
    fprintf(stderr,"[%s]&k->d[%d]:%p ; k->d+%d:%p ; &k->d[%d+1]:%p k->d+%d+1:%p val:%d\n",__FUNCTION__,i,&k->d[i],i,k->d+i,i,&k->d[i]+1,i,k->d+i+1,k->d[i]);
  size_t tt =0;
  for(int i=0;i<k->m;i++)
    if(k->d[i]==hit)
      tt++;
  fprintf(stderr,"[%s]tt:%d first:%lu last:%lu\n",__FUNCTION__,tt,k->first,k->last);
}


template<typename T>
keep<T> *keep_alloc(){
  keep<T> *ret =(keep<T>*) calloc(1,sizeof(keep<T>));
  return ret;
}

template<typename T>
void realloc(keep<T> *k,size_t newlen){
  //  fprintf(stderr,"[%s] k:%p k->m:%lu newlen:%lu\n",__FUNCTION__,k,k->m,newlen);
  kroundup32(newlen);
  k->d = (T*) realloc(k->d,sizeof(T)*newlen);
  assert(k->d!=NULL);
  memset(k->d+k->m,0,(newlen-k->m)*sizeof(T));
  k->m=newlen;
  //fprintf(stderr,"[%s] k:%p k->m:%lu newlen:%lu\n",__FUNCTION__,k,k->m,newlen);
}

template<typename T>
void keep_set(keep<T> *k,size_t pos,T val){
  //  fprintf(stderr,"[%s] first:%lu last:%lu pos:%lu\n",__FUNCTION__,k->first,k->last,pos);
    if(pos+1>k->m)
      realloc<T>(k,pos+1);
    k->d[pos]=val;
    if(k->l==0){
      k->first=pos;
      k->last=pos;
    }else{
      if(pos<k->first)
	k->first=pos;
      if(pos>k->last)
	k->last=pos;
    }
    k->l++;
    //fprintf(stderr,"[%s] first:%lu last:%lu pos:%lu\n",__FUNCTION__,k->first,k->last,pos);
}

template<typename T>
void keep_destroy(keep<T> *k){
  if(k&&k->d){//<- for some reason valgrind complains about
    free(k->d);
    k->d=NULL;
  }
  free(k);
  k=NULL;
}

template<typename T>
void keep_clear(keep<T> *k){
  memset(k->d,0,k->m*sizeof(T));
  k->m = k->first = k->last = k->l = 0;
}

template class keep<char>;


#ifdef __WITH_MAIN__
typedef char type;
int main(int argc,char **argv){
  keep<type> *k =(keep<type> *) alloc_keep<type>();

  set<type>(k,1,2);
  set<type>(k,2,1);
  set<type>(k,3,1);
  set<type>(k,3,2);
  set<type>(k,4,2);
  set<type>(k,4,1);
  set<type>(k,7,4);
  set<type>(k,8,4);

  set<type>(k,15,10);

  set<type>(k,16,10);
  destroy(k);
}
#endif
