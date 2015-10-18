#pragma once
template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};


template <typename T>
void destroy(Matrix<T> *ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret->mat[i];
  delete [] ret->mat;
  delete ret;
};


template <typename T>
void matrix_print( Matrix<T> *gls){
  for(size_t s=0;s<gls->x;s++){
    for(size_t i=0;i<gls->y;i++)
      fprintf(stderr,"\t%f",gls->mat[s][i]);
  fprintf(stderr,"\n");
  }
}




template <typename T>
Matrix<T> *alloc(size_t x,size_t y){
  Matrix<T> *ret = new Matrix<T>;
  ret->x=x;//+1;
  ret->y=y;
  ret->mat= new T*[ret->x];
  for(size_t i=0;i<ret->x;i++)
    ret->mat[i]=new T[ret->y];
  ret->x=0;
  return ret;
};



template <typename T>
void destroy(std::vector< Matrix<T> * > &gls,size_t x){
  for(size_t i=0;i<gls.size();i++)
    destroy(gls[i],x);
};


