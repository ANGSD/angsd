#include <vector>
//input are the dimension of the sfs for the two pops.
//for diploid humands sfs can be 2,4,6 etc. but we will ofcause adjoint with the zero category
#include <cmath>

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};


template <typename T>
void destroy(Matrix<T> *ret,int x){
  for(size_t i=0;i<x;i++)
    delete [] ret->mat[i];
  delete [] ret->mat;
  delete ret;
}

template <typename T>
void destroy(std::vector< Matrix<T> * > &gls,int x){
  for(size_t i=0;i<gls.size();i++)
    destroy(gls[i],x);
}



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
  ret->x=x;
  ret->y=y;
  ret->mat= new T*[x];
  for(size_t i=0;i<ret->x;i++)
    ret->mat[i]=new T[y];
  ret->x=0;
  return ret;
};


void calcCoef(int sfs1,int sfs2,double **aMat,double **bMat);
void block_coef(Matrix<float > *gl1,Matrix<float> *gl2,double *prior,double *a1,double *b1);
