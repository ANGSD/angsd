#pragma once

template <typename T>
struct Matrix
{
  size_t x;
  size_t y;
  T** mat;
  T* buffer;
};

template <typename T>
void destroy(Matrix<T> *ret, size_t x)
{
  if (ret)
  {
    if (ret->mat)
    {
      for(size_t i=0; i<x; i++)
        if (ret->mat[i])
          delete [] ret->mat[i];
      delete [] ret->mat;
    }
    if (ret->buffer)
      delete [] ret->buffer;
    delete ret;
  }
};

template <typename T>
void matrix_print(Matrix<T> *gls)
{
  for(size_t s=0;s<gls->x;s++)
  {
    int j = 0;
    for(size_t i=0; i<gls->mat[s][0]; i++)
      fprintf(stderr,"\t%f",0.);
    for(size_t i=0; i<gls->mat[s][1]; i++)
      fprintf(stderr,"\t%f",gls->mat[s][2+i]);
    for(size_t i=gls->y; i>gls->mat[s][0]+gls->mat[s][1]; i--)
      fprintf(stderr,"\t%f",0.);
    fprintf(stderr,"\n");
  }
}

template <typename T>
Matrix<T> *alloc(size_t x, size_t y){
  Matrix<T> *ret = new Matrix<T>;
  ret->x = x;
  ret->y = y;
  ret->mat = new T*[ret->x];
  for(int i=0; i<ret->x; ++i)
    ret->mat[i] = NULL;
  ret->buffer = new T[ret->y];
  ret->x = 0;
  return ret;
};

template <typename T>
void destroy(std::vector< Matrix<T> * > &gls,size_t x)
{
  for(size_t i=0; i<gls.size(); i++)
    destroy(gls[i], x);
};
