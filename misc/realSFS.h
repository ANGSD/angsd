#pragma once
#include <zlib.h>
#include <cstdio>
#include <zlib.h>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <cmath>
#include <cfloat>
#include <signal.h>
#include <cassert>
#include <pthread.h>

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include "bfgs.h"

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};


int fexists(const char* str);
size_t fsize(const char* fname);
size_t calcNsitesBanded(const char *fname,int nChr);
size_t calcNsites(const char *fname,int nChr);
Matrix<double> alloc(size_t x,size_t y);
void dalloc(Matrix<double> &ret,size_t x);
void normalize(double *tmp,int len);
size_t getTotalSystemMemory();
char *append(const char* a,const char *b);
gzFile getGz(const char *pname);
gzFile getGzBanded(const char *pname);
void readGL(gzFile fp,size_t nSites,int nChr,Matrix<double> &ret);
std::vector<int> getPosi(const char*fname);
void readSFS(const char*fname,int hint,double *ret);
void getArgs(int argc,char **argv);
double lik1(double *sfs,Matrix<double> *ret,int from,int to);
void *lik1_slave(void *p);
double lik1_master();
double lik1_bfgs(const double *sfs_1,const void *dats);
void lik1_grad_bfgs(double *sfs,double *grad,double fac,Matrix<double> *GL1,int from,int to);
void *lik1_grad_bfgs_slave(void *p);
void lik1_grad_bfgs_master(const double *sfs,double *grad);
void bfgs(double *sfs,Matrix<double> *GL1);
void emStep1(double *pre,Matrix<double> *GL1,double *post,int start,int stop);
void *emStep1_slave(void *p);
void emStep1_master(double *post);
void em1(double *sfs,Matrix<double> *GL1,double tole,int maxIter);
double lik2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,size_t start,size_t stop);
void setThreadPars(Matrix<double> *GL1,Matrix<double> *GL2,double *sfs,int nThreads);
void *lik2_slave(void *p);
double lik2_master();
void emStep2(double *pre,Matrix<double> *GL1,Matrix<double> *GL2,double *post,int start,int stop);
void *emStep2_slave(void *p);
void emStep2_master(double *post);
void handler(int s) ;
void em2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,double tole,int maxIter);
Matrix<double>  merge(Matrix<double> &pop1,Matrix<double> &pop2);
void print(Matrix<double> &mat,FILE *fp);
int main_2dsfs(int argc,char **argv);
int main_1dsfs(int argc,char **argv);
int isNewFormat(const char *fname);
int main_1dsfs_v2(const char * fname1, int chr1,long int nSites,int nThreads,char *sfsfname,double tole,int maxIter);
