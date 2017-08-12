#include "realSFS_shared.h"
#include "Matrix.hpp"


template<typename T>
void readGL(persaf *fp,size_t nSites,size_t dim,Matrix<T> *ret,int *pp, int scale2norm){
  // ret->x=nSites;
  ret->y=dim;
  size_t i;
  for(i=ret->x;SIG_COND&&i<nSites;i++){
    if(i>0 &&(i% howOften)==0  )
      fprintf(stderr,"\r\t-> Has read %fmio sites now at: %lu      ",howOften/1e6,i);
    //
    int pos;
    size_t bytes_read= iter_read(fp,ret->mat[i],sizeof(T)*dim,&pos);//bgzf_read(fp,ret->mat[i],sizeof(T)*dim);
    if(pp!=NULL)//setpos
      pp[i] =pos;//
    if(bytes_read!=0 && bytes_read<sizeof(T)*dim){
      fprintf(stderr,"\t-> Problem reading chunk from file, please check nChr is correct, will exit \n");
      exit(0);
    }
    if(bytes_read==0)
      break;

    for(size_t j=0;scale2norm&&j<dim;j++)
      ret->mat[i][j] = exp(ret->mat[i][j]);
  }
  ret->x=i;
  if(SIG_COND==0)
    exit(0);
}




//returns the number of sites read
template<typename T>
size_t readGLS(std::vector<persaf *> &adolf,size_t nSites,std::vector< Matrix<T> *> &ret,int **posi,int scale2norm){
  size_t pre=ret[0]->x;
  for(int i=0;i<adolf.size();i++){
    readGL(adolf[i],nSites,adolf[i]->nChr+1,ret[i],posi!=NULL?posi[i]:NULL,scale2norm);
    //fprintf(stderr,"adolf:%d\t%lu posi:%d tak:%d\n",i,ret[i]->x,posi[i][0],tak);
  }
     
  return ret[0]->x-pre;
}


template <typename T>
int readdata(std::vector<persaf *> &saf,std::vector<Matrix<T> *> &gls,size_t nSites,char *chooseChr,int start,int stop, int *pp,char **curChr,filt *fl,int scale2norm){
  static size_t lastread=0;
  extern int ** posiG;
  //  fprintf(stderr,"[%s] nSites:%d lastread:%d\n",__FUNCTION__,nSites,lastread);
  if(lastread==0 ){
    fprintf(stderr,"\t-> Done reading data from chromosome will prepare next chromosome\n");
    int ret = set_intersect_pos(saf,chooseChr,start,stop,curChr,fl); 
    //    fprintf(stderr,"[%s] ret:%d\n",__FUNCTION__,ret);
    //    exit(0);
    if(ret==-3)
      return -3;
  }

  lastread=readGLS(saf,nSites,gls,posiG,scale2norm);
  if(lastread>0&&saf.size()>1)
    fprintf(stderr,"\t-> [%s] lastread:%lu posi:%d\n",__FUNCTION__,lastread,posiG[0][0]);
#if 1 //<- below con be removed when we believe all is working
  if(saf.size()>1&&lastread!=0)
    for(int i=1;i<saf.size();i++){
      fprintf(stderr,"\t-> Comparing positions: %d with 0 has:%lu\n",i,gls[0]->x);
      if(memcmp(posiG[0],posiG[i],gls[0]->x*sizeof(int))!=0){
	fprintf(stderr,"SAF file is out of sync contact developer\n");
	for(int s=0;s<gls[0]->x;s++){
	  //	  fprintf(stderr,"s:%d\n",s);
	  if(posiG[0][s]!=posiG[i][s]){
	    fprintf(stderr,"Mismatch at s:%d which is i=0 vs i=%d with pos1:%d pos2:%d\n",s,i,posiG[0][s],posiG[i][s]);
	    exit(0);
	  }
	}
      }
    }
  //  fprintf(stderr,"Done checking\n");
#endif
  if(lastread==0)
    fprintf(stderr,"\t-> Only read nSites: %lu will therefore prepare next chromosome (or exit)\n",gls[0]->x);
  //fprintf(stderr,"readdata lastread:%d\n\n",lastread);
  // exit(0);
  if(pp!=NULL)
    for(int i=0;i<gls[0]->x;i++)
      pp[i] = posiG[0][i];
  if(chooseChr!=NULL&&lastread==0){
    //fprintf(stderr,"return -2\n");
    return -2;
  }
  else if(chooseChr==NULL &&lastread==0 )
    return -2;
  else
    return 1;
}
