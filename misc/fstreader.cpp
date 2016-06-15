
#include "safreader.h"
#include "fstreader.h"


int fstversion(const char *fname){
  gzFile gz=Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'",fname);
    exit(0);
  }
  char buf[8];
  gzread(gz,buf,8*sizeof(char));
  //  fprintf(stderr,"\t-> Magic nr is: \'%s\'\n",buf);
  gzclose(gz);

  if(0==strcmp(buf,"fstv1"))
    return 0;
  else {
    fprintf(stderr,"\t-> unknown magicnumber: \'%s\'\n",buf);
    exit(0);
  }
}


void destroy(myFstMap &mm){
  for(myFstMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}


void perfst_destroy(perfst *pp){
  bgzf_close(pp->fp);
  destroy(pp->mm);
  for(int i=0;i<pp->names.size();i++)
    free(pp->names[i]);
  delete pp;
}



void writefst_header(FILE *fp,perfst *pp){
  fprintf(fp,"\t-> Information from index file: nSites:%lu nChrs:%lu\n",pp->nSites,pp->mm.size());
  for(int i=0;i<pp->names.size();i++)
    fprintf(fp,"\t-> Population[%d]: %s\n",i,pp->names[i]);
  int i=0;
  for(myFstMap::const_iterator it=pp->mm.begin();it!=pp->mm.end();++it){
    dat d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\n",i++,it->first,d.nSites,(long int)d.off);
  }

}


perfst * perfst_init(char *fname){
  perfst *ret = new perfst ;
  ret->nSites =0;
  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
    exit(0);
  }
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening file:%s\n",fname);
    exit(0);
  }
  char buf[8];
  assert(fread(buf,1,8,fp)==8);
  ret->version=fstversion(fname);
  //read names
  size_t nit=0;
  assert(fread(&nit,sizeof(size_t),1,fp)==1);
  //fprintf(stderr,"nit:%lu\n",nit);
  for(int i=0;i<nit;i++){
    size_t clen;
    assert(fread(&clen,sizeof(size_t),1,fp)==1);
    //fprintf(stderr,"clen:%lu\n",clen);
    char *nam =(char*) calloc(clen+1,1);
    assert(fread(nam,sizeof(char),clen,fp)==clen);
    ret->names.push_back(nam);
  }
#if 1  
  while(fread(&clen,sizeof(size_t),1,fp)){
    char *chr = (char*)calloc(clen+1,1);
    unsigned a =(unsigned) fread(chr,1,clen,fp);
    assert(clen==a);    
    dat d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    //    exit(0);
    ret->nSites += d.nSites;
    if(1!=fread(&d.off,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s->%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    myFstMap::iterator it = ret->mm.find(chr);
    if(it==ret->mm.end())
      ret->mm[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists\n",chr);
      exit(0);
    }  
  }
#endif
  fclose(fp);
  char *tmp =(char*)calloc(strlen(fname)+100,1);//that should do it
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  // fprintf(stderr,"tmp:%s\n",tmp);
  
  char *tmp2 = (char*)calloc(strlen(fname)+100,1);//that should do it
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  fprintf(stderr,"\t-> Assuming .fst.gz file: %s\n",tmp2);
  if(ret->version!=fstversion(tmp2)){
    fprintf(stderr,"\t-> Version mismatch: %d %d\n",ret->version,fstversion(tmp2));
    return NULL;
  }
  ret->fp = bgzf_open(tmp2,"r");
  my_bgzf_seek(ret->fp,8,SEEK_SET);

  free(tmp);
  free(tmp2);
  //  writefst_header(stderr,ret);
  return ret;
 }

