#include <sys/stat.h>

#include <htslib/bgzf.h>
#include <zlib.h>
#include "header.h"
#include "psmcreader.h"



void destroy(myMap &mm){
  for(myMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}


void perpsmc_destroy(perpsmc *pp){
  bgzf_close(pp->bgzf_gls);
  bgzf_close(pp->bgzf_pos);
  destroy(pp->mm);
  
  if(pp->pos)
    delete [] pp->pos;
  if(pp->gls)
    delete [] pp->gls;

  free(pp->fname);
  delete pp;
}



void writepsmc_header(FILE *fp,perpsmc *pp){
  fprintf(fp,"\t\tInformation from index file: nSites:%lu\n",pp->nSites);
  
  int i=0;
  for(myMap::const_iterator it=pp->mm.begin();it!=pp->mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%ld\n",i++,it->first,d.nSites,(long int)d.pos,(long int)d.saf);
  }

}


int psmcversion(const char *fname){
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

  if(0==strcmp(buf,"psmcv1"))
    return 1;
  else 
    return 0;
}



perpsmc * perpsmc_init(char *fname){
  perpsmc *ret = new perpsmc ;
  ret->fname = strdup(fname);
  ret->gls =NULL;
  ret->pos = NULL;
  ret->bgzf_pos=ret->bgzf_gls=NULL;
  ret->pos = NULL;
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
  ret->version = psmcversion(fname);
  fprintf(stderr,"\t-> Version of fname: \'%s\' is:%d\n",fname,ret->version);
  if(ret->version!=1){
    fprintf(stderr,"\t-> Looks like you are trying to use a version of PSMC that does not exists\n");
    exit(0);
  }
  ret->nSites =0;
  while(fread(&clen,sizeof(size_t),1,fp)){
    char *chr = (char*)malloc(clen+1);
    assert(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';
    
    datum d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    ret->nSites += d.nSites;
    if(1!=fread(&d.pos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s->%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    if(1!=fread(&d.saf,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s->%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
  
    myMap::iterator it = ret->mm.find(chr);
    if(it==ret->mm.end())
      ret->mm[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists, psmc file needs to be sorted. (sort your -rf that you used for input)\n",chr);
      exit(0);
    }
  }
  fclose(fp);
  char *tmp =(char*)calloc(strlen(fname)+100,1);//that should do it
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  //  fprintf(stderr,"tmp:%s\n",tmp);
  
  char *tmp2 = (char*)calloc(strlen(fname)+100,1);//that should do it
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  fprintf(stderr,"\t-> Assuming .psmc.gz file: %s\n",tmp2);
  ret->bgzf_gls = bgzf_open(tmp2,"r");
  if(ret->bgzf_gls)
    my_bgzf_seek(ret->bgzf_gls,8,SEEK_SET);
  if(ret->bgzf_gls && ret->version!=psmcversion(tmp2)){
    fprintf(stderr,"\t-> Problem with mismatch of version of %s vs %s %d vs %d\n",fname,tmp2,ret->version,psmcversion(tmp2));
    exit(0);
  }

  snprintf(tmp2,strlen(fname)+100,"%spos.gz",tmp);
  fprintf(stderr,"\t-> Assuming .psmc.pos.gz: %s\n",tmp2);
  ret->bgzf_pos = bgzf_open(tmp2,"r");
  if(ret->pos)
    my_bgzf_seek(ret->bgzf_pos,8,SEEK_SET);
  if(ret->bgzf_pos&& ret->version!=psmcversion(tmp2)){
    fprintf(stderr,"Problem with mismatch of version of %s vs %s\n",fname,tmp2);
    exit(0);
  }
  //assert(ret->pos!=NULL&&ret->saf!=NULL);
  free(tmp);free(tmp2);
  
 return ret;
 }

 //chr start stop is given from commandine
 myMap::iterator iter_init(perpsmc *pp,char *chr,int start,int stop){
   assert(chr!=NULL);
   myMap::iterator it = pp->mm.find(chr);
   if(it==pp->mm.end()){
     fprintf(stderr,"\t-> [%s] Problem finding chr: %s\n",__FUNCTION__,chr);
     return it;
   }
   my_bgzf_seek(pp->bgzf_gls,it->second.saf,SEEK_SET);
   my_bgzf_seek(pp->bgzf_pos,it->second.pos,SEEK_SET);
   //fprintf(stderr,"pp->gls:%p\n",pp->gls);
   if(pp->pos)
     delete [] pp->pos;
   if(pp->gls)
     delete [] pp->gls;
   pp->pos = new int[it->second.nSites];
   my_bgzf_read(pp->bgzf_pos,pp->pos,sizeof(int)*it->second.nSites);
   pp->gls = new double[2*it->second.nSites];//<-valgrind complains about large somthing
   my_bgzf_read(pp->bgzf_gls,pp->gls,2*sizeof(double)*it->second.nSites);
   //   fprintf(stderr," end: %f %f\n",pp->gls[0],pp->gls[1]);
   pp->first=0;
   if(start!=-1)
     while(pp->first<it->second.nSites&&pp->pos[pp->first]<start)
       pp->first++;
   
   pp->last = it->second.nSites;
   if(stop!=-1&&stop<=pp->pos[pp->last-1]){
     pp->last=pp->first;
     while(pp->pos[pp->last]<stop) 
       pp->last++;
   }
   return it;
 }
