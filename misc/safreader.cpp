#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <zlib.h>
#include "safreader.h"

void destroy(myMap &mm){
  for(myMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}


void persaf_destroy(persaf *pp){
  //  fprintf(stderr,"pp:%p pp->pos:%p\n",pp,pp->pos);
  if(pp->pos)
    bgzf_close(pp->pos);
  bgzf_close(pp->saf);
  destroy(pp->mm);
  if(pp->ppos){
    delete [] pp->ppos;
  }
  keep_destroy(pp->toKeep);
  free(pp->fname);
  delete pp;
}



void writesaf_header(FILE *fp,persaf *pp){
  fprintf(fp,"\t\tInformation from index file: nChr:%lu nSites:%lu\n",pp->nChr,pp->nSites);
  
  int i=0;
  for(myMap::const_iterator it=pp->mm.begin();it!=pp->mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%ld\n",i++,it->first,d.nSites,(long int)d.pos,(long int)d.saf);
  }

}



/*
  loads first 8bytes and checks magic value

  //return 0 if it is old format
  //returns 1 if safv2
  //returns 2 if safv3

 */

int safversion(const char *fname){
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

  if(0==strcmp(buf,"safv2"))
    return 1;
  else if(0==strcmp(buf,"safv3"))
    return 2;
  else 
    return 0;
}



template <typename T>
persaf * persaf_init(char *fname){
  persaf *ret = new persaf ;
  ret->fname = strdup(fname);
  ret->pos=ret->saf=NULL;ret->toKeep=NULL;
  ret->ppos = NULL;
  ret->kind =0;
  ret->dontRead =0;
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
  ret->version = safversion(fname);
  fprintf(stderr,"\t-> Version of fname:%s is:%d\n",fname,ret->version);
  if(ret->version!=2){
    fprintf(stderr,"\t-> Looks like you are trying to use a version of realSFS that is incompatible with the old binary output from ANGSD\n\t-> Please use realSFS.old instead (or consider redoing the saf files )\n\t-> Will exit\n");
    exit(0);
  }
  if(1!=fread(&ret->nChr,sizeof(size_t),1,fp)){
    fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
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
      fprintf(stderr,"Problem with chr: %s, key already exists, saffile needs to be sorted. (sort your -rf that you used for input)\n",chr);
      exit(0);
    }
  }
  fclose(fp);
  char *tmp =(char*)calloc(strlen(fname)+100,1);//that should do it
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  //  fprintf(stderr,"tmp:%s\n",tmp);
  
  char *tmp2 = (char*)calloc(strlen(fname)+100,1);//that should do it
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  fprintf(stderr,"\t-> Assuming .saf.gz file: %s\n",tmp2);
  ret->saf = bgzf_open(tmp2,"r");
  if(!ret->saf){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n will exit\n",tmp2);
    exit(0);
  }
  if(ret->saf)
    my_bgzf_seek(ret->saf,8,SEEK_SET);
  if(ret->saf && ret->version!=safversion(tmp2)){
    fprintf(stderr,"\t-> Problem with mismatch of version of %s vs %s %d vs %d\n",fname,tmp2,ret->version,safversion(tmp2));
    exit(0);
  }

  snprintf(tmp2,strlen(fname)+100,"%spos.gz",tmp);
  fprintf(stderr,"\t-> Assuming .saf.pos.gz: %s\n",tmp2);
  ret->pos = bgzf_open(tmp2,"r");
  if(!ret->pos){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n will exit\n",tmp2);
    exit(0);
  }
  if(ret->pos)
    my_bgzf_seek(ret->pos,8,SEEK_SET);
  if(ret->pos&& ret->version!=safversion(tmp2)){
    fprintf(stderr,"Problem with mismatch of version of %s vs %s\n",fname,tmp2);
    exit(0);
  }
  //assert(ret->pos!=NULL&&ret->saf!=NULL);
  free(tmp);free(tmp2);
  
 return ret;
 }



 #if 1
 void safprint(int argc,char **argv){

   if(argc<1){
     fprintf(stderr,"Must supply afile.saf.idx [chrname]\n");
     return; 
   }

   char *bname = *argv;
   fprintf(stderr,"\t-> Assuming idxname:%s\n",bname);
   persaf *saf = persaf_init<float>(bname);
   writesaf_header(stderr,saf);

   char *chooseChr = NULL;
   if(argc>0)
     chooseChr = argv[1];
   float *flt = new float[saf->nChr+1];
   for(myMap::iterator it=saf->mm.begin();it!=saf->mm.end();++it){
     if(chooseChr!=NULL){
       it = saf->mm.find(chooseChr);
       if(it==saf->mm.end()){
	 fprintf(stderr,"\t-> Problem finding chr: %s\n",chooseChr);
	 break;
       }
     }
     my_bgzf_seek(saf->pos,it->second.pos,SEEK_SET);
     my_bgzf_seek(saf->saf,it->second.saf,SEEK_SET);
     int *ppos = new int[it->second.nSites];
     my_bgzf_read(saf->pos,ppos,sizeof(int)*it->second.nSites);
     for(size_t s=0;s<it->second.nSites;s++){
       my_bgzf_read(saf->saf,flt,sizeof(float)*(saf->nChr+1));
       fprintf(stdout,"%s\t%d",it->first,ppos[s]);
       for(size_t is=0;is<saf->nChr+1;is++)
	 fprintf(stdout,"\t%f",flt[is]);
       fprintf(stdout,"\n");
     }
     delete [] ppos;
     if(chooseChr!=NULL)
       break;
   }

   delete [] flt;
   persaf_destroy(saf);
 }
 #endif
 //chr start stop is given from commandine
 //if chr==NULL, then this function is only called once
 myMap::iterator iter_init(persaf *pp,char *chr,int start,int stop){
   //  fprintf(stderr,"kind:%d start:%d stop:%d dontread:%d\n",pp->kind,start,stop,pp->dontRead);
   pp->dontRead = 0;
   assert(chr!=NULL);
   myMap::iterator it = pp->mm.find(chr);
   if(it==pp->mm.end()){
     fprintf(stderr,"\t-> [%s] Problem finding chr: %s\n",__FUNCTION__,chr);
     return it;
   }
   my_bgzf_seek(pp->saf,it->second.saf,SEEK_SET);
   
   if(pp->toKeep==NULL)
     pp->toKeep = keep_alloc<char>();
   pp->at =-1;
   if(start==-1&&stop==-1&&pp->kind==0){
     keep_set<char>(pp->toKeep,it->second.nSites,0);
     memset(pp->toKeep->d,1,it->second.nSites);
     pp->toKeep->first = 0;
     pp->toKeep->last = it->second.nSites;
     return it;
   }
   //   fprintf(stderr,"doing pos: kind:%d nsites:%lu\n\n",pp->kind,it->second.nSites);
   my_bgzf_seek(pp->pos,it->second.pos,SEEK_SET);
   if(pp->ppos){
     delete [] pp->ppos;
   }
   pp->ppos = new int[it->second.nSites];
   my_bgzf_read(pp->pos,pp->ppos,sizeof(int)*it->second.nSites);

   keep_set<char>(pp->toKeep,it->second.nSites,0);
   keep_clear(pp->toKeep);
     
   size_t first=0;
   if(start!=-1){
     //fprintf(stderr,"ppos[%d]:%d start:%d pp<start:%d\n",first,pp->ppos[first],start,pp->ppos[first]<start);
     while(first<it->second.nSites&&pp->ppos[first]<start){
       // fprintf(stderr,"ppos[%d]:%d\n",first,pp->ppos[first]);
       first++;
     }
   }
   //   fprintf(stderr,"first:%d\n",first);
   size_t last = it->second.nSites;
   if(stop!=-1&&stop<=pp->ppos[last-1]){
     last=first;
     while(pp->ppos[last]<stop) 
       last++;
   }
   // fprintf(stderr,"last:%d\n",last);
   
   for(size_t s=first;s<last;s++)
     keep_set<char>(pp->toKeep,s,1);
   
   if(pp->kind==0){
     delete [] pp->ppos;
     pp->ppos=NULL;
     // fprintf(stderr,"yoyoyooyoyoyoyoyppos:\n");
   }
   return it;
 }

 size_t iter_read(persaf *saf, void *data, size_t length,int *pos){
   assert(data);
   //   fprintf(stderr,"[%s] kind:%d saf->ppos:%p\n",__FUNCTION__,saf->kind,saf->ppos);//exit(0);
   if(saf->dontRead==1)
     return 0;
  reread: 
   //no more to read
   //   fprintf(stderr,"saf:at:%d lst:%lu\n",saf->at,saf->toKeep->last);
   if(saf->toKeep && saf->at>=(int)saf->toKeep->last)
     return 0;
   
   size_t ret= saf->kind!=1 ? bgzf_read(saf->saf,data,length):length;
   

   if(ret==0)
     return ret;
   saf->at++;
   //   fprintf(stdout,"saf->at:%d ret:%d\n",saf->at,ret);
   if(saf->ppos&&saf->kind>0)
     *pos = saf->ppos[saf->at];
   //   fprintf(stderr,"[%s] posi:%d\n",__FUNCTION__,*pos);

   assert(ret==length);
   
   if(saf->toKeep==NULL)
     return ret;
   
   else if(saf->toKeep->d[saf->at])
     return ret;
   
  goto reread;
}
