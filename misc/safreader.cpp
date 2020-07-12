#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <zlib.h>
#include "safreader.h"

void destroy(myMap &mm)
{
  for(myMap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);
  mm.clear();
}

void persaf_destroy(persaf *pp)
{
  //  fprintf(stderr,"pp:%p pp->pos:%p\n",pp,pp->pos);
  if(pp)
  {
    if(pp->pos)
      bgzf_close(pp->pos);
    bgzf_close(pp->saf);
    destroy(pp->mm);
    if(pp->ppos)
      delete [] pp->ppos;
    keep_destroy(pp->toKeep);
    free(pp->fname);
    delete pp;
  }
}

void writesaf_header(FILE *fp, persaf *pp)
{
  if(pp->version==2)
    fprintf(fp,"\t\tInformation from index file: nChr:%lu nSites:%lu\n", pp->nChr, pp->nSites);
  if(pp->version==3)
    fprintf(fp,"\t\tInformation from index file: nChr:%lu nSites:%lu sumBand: %lu\n", pp->nChr, pp->nSites,pp->sumBand);
      
  int i=0;
  for(myMap::const_iterator it=pp->mm.begin();it!=pp->mm.end();++it)
  {
    datum d = it->second;
    if(pp->version==2)
      fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%ld\n",i++,it->first,d.nSites,(long int)d.pos,(long int)d.saf);
    if(pp->version==3)
      fprintf(fp,"\t\t%d\t%s\t%zu\t%lu\t%ld\t%ld\n",i++,it->first,d.nSites,d.sumBand,(long int)d.pos,(long int)d.saf);
  }
}

int safversion(const char *fname)
{
  /*
     Loads first 8bytes and checks magic value
     - return 0 if it is old format
     - return 1 if safv2
     - return 2 if safv3
     - return 3 if safv4
  */

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
  else if(0==strcmp(buf,"safv4"))
    return 3;
  else 
    return 0;
}

template <typename T>
persaf * persaf_init(char *fname, int verbose)
{
  persaf *ret = new persaf;
  ret->fname = strdup(fname);
  ret->pos = ret->saf = NULL;
  ret->toKeep = NULL;
  ret->ppos = NULL;
  ret->kind = 0;
  ret->dontRead = 0;
  ret->sumBand = 0;
  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr, "[persaf::persaf_init] Problem opening file: \'%s\'\n", fname);
    exit(0);
  }

  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(fp==NULL){
    fprintf(stderr, "[persaf::persaf_init] Problem opening file: \'%s\'\n", fname);
    exit(0);
  }

  char buf[8];
  assert(fread(buf,1,8,fp)==8);
  ret->version = safversion(fname);
  
  if(verbose)
    fprintf(stderr, "[persaf::persaf_init] Version of %s is %d\n", fname, ret->version);
  if(!(ret->version==2 || ret->version==3)){
    fprintf(stderr,"[persaf::persaf_init] You are using a version of realSFS that is incompatible with the old binary output from ANGSD. Please recreate the SAF files.\n");
    exit(0);
  }

  if(1!=fread(&ret->nChr,sizeof(size_t),1,fp)){
    fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
    exit(0);
  }

  ret->nSites = 0;
  while(fread(&clen,sizeof(size_t),1,fp))
  {
    char *chr = (char*)malloc(clen+1);
    assert(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';

    datum d;
    d.nSites = d.sumBand=d.pos=d.saf=0;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    ret->nSites += d.nSites;
    if(ret->version==3){
      if(1!=fread(&d.sumBand,sizeof(size_t),1,fp)){
	fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
	exit(0);
      }
      fprintf(stderr,"d.sumband: %d\n",d.sumBand);
      ret->sumBand += d.sumBand;
    }
    
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
      fprintf(stderr,"[persaf::persaf_init] Problem with chromosome %s: key already exists, SAF file needs to be sorted (e.g. sort the file -rf that was used to create the SAF)\n",chr);
      exit(0);
    }
  }

  fclose(fp);
  char *tmp =(char*) calloc(strlen(fname)+100,1);
  tmp=strncpy(tmp,fname,strlen(fname)-3);
  //  fprintf(stderr,"tmp:%s\n",tmp);

  char *tmp2 = (char*)calloc(strlen(fname)+100,1);
  snprintf(tmp2,strlen(fname)+100,"%sgz",tmp);
  if(verbose)
    fprintf(stderr,"[persaf::persaf_init] Assuming .saf.gz file is %s\n",tmp2);
  ret->saf = bgzf_open(tmp2,"r");
  if(!ret->saf){
    fprintf(stderr,"[persaf::persaf_init] Problem opening file \'%s\', will exit.\n",tmp2);
    exit(0);
  }
  if(ret->saf)
    my_bgzf_seek(ret->saf,8,SEEK_SET);
  if(ret->saf && ret->version!=safversion(tmp2)){
    fprintf(stderr,"[persaf::persaf_init] Version mismatch between %s and %s (%d and %d)\n",fname,tmp2,ret->version,safversion(tmp2));
    exit(0);
  }

  snprintf(tmp2,strlen(fname)+100,"%spos.gz",tmp);
  if(verbose)
    fprintf(stderr,"[persaf::persaf_init] Assuming .saf.pos.gz file is %s\n",tmp2);
  ret->pos = bgzf_open(tmp2,"r");
  if(!ret->pos){
    fprintf(stderr,"[persaf::persaf_init] Problem opening file \'%s\', will exit.\n",tmp2);
    exit(0);
  }
  if(ret->pos)
    my_bgzf_seek(ret->pos,8,SEEK_SET);
  if(ret->pos&& ret->version!=safversion(tmp2)){
    fprintf(stderr,"[persaf::persaf_init] Version mismatch between %s and %s (%d and %d)\n",fname,tmp2,ret->version,safversion(tmp2));
    exit(0);
  }

  //assert(ret->pos!=NULL&&ret->saf!=NULL);
  free(tmp);
  free(tmp2);

  return ret;
}

template persaf* persaf_init<float>(char *fname, int verbose);

//chr start stop is given from commandine
//if chr==NULL, then this function is only called once
myMap::iterator iter_init(persaf *pp, char *chr, int start, int stop)
{
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
  }
  return it;
}

template <typename T>
size_t iter_read(persaf *saf, T *&data, T *buffer, int *pos)
{
  assert(buffer);
  
  int band[2];
  
  if(saf->dontRead==1)
    return 0;

  reread: 

  if(saf->toKeep && saf->at >= (int)saf->toKeep->last)
    return 0;

  size_t ret = 0;

  if (saf->version==3)
  {
    ret = bgzf_read(saf->saf, band, 2*sizeof(int));
    if (ret==0)
      return ret;
    assert(ret == 2*sizeof(int));
  }
  else
  {
    band[0] = 0;
    band[1] = saf->nChr+1;
  }
  //not really sure why this was set conditional of saf->kind
  ret = saf->kind!=1 ? bgzf_read(saf->saf, buffer, band[1]*sizeof(T)) : band[1];
  if(ret==0)
    return ret;
  assert(ret==band[1]*sizeof(T) && band[0]+band[1]<=saf->nChr+1);

  saf->at++;
  //   fprintf(stdout,"saf->at:%d ret:%d\n",saf->at,ret);
  if(saf->ppos && saf->kind>0)
    *pos = saf->ppos[saf->at];
  //   fprintf(stderr,"[%s] posi:%d\n",__FUNCTION__,*pos);

  if(saf->toKeep==NULL || saf->toKeep->d[saf->at])
  {
    if(data==NULL ||band[1]>data[1]){
      delete [] data;
      data = new float[band[1]+2];  
    }
    
    data[0] = static_cast<T>(band[0]);
    data[1] = static_cast<T>(band[1]);
    memcpy(data+2, buffer, band[1]*sizeof(T));
    return ret;
  }

  goto reread;
}

template size_t iter_read(persaf *saf, float *&data, float *buffer, int *pos);
