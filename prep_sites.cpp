/*
  part of angsd:

  for standalone program compile like:
  g++ prep_sites.cpp -D__WITH_MAIN__ bgzf.o knetfile.o -lz analysisFunction.o -ggdb
*/

#ifdef __WITH_MAIN__
#define GZOPT "w6h"
#endif


#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <cassert>
#define kv_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#include "prep_sites.h"

#define BIN ".bin"
#define IDX ".idx"

template <typename T>

struct tary{
  size_t l;//contains max position with data
  size_t m;//contains max number of elements
  T *d;
};

template<typename T>
tary<T> *init(){
  tary<T> *r = new tary<T>;
  r->d = new T[1024];
  memset(r->d,0,1024*sizeof(T));
  r->m = 1024;
  r->l=0;
  return r;
}

void dalloc(tary<char> *ta){
  delete [] ta->d;
  delete ta;
}

//hint is the suggested newsize
void filt_readSites(filt*fl,char *chr,size_t hint) {
  assert(fl!=NULL);

  std::map<char*,asdf_dats,ltstr> ::iterator it = fl->offs.find(chr);
  if(it==fl->offs.end()){
    fprintf(stderr,"\n\t-> Potential problem: The filereading has reached a chromsome: \'%s\', which is not included in your \'-sites\' file.\n\t-> Please consider limiting your analysis to the chromsomes of interest \n",chr);
    fprintf(stderr,"\t-> see \'http://www.popgen.dk/angsd/index.php/Sites\' for more information\n");
    fprintf(stderr,"\t-> Program will continue reading this chromosome... \n");
    //exit(0);
    free(fl->keeps);
    free(fl->minor);
    free(fl->major);
    fl->keeps=fl->minor=fl->major=NULL;
    return;
  }

  bgzf_seek(fl->bg,it->second.offs,SEEK_SET);

  size_t nsize = std::max(fl->curLen,hint);
  nsize = std::max(nsize,it->second.len);
  if(nsize>fl->curLen) 
    fl->keeps=(char*) realloc(fl->keeps,nsize);
  memset(fl->keeps,0,nsize);
  //fprintf(stderr,"it->second.len:%lu fl->curLen:%lu fl->keeps:%p\n",it->second.len,fl->curLen,fl->keeps);
  bgzf_read(fl->bg,fl->keeps,it->second.len);

  if(fl->hasMajMin==1){
    if(nsize>fl->curLen) {
      fl->major = (char*) realloc(fl->major,nsize);
      fl->minor = (char*) realloc(fl->minor,nsize);
      memset(fl->major,0,nsize);
      memset(fl->minor,0,nsize);
    }
    bgzf_read(fl->bg,fl->major,it->second.len);
    bgzf_read(fl->bg,fl->minor,it->second.len);
  }
  fl->curNam=chr;
  fl->curLen = nsize;
}

template <typename T>
void set(tary<T> *ta,size_t pos,T val){
  if(pos>ta->m){
    size_t ndim = std::max(ta->m,pos)+1;
    kv_roundup32(ndim);
    T *tmp = new T[ndim];
    memset(tmp,0,sizeof(T)*ndim);
    for(int i=0;i<ta->m;i++)
      tmp[i] = ta->d[i];
    delete [] ta->d;
    ta->d=tmp;
    ta->m = ndim;
  }
  if((int)pos>ta->l){
    ta->l=pos+1;
  }
  ta->d[pos] =val;
}

void filt_print(FILE *fp,filt*f,char *chr){
  fprintf(stderr,"\n---------------\nfp:%p\nbg:%p\nhasMajMin:%d nchr:%lu\n",f->fp,f->bg,f->hasMajMin,f->offs.size());
  for(std::map<char*,asdf_dats,ltstr>::const_iterator it=f->offs.begin();it!=f->offs.end();++it){
    fprintf(stderr,"chr:\t\'%s\'\t",it->first);
    fprintf(stderr,"last:\t\'%lu\'\t",it->second.len);
    fprintf(stderr,"\toffs:%"PRId64"\n",it->second.offs);
  }
  fprintf(stderr,"------------\n");
  for(std::map<char*,asdf_dats,ltstr>::const_iterator it=f->offs.begin();it!=f->offs.end();++it){
    //if we have supplied a single chr
    if(chr!=NULL)
      filt_readSites(f,chr,0);
    else
      filt_readSites(f,it->first,0);
    //    fprintf(stderr,"curlen:%lu\n",f->curLen);
    for(size_t s=0;s<f->curLen;s++)
      if(f->keeps[s]){
	fprintf(fp,"%s\t%lu",f->curNam,s+1);
	if(f->hasMajMin)
	  fprintf(fp,"\t%d\t%d",f->major[s],f->minor[s]);
	
	fprintf(fp,"\n");
      }
    
    if(chr)
      break;
  }
    
  
}


//return one+two
char *append(const char *one,const char*two){
  char *ret = new char[strlen(one)+strlen(two)+1];
  strcpy(ret,one);
  strcpy(ret+strlen(ret),two);
  return ret;
}

void dalloc(filt *f){
  for(pMap::iterator it=f->offs.begin();it!=f->offs.end();++it)
    free(it->first);
  f->offs.clear();
  bgzf_close(f->bg);
  fclose(f->fp);
  free(f->keeps);
  free(f->major);
  free(f->minor);
  delete f;
  f=NULL;

}


filt *filt_read(const char *fname){
  if(!fname){
    fprintf(stderr,"\t-> Must supply filename\n");
    exit(0);
  }
  fprintf(stderr,"\t-> [%s] Reading binary representation of \'%s\'\n",__FILE__,fname);
  filt *ret = new filt;
  ret->bg =NULL;
  ret->fp =NULL;
  ret->keeps = ret->major = ret->minor=NULL;
  ret->hasMajMin =0;
  ret->curLen =0;
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);
  if(!angsd::fexists(bin_name)){
    fprintf(stderr,"\t-> Binary file doesnt exists:%s you should index:%s\n",bin_name,fname);
    fprintf(stderr,"\t-> like: \'angsd sites index %s\'\n",fname);
    exit(0);
  }
  if(!angsd::fexists(idx_name)){
    fprintf(stderr,"Binary file doesnt exists:%s you should index:%s\n",idx_name,fname);
    exit(0);
  }
  if(isNewer(fname,bin_name)||isNewer(fname,idx_name)){ 
    fprintf(stderr,"\t-> Potential problem: File: \'%s\' looks newer than files: \'%s\',\'%s\'\n\t-> Please delete %s/%s files and rerun.\n",fname,bin_name,idx_name,bin_name,idx_name);
    exit(0);
  }

  ret->fp= fopen(idx_name,"r");
  ret->bg=bgzf_open(bin_name,"r");


  while(1){
    int clen;
    if(0==fread(&clen,sizeof(int),1,ret->fp))
      break;
    char chrId[clen];
    if(clen!=(int)fread(&chrId,sizeof(char),clen,ret->fp)){
      fprintf(stderr,"[%s.%s:%d]Problem reading chunk from binary file: clen:%s\n",__FILE__,__FUNCTION__,__LINE__,chrId);
      exit(0);
    }
    asdf_dats tmp;
    if(1!=fread(&tmp.offs,sizeof(int64_t),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    if(1!=fread(&tmp.len,sizeof(size_t),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
#if 0
    fprintf(stderr,"tmp.len:%lu\n",tmp.len); 
#endif
    if(1!=fread(&ret->hasMajMin,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
#if 0
    fprintf(stderr,"tmp.hasMajMin:%d\n",ret->hasMajMin); 
#endif

    ret->offs[strdup(chrId)]=tmp;
  }
  
  fprintf(stderr,"\t-> [%s] nChr: %lu loaded from binary filter file\n",__FILE__,ret->offs.size());
  if(ret->hasMajMin==1)
    fprintf(stderr,"\t-> Filterfile contains major/minor information\n");
  
  delete [] bin_name;
  delete [] idx_name;
  return ret;
}


char *trealloc(char *ptr,size_t i,size_t j){
  fprintf(stderr,"[%s] ptr:%p i:%lu j:%lu\n",__FUNCTION__,ptr,i,j);
  char *r = (char*)calloc(j,1);
  memcpy(r,ptr,i);
  free(ptr);
  return r;
}

/*
  this function will read a file(can be gz), and depending on the number of columsn it will
  1) 2columns: assume chr pos
  2) 3columns: assume chr regionStart regionStop
  3) 4columns: assume chr pos major minor 

 */

typedef std::map<char *,int,ltstr> mmap;

//return zero if fine.
int writeDat(char *last,mmap &mm,tary<char> *keep,tary<char> *major,tary<char> *minor,BGZF *BFP,FILE *fp){
  assert(last!=NULL);
  if((major!=NULL) ^ (minor!=NULL)){
    fprintf(stderr,"major and minor should be the same\n");
    return 1;
  }
  int hasMajMin =0;
  if(major!=NULL)
    hasMajMin =1;
  fprintf(stderr,"\t-> Writing chr:\'%s\' \n",last);
  mmap::iterator it=mm.find(last);
  if(it!=mm.end()){
    return 1;
  }else
    mm[strdup(last)]=1;
  //write data and index stuff
  int64_t retVal =bgzf_tell(BFP);//now contains the offset to which we should point.
  
  //write chrname
  int clen=strlen(last)+1;
  fwrite(&clen,1,sizeof(int),fp);
  fwrite(last,clen,sizeof(char),fp);
  
  fwrite(&retVal,1,sizeof(int64_t),fp);
  fwrite(&keep->l,sizeof(size_t),1,fp);//write len of chr
  fwrite(&hasMajMin,1,sizeof(int),fp);
  bgzf_write(BFP,keep->d,keep->l);//write keep
  if(hasMajMin){
    bgzf_write(BFP,major->d,major->l);//write maj
    bgzf_write(BFP,minor->d,minor->l);//write min
  }
  return 0;
}

void filt_gen(const char *fname) {
  fprintf(stderr,"\t-> Filterfile: %s supplied will generate binary representations... \n",fname);

  gzFile gz = Z_NULL;
  if((gz = gzopen(fname,"r"))==Z_NULL){
    fprintf(stderr,"Problem opening file:%s\n",fname);
    exit(0);
  }

  char* outnames_bin = append(fname,BIN);
  char* outnames_idx = append(fname,IDX);
  
  BGZF *cfpD = bgzf_open(outnames_bin,GZOPT);
  FILE *fp=fopen(outnames_idx,"w");
  
  std::map <char*,int,ltstr> mm;//simple structure to check that input has been sorted by chr/contig
  tary<char> *keep =init<char>();
  tary<char> *major = NULL;
  tary<char> *minor = NULL;
  char *last=NULL;
  int nCols = -1;
  int hasMajMin=0;
  char buf[LENS];

  extern int SIG_COND;
  const char *delims="\t \n";
  char **parsed = new char*[4];
  
  //read a line
  while(SIG_COND && gzgets(gz,buf,LENS)) {
    if(buf[0]=='#'||buf[0]=='\n')//skip if empty or starts with #
      continue;
    
    int nRead=0;
    parsed[nRead++]=strtok(buf,delims);

    char *tok=NULL;
    while(((tok=strtok(NULL,delims)))){
      parsed[nRead++]=tok;
    }
    //it is now tokenized int int parsed
#if 0
    fprintf(stderr,"nRead:%d\n",nRead);
    for(int i=0;i<nRead;i++)
      fprintf(stderr,"c[%d] -> %s\n",i,parsed[i]);
#endif
    
    //check that we have the same number of columsn acroos all lines
    if(nCols==-1){
      nCols=nRead;
      if(nRead==4){
	hasMajMin =1;
	major= init<char>();
	minor= init<char>();
      }
    }else if(nCols!=nRead){
      fprintf(stderr,"\t-> Problem with number of columns in file. Should be 2,3,4 for entire file (%d vs %d)\n",nCols,nRead);
      SIG_COND=0;goto cleanup;
    }
    
   
    //dump
    if(last!=NULL && strcmp(last,parsed[0])){
      //situation:= we have a change of chr, so we DUMP the data
      if(writeDat(last,mm,keep,major,minor,cfpD,fp)){
	SIG_COND=0;
	goto cleanup;
      }
    }
    
    //if we start a 'new' chr
    if(last==NULL||strcmp(last,parsed[0])){
      free(last);
      last=strdup(parsed[0]);
      memset(keep->d,0,keep->m);
      keep->l=0;
      if(hasMajMin){
	memset(major->d,4,keep->m);
	memset(minor->d,4,keep->m);
	major->l=minor->l=0;
      }
	
    }

    size_t posS=atol(parsed[1]);
    assert(posS>0);
    posS--;
    set<char>(keep,posS,1);
    //fprintf(stderr,"keep->l:%lu val:%d\n",keep->l,keep->d[posS]);
    if(nCols==3){
      size_t posE=atol(parsed[2]);
      if(posS>posE){
	fprintf(stderr,"Problem parsing bedfile, end position looks before start position: %lu vs %lu\n",posS,posE);
	exit(0);
      }
      for(size_t ii=posS;ii<posE;ii++)
	set<char>(keep,ii,1);
    }
    if(nCols==4){
      //      fprintf(stderr,"This is the maj/min style\n");
      int al1 = refToInt[parsed[2][0]];
      int al2 = refToInt[parsed[3][0]];
      //      fprintf(stderr,"al1:%d al2:%d\n",al1,al2);
      if(al1==al2||al1==4||al2==4){
	fprintf(stderr,"\t-> major and minor allele should be different, and should be non \'n'/'N' chr:\'%s\' pos:%lu\n",last,posS+1);
	SIG_COND=0;
	goto cleanup;
      }
      set<char>(major,posS,al1);
      set<char>(minor,posS,al2);
    }
  }
  //now after parsing all files
  if(last!=NULL){
    if(writeDat(last,mm,keep,major,minor,cfpD,fp)){
	SIG_COND=0;
	goto cleanup;
    }
    free(last);
    last=NULL;
  }

  fprintf(stderr,"\t-> Filtering complete: Observed: %lu different chromosomes from file:%s\n",mm.size(),fname);

  for(mmap::iterator it=mm.begin();it!=mm.end();++it)
    free(it->first);


  mm.clear();
 cleanup:
  if(gz) gzclose(gz);
  if(fp) fclose(fp);
  if(cfpD) bgzf_close(cfpD);
  
  if(SIG_COND==0){
    fprintf(stderr,"\n\t-> CTRL+C was detected, we will therefore assume that the build of the binary filterfiles are incomplete.\n\t-> Will therefore delete: \'%s\',\'%s\'\n",outnames_bin,outnames_idx);
    unlink(outnames_bin);
    unlink(outnames_idx);
    exit(0);
  }
  delete [] parsed;
  delete [] outnames_idx;delete [] outnames_bin;
  dalloc(keep);
  if(major)
    dalloc(major);
  if(minor)
    dalloc(minor);
}




void filt_init(int argc,char**argv){
  if(argc==0){
    fprintf(stderr,"supply a filename with either 2, 3 or 4 columns\n");
    exit(0);
  }

  char *fname =*argv;
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);
  //  if(!aio::fexists(bin_name)||!aio::fexists(idx_name))
  filt_gen(fname);		//
  fprintf(stderr,"\t-> Generated files:\t\n\t\t'%s\'\n\t\t'%s\'\n",bin_name,idx_name);
  delete [] bin_name;
  delete [] idx_name;
}

int main_sites(int argc,char **argv){
#if 0
  for(int i=0;i<argc;i++)
    fprintf(stderr,"argv[%d]:%s\n",i,argv[i]);
#endif
  if(argc==1){
    fprintf(stderr,"sites index/print filename\n");
    return 0;
  }
  --argc;++argv;
  if(!strcasecmp(*argv,"index")){
    filt_init(--argc,++argv);
  }else if(!strcasecmp(*argv,"print")){
    filt *f = filt_read(*++argv);
    filt_print(stdout,f,NULL);
    dalloc(f);
  }

  return 0;			// 


}


#ifdef __WITH_MAIN__

int main(int argc,char**argv){
  main_sites(--argc,++argv);
}

#endif


