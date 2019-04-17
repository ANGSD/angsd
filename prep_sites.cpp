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

template<typename T>
void dalloc(tary<T> *ta){
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
    free(fl->an);
    free(fl->ac);
    free(fl->af);
    fl->keeps=fl->minor=fl->major=NULL;
    fl->af=NULL;
    fl->an=fl->ac=NULL;
    fl->curLen =0;
    return;
  }

  assert(0==bgzf_seek(fl->bg,it->second.offs,SEEK_SET));

  size_t nsize = std::max(fl->curLen,hint);//this is not the number of elements, but the last position on the reference
  nsize = std::max(nsize,it->second.len)+1;//not sure if '+1' this is needed, but it doesnt hurt...
  //fprintf(stderr,"nsize:%lu\n",nsize);
  if(nsize>fl->curLen) 
    fl->keeps=(char*) realloc(fl->keeps,nsize);
  memset(fl->keeps,0,nsize);
  //fprintf(stderr,"it->second.len:%lu fl->curLen:%lu fl->keeps:%p\n",it->second.len,fl->curLen,fl->keeps);
  assert(it->second.len==bgzf_read(fl->bg,fl->keeps,it->second.len));

  if(fl->hasExtra>0){
    if(nsize>fl->curLen) {
      fl->major = (char*) realloc(fl->major,nsize);
      fl->minor = (char*) realloc(fl->minor,nsize);
      memset(fl->major,0,nsize);
      memset(fl->minor,0,nsize);
    }
    assert(it->second.len==bgzf_read(fl->bg,fl->major,it->second.len));
    assert(it->second.len==bgzf_read(fl->bg,fl->minor,it->second.len));
  }
  if(fl->hasExtra>1){
    if(nsize>fl->curLen) {
      fl->af = (float*) realloc(fl->af,nsize*sizeof(float));
      fl->an = (int*) realloc(fl->an,nsize*sizeof(int));
      fl->ac = (int*) realloc(fl->ac,nsize*sizeof(int));
      memset(fl->af,0,nsize*sizeof(float));
      memset(fl->an,0,nsize*sizeof(int));
      memset(fl->ac,0,nsize*sizeof(int));
    }
    assert(it->second.len*sizeof(float)==bgzf_read(fl->bg,fl->af,it->second.len*sizeof(float)));
    assert(it->second.len*sizeof(int)==bgzf_read(fl->bg,fl->ac,it->second.len*sizeof(int)));
    assert(it->second.len*sizeof(int)==bgzf_read(fl->bg,fl->an,it->second.len*sizeof(int)));
  }

  fl->curNam=chr;
  fl->curLen = nsize;
}

template <typename T>
void set(tary<T> *ta,size_t pos,T val){
  if(pos>=ta->m){
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
  if((int)pos>=ta->l){
    ta->l=pos+1;
  }
  ta->d[pos] =val;
}

void filt_print(FILE *fp,filt*f,char *chr){
  fprintf(stderr,"\n---------------\nfp:%p\nbg:%p\nhasExtra:%d nchr:%lu\n",f->fp,f->bg,f->hasExtra,f->offs.size());
  for(std::map<char*,asdf_dats,ltstr>::const_iterator it=f->offs.begin();it!=f->offs.end();++it){
    fprintf(stderr,"chr:\t\'%s\'\t",it->first);
    fprintf(stderr,"last:\t\'%lu\'\t",it->second.len);
    fprintf(stderr,"\toffs:%" PRId64 "\n",it->second.offs);
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
	if(f->hasExtra>0)
	  fprintf(fp,"\t%d\t%d",f->major[s],f->minor[s]);
	if(f->hasExtra>1)
	  fprintf(fp,"\t%f\t%d\t%d",f->af[s],f->ac[s],f->an[s]);
	
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

void filt_dalloc(filt *f){
  for(pMap::iterator it=f->offs.begin();it!=f->offs.end();++it)
    free(it->first);
  f->offs.clear();
  bgzf_close(f->bg);
  fclose(f->fp);
  free(f->keeps);
  free(f->major);
  free(f->minor);
  free(f->af);
  free(f->an);
  free(f->ac);
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
  ret->af = NULL;
  ret->ac=ret->an=NULL;
  ret->hasExtra =0;
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
  if(aio::isNewer(fname,bin_name)||aio::isNewer(fname,idx_name)){ 
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
    if(1!=fread(&ret->hasExtra,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
#if 0
    fprintf(stderr,"tmp.hasMajMin:%d\n",ret->hasMajMin); 
#endif

    ret->offs[strdup(chrId)]=tmp;
  }
  
  fprintf(stderr,"\t-> [%s] nChr: %lu loaded from binary filter file\n",__FILE__,ret->offs.size());
  if(ret->hasExtra>0)
    fprintf(stderr,"\t-> Filterfile contains major/minor information\n");
  if(ret->hasExtra>1)
    fprintf(stderr,"\t-> Filterfile contains af ac an tags\n");
  
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
  4) 7columsn: assume chr pos major minor af ac an
 */

typedef std::map<char *,int,ltstr> mmap;

//return zero if fine.
int writeDat(char *last,mmap &mm,tary<char> *keep,tary<char> *major,tary<char> *minor,BGZF *BFP,FILE *fp,int doCompl,tary<float> *af,tary<int> *ac,tary<int> *an){
  assert(last!=NULL);
  if((major!=NULL) ^ (minor!=NULL)){
    fprintf(stderr,"major and minor should be the same\n");
    return 1;
  }
  int hasMajMin =0;
  if(major!=NULL)
    hasMajMin =1;
  if(af!=NULL)
    hasMajMin =2;
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
  for(int i=0;doCompl&&i<keep->l;i++)
    if(keep->d[i]==0)
      keep->d[i]=1;
    else
      keep->d[i]=0;
  

  fwrite(&keep->l,sizeof(size_t),1,fp);//write len of chr
  fwrite(&hasMajMin,1,sizeof(int),fp);
  aio::bgzf_write(BFP,keep->d,keep->l);//write keep
  if(hasMajMin>0){
    aio::bgzf_write(BFP,major->d,major->l);//write maj
    aio::bgzf_write(BFP,minor->d,minor->l);//write min
  }
  if(hasMajMin>1){
    aio::bgzf_write(BFP,af->d,af->l*sizeof(float));//write maj
    aio::bgzf_write(BFP,ac->d,ac->l*sizeof(int));//write min
    aio::bgzf_write(BFP,an->d,an->l*sizeof(int));//write min
  }
  for(int i=0;0&&i<keep->l;i++)//this is check printout loop can be removed
    if(keep->d[i])
      fprintf(stderr,"i:%d af:%f\n",i,af->d[i]);

  return 0;
}

void filt_gen(const char *fname,int posi_off,int doCompl) {
  fprintf(stderr,"\t-> Filterfile: %s supplied will generate binary representations... \n",fname);
  int printInfo = 3;
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
  tary<float> *af = NULL;
  tary<int> *ac = NULL;
  tary<int> *an = NULL;
  char *last=NULL;
  int nCols = -1;
  int hasExtra=0;
  int l = 128;
  char *buf =(char *) calloc(l,sizeof(char));

  extern int SIG_COND;
  const char *delims="\t \n";
  char **parsed = new char*[7];
  
  //read a line
  while(SIG_COND && aio::tgets(gz,&buf,&l)) {
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
    exit(0);
#endif
    
    //check that we have the same number of columsn acroos all lines
    if(nCols==-1){
      nCols=nRead;
      if(nCols==2)
	fprintf(stderr,"\t-> Input file has 2 columns, chr pos\n");
      if(nCols==3)
	fprintf(stderr,"\t-> Input file has 3 columns(bed), chr posStart posStop\n");
      if(nCols==4)
	fprintf(stderr,"\t-> Input file has 4 columns, chr pos major minor\n");
      if(nCols==7)
	fprintf(stderr,"\t-> Input file has 7 columns(vcf)(experimental) \n");
      
      if(nRead==4||nRead==7){
	hasExtra =1;
	major= init<char>();
	minor= init<char>();
      }
      if(nRead==7){
	hasExtra=2;
	af=init<float>();
	ac=init<int>();
	an=init<int>();
	
      }
    }else if(nCols!=nRead){
      fprintf(stderr,"\t-> Problem with number of columns in file. Should be 2,3,4 for entire file (%d vs %d)\n",nCols,nRead);
      SIG_COND=0;goto cleanup;
    }
    
   
    //dump
    if(last!=NULL && strcmp(last,parsed[0])){
      //situation:= we have a change of chr, so we DUMP the data
      if(writeDat(last,mm,keep,major,minor,cfpD,fp,doCompl,af,ac,an)){
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
      if(hasExtra>0){
	memset(major->d,4,keep->m);
	memset(minor->d,4,keep->m);
	major->l=minor->l=0;
      }
      if(hasExtra>1)
	assert(1);
    }
    char *position_in_sites_file = parsed[1];
    assert(atol(position_in_sites_file)>=0);
    size_t posS=atol(parsed[1]);
    assert(posS>0);
    posS--;
    posS += posi_off;

    //fprintf(stderr,"keep->l:%lu val:%d\n",keep->l,keep->d[posS]);
    if(nCols==2)
      set<char>(keep,posS,1);
    else if(nCols==3){
      set<char>(keep,posS,1);
      size_t posE=atol(parsed[2]);
      //posE += posE;
      if(posS>posE){
	fprintf(stderr,"Problem parsing bedfile, end position looks before start position: %lu vs %lu\n",posS,posE);
	exit(0);
      }
      for(size_t ii=posS;ii<posE;ii++)
	set<char>(keep,ii,1);
    }else if(nCols==4||nCols==7){
      //      fprintf(stderr,"This is the maj/min style\n");
      if(strlen(parsed[2])>1||strlen(parsed[3])>1){
	if(printInfo>0){
	  fprintf(stderr,"\t-> major and major are not allowed to be insertions chr:%s posS+1:%d will only print this msg three times\n",last,(int)posS+1);
	  printInfo--;
	}
	continue;
      }
	
      int al1 = refToInt[parsed[2][0]];
      int al2 = refToInt[parsed[3][0]];
      //      fprintf(stderr,"al1:%d al2:%d\n",al1,al2);
      if(al1==al2||al1==4||al2==4){
	fprintf(stderr,"\t-> major and minor allele should be different, and should be non \'n'/'N' chr:\'%s\' pos:%lu\n",last,posS+1);
	SIG_COND=0;
	goto cleanup;
      }
      set<char>(keep,posS,1);
      set<char>(major,posS,al1);
      set<char>(minor,posS,al2);
      if(nCols==7){
	//fprintf(stderr,"This is the vcf af ac an tag style\n");
	
	float freq = atof(parsed[4]);
	int acl = atoi(parsed[5]);
	int anl = atoi(parsed[6]);
	//      fprintf(stderr,"freq:%f acl:%d anl:%d\n",freq,acl,anl);
	set<float>(af,posS,freq);
	set<int>(ac,posS,acl);
	set<int>(an,posS,anl);
      }
    }


  }
  //now after parsing all files
  if(last!=NULL){
    if(writeDat(last,mm,keep,major,minor,cfpD,fp,doCompl,af,ac,an)){
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
  if(af)
    dalloc(af);
  if(ac)
    dalloc(ac);
  if(an)
    dalloc(an);
  free(buf);
}




void filt_init(int argc,char**argv){
  if(argc==0){
    fprintf(stderr,"supply a filename with either 2, 3 or 4 columns\n");
    exit(0);
  }

  char *fname = *argv;
  int posi_offs = 0;
  int doCompl =0;
  argv++;argc--;
  while(*argv){
    if(strcmp(*argv,"-o")==0)
      posi_offs =  atoi(*(++argv));
    else if(strcmp(*argv,"-compl")==0){
      doCompl =  atoi(*(++argv));
    }else{
      fprintf(stderr,"\t-> Problem parsing option: \'%s\' will exit. -o -compl is allowed\n",*argv);
      exit(0);
    }
    ++argv;
  }
  fprintf(stderr,"\t-> Indexing %s and will add \'%d\' to pos column\n",fname,posi_offs);
  if(doCompl)
    fprintf(stderr,"\t-> Your are using -compl , consider adding a dummy last position\n");
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);
  //  if(!aio::fexists(bin_name)||!aio::fexists(idx_name))
  filt_gen(fname,posi_offs,doCompl);		//
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
    fprintf(stderr,"\tsites print filename\t\tPrint index file\n\tsites index filename [-r offset -compl doCompl]\tgenerate binary index file\n");
    return 0;
  }
  --argc;++argv;
  if(!strcasecmp(*argv,"index")){
    filt_init(--argc,++argv);
  }else if(!strcasecmp(*argv,"print")){
    filt *f = filt_read(*++argv);
    filt_print(stdout,f,NULL);
    filt_dalloc(f);
  }else
    fprintf(stderr,"Unknown option: \'%s\'\n",*argv);

  return 0;
}


#ifdef __WITH_MAIN__

int main(int argc,char**argv){
  main_sites(--argc,++argv);
}

#endif


