#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "khash.h"
#include "kvec.h"
#include <libgen.h>
#define LENS  4096

const char * BIN= ".counts.bin";
const char * IDX= ".counts.idx";

//does a file exists?
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

typedef struct{
  int nSites;
  int64_t fpos;
}datum;

KHASH_MAP_INIT_STR(s, datum)

void des(khash_t(s) *h){

  for(khiter_t iter=kh_begin(h);iter!=kh_end(h);++iter){
    if(kh_exist(h,iter)){
      free((char*)kh_key(h,iter));
    }
  }
  
  kh_destroy(s,h);
  h=NULL;
}


//return one+two
char *append(char *one,const char*two){
  int len=strlen(one)+strlen(two)+1;

  char *ret =(char*) malloc(strlen(one)+strlen(two)+1);
  strcpy(ret,one);
  strcpy(ret+strlen(ret),two);
  return ret;
}




void writemap(FILE *fp,khash_t(s) *h){
  fprintf(fp,"\t\tInformation from index file:\n");
  int i=0;
  khiter_t iter;// = kh_put(s,h,chr,&ret);
  for(iter=kh_begin(h);iter!=kh_end(h);++iter)
    if (kh_exist(h, iter)){
      datum d = kh_value(h,iter);
      fprintf(fp,"\t%d\t%s\t%d\t%ld\n",i++,kh_key(h,iter),d.nSites,(long int)d.fpos);

    }
  //  exit(0);
}


khash_t(s) *getMap(const char *fname){
  int clen;
  if(!fexists(fname)){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }
  khash_t(s) *h = NULL;
  h =  kh_init(s);
  FILE *fp = fopen(fname,"r");
  while(fread(&clen,sizeof(int),1,fp)) {
  
    char *chr =(char*) malloc(clen+1);
    int hasRead =fread(chr,1,clen,fp); 
  
    if(clen!=hasRead){
      fprintf(stderr,"Problem with reading chr: clen:%d hasread:%d\n",clen,hasRead);
      exit(0);
    }
    chr[clen] = '\0';
    //    fprintf(stderr,"chr is:%s\n",chr);
    datum d;
    if(1!=fread(&d.nSites,sizeof(int),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    //fprintf(stderr,"len:%d\n",d.nSites);
    if(1!=fread(&d.fpos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    
    int ret;
    khiter_t iter = kh_put(s,h,chr,&ret);
    if(ret!=1){
      fprintf(stderr,"\t-> Duplicate key: chr:%s",chr);
      exit(0);
    }
    kh_value(h,iter) =d;
  }
  fclose(fp);
  return h;
}


void stats(int arc,char**argv){
  if(arc==0){
    fprintf(stderr,"./smartCount stats filename.bam\n");
    return;
  }
  char *base = *argv;
  char* outnames_bin = append(base,BIN);
  char* outnames_idx = append(base,IDX);
  fprintf(stderr,"Assuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);

  int clen;
  if(!fexists(outnames_idx)){
    fprintf(stderr,"Problem opening file: %s\n",outnames_idx);
    exit(0);
  }
  khash_t(s) *h = NULL;
  h =  kh_init(s);
  FILE *fp = fopen(outnames_idx,"r");
  while(fread(&clen,sizeof(int),1,fp)) {
  
    char *chr =(char*) malloc(clen+1);
    int hasRead =fread(chr,1,clen,fp); 
  
    if(clen!=hasRead){
      fprintf(stderr,"Problem with reading chr: clen:%d hasread:%d\n",clen,hasRead);
      exit(0);
    }
    chr[clen] = '\0';
    fprintf(stderr,"chrname is \'%s\'\t",chr);
    datum d;
    if(1!=fread(&d.nSites,sizeof(int),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,outnames_idx);
      exit(0);
    }
    fprintf(stderr," length is %d\n",d.nSites);
    if(1!=fread(&d.fpos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,outnames_idx);
      exit(0);
    }
    
    int ret;
    khiter_t iter = kh_put(s,h,chr,&ret);
    if(ret!=1){
      fprintf(stderr,"\t-> Duplicate key: chr:%s",chr);
      exit(0);
    }
    kh_value(h,iter) =d;
  }
  fclose(fp);
  //return h;
}

typedef struct{
  char *chr;
  int l;
  int m;
  unsigned char **counts;
}perChr;

void dalloc(perChr *pc){
  free(pc->chr);
  for(int i=0;i<4;i++)
    free(pc->counts[i]);
  free(pc->counts);
  free(pc);pc=NULL;
}

perChr *alloc(){
  perChr *r =(perChr*) malloc(sizeof(perChr));
  r->l=r->m=0;
  r->counts =(unsigned char**) malloc(4*sizeof(unsigned char*));
  for(int i=0;i<4;i++)
    r->counts[i] = NULL;
  return r;
}


void getPerChr(BGZF *fp,perChr *ret){
  int clen;
  if(bgzf_read(fp,&clen,sizeof(int))==0){
    ret->l =0;
    return;
  }

  ret->chr =(char*) malloc(clen+1);
  bgzf_read(fp,ret->chr,clen);
  ret->chr[clen] = '\0';

  int newsize;
  bgzf_read(fp,&newsize,sizeof(int));
  if(newsize>ret->m){
    for(int i=0;i<4;i++)
      ret->counts[i]=(unsigned char *)realloc(ret->counts[i],newsize);
    ret->m=newsize;
  }
  ret->l=newsize;
  for(int i=0;i<4;i++){
    bgzf_read(fp,ret->counts[i],ret->l*sizeof(char));  
  }
}



void print_main(perChr *pc,FILE *fp) {
  for(int i=0;i<pc->l;i++)
    fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\n",pc->chr,i+1,pc->counts[0][i],pc->counts[1][i],pc->counts[2][i],pc->counts[3][i]);
}

int print(int argc, char**argv){
  if(argc==0){
    fprintf(stderr,"print FILE [-r chrName]\n");
    exit(0);
  }

  char *base = *argv;
  char* outnames_bin = append(base,BIN);
  char* outnames_idx = append(base,IDX);
  fprintf(stderr,"Assuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);
  
  khash_t(s) *h = getMap(outnames_idx);
  BGZF *fp = bgzf_open(outnames_bin,"r");

  --argc;++argv;
  //  fprintf(stderr,"argc=%d\n",argc);
  int argP =0;
  char *chr=NULL;

  while(argP<argc){
    //   fprintf(stderr,"args=%s\n",argv[argP]);
    if(argP==argc){
      fprintf(stderr,"incomplete arguments list\n");
      exit(0);
    }
    if(strcmp("-r",argv[argP])==0)
      chr = argv[argP+1];
    else {
      fprintf(stderr,"Unknown argument:%s\n",argv[argP]);
      exit(0);

    }
    argP +=2;
  }
  
  
  if(chr!=NULL){

    khiter_t iter=kh_get(s,h,chr);
    if(iter==kh_end(h)){
      fprintf(stderr,"Problem finding chr: \'%s\' in index\n",chr);
      exit(0);

    }
    datum d = kh_value(h,iter);

    bgzf_seek(fp,d.fpos,SEEK_SET);
  }

  perChr *pc = alloc();
  while(1){
    getPerChr(fp,pc);

    if(pc->l==0)
      break;
    fprintf(stderr,"pc.chr=%s pc.nSites=%d \n",pc->chr,pc->l);
    print_main(pc,stdout);
    if(chr!=NULL)
      break;

  }
  
  dalloc(pc);
  free(outnames_bin);
  free(outnames_idx);
  bgzf_close(fp);
  des(h);
  return 0;
}
int extract(int argc, char**argv,FILE *fpout) {
  if(argc!=2){
    fprintf(stderr,"./smartCount extract filename.bam sites.txt (sites.txt is chr pos) \n");
    return 0;
  }
  char *base = *argv;
  char* outnames_bin = append(base,BIN);
  char* outnames_idx = append(base,IDX);
  fprintf(stderr,"Assuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);
  
  khash_t(s) *h = getMap(outnames_idx);
  writemap(stderr,h);
  BGZF *fp = bgzf_open(outnames_bin,"r");

  --argc;++argv;
  
  //  fprintf(stderr,"adfasdf:%s\n",*argv);
  
  gzFile gz=Z_NULL;
  gz=gzopen(*argv,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file:%s\n",*argv);
    exit(0);
  }

  perChr *pc = alloc();


  char *buf =(char*) malloc(LENS);
  
  char *last =NULL;
  while(gzgets(gz,buf,LENS)){
    char *chr = strtok(buf,"\n\t ");
    int inflate =0;
    if(last==NULL||strcmp(chr,last)!=0){
      //load chr
      khiter_t iter=kh_get(s,h,chr);
      if(iter==kh_end(h)){
	fprintf(stderr,"Problem finding chr: %s in index\n",chr);
	//	exit(0);
	inflate =1;
      }else{
	datum d = kh_value(h,iter);
	bgzf_seek(fp,d.fpos,SEEK_SET);
	getPerChr(fp,pc);
	free(last);
	last=strdup(chr);//<- small leak, only happens when we change chr
      }      
    }
    int p=atoi(strtok(NULL,"\t\n "))-1;
    if(inflate==0&&p>pc->l){
      fprintf(stderr,"position to extract is after end of chromosome:\n\n");fflush(stderr);
      exit(0);
    }

    fprintf(fpout,"%s\t%d",chr,p+1);
    if(inflate==0)
      for(int i=0;i<4;i++)
	fprintf(fpout,"\t%d",(int)pc->counts[i][p]);
    else
      fprintf(fpout,"\t0\t0\t0\t0");
    fprintf(fpout,"\n");
  }
  free(last);
  dalloc(pc);
  gzclose(gz);
  bgzf_close(fp);
  free(outnames_bin);
  free(outnames_idx);
  free(buf);
  des(h);
  return 0;
}

typedef struct{
  char *fname;
  BGZF *fp;
  khash_t(s) *hash;
}finfo;


char sample(int a,int b,int c,int d){
  double tmp[4] ={a,b,c,d};
  for(int i=1;i<4;i++)
    tmp[i] += tmp[i-1];
  double ts = 0;
  for(int i=0;i<4;i++)
    ts += tmp[i];
  
  for(int i=0;i<4;i++){
    tmp[i] /= ts;
    // fprintf(stderr,"%d=%f\n",i,tmp[i]);
  }
  double dr = drand48();
  for(int i=0;i<3;i++)
    if(dr>=tmp[i]&&dr<tmp[i+1])
      return i;
  return 0;
}

int doAnal(perChr *pc,int nInd){
  assert(nInd>1);
  //  fprintf(stderr,"[%s] chr:%s\n",__FUNCTION__,pc[0].chr);
  //validate same reference lengths
  for(int i=1;i<nInd;i++)
    if(pc[0].l!=pc[i].l){
      fprintf(stderr,"Problem with size of lengths of ref chrs: %d vs %d\n",pc[0].l,pc[i].l);
      exit(0);
    }
  
  //loop over sites
  for(int s=0;s<pc[0].l;s++){
    int hasData = 0; //how many samples has data
    int counts[4] = {0,0,0,0};//counts number of A,C,G,T
    for(int i=0;i<nInd;i++){
      int iDepth =0;//individual depth
      for(int o=0;o<4;o++){
	iDepth += pc[i].counts[o][s];
	counts[o] += pc[i].counts[o][s];
	
      }
      if(iDepth>0)//if individual depth >0 then increment effective size 'hasData'
	hasData++;
    }
    int nobs=0;//check number of different alleles observed
    for(int i=0;i<4;i++)
      if(counts[i])
	nobs++;

    if(hasData==nInd&&nobs==2) {//only do stuff or diallelic sites AND we have data for alle samples

      char res[nInd];//make vector where we sample from 0,1,2,3
      for(int i=0;i<nInd;i++){
	unsigned char **tmp = pc[i].counts;
	res[i] = sample(tmp[0][s],tmp[1][s],tmp[2][s],tmp[3][s]);
      }
      
      //count the number of A,C,G,T
      int counts[4]={0,0,0,0};
      for(int i=0;i<nInd;i++)
	counts[res[i]]++;
      
      int nobs =0;
      int first=-1;
      for(int i=0;i<4;i++){
	if(counts[i])
	  nobs++;
	if(first==-1&& counts[i])
	  first=i;
      }
      
      if(nobs!=2)
	continue;

#if 1
      //      fprintf(stdout,"%s\t%d",pc[0].chr,s+1);
      for(int i=0;i<nInd;i++)
	if(res[i]==first)
	  fprintf(stdout,"\t1,0");
	else
	  fprintf(stdout,"\t0,1");
      fprintf(stdout,"\n");
#endif
    }
  }
  
  return 0;
}



int treemixer(int argc,char **argv){
  if(argc==0){
    fprintf(stderr,"./smartCounts treemixer -bam FILE.list -bin/-idx FILE.list -r chrName\n");
    fprintf(stderr,"FILE.list should be either a list of bam or a list of bins/idxs");
    
  }
  char *flis = NULL;
  char *reg = NULL;
  int type =0;
  while(*argv){
    if(!strcmp("-r",*argv))
      reg=strdup(*++argv);
    else if(!strcmp("-bin",*argv)||!strcmp("-idx",*argv)){
      flis=strdup(*++argv);
    }
    else if(!strcmp("-bam",*argv)){
      flis=strdup(*++argv);
      type++;
    }
    ++argv;
  }
  fprintf(stderr,"reg:%s flis:%s\n",reg,flis );
  if(!flis){
    fprintf(stderr,"Must supply -bin/idx or bam\n");
    return 0;
  }


  kvec_t(finfo) alls;
  kv_init(alls);
  FILE *fp=NULL;
  //  fprintf(stderr,"argc:%d argv[0]=%s\n",argc,argv[0]);
  if(NULL==((fp=fopen(flis,"r")))){
    fprintf(stderr,"Problem opening filelise: %s\n",flis);
    exit(0);
  }
  char *buf = malloc(LENS);
  while(fgets(buf,LENS,fp)){
    char *base = strtok(buf,"\n\t\r ");
    

    finfo tmp;
    tmp.fname = NULL;
    tmp.fp = NULL;
    tmp.hash = NULL;
    char* outnames_bin=NULL;
    char *outnames_idx=NULL;

    if(type==1){
      outnames_bin = append(base,BIN);
      outnames_idx = append(base,IDX);
    }else{
      outnames_bin = strdup(base);
      outnames_idx = strdup(base);
      strcpy(outnames_bin+strlen(outnames_bin)-3,"bin");
      strcpy(outnames_idx+strlen(outnames_idx)-3,"idx");
    }
    // fprintf(stderr,"obin:%s oidX:%s\n",outnames_bin,outnames_idx);
    tmp.fp=bgzf_open(outnames_bin,"r");
    tmp.hash = getMap(outnames_idx);
    tmp.fname = strdup(base);
    fprintf(stdout,"%s\t",basename(tmp.fname));
    kv_push(finfo,alls,tmp);
  }
  fprintf(stdout,"\n");
  fprintf(stderr,"bamlist contains: %zu files\n",kv_size(alls));
#if 0// this prints out the assoc array
  for(uint i =0;i<kv_size(alls);i++){
    finfo tmp = kv_A(alls,i);
    fprintf(stderr,"%u %s %p %p\n",i,tmp.fname,tmp.fp,tmp.hash);
  }
#endif
  
  --argc;++argv;

  char *chr = reg;

  if(chr!=NULL){
    for(uint i =0;i<kv_size(alls);i++) {
      finfo tmp = kv_A(alls,i);
      khiter_t iter=kh_get(s,tmp.hash,chr);
      if(iter==kh_end(tmp.hash)){
	fprintf(stderr,"Problem finding chr: \'%s\' in index\n",chr);
	exit(0);
	
      }
      datum d = kh_value(tmp.hash,iter);
      bgzf_seek(tmp.fp,d.fpos,SEEK_SET);
    }
  }
  perChr *pc = malloc(sizeof(perChr)*kv_size(alls));
  for(uint i =0;i<kv_size(alls);i++){
    pc[i].l=pc[i].m=0;
    pc[i].counts =(unsigned char**) malloc(4*sizeof(unsigned char*));
    for(int ii=0;ii<4;ii++)
      pc[i].counts[ii] = NULL;
  }
  
  while(1){

    //read data for a chr
    for(uint i =0;i<kv_size(alls);i++) {
      finfo tmp = kv_A(alls,i);
      getPerChr(tmp.fp,&pc[i]);
      //    fprintf(stderr,"pc.chr=%s pc.nSites=%d firstpos=%d lastpos=%d\n",pc->chr,pc->nSites,pc->posi[0],pc->posi[pc->nSites-1]);
    }
    if(pc[0].l==0)
      break;
    
    doAnal(pc,kv_size(alls));
    if(chr!=NULL)
      break;
    
    
  }
  return 0;
}


int main(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"./smartCount [stats print extract treemixer]\n");
    return 0;
  }
  //  fprintf(stderr,"zlibversion=%s zlibversion=%s file:%s\n",ZLIB_VERSION,zlib_version,__FILE__);
  --argc,++argv;
  if(!strcmp(*argv,"stats")){
    stats(--argc,++argv);
  }else if(strcmp(*argv,"print")==0){
    print(--argc,++argv);
  }else if(strcmp(*argv,"extract")==0){
    extract(--argc,++argv,stdout);
  }else if(strcmp(*argv,"treemixer")==0){
    //    fprintf(stderr,"treemixer talks\n");
    treemixer(--argc,++argv);
  }else{
    fprintf(stderr,"Unknown argument: \'%s' please supply:\n",*argv);
    fprintf(stderr,"./smartCount [print|extract]\n");
  }
  return 0;
}

