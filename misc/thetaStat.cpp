#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <zlib.h>
#include <sys/stat.h>
#include "../htslib/bgzf.h"
#include "kstring.h"
#include "stats.cpp"
#include <cassert>

#define VERSION "0.01"


const char * BIN= ".bin";
const char * IDX= ".idx";
const char * RES=  ".pestPG";
//does a file exists?
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}


struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef struct{
  size_t nSites;
  int64_t fpos;
}datum;

typedef std::map<char*,datum,ltstr> myMap;



//return one+two
char *append(char *one,const char*two){
  char *ret = new char[strlen(one)+strlen(two)+1];
  strcpy(ret,one);
  strcpy(ret+strlen(ret),two);
  return ret;
}


#define LENS  4096

typedef struct{
  int posi;
  float *vals;
}the_t;


int64_t writeAll(std::vector<the_t> &thetas, char *chr,BGZF *fp){
  fprintf(stderr,"\tWriting: chr:%s with nSites:%zu\n",chr,thetas.size());
  int64_t retVal =bgzf_tell(fp); 
  size_t clen=strlen(chr);
  bgzf_write(fp,&clen,sizeof(size_t));//write len of chr
  bgzf_write(fp,chr,clen);//write chr
  size_t vLen = thetas.size();
  bgzf_write(fp,&vLen,sizeof(size_t));//write len of positions;
  int *posi = new int[thetas.size()];
  static float **the = new float*[5];
  for(int i=0;i<5;i++)
    the[i] = new float[thetas.size()];
  for(size_t i=0;i<thetas.size();i++){
    posi[i] =thetas[i].posi;
    for(int j=0;j<5;j++)
      the[j][i] = thetas[i].vals[j];
    delete [] thetas[i].vals;
  }
  bgzf_write(fp,posi,sizeof(int)*thetas.size());
  for(int j=0;j<5;j++){
    bgzf_write(fp,the[j],sizeof(float)*thetas.size());
    delete [] the[j];
  }
  
  delete [] posi;
  fprintf(stderr,"\tDone writing: %s\n",chr);
  return retVal;
}


void write_index(size_t nSites,char*chr,FILE *fp,int64_t fpos){
  //write chr
  size_t clen=strlen(chr);
  fwrite(&clen,sizeof(size_t),1,fp);
  fwrite(chr,1,clen,fp);
  //write nSites;
  fwrite(&nSites,sizeof(size_t),1,fp);
  //write fpos
  fwrite(&fpos,sizeof(int64_t),1,fp);
  fflush(fp);
}

void make_bed(int argc,char **argv){
  //  fprintf(stderr,"[%s] \n",__FUNCTION__);
  if(argc==0){
    fprintf(stderr,"make_bed FILE.theta.gz [OUTNAMES] (if OUTNAMES is supplied, this will be used as prefix \n");
    exit(0);
  }
  
  if(!fexists(argv[0])){
    fprintf(stderr,"Problem opening file: %s\n",argv[0]);
    exit(0);
  }
  char *base = argv[0];
  if(argc==2)
    base = argv[1];

  char* outnames_bin = append(base,BIN);
  char* outnames_idx = append(base,IDX);
    
  const char *delims = "\t \n";
  gzFile gfp = gzopen(argv[0],"r");
  char *buf = new char[LENS];
  BGZF *cfpD = bgzf_open(outnames_bin,"w9");
  FILE *fp =fopen(outnames_idx,"w");
  
  std::vector<the_t> vec;
  char *lastChr = NULL;
  
  
  while(gzgets(gfp,buf,LENS)){
    char *chr = strtok(buf,delims);
    if(chr[0]=='#')
      continue;
    int posi=atoi(strtok(NULL,delims)) ;
    
    if(lastChr==NULL){
      lastChr = strdup(chr);
    }else if(strcmp(lastChr,chr)!=0){
      int64_t id=writeAll(vec,lastChr,cfpD);//write data
      write_index(vec.size(),lastChr,fp,id);//write index;
      
      vec.clear();
      free(lastChr);
      lastChr=strdup(chr);
    }
    the_t t;
    t.posi =posi;
    float *the =new float[5];
    for(int i=0;i<5;i++)
      the[i] = atof(strtok(NULL,delims)) ;
    t.vals = the;
    vec.push_back(t);
#if 0
    fprintf(stderr,"%s %d ",chr,posi);
    for(int i=0;i<5;i++)
      fprintf(stderr," %f",the[i]);
    fprintf(stderr,"\n");
#endif
  }
  int64_t id=writeAll(vec,lastChr,cfpD);//write data
  write_index(vec.size(),lastChr,fp,id);//write index;
  vec.clear();
  free(lastChr);
  
  fprintf(stderr,"\tHas dumped files:\n\t\t'%s\'\n\t\t\'%s\'\n",outnames_bin,outnames_idx);
  bgzf_close(cfpD);
  fclose(fp);
  
  gzclose(gfp);
  delete [] buf;
  delete [] outnames_bin; delete [] outnames_idx;
}

void writemap(FILE *fp,const myMap &mm){
  fprintf(fp,"\t\tInformation from index file:\n");
  int i=0;
  for(myMap::const_iterator it=mm.begin();it!=mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\n",i++,it->first,d.nSites,(long int)d.fpos);

  }

}


myMap getMap(const char *fname){
  myMap ret;
  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(0);
  }
  FILE *fp = fopen(fname,"r");
  while(fread(&clen,sizeof(size_t),1,fp)){
    char *chr = new char[clen+1];
    assert(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';
    
    datum d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    if(1!=fread(&d.fpos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(0);
    }
    myMap::iterator it = ret.find(chr);
    if(it==ret.end())
      ret[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists\n",chr);
      exit(0);
    }
  }

  return ret;
}

typedef struct{
  char *chr;
  size_t nSites;
  int *posi;
  float *tW;
  float *tP;
  float *tF;
  float *tH;
  float *tL;
}perChr;

void dalloc(perChr &pc){
  delete [] pc.chr;
  delete [] pc.posi;
  delete [] pc.tW;
  delete [] pc.tP;
  delete [] pc.tF;
  delete [] pc.tH;
  delete [] pc.tL;

}



perChr getPerChr(BGZF *fp){
  perChr ret;
  ret.nSites =0;
  ret.posi=NULL;
  ret.tW=ret.tP=ret.tF=ret.tH=ret.tL=NULL;
  size_t clen;
  
  if(bgzf_read(fp,&clen,sizeof(size_t))==0)
    return ret;
  ret.chr = new char[clen+1];
  bgzf_read(fp,ret.chr,clen);
  ret.chr[clen] = '\0';
  bgzf_read(fp,&ret.nSites,sizeof(size_t));
  ret.posi = new int[ret.nSites];
  ret.tW = new float[ret.nSites];
  ret.tP = new float[ret.nSites];
  ret.tF = new float[ret.nSites];
  ret.tH = new float[ret.nSites];
  ret.tL = new float[ret.nSites];
  
  //read positions and thetas
  bgzf_read(fp,ret.posi,ret.nSites*sizeof(int));
  bgzf_read(fp,ret.tW,ret.nSites*sizeof(float));
  bgzf_read(fp,ret.tP,ret.nSites*sizeof(float));
  bgzf_read(fp,ret.tF,ret.nSites*sizeof(float));
  bgzf_read(fp,ret.tH,ret.nSites*sizeof(float));
  bgzf_read(fp,ret.tL,ret.nSites*sizeof(float));
  
  //make thetas into normal space
  for(size_t i=0;i<ret.nSites;i++){
    ret.tW[i] = exp(ret.tW[i]);
    ret.tP[i] = exp(ret.tP[i]);
    ret.tF[i] = exp(ret.tF[i]);
    ret.tH[i] = exp(ret.tH[i]);
    ret.tL[i] = exp(ret.tL[i]);
  }


  return ret;
}

double slice(int begI, int endI,float *ary){
  double ret=0;
  for(int i=begI;i<endI;i++)
    ret += ary[i];
  return ret;
}


void calc_stat(int begI,int endI,perChr &pc,kstring_t &str,int nChr){
  double tW = slice(begI,endI,pc.tW); 
  double tP = slice(begI,endI,pc.tP); 
  double tF = slice(begI,endI,pc.tF); 
  double tH = slice(begI,endI,pc.tH); 
  double tL = slice(begI,endI,pc.tL); 
  
  ksprintf(&str,"%f\t%f\t%f\t%f\t%f\t",tW,tP,tF,tH,tL);
  
  double tajima=tajd(nChr,tW,tP);
  double fuf = fulif(nChr,tW,tF,tP);
  double fud = fulid(nChr,tW,tF);
  double fayh_val = fayh(nChr,tW,tH,tP);
  double zeng_val = zenge(nChr,tW,tL);
  //  fprintf(stderr,"%f\n",tajima);
  ksprintf(&str,"%f\t%f\t%f\t%f\t%f",tajima,fuf,fud,fayh_val,zeng_val);
}


//type=0 regular=old, type=1 start isfirst pos,type=2 is from pos0.
kstring_t do_stat_main(perChr &pc,int step,int win,int nChr,int type){
  int pS,pE;//physical start,physical end
  int begI,endI;//position in array for pS, pE;
  
  kstring_t str;
  str.l=str.m=0;
  str.s=NULL;

  if(step==0&&win==0){
    //assuming whole chromosome as window
    begI=0;
    endI=pc.nSites-1;
    pS=0;
    pE=pc.posi[endI];
    ksprintf(&str,"(%d,%d)(%d,%d)(%d,%d)\t%s\t%d\t",begI,endI,pc.posi[begI],pc.posi[endI],pS,pE,pc.chr,pS+(pE-pS)/2);
    calc_stat(begI,endI,pc,str,nChr);
    ksprintf(&str,"\t%d\n",endI-begI);
    return str;

  }


  if(type==0)
    pS = (pc.posi[0]/step)*step +step;
  else if(type==1)
    pS = pc.posi[0];
  else if(type==2)
    pS = 1;
  pE = pS+win;
  begI=endI=0;
  
  
  if(pE>pc.posi[pc.nSites-1]){
    fprintf(stderr,"end of dataset is before end of window: end of window:%d last position in chr:%d\n",pE,pc.posi[pc.nSites-1]);
    return str;
  }

  while(pc.posi[begI]<pS) begI++;
  
  endI=begI;
  while(pc.posi[endI]<pE) endI++;

  while(1){
    //    fprintf(estpgF,"(%d,%d)(%d,%d)\t%s\t%d\t",begI,endI,the[begI].pos,the[endI].pos,chrCur,the[begI].pos+(the[endI].pos-the[begI].pos)/2);
    ksprintf(&str,"(%d,%d)(%d,%d)(%d,%d)\t%s\t%d\t",begI,endI,pc.posi[begI],pc.posi[endI],pS,pE,pc.chr,pS+(pE-pS)/2);
    calc_stat(begI,endI,pc,str,nChr);
    ksprintf(&str,"\t%d\n",endI-begI);

    //str.l=0;
    pS += step;
    pE =pS+win;
    if(pE>pc.posi[pc.nSites-1])
      break;

    while(pc.posi[begI]<pS) begI++;
    while(pc.posi[endI]<pE) endI++;
    
  }
  
  return str;

}





int do_stat(int argc, char**argv){
  if(argc==0){
    fprintf(stderr,"do_stat FILE -win -step -nChr [-r chrName -type [0,1,2]]\n");
    exit(0);
  }
  char *base = *argv;
  char* outnames_bin = append(base,BIN);
  char* outnames_idx = append(base,IDX);
  fprintf(stderr,"\tAssuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);
  
  myMap mm = getMap(outnames_idx);
  writemap(stderr,mm);
  BGZF *fp = bgzf_open(outnames_bin,"r");

  --argc;++argv;
  //  fprintf(stderr,"argc=%d\n",argc);
  int argP =0;
  char *chr=NULL;
  char *outnames = NULL;
  int nChr =0;
  int win =0;
  int step =0;
  int type =0;
  while(argP<argc){
    //   fprintf(stderr,"args=%s\n",argv[argP]);
    if(argP==argc){
      fprintf(stderr,"incomplete arguments list\n");
      exit(0);
    }
    if(strcmp("-r",argv[argP])==0)
      chr = argv[argP+1];
    else if(strcmp("-outnames",argv[argP])==0)
      outnames = argv[argP+1];
    else if(strcmp("-step",argv[argP])==0)
      step = atoi(argv[argP+1]);
    else if(strcmp("-win",argv[argP])==0)
      win = atoi(argv[argP+1]);
    else if(strcmp("-nChr",argv[argP])==0)
      nChr = atoi(argv[argP+1]);
    else if(strcmp("-type",argv[argP])==0)
      type = atoi(argv[argP+1]);
    
    else {
      fprintf(stderr,"Unknown argument:%s\n",argv[argP]);
      exit(0);
    }
    argP +=2;
  }

  fprintf(stderr,"\t -r=%s outnames=%s step: %d win: %d nChr:%d\n",chr,outnames,step,win,nChr);
  if(nChr==0){
    fprintf(stderr,"nChr must be different from zero\n");
    exit(0);
  }
  if(win==0||step==0){
    fprintf(stderr,"\tWinsize equals zero or step size equals zero. Will use entire chromosome as window\n");
    win=step=0;
  }  
  
  if(chr!=NULL){  
    myMap::iterator it = mm.find(chr);
    if(it==mm.end()){
      fprintf(stderr,"\tProblem finding chr: %s in index\n",chr);
      exit(0);
    }
    datum d = it->second;
    bgzf_seek(fp,d.fpos,SEEK_SET);
  }
  if(outnames==NULL)
    outnames = base;

  char *resname = append(outnames,RES);
  FILE *fpres = fopen(resname,"w");
  fprintf(fpres,"## thetaStat VERSION: %s build:(%s,%s)\n",VERSION,__DATE__,__TIME__);
  fprintf(fpres,"#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)\t");
  fprintf(fpres,"Chr\tWinCenter\t");
  fprintf(fpres,"tW\ttP\ttF\ttH\ttL\t");
  fprintf(fpres,"Tajima\tfuf\tfud\tfayh\tzeng\tnSites\n");
  while(1){
    perChr pc = getPerChr(fp);
    if(pc.nSites==0)
      break;
    fprintf(stderr,"\tpc.chr=%s pc.nSites=%zu firstpos=%d lastpos=%d\n",pc.chr,pc.nSites,pc.posi[0],pc.posi[pc.nSites-1]);
    kstring_t str = do_stat_main(pc,step,win,nChr,type);
    fwrite(str.s,1,str.l,fpres);//should clean up str, doesn't matter for this program;
    fflush(fpres);
    if(chr!=NULL)
      break;
    dalloc(pc);
  }
  fclose(fpres);
  fprintf(stderr,"\tDumping file: \"%s\"\n",resname);
  return 0;
}

void print_main(perChr &pc,FILE *fp){
  for(size_t i=0;i<pc.nSites;i++)
    fprintf(fp,"%s\t%d\t%f\t%f\t%f\t%f\t%f\n",pc.chr,pc.posi[i],log(pc.tW[i]),log(pc.tP[i]),log(pc.tF[i]),log(pc.tH[i]),log(pc.tL[i]));
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
  
  myMap mm = getMap(outnames_idx);
  writemap(stderr,mm);
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
    myMap::iterator it = mm.find(chr);
    if(it==mm.end()){
      fprintf(stderr,"Problem finding chr: %s in index\n",chr);
      exit(0);
    }
    datum d = it->second;
    bgzf_seek(fp,d.fpos,SEEK_SET);
  }

  while(1){
    perChr pc = getPerChr(fp);
    if(pc.nSites==0)
      break;
    fprintf(stderr,"pc.chr=%s pc.nSites=%zu firstpos=%d lastpos=%d\n",pc.chr,pc.nSites,pc.posi[0],pc.posi[pc.nSites-1]);
    print_main(pc,stdout);
    if(chr!=NULL)
      break;
    dalloc(pc);
  }

  return 0;
}
void fun(float a,float b){
  if(fabs(exp(a)-b)>1e-6){
    fprintf(stderr,"\tDifference between binary dump and raw text is rather big:%f vs %f\n",a,b);
    exit(0);
  }

}

int val_bed(int argc, char**argv){
  if(argc!=1){
    fprintf(stderr,"val_bed FILE.gz \n");
    exit(0);
  }
  char *base = *argv;
  char* outnames_bin = append(base,BIN);
  char* outnames_gz = base;
  fprintf(stderr,"Assuming binfile:%s and gzfile:%s\n",outnames_bin,outnames_gz);
  
  
  BGZF *fp = bgzf_open(outnames_bin,"r");
  gzFile gz =gzopen(outnames_gz,"r");
  char buf[4096];
  gzgets(gz,buf,4096);
  while(1){
    perChr pc = getPerChr(fp);
    if(pc.nSites==0)
      break;
    fprintf(stderr,"pc.chr=%s pc.nSites=%zu firstpos=%d lastpos=%d\n",pc.chr,pc.nSites,pc.posi[0],pc.posi[pc.nSites-1]);
    for(size_t i=0;i<pc.nSites;i++){
      gzgets(gz,buf,4096);
      char *chr = strtok(buf,"\n\t ");
      if(strcmp(chr,pc.chr)!=0){
	fprintf(stderr,"Problem with nonmatching chromosome: \'%s\' vs \'%s\'\n",chr,pc.chr);
	exit(0);
      }
      int posi =atoi(strtok(NULL,"\t\n "));
      if(posi!=pc.posi[i]){
	fprintf(stderr,"Problem with nonmatching position\n");
	exit(0);
      }
      float tW = atof(strtok(NULL,"\t\n "));
      float tP = atof(strtok(NULL,"\t\n "));
      float tF = atof(strtok(NULL,"\t\n "));
      float tH = atof(strtok(NULL,"\t\n "));
      float tL = atof(strtok(NULL,"\t\n "));
      fun(tW,pc.tW[i]);
      fun(tP,pc.tP[i]);
      fun(tF,pc.tF[i]);
      fun(tH,pc.tH[i]);
      fun(tL,pc.tL[i]);
    }
    fprintf(stderr,"FILE: %s chr: %s OK\n",base,pc.chr); 
    dalloc(pc);
  }

  fprintf(stderr,"ALL OK: %s\n",base);
  return 0;
}



int main(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"\t\'thetaStat\', a program to do neutrality test statistics using thetas.gz output from angsd\n");
    fprintf(stderr,"\tSYNOPSIS:\n\t\t./thetaStat [make_bed|do_stat|validate_bed|print]\n");
    fprintf(stderr,"\tEXAMPLE:\n\t\t1) \'./thetaStat make_bed N00200.thetas.gz\'\n");
    fprintf(stderr,"\t\t2) \'./thetaStat do_stat N00200.thetas.gz -nChr 16 -win 40000 -step 10000\'\n");
    fprintf(stderr,"\tOUTPUT IS THEN CALLED  \'N00200.thetas.gz.pestPG\n");
    fprintf(stderr,"\n\tYOU CAN TRY WITH DIFFERENT WINDOWSIZE LIKE:\n\t \'./thetaStat do_stat N00200.thetas.gz -nChr 16 -win 20000 -step 10000\' \n");
    return 0;
  }
  //  fprintf(stderr,"zlibversion=%s zlibversion=%s file:%s\n",ZLIB_VERSION,zlib_version,__FILE__);
  --argc,++argv;
  if(strcmp(*argv,"make_bed")==0){
    make_bed(--argc,++argv);
  }else if(strcmp(*argv,"do_stat")==0){
    do_stat(--argc,++argv);
  }else if(strcmp(*argv,"validate_bed")==0){
    val_bed(--argc,++argv);
  }
  else if(strcmp(*argv,"print")==0){
    print(--argc,++argv);
  }else{
    fprintf(stderr,"Unknown argument: \'%s' please supply:\n",*argv);
    fprintf(stderr,"./thetaStat [make_bed|do_stat|validate_bed|print]\n");
  }
}

