#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <string>
#include <zlib.h>
#include <sys/stat.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include "stats.cpp"
#include "misc.h"


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
  int nChr;
}datum;

typedef std::map<char*,datum,ltstr> myMap;
//ssize_t bgzf_read(BGZF *fp, void *data, size_t length) HTS_RESULT_USED;
void my_bgzf_read(BGZF *fp, void *data, size_t length){
  ASSERT(length==bgzf_read(fp,data,length));
}

//checks that the index and data files start with the magic number from realSFS saf2theta,
//exits with a message if a file is missing, empty, truncated or from an old version
int isok(char *bgzf_name,char *fp_name){
  const char *rerun = "\t-> Was 'realSFS saf2theta' interrupted (eg. time or memory limit, or a full disk)? Then rerun it\n";
  const char *old = "\t-> It looks like input files have been generated with an older version of ANGSD. Either use older version of angsd <0.915 or rerun main angsd with later version of angsd >0.917\n";
  char buf[8];

  FILE *fp = NULL;
  fp=fopen(fp_name,"rb");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",fp_name);
    exit(1);
  }
  size_t nread = fread(buf,sizeof(char),8,fp);
  fclose(fp);
  if(nread!=8){
    fprintf(stderr,"\t-> The index file \'%s\' is empty or truncated\n%s",fp_name,rerun);
    exit(1);
  }
  if(memcmp(buf,"thetav2",8)!=0){
    fprintf(stderr,"\t-> Unknown magic number in \'%s\'\n%s",fp_name,old);
    exit(1);
  }

  BGZF *bgzf=NULL;
  bgzf=bgzf_open(bgzf_name,"rb");
  if(bgzf==NULL){
    fprintf(stderr,"\t-> Problems opening file: \'%s\'\n",bgzf_name);
    exit(1);
  }
  ssize_t nbg = bgzf_read(bgzf,buf,8*sizeof(char));
  bgzf_close(bgzf);
  if(nbg<0){
    fprintf(stderr,"\t-> Problem reading \'%s\', the file is damaged or not a bgzf/gzip file\n%s",bgzf_name,rerun);
    exit(1);
  }
  if(nbg!=8){
    fprintf(stderr,"\t-> The file \'%s\' is empty or truncated (%zd of the first 8 bytes could be read)\n%s",bgzf_name,nbg,rerun);
    exit(1);
  }
  if(memcmp(buf,"thetav2",8)!=0){
    fprintf(stderr,"\t-> Unknown magic number in \'%s\'\n%s",bgzf_name,old);
    exit(1);
  }
  return 1;
}



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


void writemap(FILE *fp,const myMap &mm){
  fprintf(fp,"\t\tInformation from index file:\n");
  int i=0;
  for(myMap::const_iterator it=mm.begin();it!=mm.end();++it){
    datum d = it->second;
    fprintf(fp,"\t\t%d\t%s\t%zu\t%ld\t%d\n",i++,it->first,d.nSites,(long int)d.fpos,d.nChr);

  }

}


myMap getMap(const char *fname){
  myMap ret;
  size_t clen;
  if(!fexists(fname)){
    fprintf(stderr,"Problem opening file: %s\n",fname);
    exit(1);
  }
  FILE *fp = fopen(fname,"r");
  char tmp[8];
  ASSERT(8==fread(tmp,1,8,fp));
  while(fread(&clen,sizeof(size_t),1,fp)){
    char *chr = new char[clen+1];
    ASSERT(clen==fread(chr,1,clen,fp));
    chr[clen] = '\0';

    datum d;
    if(1!=fread(&d.nSites,sizeof(size_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(1);
    }

    if(1!=fread(&d.nChr,sizeof(int),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(1);
    }

    if(1!=fread(&d.fpos,sizeof(int64_t),1,fp)){
      fprintf(stderr,"[%s.%s():%d] Problem reading data: %s \n",__FILE__,__FUNCTION__,__LINE__,fname);
      exit(1);
    }

    myMap::iterator it = ret.find(chr);
    if(it==ret.end())
      ret[chr] =d ;
    else{
      fprintf(stderr,"Problem with chr: %s, key already exists\n",chr);
      exit(1);
    }
  }
  fclose(fp);
  return ret;
}

void deleteMyMap(myMap &mm){
  for(myMap::iterator it=mm.begin();it!=mm.end();it++){
    delete [] it->first;
  }

}


typedef struct{
  char *chr;
  int nChr;
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
  my_bgzf_read(fp,ret.chr,clen);
  ret.chr[clen] = '\0';
  my_bgzf_read(fp,&ret.nSites,sizeof(size_t));
  my_bgzf_read(fp,&ret.nChr,sizeof(int));
  ret.posi = new int[ret.nSites];
  ret.tW = new float[ret.nSites];
  ret.tP = new float[ret.nSites];
  ret.tF = new float[ret.nSites];
  ret.tH = new float[ret.nSites];
  ret.tL = new float[ret.nSites];
  
  //read positions and thetas
  my_bgzf_read(fp,ret.posi,ret.nSites*sizeof(int));
  my_bgzf_read(fp,ret.tW,ret.nSites*sizeof(float));
  my_bgzf_read(fp,ret.tP,ret.nSites*sizeof(float));
  my_bgzf_read(fp,ret.tF,ret.nSites*sizeof(float));
  my_bgzf_read(fp,ret.tH,ret.nSites*sizeof(float));
  my_bgzf_read(fp,ret.tL,ret.nSites*sizeof(float));

  
  for(int i=0;i<ret.nSites;i++){
    ret.posi[i] = ret.posi[i]+1; //Old implemenation assummed positions was one indexed
  }

  
  //make thetas into normal space
  for(size_t i=0;i<ret.nSites;i++){
    if(1){
      if(std::isnan(ret.tW[i]))
	ret.tW[i] = log(0);
      if(std::isnan(ret.tP[i]))
	ret.tP[i] = log(0);
      if(std::isnan(ret.tF[i]))
	ret.tF[i] = log(0);
      if(std::isnan(ret.tH[i]))
	ret.tH[i] = log(0);
      if(std::isnan(ret.tL[i]))
	ret.tL[i] = log(0);
    }
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


char *idxToGz(char *one){
  char *ret = strdup(one);
  strcpy(ret+strlen(ret)-3,"gz\0");
  //  fprintf(stderr,"ret:%s\n",ret);
  return ret;
}




int do_stat(int argc, char**argv){
  if(argc==0){
    fprintf(stderr,"\n\t./thetaStat do_stat .thetas.idx [-win INT -step INT -type [0,1,2] -r chrName -nChr INT -outnames outputprefix]\n");
    fprintf(stderr,"\n\tSums the thetas in windows and calculates Tajima's D, Fu & Li's F and D, Fay's H and Zeng's E.\n");
    fprintf(stderr,"\tWrites one line per window to outputprefix.pestPG\n");
    fprintf(stderr,"\n\t-win INT\tWindow size in bases. Without -win and -step there is one window per chromosome\n");
    fprintf(stderr,"\t-step INT\tStep size in bases\n");
    fprintf(stderr,"\t-type INT\tWhere the first window starts (default 0)\n");
    fprintf(stderr,"\t\t\t0: the first multiple of -step after the first position with data\n");
    fprintf(stderr,"\t\t\t1: the first position with data\n");
    fprintf(stderr,"\t\t\t2: position 1\n");
    fprintf(stderr,"\t-r chrName\tOnly use this chromosome\n");
    fprintf(stderr,"\t-nChr INT\tNumber of chromosomes in the sample (2 x individuals for diploids).\n\t\t\tDefault is the number stored in the thetas file\n");
    fprintf(stderr,"\t-outnames STR\tOutput prefix, default is the input filename (angsdput.thetas.idx.pestPG)\n");
    fprintf(stderr,"\n\tExamples:\n\t1)./thetaStat do_stat angsdput.thetas.idx\n");
    fprintf(stderr,"\t2)./thetaStat do_stat angsdput.thetas.idx -win 5000 -step 1000\n");
    fprintf(stderr,"\t3)./thetaStat do_stat angsdput.thetas.idx -win 5000 -step 1000 -r chr1\n");
    fprintf(stderr,"\t4)./thetaStat do_stat angsdput.thetas.idx -win 5000 -step 1000 -r chr1 -nChr 20 -outnames newoutputname\n\n");
    return 0;
  }
  char *base = *argv;
  char* outnames_bin = idxToGz(base);
  char* outnames_idx = base;
  fprintf(stderr,"\tAssuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);
  if(isok(outnames_bin,outnames_idx)==0)
    return 0;

  myMap mm = getMap(outnames_idx);
  writemap(stderr,mm);
  BGZF *fp = bgzf_open(outnames_bin,"r");
  char tmp[8];
  my_bgzf_read(fp,tmp,8);

  --argc;++argv;
  //  fprintf(stderr,"argc=%d\n",argc);
  int argP =0;
  char *chr=NULL;
  char *outnames = NULL;
  int nChr =-1;
  int win =0;
  int step =0;
  int type =0;
  while(argP<argc){
    //   fprintf(stderr,"args=%s\n",argv[argP]);
    if(argP+1>=argc){
      fprintf(stderr,"\t-> %s needs an argument\n",argv[argP]);
      exit(1);
    }
    if(strcmp("-r",argv[argP])==0)
      chr = argv[argP+1];
    else if(strcmp("-outnames",argv[argP])==0)
      outnames = argv[argP+1];
    else if(strcmp("-step",argv[argP])==0)
      step = atoi(argv[argP+1]);
    else if(strcmp("-win",argv[argP])==0)
      win = atoi(argv[argP+1]);
    else if(strcmp("-type",argv[argP])==0)
      type = atoi(argv[argP+1]);
    else if(strcmp("-nChr",argv[argP])==0)
      nChr = atoi(argv[argP+1]);
    
    else {
      fprintf(stderr,"Unknown argument:%s\n",argv[argP]);
      exit(1);
    }
    argP +=2;
  }

  fprintf(stderr,"\t -r=%s outnames=%s step: %d win: %d\n",chr,outnames,step,win);

  if(win==0||step==0){
    fprintf(stderr,"\tWinsize equals zero or step size equals zero. Will use entire chromosome as window\n");
    win=step=0;
  }  
  
  if(chr!=NULL){  
    myMap::iterator it = mm.find(chr);
    if(it==mm.end()){
      fprintf(stderr,"\tProblem finding chr: %s in index\n",chr);
      exit(1);
    }
    datum d = it->second;
    ASSERT(bgzf_seek(fp,d.fpos,SEEK_SET)==0);
  }
  if(outnames==NULL)
    outnames = base;

  char *resname = append(outnames,RES);
  FILE *fpres = fopen(resname,"w");
  //fprintf(fpres,"## thetaStat VERSION: %s build:(%s,%s)\n",VERSION,__DATE__,__TIME__);
  fprintf(fpres,"#(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)\t");
  fprintf(fpres,"Chr\tWinCenter\t");
  fprintf(fpres,"tW\ttP\ttF\ttH\ttL\t");
  fprintf(fpres,"Tajima\tfuf\tfud\tfayh\tzeng\tnSites\n");
  while(1){
    perChr pc = getPerChr(fp);
    if(pc.nSites==0)
      break;
    fprintf(stderr,"\tpc.chr=%s pc.nSites=%zu firstpos=%d lastpos=%d\n",pc.chr,pc.nSites,pc.posi[0],pc.posi[pc.nSites-1]);
    if(nChr == -1)
      nChr = pc.nChr;
    kstring_t str = do_stat_main(pc,step,win,nChr,type);
    fwrite(str.s,1,str.l,fpres);//should clean up str, doesn't matter for this program;
    free(str.s);
    fflush(fpres);
    if(chr!=NULL)
      break;
    dalloc(pc);
  }
  free(outnames_bin);
  fclose(fpres);
  fprintf(stderr,"\tDumping file: \"%s\"\n",resname);
  bgzf_close(fp);
  deleteMyMap(mm);
  delete [] resname;
  return 0;
}

void print_main(perChr &pc,FILE *fp){
  for(size_t i=0;i<pc.nSites;i++)
    fprintf(fp,"%s\t%d\t%f\t%f\t%f\t%f\t%f\n",pc.chr,pc.posi[i],log(pc.tW[i]),log(pc.tP[i]),log(pc.tF[i]),log(pc.tH[i]),log(pc.tL[i]));
}


int print(int argc, char**argv){
  if(argc==0){
    fprintf(stderr,"\n\t./thetaStat print angsdput.thetas.idx [-r chrName]\n");
    fprintf(stderr,"\n\tPrints the thetas per site to stdout: chromosome, position and the log of the\n");
    fprintf(stderr,"\tWatterson, pairwise, singleton, H and L thetas\n");
    fprintf(stderr,"\n\t-r chrName\tOnly print this chromosome\n");
    fprintf(stderr,"\n\tExamples:\n\t1)./thetaStat print angsdput.thetas.idx\n");
    fprintf(stderr,"\t2)./thetaStat print angsdput.thetas.idx -r chr2\n\n");
    return 0;
  }
  char *base = *argv;
  char* outnames_bin = idxToGz(base);
  char* outnames_idx = base;
  if(isok(outnames_bin,outnames_idx)==0)
    return 0;
  fprintf(stderr,"Assuming binfile:%s and indexfile:%s\n",outnames_bin,outnames_idx);
  
  myMap mm = getMap(outnames_idx);
  writemap(stderr,mm);
  BGZF *fp = bgzf_open(outnames_bin,"r");
  char tmp[8];
  my_bgzf_read(fp,tmp,8);
  --argc;++argv;
  //  fprintf(stderr,"argc=%d\n",argc);
  int argP =0;
  char *chr=NULL;

  while(argP<argc){
    //   fprintf(stderr,"args=%s\n",argv[argP]);
    if(argP+1>=argc){
      fprintf(stderr,"\t-> %s needs an argument\n",argv[argP]);
      exit(1);
    }
    if(strcmp("-r",argv[argP])==0)
      chr = argv[argP+1];
    else {
      fprintf(stderr,"Unknown argument:%s\n",argv[argP]);
      exit(1);

    }
    argP +=2;
  }
  
  
  if(chr!=NULL){  
    myMap::iterator it = mm.find(chr);
    if(it==mm.end()){
      fprintf(stderr,"Problem finding chr: %s in index\n",chr);
      exit(1);
    }
    datum d = it->second;
    ASSERT(0==bgzf_seek(fp,d.fpos,SEEK_SET));
  }
  fprintf(stdout,"#Chromo\tPos\tWatterson\tPairwise\tthetaSingleton\tthetaH\tthetaL\n");
  while(1){
    perChr pc = getPerChr(fp);
    if(pc.nSites==0)
      break;
    fprintf(stderr,"pc.chr=%s pc.nSites=%zu pc.nChr=%d firstpos=%d lastpos=%d\n",pc.chr,pc.nSites,pc.nChr,pc.posi[0],pc.posi[pc.nSites-1]);
    print_main(pc,stdout);
    if(chr!=NULL)
      break;
    dalloc(pc);
  }
  deleteMyMap(mm);
  free(outnames_bin);
  bgzf_close(fp);
  return 0;
}


//reads exactly length bytes or exits, a damaged or truncated file should not give a partial output
void cat_read(BGZF *fp,void *data,size_t length,const char *fname){
  ssize_t nread = bgzf_read(fp,data,length);
  if(nread<0||(size_t)nread!=length){
    fprintf(stderr,"\t-> Problem reading \'%s\', the file is damaged or truncated\n",fname);
    exit(1);
  }
}

void cat_write(BGZF *fp,const void *data,size_t length,const char *fname){
  if(bgzf_write(fp,data,length)!=(ssize_t)length){
    fprintf(stderr,"\t-> Problem writing \'%s\'\n",fname);
    exit(1);
  }
}

//one filename per line, empty lines are skipped
std::vector<char*> readFilelist(const char *fname){
  std::vector<char*> ret;
  FILE *fp = fopen(fname,"r");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem opening file: \'%s\'\n",fname);
    exit(1);
  }
  char buf[LENS];
  while(fgets(buf,LENS,fp)){
    char *tok = strtok(buf," \t\r\n");
    if(tok)
      ret.push_back(strdup(tok));
  }
  fclose(fp);
  return ret;
}

//concatenates .thetas.idx/.thetas.gz files from realSFS saf2theta, the files must have different chromosomes
int theta_cat(int argc,char **argv){
  if(argc==0){
    fprintf(stderr,"\n\t./thetaStat cat file1.thetas.idx file2.thetas.idx [...] -outnames outputprefix\n");
    fprintf(stderr,"\t./thetaStat cat -b file.list -outnames outputprefix\n");
    fprintf(stderr,"\n\tMerges .thetas.idx/.thetas.gz files from 'realSFS saf2theta' into one, in the order given.\n");
    fprintf(stderr,"\tWrites outputprefix.thetas.idx and outputprefix.thetas.gz\n");
    fprintf(stderr,"\n\tA chromosome can only be in one of the files, and the files must be made from the same\n");
    fprintf(stderr,"\tnumber of sampled chromosomes (nChr, 2 x individuals for diploids)\n");
    fprintf(stderr,"\n\t-b file.list\tA file with one .thetas.idx filename per line\n");
    fprintf(stderr,"\t-outnames STR\tOutput prefix\n\n");
    return 0;
  }
  char *outnames = NULL;
  std::vector<char*> files;
  for(int i=0;i<argc;i++){
    if(strcmp("-outnames",argv[i])==0||strcmp("-b",argv[i])==0){
      if(i+1==argc){
	fprintf(stderr,"\t-> %s needs an argument\n",argv[i]);
	exit(1);
      }
      if(strcmp("-outnames",argv[i])==0)
	outnames = argv[++i];
      else{
	std::vector<char*> tmp = readFilelist(argv[++i]);
	files.insert(files.end(),tmp.begin(),tmp.end());
      }
    }else if(argv[i][0]=='-'){
      fprintf(stderr,"\t-> Unknown argument: %s\n",argv[i]);
      exit(1);
    }else
      files.push_back(strdup(argv[i]));
  }
  if(outnames==NULL){
    fprintf(stderr,"\t-> Must supply -outnames\n");
    exit(1);
  }
  if(files.size()==0){
    fprintf(stderr,"\t-> Must supply .thetas.idx files\n");
    exit(1);
  }
  char *out_idx_name = append(outnames,".thetas.idx");
  char *out_bin_name = append(outnames,".thetas.gz");
  for(size_t f=0;f<files.size();f++){
    size_t l = strlen(files[f]);
    if(l<4||strcmp(files[f]+l-4,IDX)!=0){
      fprintf(stderr,"\t-> \'%s\' is not a .thetas.idx file\n",files[f]);
      exit(1);
    }
    if(strcmp(files[f],out_idx_name)==0){
      fprintf(stderr,"\t-> The output \'%s\' is also an input file\n",out_idx_name);
      exit(1);
    }
    char *bin = idxToGz(files[f]);
    isok(bin,files[f]);
    free(bin);
  }

  FILE *out_idx = fopen(out_idx_name,"wb");
  BGZF *out_bin = bgzf_open(out_bin_name,"w");
  if(out_idx==NULL||out_bin==NULL){
    fprintf(stderr,"\t-> Problem opening output files: \'%s\' \'%s\'\n",out_idx_name,out_bin_name);
    exit(1);
  }
  char magic[8] = "thetav2";
  cat_write(out_bin,magic,8,out_bin_name);
  if(fwrite(magic,1,8,out_idx)!=8){
    fprintf(stderr,"\t-> Problem writing \'%s\'\n",out_idx_name);
    exit(1);
  }

  std::map<std::string,char*> seen;//chromosome -> file it came from
  int nChrAll = -1;
  size_t nChromo = 0;
  size_t bufLen = 1<<20;
  char *buf = new char[bufLen];
  for(size_t f=0;f<files.size();f++){
    char *bin_name = idxToGz(files[f]);
    BGZF *fp = bgzf_open(bin_name,"r");
    cat_read(fp,magic,8,bin_name);
    fprintf(stderr,"\t-> Adding: %s\n",bin_name);
    while(1){
      size_t clen;
      ssize_t nread = bgzf_read(fp,&clen,sizeof(size_t));
      if(nread==0)
	break;
      if(nread!=sizeof(size_t)||clen==0||clen>65536){
	fprintf(stderr,"\t-> Problem reading \'%s\', the file is damaged or truncated\n",bin_name);
	exit(1);
      }
      char *chr = new char[clen+1];
      cat_read(fp,chr,clen,bin_name);
      chr[clen] = '\0';
      size_t nSites;
      int nChr;
      cat_read(fp,&nSites,sizeof(size_t),bin_name);
      cat_read(fp,&nChr,sizeof(int),bin_name);

      std::map<std::string,char*>::iterator it = seen.find(chr);
      if(it!=seen.end()){
	fprintf(stderr,"\t-> Chromosome \'%s\' is in both \'%s\' and \'%s\', the files must have different chromosomes\n",chr,it->second,files[f]);
	exit(1);
      }
      seen[chr] = files[f];
      if(nChrAll==-1)
	nChrAll = nChr;
      else if(nChrAll!=nChr){
	fprintf(stderr,"\t-> Chromosome \'%s\' in \'%s\' is made from %d sampled chromosomes (nChr), the earlier files from %d\n",chr,files[f],nChr,nChrAll);
	exit(1);
      }

      int64_t offs = bgzf_tell(out_bin);
      cat_write(out_bin,&clen,sizeof(size_t),out_bin_name);
      cat_write(out_bin,chr,clen,out_bin_name);
      cat_write(out_bin,&nSites,sizeof(size_t),out_bin_name);
      cat_write(out_bin,&nChr,sizeof(int),out_bin_name);
      //positions and the 5 thetas are copied as they are
      size_t left = nSites*(sizeof(int)+5*sizeof(float));
      while(left>0){
	size_t n = left<bufLen?left:bufLen;
	cat_read(fp,buf,n,bin_name);
	cat_write(out_bin,buf,n,out_bin_name);
	left -= n;
      }

      if(fwrite(&clen,sizeof(size_t),1,out_idx)!=1||fwrite(chr,1,clen,out_idx)!=clen||
	 fwrite(&nSites,sizeof(size_t),1,out_idx)!=1||fwrite(&nChr,sizeof(int),1,out_idx)!=1||
	 fwrite(&offs,sizeof(int64_t),1,out_idx)!=1){
	fprintf(stderr,"\t-> Problem writing \'%s\'\n",out_idx_name);
	exit(1);
      }
      nChromo++;
      delete [] chr;
    }
    bgzf_close(fp);
    free(bin_name);
  }
  delete [] buf;
  if(bgzf_close(out_bin)!=0||fclose(out_idx)!=0){
    fprintf(stderr,"\t-> Problem closing output files: \'%s\' \'%s\'\n",out_idx_name,out_bin_name);
    exit(1);
  }
  fprintf(stderr,"\t-> Wrote %zu chromosomes from %zu files to: \'%s\' and \'%s\'\n",nChromo,files.size(),out_idx_name,out_bin_name);
  for(size_t f=0;f<files.size();f++)
    free(files[f]);
  delete [] out_idx_name;
  delete [] out_bin_name;
  return 0;
}


int main(int argc,char **argv){
  if(argc==1){
    fprintf(stderr,"\n\t\'./thetaStat\': a program to do neutrality test statistics using thetas.idx output from angsd\n");
    fprintf(stderr,"\n\tdo_stat\tthetas and neutrality test statistics in windows or per chromosome\n");
    fprintf(stderr,"\tprint\tprint the thetas per site\n");
    fprintf(stderr,"\tcat\tmerge thetas files with different chromosomes\n");
    fprintf(stderr,"\n\tExamples:\n\t1) ./thetaStat print angsdput.thetas.idx\n");
    fprintf(stderr,"\t2) ./thetaStat do_stat angsdput.thetas.idx -win 5000 -step 1000\n");
    fprintf(stderr,"\t3) ./thetaStat cat chr1.thetas.idx chr2.thetas.idx -outnames merged\n");
    fprintf(stderr,"\n\tType \'./thetaStat do_stat\', './thetaStat print' or './thetaStat cat' for more information\n\n");
    return 0;
  }
  for(int i=0;i<argc;i++)
    fprintf(stderr," %s",argv[i]);
  fprintf(stderr,"\n");

  --argc,++argv;
  if(strcmp(*argv,"do_stat")==0)
    return do_stat(--argc,++argv);
  else if(strcmp(*argv,"print")==0)
    return print(--argc,++argv);
  else if(strcmp(*argv,"cat")==0)
    return theta_cat(--argc,++argv);

  fprintf(stderr,"Unknown argument: \'%s' please supply:\n",*argv);
  fprintf(stderr,"./thetaStat [do_stat||print||cat]\n");
  return 1;
}

