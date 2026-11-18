#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <zlib.h>
#include <sys/stat.h>
#include <vector>
#include <map>
#include <libgen.h>


typedef unsigned char uchar;

typedef struct{
  unsigned int rel_pos;
  uchar A:2;
  uchar C:2;
  uchar G:2;
  uchar T:2;
}counts;


typedef struct{
  char *chr;
  int pos;
  char al1;
  char al2;
}datum;

typedef std::map<int, char *> bMap;
typedef std::map<int, datum> bMap2;

int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

//next whitespace separated token of a line, exits with an error if the line has too few columns
char *nextTok(char *str,const char *fname,int lineno){
  char *tok = strtok(str,"\t \n\r");
  if(tok==NULL){
    fprintf(stderr,"\t-> Line %d in file \'%s\' has too few columns\n",lineno,fname);
    exit(1);
  }
  return tok;
}

//true if the line is empty or only whitespace
int isBlank(const char *buf){
  for(;*buf;buf++)
    if(!isspace((unsigned char)*buf))
      return 0;
  return 1;
}

//the sites file is indexed by line number (rel_pos), so every line must be a site
bMap readvcf(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(1);
  }
  int at=0;
  bMap bm;
  char buf[1024];
  while(gzgets(gz,buf,1024)){
    char *line = strtok(buf,"\n\r");
    if(line==NULL){
      fprintf(stderr,"\t-> Line %d in file \'%s\' is empty\n",at+1,fname);
      exit(1);
    }
    bm[at++]= strdup(line);
  }
  gzclose(gz);
  fprintf(stderr,"\t-> read nsites from vcf:%lu\n",bm.size());
  return bm;
}

bMap2 readvcf2(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(1);
  }
  int at=0;
  bMap2 bm;
  char buf[1024];
  while(gzgets(gz,buf,1024)){
    datum d;
    d.chr = strdup(nextTok(buf,fname,at+1));
    d.pos = atoi(nextTok(NULL,fname,at+1));
    d.al1 = toupper(nextTok(NULL,fname,at+1)[0]);
    d.al2 = toupper(nextTok(NULL,fname,at+1)[0]);
    bm[at++]= d;
  }
  gzclose(gz);
  fprintf(stderr,"\t-> read nsites from vcf:%lu\n",bm.size());
  return bm;
}


void printit(counts &cnt,FILE *fp){
  fprintf(fp,"%u\t%d\t%d\t%d\t%d\n",cnt.rel_pos,cnt.A,cnt.C,cnt.G,cnt.T);
}

FILE *myfopen(const char *file,const char *mode){
  FILE *fp = NULL;
  fp = fopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(1);
  }
  return fp;
}

gzFile mygzopen(const char *file,const char *mode){
  gzFile fp = Z_NULL;
  fp = gzopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(1);
  }
  return fp;
}

//read one counts record. Returns 1 on success and 0 at end of file, exits on errors or truncated records
int readRecord(gzFile gz,counts &cnts,const char *fname){
  int ret = gzread(gz,&cnts,sizeof(counts));
  if(ret==0)
    return 0;
  if(ret!=(int)sizeof(counts)){
    fprintf(stderr,"\t-> Problem reading \'%s\' (got %d of %lu bytes), file truncated or corrupt?\n",fname,ret,(unsigned long)sizeof(counts));
    exit(1);
  }
  return 1;
}

typedef std::map<int,int*> smap;


//for every site, sample one read from the individual and count its allele
void counter(smap &sm,char *fname){
  gzFile dat = mygzopen(fname,"rb");
  counts cnts;
  size_t nEmpty = 0;
  while(readRecord(dat,cnts,fname)){
    int tmp[4] = {cnts.A,cnts.C,cnts.G,cnts.T};
    int depth = tmp[0]+tmp[1]+tmp[2]+tmp[3];
    if(depth==0){//older angsd versions wrote sites without data
      nEmpty++;
      continue;
    }
    smap::iterator it=sm.find(cnts.rel_pos);

    if(it==sm.end()){
      int *ary = new int[4];
      ary[0]=ary[1]=ary[2]=ary[3]=0;
      sm[cnts.rel_pos] = ary;
      it = sm.find(cnts.rel_pos);
    }
    //pick a read at random, so each base is chosen proportional to its depth
    int which = lrand48() % depth;
    int i=0;
    while(which>=tmp[i]){
      which -= tmp[i];
      i++;
    }
    it->second[i]++;
  }
  gzclose(dat);
  if(nEmpty>0)
    fprintf(stderr,"\t-> Skipped %lu records without data in: %s\n",(unsigned long)nEmpty,fname);
}


//number of files in the list, or the number belonging to popname if it is given
int nlines(char *fname,char *popname){
  gzFile gz = mygzopen(fname,"rb");
  char buf[4096];

  int ret =0;
  int lineno =0;
  while(gzgets(gz,buf,4096)){
    lineno++;
    if(isBlank(buf))
      continue;
    char *file = nextTok(buf,fname,lineno);
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(1);
    }
    if(popname==NULL)
      ret++;
    else if(0==strcasecmp(popname,nextTok(NULL,fname,lineno))){
      ret++;
    }
  }
  gzclose(gz);
  return ret;
}


void merge1(char *flist,bMap &bm, char *pop,char *odir){
  fprintf(stderr,"\t-> printing merge1\n");
  smap sm;
  int nind = nlines(flist,pop);
  fprintf(stderr,"\t-> Number of individuals in pop:%s is:%d\n",pop,nind);
  gzFile gz = mygzopen(flist,"rb");
  char buf[4096];

  char oname[4096];
  if(odir==NULL)
    snprintf(oname,4096,"%s.%d.txt",pop,nind);
  else
    snprintf(oname,4096,"%s/%s.%d.txt",odir,pop,nind);
  FILE *fp = myfopen(oname,"wb");

  int lineno =0;
  while(gzgets(gz,buf,4096)){
    lineno++;
    if(isBlank(buf))
      continue;
    char *file = nextTok(buf,flist,lineno);
    char *pop1 = nextTok(NULL,flist,lineno);
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(1);
    }

    if(0==strcasecmp(pop,pop1)){
      counter(sm,file);
    }
  }
  gzclose(gz);
  for(smap::iterator it=sm.begin();it!=sm.end();it++){
    int *ary = it->second;
    if(bm.size()!=0){
      bMap::iterator bit = bm.find(it->first);
      if(bit==bm.end()){
	fprintf(stderr,"\t-> Site index %d is not in the -info file\n",it->first);
	exit(1);
      }
      fprintf(fp,"%s\t%d\t%d\t%d\t%d\n",bit->second,ary[0],ary[1],ary[2],ary[3]);
    }else
      fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",it->first,ary[0],ary[1],ary[2],ary[3]);
  }
  fclose(fp);
  fprintf(stderr,"\t-> end merge1\n");
}

int get(char c){
  switch(toupper(c)){
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  }
  fprintf(stderr,"\t-> Unknown base in -info file: \'%c\'\n",c);
  exit(1);
}

void helper(char *fname,std::map<int,double *> &freq,std::map<int,int*> &nHap,bMap2 &bm2,int total,int at){
  // fprintf(stderr,"fname:%s\n",fname);
  gzFile gz = mygzopen(fname,"rb");
  char buf[4096];
  int lineno =0;
  while(gzgets(gz,buf,4096)){
    lineno++;
    if(isBlank(buf))
      continue;
    int index = atoi(nextTok(buf,fname,lineno));

    int cnts[4];
    for(int i=0;i<4;i++)
      cnts[i] = atoi(nextTok(NULL,fname,lineno));
    //    fprintf(stderr,"%d %d %d %d\n",cnts[0],cnts[1],cnts[2],cnts[3]);
    bMap2::iterator bit = bm2.find(index);
    if(bit==bm2.end()){
      fprintf(stderr,"problem finding:%d\n",index);
      exit(1);
    }
    std::map<int, double *>::iterator fit = freq.find(index);
    if(fit==freq.end()){
      double *tmp = new double[total];
      int *tmp2 = new int[total];
      for(int i=0;i<total;i++){
	tmp[i] = 0;//no data for this population, nhap is 0
	tmp2[i] = 0;
      }
      freq[index] = tmp;
      nHap[index] = tmp2;
      fit = freq.find(index);
    }
    int tot = cnts[get(bit->second.al1)]+cnts[get(bit->second.al2)];
    if(tot>0)
      fit->second[at] = ((double) cnts[get(bit->second.al2)])/tot;
    nHap[index][at] = tot;
  }
  gzclose(gz);
}

void merge2(char *flist,bMap2 &bm2,char *odir){
  fprintf(stderr,"\t-> printing merge2 flist:%s\n",flist);
  int npops = nlines(flist,NULL);
  fprintf(stderr,"\t-> npops: %d\n",npops);
  std::vector<int> nind;
  std::vector<char*> popnames;
  gzFile gz = mygzopen(flist,"rb");
  char buf[4096];
  std::map<int, double *> freqs;
  std::map<int, int *> nHap;
  int at =0;
  int lineno =0;
  while(gzgets(gz,buf,4096)){
    lineno++;
    if(isBlank(buf))
      continue;
    char *file = nextTok(buf,flist,lineno);
    char *dup = strdup(file);
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(1);
    }
    //files from -mergeSingle are named POP.NIND.txt
    char *bname=basename(file);
    char *pop1 = strtok(bname,".\n");
    char *nindStr = pop1?strtok(NULL,".\n"):NULL;
    if(pop1==NULL||nindStr==NULL){
      fprintf(stderr,"\t-> File \'%s\' is not named POP.NIND.txt as written by -mergeSingle\n",dup);
      exit(1);
    }
    nind.push_back(atoi(nindStr));
    popnames.push_back(strdup(pop1));
    helper(dup,freqs,nHap,bm2,npops,at++);
    free(dup);
  }
  gzclose(gz);
  //  fprintf(stderr,"number of popnames:%lu\n",popnames.size());
  fprintf(stderr,"\t-> end merge2\n");

  char oname[4096];
  if(odir==NULL)
    snprintf(oname,4096,"freqs.txt");
  else
    snprintf(oname,4096,"%s/freqs.txt",odir);
  FILE *fp = myfopen(oname,"wb");
  fprintf(fp,"#chr\tpos");
  for(size_t i=0;i<popnames.size();i++)
    fprintf(fp,"\t%s",popnames[i]);
  fprintf(fp,"\n");

  for(std::map<int, double *>::iterator fit=freqs.begin();fit!=freqs.end();fit++){
    //    fprintf(fp,"%d",fit->first);
    bMap2::iterator it = bm2.find(fit->first);
    fprintf(fp,"%s\t%d",it->second.chr,it->second.pos);
    double *ff=fit->second;
    for(int i=0;i<npops;i++)
      fprintf(fp,"\t%f",ff[i]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  if(odir==NULL)
    snprintf(oname,4096,"nhap.txt");
  else
    snprintf(oname,4096,"%s/nhap.txt",odir);
  fp = myfopen(oname,"wb");
  fprintf(fp,"#chr\tpos");
  for(size_t i=0;i<popnames.size();i++)
    fprintf(fp,"\t%s",popnames[i]);
  fprintf(fp,"\n");
  for(std::map<int, int *>::iterator hit=nHap.begin();hit!=nHap.end();hit++){
    bMap2::iterator it = bm2.find(hit->first);
    //    fprintf(fp,"%d",hit->first);
    fprintf(fp,"%s\t%d",it->second.chr,it->second.pos);
    int *ff=hit->second;
    for(int i=0;i<npops;i++)
      fprintf(fp,"\t%d",ff[i]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  if(odir==NULL)
    snprintf(oname,4096,"siteinfo.txt");
  else
    snprintf(oname,4096,"%s/siteinfo.txt",odir);

  fp = myfopen(oname,"wb");
  fprintf(fp,"#chromo\tposition\ttype\n");
  for(std::map<int,double *>::iterator fit=freqs.begin();fit!=freqs.end();fit++){
    bMap2::iterator it = bm2.find(fit->first);
    fprintf(fp,"%s\t%d\t%c%c\n",it->second.chr,it->second.pos,it->second.al1,it->second.al2);
  }
  fclose(fp);
  if(odir==NULL)
    snprintf(oname,4096,"nind.txt");
  else
    snprintf(oname,4096,"%s/nind.txt",odir);

  fp = myfopen(oname,"wb");
  for(size_t i=0;i<popnames.size();i++)
    fprintf(fp,"%s\t%d\n",popnames[i],nind[i]);
  fclose(fp);
}

void usage(){
  fprintf(stderr,"Examples:\ninputfiles are)\n ./angsd  -doscounts 1 -i file.bam -vcfname all.ceu.sites.gz -docounts 1 -sites all.ceu.sites.gz\n\n1) ./scounts -seed 10 -list scounts.list -mergeSingle -pop CEU [-info all.ceu.sites.gz]\n  (output are CEU.NSAMPLES.txt)(intputfiles are filelist of the angsd outputfiles)\n  (with -info the first column is the site from the info file instead of its index)\n2)./scounts -list flist -mergeMulti -info all.ceu.sites.gz\n (inputfiles are from step1 above)(outputfiles are: \'siteinfo.txt\',\'nind.txt\',\'nhap.txt\',\'freqs.txt\')\n\nPrinting textversion of binary files)\n\ta)./scounts print file.scounts.gz\n");
}

int main(int argc,char **argv){
  long seed =0;
  fprintf(stderr,"\t-> NOTE: scounts is deprecated and will be removed in a future version of ANGSD\n");

  if(argc==1){
    usage();
    return 1;
  }
  if(0==strcmp(argv[1],"print")){
    if(argc<3){
      fprintf(stderr,"\t-> print requires a .scounts.gz file\n");
      return 1;
    }
    fprintf(stderr,"printing\n");
    gzFile dat = mygzopen(argv[2],"rb");
    counts cnts;
    while(readRecord(dat,cnts,argv[2]))
      printit(cnts,stdout);
    gzclose(dat);
    fprintf(stderr,"done printing\n");
    return 0;
  }
  char *flist = NULL;
  char *info = NULL;
  char *pop = NULL;
  int mergeSingle =0;
  int mergeMulti =0;
  char *odir = NULL;
  ++argv;
  while(*argv){
    // fprintf(stderr,"argv:%s\n",*argv);
    int takesValue = !strcmp(*argv,"-seed")||!strcmp(*argv,"-list")||!strcmp(*argv,"-info")||!strcmp(*argv,"-odir")||!strcmp(*argv,"-pop");
    if(takesValue&&argv[1]==NULL){
      fprintf(stderr,"\t-> Argument \'%s\' requires a value\n",*argv);
      return 1;
    }
    if(!strcmp(*argv,"-seed"))
      seed = atol(*(++argv));
    else if(!strcmp(*argv,"-list"))
      flist = *(++argv);
    else if(!strcmp(*argv,"-info"))
      info = *(++argv);
    else if(!strcmp(*argv,"-mergeSingle"))
      mergeSingle = 1;
    else if(!strcmp(*argv,"-mergeMulti"))
      mergeMulti = 1;
    else if(!strcmp(*argv,"-odir"))
      odir = *(++argv);
    else if(!strcmp(*argv,"-pop"))
       pop = *(++argv);
    else {
      fprintf(stderr,"\t-> unknown argument: %s\n",*argv);
      return 1;
    }
    argv++;
  }
  fprintf(stderr,"seed:%ld filelist:%s info:%s mergeSingle:%d mergeMulti:%d pop:%s odir:%s\n",seed,flist,info,mergeSingle,mergeMulti,pop,odir);

  if(!mergeSingle&&!mergeMulti){
    fprintf(stderr,"\t-> Must supply -mergeSingle or -mergeMulti\n");
    usage();
    return 1;
  }
  if(flist==NULL){
    fprintf(stderr,"\t-> Must supply -list\n");
    return 1;
  }
  if(mergeSingle&&pop==NULL){
    fprintf(stderr,"\t-> Must supply -pop for merging single populations\n");
    return 1;
  }
  if(mergeMulti&&pop){
    fprintf(stderr,"\t-> shouldnt supply -pop for merging population files\n");
    return 1;
  }
 if(mergeMulti&&info==NULL){
    fprintf(stderr,"\t-> need to supply -info when using mergeMulti\n");
    return 1;
  }
  srand48(seed);
  bMap bm;
  bMap2 bm2;

  if(info&&mergeSingle)
    bm = readvcf(info);
  if(info&&mergeMulti)
    bm2 = readvcf2(info);

  if(mergeSingle)
    merge1(flist,bm,pop,odir);
  if(mergeMulti)
    merge2(flist,bm2,odir);

  return 0;
}
