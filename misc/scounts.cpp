#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cassert>
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


bMap readvcf(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(0);
  }
  int at=0;
  bMap bm;
  char buf[1024];
  while(gzgets(gz,buf,1024)){
    char *line = strtok(buf,"\n");
    bm[at++]= strdup(line);
  }
  fprintf(stderr,"\t-> read nsites from vcf:%lu\n",bm.size());
  return bm;
}

bMap2 readvcf2(const char *fname){
  gzFile gz = Z_NULL;
  gz=gzopen(fname,"rb");
  if(gz==Z_NULL){
    fprintf(stderr,"\t-> Problem opening file: %s\n",fname);
    exit(0);
  }
  int at=0;
  bMap2 bm;
  char buf[1024];
  while(gzgets(gz,buf,1024)){
    datum d;
    d.chr = strdup(strtok(buf,"\t \n"));
    d.pos = atoi(strtok(NULL,"\t \n"));
    d.al1 = strtok(NULL,"\t \n")[0];
    d.al2 = strtok(NULL,"\t \n")[0];
    bm[at++]= d;
  }
  fprintf(stderr,"\t-> read nsites from vcf:%lu\n",bm.size());
  return bm;
}




void set(counts &cnt,int type,int val){
  if(type==0)
    cnt.A=val;
  if(type==1)
    cnt.C=val;
  if(type==2)
    cnt.G=val;
  if(type==3)
    cnt.T=val;
}

void printit(counts &cnt,FILE *fp){
  fprintf(stdout,"%u\t%d\t%d\t%d\t%d\n",cnt.rel_pos,cnt.A,cnt.C,cnt.G,cnt.T);
}

FILE *myfopen(const char *file,char *mode){
  FILE *fp = NULL;
  fp = fopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(0);
  }
  return fp;
}

gzFile mygzopen(const char *file,char *mode){
  gzFile fp = Z_NULL;
  fp = gzopen(file,mode);
  if(fp==NULL){
    fprintf(stderr,"problem opening file: %s\n",file);
    exit(0);
  }
  return fp;
}

typedef std::map<int,int*> smap;


void counter(smap &sm,char *fname){
  gzFile dat = mygzopen(fname,(char*)"rb");
  counts cnts;
  while(gzread(dat,&cnts,sizeof(counts))){
    smap::iterator it=sm.find(cnts.rel_pos);
    
    if(it==sm.end()){
      int *tmp = new int[4];
      tmp[0]=tmp[1]=tmp[2]=tmp[3]=0;
      sm[cnts.rel_pos] = tmp;
      it = sm.find(cnts.rel_pos);
    }
    int *total = it->second;
    int tmp[4] = {cnts.A,cnts.C,cnts.G,cnts.T};
    int tmp2[4] = {0,0,0,0};
    int at=0;
    for(int i=0;i<4;i++)
      if(tmp[i]>0)
	tmp2[at++] =i;
    //    fprintf(stderr,"informative sites:%d\n",at);
    int which = lrand48() % at;
    it->second[tmp2[which]]++;
  }
}


int nlines(char *fname,char *popname){
  gzFile gz = mygzopen(fname,(char*)"rb");
  char buf[4096];

  int ret =0;
  while(gzgets(gz,buf,4096)){
    char *file = strtok(buf," \t\n");
    char *pop1 = strtok(NULL," \t\n");
    //  fprintf(stderr,"\t-> file:%s pop1:%s\n",file,pop1);
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(0);
    }
    if(popname==NULL)
      ret++;
    else if(0==strcasecmp(popname,pop1)){
      ret++;
    }
  }
  return ret;
}


void merge1(char *flist,bMap &bm, char *pop,char *odir){
  fprintf(stderr,"\t-> printing merge1\n");
  smap sm;
  int nind = nlines(flist,pop);
  fprintf(stderr,"\t-> Number of individuals in pop:%s is:%d\n",pop,nind);
  gzFile gz = mygzopen(flist,(char*)"rb");
  char buf[4096];

  char oname[1024];
  if(odir==NULL)
    sprintf(oname,"%s.%d.txt",pop,nind);
  else
    sprintf(oname,"%s/%s.%d.txt",odir,pop,nind);
  FILE *fp = myfopen(oname,(char*)"wb");
  
  while(gzgets(gz,buf,4096)){
    char *file = strtok(buf," \t\n");
    char *pop1 = strtok(NULL," \t\n");
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(0);
    }
    
    if(0==strcasecmp(pop,pop1)){
      counter(sm,file);
    }
  }
  for(smap::iterator it=sm.begin();it!=sm.end();it++){
    int *ary = it->second;
    if(bm.size()!=0)
      fprintf(fp,"%s\t%d\t%d\t%d\t%d\n",bm[it->first],ary[0],ary[1],ary[2],ary[3]);
    else
      fprintf(fp,"%d\t%d\t%d\t%d\t%d\n",it->first,ary[0],ary[1],ary[2],ary[3]);
  }
  fprintf(stderr,"\t-> end merge1\n");
}

int get(char c){
  if(c=='A')
    return 0;
  if(c=='C')
    return 1;
  if(c=='G')
    return 2;
  if(c=='T')
    return 3;
  fprintf(stderr,"unknown base:%c\n",c);
  return 0;
}

void helper(char *fname,std::map<int,double *> &freq,std::map<int,int*> &nHap,bMap2 &bm2,int total,int at){
  // fprintf(stderr,"fname:%s\n",fname);
  gzFile gz = mygzopen(fname,(char*)"rb");
  char buf[4096];
  while(gzgets(gz,buf,4096)){
    int index = atoi(strtok(buf,"\t\n "));

    int cnts[4];
    for(int i=0;i<4;i++)
      cnts[i] = atoi(strtok(NULL,"\t\n "));
    //    fprintf(stderr,"%d %d %d %d\n",cnts[0],cnts[1],cnts[2],cnts[3]);
    std::map<int, double *>::iterator fit = freq.find(index);
    if(fit==freq.end()){
      double *tmp = new double[total];
      int *tmp2 = new int[total];
      for(int i=0;i<total;i++){
	tmp[i] = 0;
	tmp2[i] = 0;
      }
      freq[index] = tmp;
      nHap[index] = tmp2;
    }
    bMap2::iterator bit = bm2.find(index);
    if(bit==bm2.end()){
      fprintf(stderr,"problem finding:%d\n",index);
      exit(0);

    }
    double tot = cnts[get(bit->second.al1)]+cnts[get(bit->second.al2)];
    double f =((double) cnts[get(bit->second.al2)])/tot ;
    fit = freq.find(index);
    double *ff = fit->second; ff[at] = f;
    std::map<int, int *>::iterator nit = nHap.find(index);
    int *nh = nit->second; nh[at] = tot;
  }
  
}

void merge2(char *flist,bMap2 &bm2,char *odir){
  fprintf(stderr,"\t-> printing merge2 flist:%s\n",flist);
  int npops = nlines(flist,NULL);
  fprintf(stderr,"\t-> npops: %d\n",npops);
  std::vector<int> nind;
  std::vector<char*> popnames;
  smap sm;
  gzFile gz = mygzopen(flist,(char*)"rb");
  char buf[4096];
  std::map<int, double *> freqs;
  std::map<int, int *> nHap;
  int at =0;
  while(gzgets(gz,buf,4096)){
    char *file = strtok(buf," \t\n");
    char *dup = strdup(file);
    if(!fexists(file)){
      fprintf(stderr,"file: %s doesnt exists\n",file);
      exit(0);
    }
    char *bname=basename(file);
    //    fprintf(stderr,"tmpdup:%s\n",bname);
    char *pop1 = strtok(bname,".\n");
    int nind_pop = atoi(strtok(NULL,".\n"));
    //    fprintf(stderr,"\npop: %s has:%d\n",pop1,nind_pop);
    nind.push_back(nind_pop);
    popnames.push_back(strdup(pop1));
    helper(dup,freqs,nHap,bm2,npops,at++);
  }
  //  fprintf(stderr,"number of popnames:%lu\n",popnames.size());
  fprintf(stderr,"\t-> end merge2\n");

  char oname[1024];
  if(odir==NULL)
    sprintf(oname,"freqs.txt");
  else
    sprintf(oname,"%s/freqs.txt",odir);
  FILE *fp = myfopen(oname,(char*)"wb");
  fprintf(fp,"#chr\tpos");
  for(int i=0;i<popnames.size();i++)
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
    sprintf(oname,"nhap.txt");
  else
    sprintf(oname,"%s/nhap.txt",odir);
  fp = myfopen(oname,(char*)"wb");
  fprintf(fp,"#chr\tpos");
  for(int i=0;i<popnames.size();i++)
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
    sprintf(oname,"siteinfo.txt");
  else
    sprintf(oname,"%s/siteinfo.txt",odir);

  fp = myfopen(oname,(char*)"wb");
  fprintf(fp,"#chromo\tposition\ttype\n");
  for(std::map<int,double *>::iterator fit=freqs.begin();fit!=freqs.end();fit++){
    bMap2::iterator it = bm2.find(fit->first); 
    fprintf(fp,"%s\t%d\t%c%c\n",it->second.chr,it->second.pos,it->second.al1,it->second.al2);
  }
  fclose(fp);
  if(odir==NULL)
    sprintf(oname,"nind.txt");
  else
    sprintf(oname,"%s/nind.txt",odir);

  fp = myfopen(oname,(char*)"wb");
  for(int i=0;i<popnames.size();i++)
    fprintf(fp,"%s\t%d\n",popnames[i],nind[i]);
  fclose(fp);
}


int main(int argc,char **argv){
  long seed =0;
  unsigned nsites = 70e6;
  
  if(argc==1){
    fprintf(stderr,"Examples:\ninputfiles are)\n ./angsd  -doscounts 1 -i file.bam -vcfname all.ceu.sites.gz -docounts 1 -sites all.ceu.sites.gz\n\n1) ./scounts -seed 10 -list scounts.list -mergeSingle -pop CEU \n  (output are CEU.NSAMPLES.txt)(intputfiles are filelist of the angsd outputfiles)\n2)./scounts -list flist -mergeMulti -info all.ceu.sites.gz\n (inputfiles are from step1 above)(outputfiles are: \'siteinfo.txt\',\'nind.txt\',\'nhap.txt\',\'freqs.txt\')\n\nPrinting textversion of binary files)\n\ta)./scounts print CEU.16.txt\n");
    return 0;
  }
  if(0==strcmp(argv[1],"print")){
    fprintf(stderr,"printing\n");
    gzFile dat = mygzopen(argv[2],(char*)"rb");
    counts cnts;
    while(gzread(dat,&cnts,sizeof(counts)))
      printit(cnts,stdout);
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
      odir = strdup(*(++argv));
    else if(!strcmp(*argv,"-pop"))
       pop = *(++argv);
    else {
      fprintf(stderr,"\t-> unknown argument: %s\n",*argv);
      return 0;
    }
    argv++;
  }
  fprintf(stderr,"seed:%ld filelist:%s info:%s mergeSingle:%d mergeMulti:%d pop:%s odir:%s\n",seed,flist,info,mergeSingle,mergeMulti,pop,odir);
  if(mergeSingle&&info)
    fprintf(stderr,"\t-> remove -info vcffile if you want to merge population files\n");
  
  if(mergeSingle&&pop==NULL){
    fprintf(stderr,"\t-> Must supply -pop for merging single populations\n");
    return 0;
  }
  if(mergeMulti&&pop){
    fprintf(stderr,"\t-> shouldnt supply -pop for merging population files\n");
    return 0;
  }
 if(mergeMulti&&info==NULL){
    fprintf(stderr,"\t-> need to supply -info when using mergeMulti\n");
    return 0;
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
