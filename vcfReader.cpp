#include <cassert>
#include "analysisFunction.h"
#include "vcfReader.h"
#include <cfloat>

vcfReader::~vcfReader(){
  if(buf!=NULL)
    free(buf);
};

vcfReader::vcfReader(int &nInd_a,gzFile gz_a,int bytesPerLine,const aMap *revMap_a){

  curChr=-1;
  nInd = nInd_a;
  gz=gz_a;
  len=bytesPerLine;
  buf=NULL;
  buf=(char*)malloc(len);
  revMap=revMap_a;
  while(gzgets(gz,buf,len)){
    saveptr=buf;
    if(!strncmp(saveptr,"#CHROM",5)) {
      int i=0;
      char *tok=NULL;
      while(((tok=angsd::strpop(&saveptr,'\t')))[0]!='\0'){
	i++;
      }
      
      if(i-9!=nInd){
	fprintf(stderr,"\t-> Looks like wrong -nInd has supplied compared to header of vcf file: %d vs %d\n",nInd,i-9);
	exit(0);
      }
      break;
    }
    

  }
}

//return value is pos. If this is zero means dont use site
int vcfReader::parseline(double **lk,double **gp,char &major,char &minor){
  int pos = atoi(angsd::strpop(&saveptr,'\t'));
  angsd::strpop(&saveptr,'\t');//ID
  char *ref = angsd::strpop(&saveptr,'\t');
  char *alt = angsd::strpop(&saveptr,'\t');
  char *qual = angsd::strpop(&saveptr,'\t');
  char *filter =angsd::strpop(&saveptr,'\t');
  char *info = angsd::strpop(&saveptr,'\t');
  char *format= angsd::strpop(&saveptr,'\t');
  
  // fprintf(stderr,"\npos:%d ref:%s alt:%s qual:%s flt:%s info:%s format:%s strlen(ref):%zu strlen(alt):%zu\n",pos,ref,alt,qual,filter,info,format,strlen(ref),strlen(alt)); exit(0);
 
 if(strlen(ref)!=1||strlen(alt)!=1){
    fprintf(stderr,"skipping site:%d multiref/multialt (indel)\n",pos+1);
    return 0;
  }
  if(strcmp(filter,"PASS")){
    fprintf(stderr,"skipping site:%d (non PASS)\n",pos+1);
    return 0;
    
  }
    
  if(strlen(ref)!=1||strlen(alt)!=1){
    fprintf(stderr,"skipping site:%d (indel)\n",pos+1);
    return 0;
  }
  ref[0]=refToInt[ref[0]];
  alt[0]=refToInt[alt[0]];
  if(ref[0]==4|| alt[0]==4){
    fprintf(stderr,"REF/ALT is 'N' will discard site: %d \n",pos+1);
    return 0;
  }
  if(ref[0]==4||alt[0]==4){
    fprintf(stderr,"skipping site:%d (ref and alt must be A,C,G,T not N/n)\n",pos);
    return 0;

  }
  major=ref[0];
  minor=alt[0];
  char pick[2]={-1,-1};//pick[0]=GL,pick[1]=GP;
  int i=0;
  
  while(1){
    if(format[0]=='G'){
      if(format[1]=='L')
	pick[0] = i;
      if(format[1]=='P'){
	pick[1] = i;
      }
    }
    i++;
    format=strchr(format,':');
    if(format==NULL)
      break;
    format++;
  }

  if(pick[0]+pick[1]==-2){
    fprintf(stderr,"skipping site:%d (no GT/GP tag)\n",pos);
    return 0;
  }
  //  fprintf(stderr,"pick[0]:%d pick[1]:%d\n",pick[0],pick[1]);
  *lk = new double[nInd*10];
  *gp = new double[nInd*3];

  for(int i=0;1&&i<nInd*10;i++)
    lk[0][i]=-DBL_MAX;
  for(int i=0;i<3*nInd;i++)
    gp[0][i]=0;
 
  //now parse it
  int n=0;
  
  char *pi=NULL;
  i=0;
  while(((pi=angsd::strpop(&saveptr,'\t')))[0]!='\0') {//loop over inds
    char *tok=NULL;//loop over GT:PL:GP fields
    int p=0;
    while(((tok = angsd::strpop(&pi,':')))[0]!='\0') {
      if(p==pick[0]){
	//	fprintf(stderr,"GL tag not implemented properly yet:%s %lu\n",tok,strlen(tok));
	
	double tre[3];
	for(int j=0;j<3;j++){
	  tre[j] = atof(angsd::strpop(&tok,','));
	  //	  fprintf(stdout,"j[%d]:%f\n",j,tre[j]);
	}
	lk[0][i*10+angsd::majorminor[ref[0]][ref[0]]]=tre[0]/M_LOG10E;
	lk[0][i*10+angsd::majorminor[ref[0]][alt[0]]]=tre[1]/M_LOG10E;
	lk[0][i*10+angsd::majorminor[alt[0]][alt[0]]]=tre[2]/M_LOG10E;
		 //exit(0);

      }
      if(p==pick[1]){
	int pp=0;
	char *val=NULL;
	while(((val=angsd::strpop(&tok,',')))[0]!='\0'){
	  (*gp)[i*3+pp++]=atof(val);
	}
	assert(pp==3);
      }

      p++;
    }
    i++;
  }
 
  return pos;
}


funkyPars *vcfReader::fetch(int chunkSize){
  static int eof=0;
  if(eof)
    return NULL;
  funkyPars *r = allocFunkyPars();  
  r->likes=NULL;
  r->likes=new double*[chunkSize];
  r->post=new double*[chunkSize];
  r->posi=new int[chunkSize];
  r->major = new char[chunkSize];
  r->minor = new char[chunkSize];
  
  memset(r->major,0,chunkSize);
  memset(r->minor,0,chunkSize);
  r->refId = curChr;
  
  static int changed =0;
  if(changed){
    fprintf(stderr,"inchange\n");
    double *lk = new double [10*nInd];
    int p = parseline(r->likes,r->post,r->major[0],r->minor[0]);
    if(p>0)
      r->posi[0]=p-1;
    
    
  }

  int i=changed>0?1:0;
  changed=0;
  int cnt=i;//counter for the in-array position.

  for(;i<chunkSize;i++){
    if(Z_NULL==gzgets(gz,buf,len)){
      fprintf(stderr,"\t-> Done reading vcffile\n");
      eof=1;
      break;
    }
    buf[strlen(buf)-1] = '\0';
    saveptr=buf;
  
    // double *lk = new double [10*nInd];
    //    double *lk =NULL;
    //double *gp = new double [3*nInd];

    char *tmp = strchr(saveptr,'\t');
    tmp[0]='\0';
    char *chr = saveptr;
    saveptr=tmp+1;

    aMap::const_iterator it=revMap->find(chr);
    if(it==revMap->end()){
      fprintf(stderr,"Problem finding chromosome: %s in fai file\n",chr);
      exit(0);
    }
    if(curChr==-1){
      //init chr
      r->refId=curChr=it->second;
    }
    if(curChr!=it->second){
      changed =1;
      curChr=it->second;
      break;
    }
    int p = parseline(r->likes+cnt,r->post+cnt,r->major[cnt],r->minor[cnt]);
    if(p==0)
      continue;
    else{
      r->posi[cnt]=p-1;
      cnt++;//only update cnt
    }
  }
  
  r->nInd=nInd;
  r->numSites=cnt;
  //  fprintf(stderr,"\t-> r->nind:%d r->numSites:%d r->refId:%d\n",r->nInd,r->numSites,r->refId);
  assert(r->refId!=-1);
  return r;
}

