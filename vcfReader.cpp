#include <cassert>
#include "vcfReader.h"


vcfReader::~vcfReader(){

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
    if(!strncmp(buf,"#CHROM",5)) {
      int i=1;
      tmp=buf;
      while(((tmp=strchr(tmp,'\t')))){
	i++;
	tmp++;
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
int vcfReader::parseline(double *lk,double *gp,char &major,char &minor){
  tmp = strchr(saveptr,'\t');
  tmp[0] = '\0';
  int pos = atoi(saveptr);
  saveptr = tmp+1;
  saveptr = strchr(saveptr,'\t')+1;//ID
  char *ref = saveptr;
  tmp = strchr(saveptr,'\t'); tmp[0]='\0';
  char *alt = saveptr = tmp+1;
  tmp =  strchr(saveptr,'\t'); tmp[0]='\0';
  char *qual = saveptr=tmp+1;
  tmp =  strchr(saveptr,'\t'); tmp[0]='\0';
  char *filter =saveptr=tmp+1;
  tmp =  strchr(saveptr,'\t'); tmp[0]='\0';
  char *info =saveptr=tmp+1;
  tmp =  strchr(saveptr,'\t'); tmp[0]='\0';
  char *format= saveptr=tmp+1;
  tmp =  strchr(saveptr,'\t'); tmp[0]='\0';
  //fprintf(stderr,"\npos:%d ref:%s alt:%s qual:%s flt:%s info:%s format:%s strlen(ref):%zu strlen(alt):%zu\n",pos,ref,alt,qual,filter,info,format,strlen(ref),strlen(alt));
  saveptr=tmp+1;
  if(strlen(ref)!=1||strlen(alt)!=1){
    fprintf(stderr,"skipping site:%d (indel)\n",pos);
    return 0;
  }
  if(strcmp(filter,"PASS")){
    fprintf(stderr,"skipping site:%d (non PASS)\n",pos);
    return 0;
    
  }
    
  if(strlen(ref)!=1||strlen(alt)!=1){
    fprintf(stderr,"skipping site:%d (indel)\n",pos);
    return 0;
  }
  ref[0]=refToInt[ref[0]];
  alt[0]=refToInt[alt[0]];
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
  fprintf(stderr,"pick[0]:%d pick[1]:%d\n",pick[0],pick[1]);
  for(int i=0;i<nInd*10;i++)
    lk[i]=-0.0;
  for(int i=0;i<3*nInd;i++)
    gp[i]=0;
 
  //now parse it
  //  fprintf(stderr,"\nsaveptr:%s\n",saveptr);exit(0);
  for(int i=0;i<nInd;i++){
    int n=0;
    tmp=saveptr;
    saveptr=strchr(saveptr,'\t');
    //saveptr now contains the head to the tail.
    for(int n=0;n<std::max(pick[0],pick[1]);n++){
      
      if(n==pick[0])//gls
	fprintf(stderr,"Need to check GL tag filereading\n");
	snprintf(tmp,strlen(tmp),"%f,%f,%f",				\
		 lk[i*10+angsd::majorminor[ref[0]][ref[0]]]/M_LOG10E,	\
		 lk[i*10+angsd::majorminor[ref[0]][alt[0]]]/M_LOG10E,	\
		 lk[i*10+angsd::majorminor[alt[0]][alt[0]]]/M_LOG10E);
      
      if(n==pick[1]){//gps
	fprintf(stderr,"hit\n");
	int jj=0;
	for(char *ttt = strtok(NULL,",");ttt!=NULL;ttt=strtok(NULL,"")){
	  fprintf(stderr,"\t-> ttt:%s ii:%d \n",ttt,i*3+jj);
	  gp[i*3+jj] = atof(ttt);
	  
	}
	//	fprintf(stderr,"%f %f %f\n",gp[i*3],gp[i*3+1],gp[i*3+2]);
      }
    }
    exit(0);
  }
  
  return pos;
}


funkyPars *vcfReader::fetch(int chunkSize){
  static int eof=0;
  if(eof)
    return NULL;
  fprintf(stderr,"CAlling fecth:%d len:%d\n",chunkSize,len);
  funkyPars *r = allocFunkyPars();  

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
    double *gp = new double [3*nInd];
    int p = parseline(lk,gp,r->major[0],r->minor[0]);
    if(p==0){
      delete [] lk;
      delete [] gp;
    }else{
      r->posi[0]=p-1;
      r->post[0]=gp;
      r->likes[0]=lk;
    }
    
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
  
    double *lk = new double [10*nInd];
    double *gp = new double [3*nInd];

    tmp = strchr(saveptr,'\t');
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
    int p = parseline(lk,gp,r->major[cnt],r->minor[cnt]);
    if(p==0)
      continue;
    else{
      r->posi[cnt]=p-1;
      r->likes[cnt] = lk;
      r->post[cnt] = gp;
      cnt++;//only update cnt
    }
  }
  
  r->nInd=nInd;
  r->numSites=cnt;
  fprintf(stderr,"r->nind:%d r->numSites:%d\n",r->nInd,r->numSites);
  
  return r;
}

