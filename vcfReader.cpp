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
  while(gzgets(gz,buf,len)){
    if(!strncmp(buf,"#CHROM",5)){
      int i=0;
      char *tmp=NULL;
      for(tmp=strtok(buf,"\t\n ");tmp!=NULL;tmp=strtok(NULL,"\n\t "))
	i++;
      if(i-9!=nInd){
	fprintf(stderr,"\t-> Looks like wrong -nInd has supplied compared to header of vcf file: %d vs %d\n",nInd,i-9);
	exit(0);
      }
    }
    

  }
  revMap=revMap_a;
}

//return value is pos. If this is zero means dont use site
int vcfReader::parseline(double *lk,double *gp,char &major,char &minor){
  int pos = atoi(strtok_r(NULL,"\n\t ",&saveptr));
  strtok_r(NULL,"\n\t ",&saveptr);//ID
  char *ref = strtok_r(NULL,"\n\t ",&saveptr);//REF
  char *alt = strtok_r(NULL,"\n\t ",&saveptr);//ALT
  strtok_r(NULL,"\n\t ",&saveptr);//QUAL
  char *filter = strtok_r(NULL,"\n\t ",&saveptr);//FILTER
  strtok_r(NULL,"\n\t ",&saveptr);//INFO
  char *format= strtok_r(NULL,"\n\t ",&saveptr);//FORMAT
  fprintf(stderr,"pos:%d ref:%s alt:%s flt:%s format:%s strlen(ref):%zu strlen(alt):%zu\n",pos,ref,alt,filter,format,strlen(ref),strlen(alt));
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
  char pick[2]={-1,-1};//pick[0]=GL,pick[1]=GP;

  int i=0;
  char *tmp=NULL;
  for(tmp=strtok(format,":");tmp!=NULL;tmp=strtok(NULL,":")){
    fprintf(stderr,"tmp:%s\n",tmp);
    if(strcmp(tmp,"GL"))
      pick[0]=i;
    else if(!strcmp(tmp,"GP"))
      pick[1]=i;
    i++;
  }
  fprintf(stderr,"pick[0]:%d pick[1]:%d\n",pick[0],pick[1]);
  if(pick[0]+pick[1]==-2){
    fprintf(stderr,"skipping site:%d (no GT/GP tag)\n",pos);
    return 0;
  }
  for(int i=0;i<nInd*10;i++)
    lk[i]=-0.0;
  for(int i=0;i<3*nInd;i++)
    gp[i]=0;
  
  //now parse it
  for(int i=0;i<nInd;i++){
    char *tmp = strtok_r(NULL,"\n\t ",&saveptr);//contains something like: \'0|0:0:1,0,0\'
    int n=0;
    for(tmp=strtok(tmp,":");tmp!=NULL;tmp=strtok(NULL,":")){
      if(n==pick[0])//gls
	snprintf(tmp,strlen(tmp),"%f,%f,%f",				\
		 lk[i*10+angsd::majorminor[ref[0]][ref[0]]]/M_LOG10E,	\
		 lk[i*10+angsd::majorminor[ref[0]][alt[0]]]/M_LOG10E,	\
		 lk[i*10+angsd::majorminor[alt[0]][alt[0]]]/M_LOG10E);
      
      if(n==pick[1])//gps
	snprintf(tmp,strlen(tmp),"%f,%f,%f",gp[i*3],gp[i*3+1],gp[i*3+2]);
	      
      n++;
    }
  }
  
  return pos;
}


funkyPars *vcfReader::fetch(int chunkSize){
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
    double *lk = new double [10*nInd];
    double *gp = new double [3*nInd];
    int p = parseline(lk,gp,r->major[0],r->minor[0]);
    if(p==0){
      delete [] lk;
      delete [] gp;
    }else{
      r->posi[0]=p;
      r->post[0]=gp;
      r->likes[0]=lk;
    }
    
  }

  int i=changed>0?1:0;
  changed=0;
  int cnt=i;//counter for the in-array position.
  for(;i<chunkSize;i++){
    if(gzgets(gz,buf,len)==Z_NULL)
      break;
    
    double *lk = new double [10*nInd];
    double *gp = new double [3*nInd];

    char *chr = strtok_r(buf,"\n\t ",&saveptr);

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
    int p = parseline(lk,gp,r->major[i],r->minor[i]);
    if(p==0)
      continue;
    else{
      r->posi[cnt]=p;
      //      fprintf(stderr,"posi[i]:%d\n",r->posi[cnt]);
      r->likes[cnt] = lk;
      r->post[cnt] = gp;
      cnt++;//only update cnt
    }
  }
  if(cnt==0){
    //cleanup
    return NULL;
  }
  
  r->nInd=nInd;
  
  r->numSites=i;
 
  
  return r;
}

