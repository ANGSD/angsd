#include <cassert>
#include "analysisFunction.h"
#include "vcfReader.h"
#include <cfloat>

vcfReader::~vcfReader(){
  if(buf!=NULL)
    free(buf);
};

vcfReader::vcfReader(int &nInd_a,gzFile gz_a,const aMap *revMap_a){

  curChr=-1;
  nInd = nInd_a;
  gz=gz_a;
  len=128;
  buf=NULL;
  buf=original=saveptr=(char*)calloc(len,sizeof(char));
  revMap=revMap_a;
  while(aio::tgets(gz,&buf,&len)){
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
  fprintf(stderr,"\t-> BETA when using vcf-gl make sure that you use -doMajorMinor 1");
  original=saveptr=buf;
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
   fprintf(stderr,"skipping site:%d multiref/multialt (indel) ref:%s alt:%s\n",pos,ref,alt);
   return 0;
  }
  if(0&&strcmp(filter,"PASS")){
    fprintf(stderr,"skipping site:%d (non PASS)\n",pos);
    return 0;
    
  }
    
  if(strlen(ref)!=1||strlen(alt)!=1){
    fprintf(stderr,"skipping site:%d (indel)\n",pos);
    return 0;
  }

  ref[0]=refToInt[ref[0]];
  alt[0]=refToInt[alt[0]];
  if(ref[0]==4|| alt[0]==4){
    fprintf(stderr,"REF/ALT is 'N' will discard site: %d \n",pos);
    return 0;
  }
  if(ref[0]==4||alt[0]==4){
    fprintf(stderr,"skipping site:%d (ref and alt must be A,C,G,T not N/n)\n",pos);
    return 0;

  }

  major=ref[0];
  minor=alt[0];
  char pick[3]={-1,-1,-1};//pick[0]=GL,pick[1]=GP,pick[2]=PL;
  int i=0;

  while(1){
    if(format[0]=='G'){
      if(format[1]=='L')
	pick[0] = i;
      if(format[1]=='P'){
	pick[1] = i;
      }
    }
    if(format[0]=='P'&&format[1]=='L')
	pick[2] = i;

    i++;
    format=strchr(format,':');
    if(format==NULL)
      break;
    format++;
  }
  //  fprintf(stderr,"pick[0]:%d pick[1]:%d pick[2]:%d\n",pick[0],pick[1],pick[2]);
  if(pick[0]+pick[1]+pick[2]==-3){
    fprintf(stderr,"skipping site:%d (no GT/GP/PL tag)\n",pos);
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
      }if(p==pick[2]){
	//PL tag
	
	double tre[3];
	for(int j=0;j<3;j++){
	  tre[j] = log(pow(10,-atof(angsd::strpop(&tok,','))/10.0));//hack
	  // fprintf(stdout,"j[%d]:%f\n",j,tre[j]);
	}
	lk[0][i*10+angsd::majorminor[ref[0]][ref[0]]]=tre[0]/M_LOG10E;
	lk[0][i*10+angsd::majorminor[ref[0]][alt[0]]]=tre[1]/M_LOG10E;
	lk[0][i*10+angsd::majorminor[alt[0]][alt[0]]]=tre[2]/M_LOG10E;


      }

      p++;
    }
    i++;
  }
 
  return pos;
}


funkyPars *vcfReader::fetch(int chunkSize){
  //  fprintf(stderr,"fetch:%d curChr:%d\n\n",chunkSize,curChr);
  static int eof=0;
  if(eof)
    return NULL;
  funkyPars *r = allocFunkyPars();  
  //  r->likes=r->post=NULL;
  r->likes=new double*[chunkSize];
  r->post=new double*[chunkSize];
  for(int i=0;i<chunkSize;i++){
    memset(r->likes,0,chunkSize*sizeof(double*));
    memset(r->post,0,chunkSize*sizeof(double*));
  }
  r->posi=new int[chunkSize];
  r->major = new char[chunkSize];
  r->minor = new char[chunkSize];
  
  memset(r->major,0,chunkSize);
  memset(r->minor,0,chunkSize);
  //  fprintf(stderr,"curChr:%d\n",curChr);
  r->refId = curChr;
  
  static int changed =0;
 reread:
  int i=0;
  if(changed){
    //    fprintf(stderr,"inchange\n");
    int p = parseline(r->likes,r->post,r->major[0],r->minor[0]);
    if(p>0)
      r->posi[i++]=p-1;
  }

  changed=0;
  int cnt=i;//counter for the in-array position.
  buf=original;
  for(;i<chunkSize;i++){
    if(0==aio::tgets(gz,&buf,&len)){
      fprintf(stderr,"\t-> Done reading vcffile\n");fflush(stderr);
      eof=1;
      break;
    }
    if(buf!=original)
      original=buf;
    //    fprintf(stderr,"in for loop2\n");
    buf[strlen(buf)-1] = '\0';
    saveptr=buf;
  
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
      //      fprintf(stderr,"breaking for change of chr\n");
      break;
    }
    //  fprintf(stderr,"in for loop3\n");
    int p = parseline(r->likes+cnt,r->post+cnt,r->major[cnt],r->minor[cnt]);
    //  fprintf(stderr,"in for loop4\n");
    buf=original;
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
  if(r->numSites==0&&eof!=1)
    goto reread;
  //  fprintf(stderr,"\t-> r->nind:%d r->numSites:%d r->refId:%d\n",r->nInd,r->numSites,r->refId);
  return r;
}

