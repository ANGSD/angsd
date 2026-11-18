#include <climits>
#include <cstdarg>
#include "analysisFunction.h"
#include "glfReader_text.h"
#include "aio.h"

static const char *delims = "\t \r\n";

static void lineError(int lineNo,const char *fmt,...){
  fprintf(stderr,"\t-> Problem in text glf file at line %d: ",lineNo);
  va_list ap;
  va_start(ap,fmt);
  vfprintf(stderr,fmt,ap);
  va_end(ap);
  fprintf(stderr,". Will exit\n");
  exit(1);
}

glfReader_text::~glfReader_text(){
  free(buffer);
  free(buf);
}


glfReader_text::glfReader_text(int nInd_a,gzFile gz_a,const aMap *revMap_a){
  fprintf(stderr,"\t-> NB text glf files are still in the developing stage, please report strange behaviour to thorfinn@binf.ku.dk\n\t-> NB So far only tested with -dosaf -domaf -domajorminor\n");
  gz=gz_a;
  revMap=revMap_a;
  nInd = nInd_a;
  lineNo=0;
  changed=0;
  l=128;
  buffer =(char *) calloc(l,sizeof(char));
  bufSize=l;
  buf =(char *) calloc(bufSize,sizeof(char));
}


//parse the line in buffer. Returns 0 for an empty line.
int glfReader_text::parseLine(int *refId,int *posi,double *likes){
  //tokenize a copy, buffer is kept for when the line starts a new chromosome
  if(bufSize<l){
    buf=(char*)realloc(buf,l);
    bufSize=l;
  }
  strcpy(buf,buffer);

  fields.clear();
  char *save=NULL;
  for(char *tok=strtok_r(buf,delims,&save);tok;tok=strtok_r(NULL,delims,&save))
    fields.push_back(tok);
  if(fields.size()==0)
    return 0;

  int nLikes=10*nInd;
  char *chr=fields[0];
  char *pos=NULL;
  int first;//index of the first likelihood
  if((int)fields.size()==nLikes+2){//chr pos likes
    pos=fields[1];
    first=2;
  }else if((int)fields.size()==nLikes+1){//chr_pos likes; chr can contain '_'
    char *us=strrchr(chr,'_');
    if(us==NULL)
      lineError(lineNo,"first column \'%s\' is not chr_pos",chr);
    *us='\0';
    pos=us+1;
    first=1;
  }else{
    lineError(lineNo,"found %zu columns, expected chr, pos and 10 likelihoods for each of the %d individuals (%d columns), check -nInd",fields.size(),nInd,nLikes+2);
    return 0;//not reached
  }

  char *end=NULL;
  long p=strtol(pos,&end,10);
  if(end==pos||*end!='\0'||p<1||p>INT_MAX)
    lineError(lineNo,"invalid position \'%s\'",pos);
  aMap::const_iterator it = revMap->find(chr);
  if(it==revMap->end())
    lineError(lineNo,"chromosome \'%s\' is not in the fai file",chr);
  *refId=it->second;
  *posi=p-1;

  for(int i=0;i<nLikes;i++){
    char *tok=fields[first+i];
    likes[i]=strtod(tok,&end);
    if(end==tok||*end!='\0')
      lineError(lineNo,"likelihood \'%s\' is not a number",tok);
  }
  return 1;
}


funkyPars *glfReader_text::fetch(int chunksize){
  funkyPars * myfunky =funkyPars_init();
  myfunky->posi = new int[chunksize];
  double **likes = new double*[chunksize];
  for(int s=0;s<chunksize;s++)
    likes[s] = new double[nInd*10];

  int nSites=0;
  //buffer still holds a line from the previous call that started a new chromosome,
  //so parse it before reading more
  int haveLine = changed;
  changed=0;
  while(haveLine||aio::tgets(gz,&buffer,&l)){
    if(!haveLine)
      lineNo++;
    haveLine=0;

    int refId=0;
    int posi=0;
    if(!parseLine(&refId,&posi,likes[nSites]))
      continue;//empty line
    if(nSites>0&&refId!=myfunky->refId){
      //new chromosome: keep the line in buffer for the next call
      changed=1;
      break;
    }
    myfunky->refId=refId;
    myfunky->posi[nSites]=posi;
    nSites++;
    if(nSites>=chunksize)
      break;
  }

  for(int s=nSites;s<chunksize;s++)
    delete[] likes[s];
  myfunky->nInd=nInd;
  myfunky->likes=likes;
  myfunky->numSites = nSites;

  if(nSites==0){
    funkyPars_destroy(myfunky);
    return(NULL);
  }
  return(myfunky);
}
