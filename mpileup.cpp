#include <ctype.h>
#include <climits>
#include <cstdarg>
#include "pooled_alloc.h"
#include "mpileup.h"
#include "mUpPile.h"
#include "aio.h"

extern tpool_alloc_t *tnodes;

static void lineError(int lineNo,const char *fmt,...){
  fprintf(stderr,"\t-> Problem in pileup file at line %d: ",lineNo);
  va_list ap;
  va_start(ap,fmt);
  vfprintf(stderr,fmt,ap);
  va_end(ap);
  fprintf(stderr,". Will exit\n");
  exit(1);
}

mpileup::mpileup(int nInd_a,gzFile gz_a,const aMap* revMap_a,int minQ_a){
  fprintf(stderr,"\t-> You are using -pileup, this means:\n");
  fprintf(stderr,"\t-> 1) Internal positions both from front and back is coded to 255\n");
  fprintf(stderr,"\t-> 2) All mapping qualities (mapQ) are set to 30\n");
  fprintf(stderr,"\t-> 3) Program will not represent insertions, use raw BAM/CRAM for that\n");

  nInd = nInd_a;
  gz = gz_a;
  revMap = revMap_a;
  minQ = minQ_a;
  lineNo = 0;
  changed = 0;
  nIndelWarn = 0;
  l = 128;
  buffer = (char*) calloc(l,sizeof(char));
  bufSize = l;
  buf = (char*) calloc(bufSize,sizeof(char));
}

mpileup::~mpileup(){
  free(buffer);
  free(buf);
}

//split a copy of buffer into fields. Fields are separated by tabs (empty fields are kept),
//or by spaces if the line has no tabs. Returns the number of fields, 0 for an empty line.
int mpileup::splitLine(){
  if(bufSize<l){
    buf=(char*)realloc(buf,l);
    bufSize=l;
  }
  strcpy(buf,buffer);
  size_t len=strlen(buf);
  while(len>0&&(buf[len-1]=='\n'||buf[len-1]=='\r'))
    buf[--len]='\0';

  fields.clear();
  if(strchr(buf,'\t')){
    char *p=buf;
    fields.push_back(p);
    while((p=strchr(p,'\t'))){
      *p++='\0';
      fields.push_back(p);
    }
  }else{
    char *save=NULL;
    for(char *tok=strtok_r(buf," ",&save);tok;tok=strtok_r(NULL," ",&save))
      fields.push_back(tok);
  }
  if(fields.size()==1&&fields[0][0]=='\0')
    fields.clear();
  return fields.size();
}

//parse depth, bases and quals for each individual, fields start at index 3
tNode **mpileup::parseNd(char ref){
  tNode **ret =(tNode**) calloc(nInd,sizeof(tNode*));

  for(int i=0;i<nInd;i++) {
    char *depthStr = fields[3+3*i];
    char *bases = fields[4+3*i];
    char *quals = fields[5+3*i];
    char *end=NULL;
    long depth = strtol(depthStr,&end,10);
    if(end==depthStr||*end!='\0'||depth<0||depth>INT_MAX)
      lineError(lineNo,"invalid depth \'%s\' for individual %d",depthStr,i);
    ret[i] = initNodeT(depth);
    ret[i]->l = depth;
    if(depth==0)//bases and quals are '*' or empty
      continue;

    //one entry in bases (and one quality) per read
    int at=0;
    for(char *b=bases;*b;){
      if(*b=='$'){//end of read
	b++;
	continue;
      }
      if(*b=='^'){//start of read, followed by its mapping quality
	if(b[1]=='\0')
	  lineError(lineNo,"'^' at the end of the bases for individual %d",i);
	b+=2;
	continue;
      }
      if(*b=='+'||*b=='-'){//indel: +3ACG or -3ACG, belongs to the previous base
	if(nIndelWarn++==0){
	  fprintf(stderr,"\t-> Problem with indels in -pileup: please use raw binary for bam files and NOT pileup files for inserts:%s\n",bases);
	  fprintf(stderr,"\t-> Parsing will only use non-indel information. This msg will appear only once\n");
	}
	char *e=NULL;
	long opLen=strtol(b+1,&e,10);
	if(e==b+1||opLen<0||(long)strlen(e)<opLen)
	  lineError(lineNo,"malformed indel in the bases for individual %d",i);
	b=e+opLen;
	continue;
      }
      char base;
      if(isalpha(*b))
	base=*b;
      else if(*b=='.'||*b==',')
	base=ref;
      else if(*b=='*'||*b=='#'||*b=='>'||*b=='<')//deletion or reference skip
	base='N';
      else{
	lineError(lineNo,"unexpected character \'%c\' in the bases for individual %d",*b,i);
	return NULL;//not reached
      }
      if(at>=depth)
	lineError(lineNo,"more bases than the depth %ld for individual %d",depth,i);
      ret[i]->seq[at++]=base;
      b++;
    }
    if(at!=depth)
      lineError(lineNo,"found %d bases, expected the depth %ld for individual %d",at,depth,i);
    if((long)strlen(quals)!=depth)
      lineError(lineNo,"found %zu base qualities, expected the depth %ld for individual %d",strlen(quals),depth,i);
    for(int j=0;j<depth;j++){
      ret[i]->qs[j] = quals[j]-33;
      ret[i]->posi[j] = 255;
      ret[i]->isop[j] = 255;
      ret[i]->mapQ[j] = 30;
      if(ret[i]->qs[j]<minQ)
	ret[i]->seq[j] = 'n';
    }
  }
  return ret;
}

funkyPars *mpileup::fetch(int chunksize){
  funkyPars * myfunky =funkyPars_init();
  myfunky->posi = new int[chunksize];
  myfunky->ref = new char[chunksize];
  myfunky->chk = new chunkyT;
  myfunky->chk->nd = new tNode **[chunksize];
  myfunky->chk->refPos = NULL;
  int nSites=0;

  //buffer still holds a line from the previous call that started a new chromosome,
  //so parse it before reading more
  int haveLine = changed;
  changed=0;
  while(haveLine||aio::tgets(gz,&buffer,&l)){
    if(!haveLine)
      lineNo++;
    haveLine=0;

    int nFields = splitLine();
    if(nFields==0)
      continue;//empty line
    if(nFields!=3+3*nInd)
      lineError(lineNo,"found %d columns, expected chr, pos, ref and depth, bases, quals for each of the %d individuals (%d columns), check -nInd",nFields,nInd,3+3*nInd);

    aMap::const_iterator it=revMap->find(fields[0]);
    if(it==revMap->end())
      lineError(lineNo,"chromosome \'%s\' is not in the fai file",fields[0]);
    if(nSites>0&&it->second!=myfunky->refId){
      //new chromosome: keep the line in buffer for the next call
      changed=1;
      break;
    }
    myfunky->refId = it->second;

    char *end=NULL;
    long pos=strtol(fields[1],&end,10);
    if(end==fields[1]||*end!='\0'||pos<1||pos>INT_MAX)
      lineError(lineNo,"invalid position \'%s\'",fields[1]);
    if(fields[2][0]=='\0')
      lineError(lineNo,"empty reference base");
    char ref = fields[2][0];

    myfunky->posi[nSites] = pos-1;
    myfunky->ref[nSites] = refToInt[(unsigned char)ref];
    myfunky->chk->nd[nSites++] = parseNd(ref);
    if(nSites>=chunksize)
      break;
  }
  myfunky->nInd=myfunky->chk->nSamples= nInd;
  myfunky->numSites = myfunky->chk->nSites=nSites;
  if(nSites==0){
    funkyPars_destroy(myfunky);
    return(NULL);
  }
  return(myfunky);
}
