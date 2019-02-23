#include <vector>
#include <cstring>
#include <cstdlib>
#include <ctype.h>
#include <pthread.h>
#include <cassert>
#include <htslib/hts.h>
#include "mUpPile.h"
#include "abcGetFasta.h"
#include "analysisFunction.h"
#include "makeReadPool.h"
#include "pop1_read.h"
#define USE_MALLOC_WRAPPERS
#include "malloc_wrap.h"


#include "pooled_alloc.h"
void *tail=NULL;//<- this will be point to the pool->free
void *head=NULL;//<- this will be the last node adjoint
tpool_alloc_t *tnodes = NULL;
size_t currentnodes=0;
size_t freenodes =0;
pthread_mutex_t slist_mutex = PTHREAD_MUTEX_INITIALIZER;
size_t when_to_flush = 5000;


pthread_mutex_t mUpPile_mutex = PTHREAD_MUTEX_INITIALIZER;

//slist *sl = NULL;

extern int SIG_COND;
extern int minQ;
extern int trim;
extern int trim5;
extern int trim3;
#define bam_nt16_rev_table seq_nt16_str
#define __NEW__
//#define MAX_SEQ_LEN 200 //this is used for getting some secure value for when a site is securely covered for sample
/*
  When merging different nodes, we are using a greedy but fast approach,
  however if we have a large genomic region with no data we are are allocating huge chunks.
  So the cutoff below will make a fallback to a slower but memsecure version
 */
#define BUG_THRES 1000000 //<- we allow 1mio sites to be allocated,
#define UPPILE_LEN 4

extern int MAX_SEQ_LEN;

/*
  node for uppile for a single individual for a single site
  this is used for textoutput emulating samtools mpileup
 */

typedef struct{
  int len;//length of seq,qs,pos
  int maxLen;//maxlength of seq,qs.pos, maybe do realloc
  kstring_t seq;
  kstring_t qs;
  kstring_t pos;
  int depth;
  int refPos;
}node;

template<typename T>
T getmax(const T *ary,size_t len){
  assert(len>0&&ary!=NULL);
      
  T high = ary[0];
  for(size_t i=1;i<len;i++){
    if(ary[i]>high)
      high=ary[i];
  }
  return high;
}

void dalloc_node(node &n){
  free(n.seq.s);
  free(n.qs.s);
  free(n.pos.s);
}

void tnode_destroy1(tNode *n){
  free(n->posi);
  free(n->isop);
  free(n->mapQ);
  free(n->seq);
  free(n->qs);
}

//this will only be called from a threadsafe context
void tnode_destroy(tNode *n){
  if(n==NULL)
    return;

  if(tail==NULL&&head==NULL)
    head = tail = n;
  else{
    *(void **)n = head;
    head=n;
  }
  freenodes++;
  /*
  tnode_destroy1(n);
  free(n);
  n=NULL;
  */
}

typedef struct{
  char *refName;
  int regStart;
  int regStop;
  int nSites;
  int nSamples;
  node **nd;//nd[site][ind]
  int *refPos;//length is nSites
  int first;
  int last;
  int start;
  int length;
  int refId;
}chunky;



void printChunky2(const chunky* chk,FILE *fp,char *refStr,abcGetFasta *gf) {
  //  fprintf(stderr,"[%s] nsites=%d region=(%d,%d) itrReg=(%d,%d)\n",__FUNCTION__,chk->nSites,chk->refPos[0],chk->refPos[chk->nSites-1],chk->regStart,chk->regStop);
  if(chk->refPos[0]>chk->regStop){
    fprintf(stderr,"\t->Problems with stuff\n");
    exit(0);
  }
  int refId = chk->refId;
  for(int s=0;s<chk->nSites;s++) {
    if(chk->refPos[s]<chk->regStart || chk->refPos[s]>chk->regStop-1 ){
      for(int i=0;i<chk->nSamples;i++)
	dalloc_node(chk->nd[s][i]);
      delete [] chk->nd[s];
      continue;
    }
    fprintf(fp,"%s\t%d",refStr,chk->refPos[s]+1);     
    if(gf->ref!=NULL){
      if(refId!=gf->ref->curChr)
	gf->loadChr(gf->ref,refStr,refId);
      if(gf->ref->seqs!=NULL){
	fprintf(fp,"\t%c",gf->ref->seqs[chk->refPos[s]]);
      }
    }
    for(int i=0;i<chk->nSamples;i++) {

      //      fprintf(stderr,"seqlen[%d,%d]=%lu\t",s,i,chk->nd[s][i].seq->l);
      if(chk->nd[s][i].seq.l!=0){
	fprintf(fp,"\t%d\t",chk->nd[s][i].depth);
	for(size_t l=0;l<chk->nd[s][i].seq.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].seq.s[l]);
	fprintf(fp,"\t");
		
	for(size_t l=0;l<chk->nd[s][i].qs.l;l++)
	  fprintf(fp,"%c",chk->nd[s][i].qs.s[l]);
	//	fprintf(fp,"\t");

      }else
	fprintf(fp,"\t0\t*\t*");
      dalloc_node(chk->nd[s][i]);
    }
    //    fprintf(stderr,"\n");
    fprintf(fp,"\n");
    delete [] chk->nd[s];
  }
  delete [] chk->nd;
  delete [] chk->refPos;
  delete chk;
} 


typedef struct{
  int l;//number of nodes in nodes
  int m; //possible number of of nodes
  int first;//simply a value which is equal to nodes[0].refPos;
  int last;
  node *nds;//this length can maximum be the maxlenght of a read.NOTANYMORE

}nodePool;



nodePool allocNodePool(int l){
  nodePool np;
  np.l=0;
  np.m = l;
  kroundup32(np.m);

  np.first=np.last=-1;
  np.nds = new node[np.m];
  
  return np;
}

nodePoolT allocNodePoolT(int l){
  nodePoolT np;
  np.l=0;
  np.m = l;
  kroundup32(np.m);

  np.first=np.last=-1;
  np.nds =(tNode**) calloc(np.m,sizeof(tNode*));
  
  return np;
}



void dalloc_nodePool (nodePool& np){
  for(int i=0;i<np.l;i++)
    dalloc_node(np.nds[i]);
  
}

void cleanUpChunkyT(chunkyT *chk){
  pthread_mutex_lock(&slist_mutex);//just to make sure, its okey to clean up
  for(int s=0;s<chk->nSites;s++) {
    for(int i=0;i<chk->nSamples;i++) {
      if(chk->nd[s][i]==NULL)
	continue;
      if(chk->nd[s][i]&&chk->nd[s][i]->l2!=0){
	for(int j=0;j<chk->nd[s][i]->l2;j++){
	  tnode_destroy(chk->nd[s][i]->insert[j]);
	}
	free(chk->nd[s][i]->insert);

      }
      tnode_destroy(chk->nd[s][i]);
    }
    free(chk->nd[s]);
  }
  delete [] chk->refPos;
  delete [] chk->nd;
  delete chk;
  pthread_mutex_unlock(&slist_mutex);//just to make sure, its okey to clean up
}




node initNode(int l,int refPos,int arypos){
  node nd;
  nd.maxLen = l;
  nd.seq.s=NULL; nd.seq.l=nd.seq.m=0;
  nd.qs.s =NULL; nd.qs.l=nd.qs.m=0;
  nd.pos.s=NULL; nd.pos.l=nd.pos.m=0;
  nd.refPos = refPos;//not used anymore
  nd.depth = 0;
  return nd;
}

void tnode_realloc(tNode *d,int newsize){

  kroundup32(newsize);
  if(newsize<=d->m)
    fprintf(stderr,"[%s] problems newsize should be bigger than oldesize\n",__FUNCTION__);
  d->m=newsize;

  d->seq=(char*)realloc(d->seq,d->m*sizeof(char));
  d->qs =(char*) realloc(d->qs,d->m*sizeof(char));
  d->posi = (suint*)realloc(d->posi,d->m*sizeof(suint));
  d->isop = (suint*)realloc(d->isop,d->m*sizeof(suint));
  d->mapQ = (unsigned char*)realloc(d->mapQ,d->m*sizeof(unsigned char));
 
}


//this is called from a threadsafe context

void flush_queue(){
  //  fprintf(stderr,"[%s]IN currentnodes:%lu freenodes:%lu when_to_flush:%lu\n",__FUNCTION__,currentnodes,freenodes,when_to_flush);
  pthread_mutex_lock(&slist_mutex);
  if(freenodes>0){
    if(tail!=NULL&&head!=NULL){
      //  fprintf(stderr,"PRE: infree:%lu\n",tpool_infree(tnodes));
      *(void **)tail = tnodes->free;
      tnodes->free = head;
      head=tail=NULL;
      //fprintf(stderr,"POST: infree:%lu\n",tpool_infree(tnodes));
    }
    currentnodes -=freenodes;
    freenodes=0;
  }
  pthread_mutex_unlock(&slist_mutex);
  // fprintf(stderr,"[%s]OUT currentnodes:%lu freenodes:%lu when_to_flush:%lu\n",__FUNCTION__,currentnodes,freenodes,when_to_flush);
}

//int staaa =0;
tNode *initNodeT(int l){
  if(l<UPPILE_LEN)
    l=UPPILE_LEN;
  tNode *d=(tNode*)tpool_alloc(tnodes);
  //  tNode *d=(tNode*)calloc(1,sizeof(tNode));


  if(d->m==0){
    d->m = l;
    kroundup32(d->m);
    d->seq=(char *)malloc(d->m);
    d->qs=(char *)malloc(d->m);
    d->posi=(suint *)malloc(sizeof(suint)*d->m);
    d->isop=(suint *)malloc(sizeof(suint)*d->m);
    d->mapQ=(unsigned char *)malloc(d->m);
    d->m2=0;
  }else if(l>d->m){
    tnode_realloc(d,l);
  }
  d->l = d->l2 =d->m2 =0;
  //d->l = d->l2 = d->m2 = 0;
  d->refPos= -999;
  d->insert = NULL;
  currentnodes++;
  if(currentnodes>when_to_flush &&(currentnodes % 5000)==0)
    flush_queue();
  return d;
}


/*
  What does this do?
  returns the number of basepairs covered by np and sgl
*/

int coverage_in_bp(nodePool *np, readPool *sgl){
#if 0
  for(int i=0;0&&i<sgl->readIDstop;i++)
    fprintf(stderr,"sglrange[%d]\t %d\t%d\n",i,sgl->first[i],sgl->last[i]);
#endif
  int sumReg =0;
  int start=0;
  int first = np->first;
  int last = np->last +1;//because we are comparing with calc_end, which is noncontain

  if(np->l==0&&sgl->readIDstop!=0){//LAST PART OF CONDTIONAL WAS A MINDBUGGING BUG

    first = sgl->first[0];
    last = sgl->last[0];
    start++;
  }
  if(np->l==0&&sgl->readIDstop==0)//JESUS CHRIST
    return 0;
  for(int i=start;i<sgl->readIDstop;i++) {
    if(sgl->first[i]>last){
      sumReg += last-first;
      first =sgl->first[i];
    }
    last = std::max(sgl->last[i],last);
  }
  sumReg += last-first;
  return sumReg;
}

int coverage_in_bpT(nodePoolT *np, readPool *sgl){
#if 0
  for(int i=0;0&&i<sgl->readIDstop;i++)
    fprintf(stderr,"sglrange[%d]\t %d\t%d\n",i,sgl->first[i],sgl->last[i]);
#endif
  int sumReg =0;
  int start=0;
  int first = np->first;
  int last = np->last +1;//because we are comparing with calc_end, which is noncontain
  
  if(np->l==0&&sgl->readIDstop!=0){//LAST PART OF CONDTIONAL WAS A MINDBUGGING BUG

    first = sgl->first[0];
    last = sgl->last[0];
    start++;
  }
  if(np->l==0&&sgl->readIDstop==0)//JESUS CHRIST
    return 0;
  for(int i=start;i<sgl->readIDstop;i++) {
    if(sgl->first[i]>last){
      sumReg += last-first;
      first =sgl->first[i];
    }
    last = std::max(sgl->last[i],last);
  }
  sumReg += last-first;
  return sumReg;
}

nodePool mkNodes_one_sampleb(readPool *sgl,nodePool *np,abcGetFasta *gf) {
  int regionLen = coverage_in_bp(np,sgl);//true covered regions
  nodePool dn;
  dn.l =0;
  if(regionLen==0)
    return dn;
  int lastSecure = sgl->lowestStart;
  
  dn = allocNodePool(regionLen);
  nodePool *ret = &dn;

  node *nds = ret->nds;


  int offs = np->first;//first position from buffered
  int last = np->last+1;//because we are comparing with calc_end

  //plug in the old buffered nodes
  for(int i=0;i<np->l;i++)
    nds[np->nds[i].refPos-offs] = np->nds[i];

  //initialize new nodes for the rest of the region
  for(int i=np->l;i<regionLen;i++)
    nds[i] = initNode(UPPILE_LEN,-1,-1);

  if(np->l==0){//if we didn't have have anybuffered then
    offs = sgl->first[0];
    last = sgl->last[0];
  }

  //parse all reads
  int r;
  for( r=0;r<sgl->readIDstop;r++) {
    bam1_t *rd=sgl->reads[r];

    //    fprintf(stderr,"r=%d\tpos=%d\n",r,rd.pos);
    if(sgl->first[r] > last){
      int diffs = (rd->core.pos-last);
      offs = offs + diffs;
      last = sgl->last[r];
    }else
      last = std::max(sgl->last[r],last);

    char *seq =(char *) bam_get_seq(rd);
    char *quals =(char *) bam_get_qual(rd);
    int nCig = rd->core.n_cigar;

    uint32_t *cigs = bam_get_cigar(rd);
    int seq_pos =0; //position within sequence
    int wpos = rd->core.pos-offs;//this value is the current position assocatied with the positions at seq_pos
    node *tmpNode =NULL; //this is a pointer to the last node beeing modified
    int hasInfo =0;//used when first part is a insertion or deletion
    int hasPrintedMaq =0;
    //loop through read by looping through the cigar string
    for(int i=0;i<nCig;i++) {
      int opCode = cigs[i]&BAM_CIGAR_MASK; //what to do
      int opLen = cigs[i]>>BAM_CIGAR_SHIFT; //length of what to do

      if(opCode==BAM_CINS||opCode==BAM_CDEL){//handle insertions and deletions
	if(i==0){ //skip indels if beginning of a read, print mapQ
	  tmpNode = &nds[wpos];
	  kputc('^', &tmpNode->seq);
	  if(rd->core.qual!=255)
	    kputc(rd->core.qual+33, &tmpNode->seq);
	  else
	    kputc('~', &tmpNode->seq);
	  hasPrintedMaq =1;
	  if(opCode==BAM_CINS){
	    seq_pos += opLen;
	    continue;
	  } 
	}
	if(opCode==BAM_CINS&&hasInfo==0){//when should this happen?
	  seq_pos += opLen;	  
	  continue;
	}
	if(i!=0){
	  wpos--; //insertion/deletion is bound to the last position of the read
	  tmpNode = &nds[wpos];
	  kputc(opCode&BAM_CINS?'+':'-',&tmpNode->seq);
	  kputw(opLen,&tmpNode->seq);
	  hasInfo++;
	}
	if(opCode==BAM_CINS){
	  for(int ii=0;ii<opLen;ii++){
	    char c = bam_nt16_rev_table[bam_seqi(seq, seq_pos)];
	    kputc(bam_is_rev(rd)? tolower(c) : toupper(c), &tmpNode->seq);
	    kputw(seq_pos+1,&tmpNode->pos);
	    kputc(',',&tmpNode->pos);
	    seq_pos++;
	  }
	  wpos++;
	}else {//this is the deletion part
	  if(i!=0){
	    if(gf->ref==NULL)
	      for(int ii=0;ii<opLen;ii++)
		kputc(bam_is_rev(rd)? tolower('N') : toupper('N'),&tmpNode->seq);
	    else
	      for(int ii=0;ii<opLen;ii++)
		kputc(bam_is_rev(rd)? tolower(gf->ref->seqs[offs+wpos+ii+1]) : toupper(gf->ref->seqs[offs+wpos+ii+1]),&tmpNode->seq);
	    wpos++;//write '*' from the next position and opLen more
	  }
	  for(int fix=wpos;wpos<fix+opLen;wpos++){

	    tmpNode = &nds[wpos];
	    tmpNode->refPos=wpos+offs;
	    tmpNode->depth ++;
	    kputc('*',&tmpNode->seq);
	    kputw(seq_pos+1, &tmpNode->pos);
	    kputc(quals[seq_pos]+33, &tmpNode->qs);
	  }
	}

      }else if(opCode==BAM_CSOFT_CLIP){
	//occurs only at the ends of the read
	if(seq_pos == 0){
	  //then we are at beginning of read and need to write mapQ
	  tmpNode = &nds[wpos];
	  tmpNode->refPos=wpos+offs;
	  kputc('^', &tmpNode->seq);
	  if(rd->core.qual!=255)
	    kputc(rd->core.qual+33, &tmpNode->seq);
	  else
	    kputc('~', &tmpNode->seq);
	  seq_pos += opLen;
	  //	  wpos -= opLen;
	}else//we are at the end of read, then break CIGAR loop
	  break;
      }else if(opCode==BAM_CMATCH||opCode==BAM_CEQUAL||opCode==BAM_CDIFF) {
	hasInfo++;
	for(int fix=wpos ;wpos<(fix+opLen) ;wpos++) {
	  tmpNode =  &nds[wpos];
	  tmpNode->refPos=wpos+offs;
	  tmpNode->depth++;
	  if(seq_pos==0 &&hasPrintedMaq==0){
	    kputc('^', &tmpNode->seq);
	    if(rd->core.qual!=255)
	      kputc(rd->core.qual+33, &tmpNode->seq);
	    else
	      kputc('~', &tmpNode->seq);
	  }
	  char c = bam_nt16_rev_table[bam_seqi(seq, seq_pos)];

	  if(gf->ref==NULL ||gf->ref->chrLen<wpos+offs)//prints the oberved allele
	    kputc(bam_is_rev(rd)? tolower(c) : toupper(c), &tmpNode->seq);
	  else{
	    if(refToInt[c]==refToInt[gf->ref->seqs[wpos+offs]])
	      kputc(bam_is_rev(rd)? ',' : '.', &tmpNode->seq);
	    else
	      kputc(bam_is_rev(rd)? tolower(c) : toupper(c), &tmpNode->seq);
	  }
	  kputc(quals[seq_pos]+33, &tmpNode->qs);
	  kputw(seq_pos+1, &tmpNode->pos);
	  kputc(',',&tmpNode->pos);
	  tmpNode->len++;
	  seq_pos++;
	}
      }else if(opCode==BAM_CREF_SKIP) {
	for(int fix=wpos;wpos<fix+opLen;wpos++){
	  tmpNode = &nds[wpos];
	  tmpNode->refPos=wpos+offs;
	  tmpNode->depth ++;
	  bam_is_rev(rd)?kputc('<',&tmpNode->seq):kputc('>',&tmpNode->seq);
	  kputw(seq_pos+1, &tmpNode->pos);
	  kputc(quals[seq_pos]+33, &tmpNode->qs);
	}
      }else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP) {
	//dont care
      }else{
	fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
      }
    }
    //after end of read/parsing CIGAR always put the endline char
    //  fprintf(stderr,"printing endpileup for pos=%d\n",rd->pos);
    kputc('$', &tmpNode->seq);
  }

  //plug the reads back up //FIXME maybe do list type instead

  int miss= sgl->l-r;
  for(int i=0;i<(sgl->l-r);i++){
    sgl->first[i] =sgl->first[i+r];
    sgl->last[i] =sgl->last[i+r];
    std::swap(sgl->reads[i],sgl->reads[i+r]);
  }
  sgl->l = miss;
  //copy the part not meant for printing in this round into the buffer



  int lastSecureIndex =regionLen;
  //fprintf(stderr,"lastSecureIndex=%d\tregionLen=%d\t ret->nds.refpos=%d\n",lastSecureIndex,regionLen,ret->nds[regionLen-1].refPos);
  int tailPos = ret->nds[regionLen-1].refPos;
  if(tailPos>lastSecure)
    lastSecureIndex = regionLen-tailPos+lastSecure;
  ret->l = lastSecureIndex;
  
  if(regionLen-lastSecureIndex+4>np->m){
    delete [] np->nds;
    np->m=regionLen;
    kroundup32(np->m);
    np->nds = new node[np->m];
  }
  np->l=0;
  assert(regionLen-lastSecureIndex+4<=np->m);
  for(int i=lastSecureIndex;i<regionLen;i++)
    np->nds[np->l++] = nds[i];


  if(np->l!=0){
    np->first = np->nds[0].refPos;
    np->last = np->nds[np->l-1].refPos;
  }

  if(ret->l!=0){
    ret->first = ret->nds[0].refPos;
    ret->last = ret->nds[ret->l-1].refPos;
  }
  return dn;
}



nodePoolT mkNodes_one_sampleTb(readPool *sgl,nodePoolT *np,int refID) {
  int regionLen = coverage_in_bpT(np,sgl);//true covered regions
  nodePoolT dn;
  dn.l =0;
  if(regionLen==0)
    return dn;
  int lastSecure = sgl->lowestStart;
  
  dn = allocNodePoolT(regionLen);
  int offs = np->first;//first position from buffered
  int last = np->last+1;//because we are comparing with calc_end

  //plug in the old buffered nodes
  for(int i=0;i<np->l;i++)
    dn.nds[np->nds[i]->refPos-offs] = np->nds[i];

  if(np->l==0){//if we didn't have have anybuffered then
    offs = sgl->first[0];
    last = sgl->last[0];
  }

  //parse all reads
  int r;
  
  for( r=0;r<sgl->readIDstop;r++) {

    bam1_t *rd = sgl->reads[r];
    if(rd->core.tid!=refID){
      fprintf(stderr,"ReferenceID:(r:%d) for read:%s is not: %d but is:%d\n",r,bam_get_qname(rd),refID,rd->core.tid);
      exit(0);
    }
    int mapQ = rd->core.qual;
    if(mapQ>=255) mapQ = 20;

    if(sgl->first[r] > last){
      int diffs = (rd->core.pos-last);
      offs = offs + diffs;
      last = sgl->last[r];
    }else
      last = std::max(sgl->last[r],last);

    char *seq =(char *) bam_get_seq(rd);
    char *quals =(char *) bam_get_qual(rd);
    int nCig = rd->core.n_cigar;

    uint32_t *cigs = bam_get_cigar(rd);
    int seq_pos =0; //position within sequence
    int wpos = rd->core.pos-offs;//this value is the current position assocatied with the positions at seq_pos
    tNode *tmpNode =NULL; //this is a pointer to the last node beeing modified
    int hasInfo = 0;
    //loop through read by looping through the cigar string
    for(int i=0;i<nCig;i++) {
      int opCode = cigs[i]&BAM_CIGAR_MASK; //what to do
      int opLen = cigs[i]>>BAM_CIGAR_SHIFT; //length of what to do
      //fprintf(stderr,"opCode=%d opLen=%d seqPos=%d wpos=%d\n",opCode,opLen,seq_pos,wpos);
      if(opCode==BAM_CINS||opCode==BAM_CDEL){//handle insertions and deletions
	if(opCode==BAM_CINS&&i==0){ //skip indels if beginning of a read
	  seq_pos += opLen;	  
	  continue;
	}
	if(opCode==BAM_CINS&&hasInfo==0){
	  seq_pos += opLen;	  
	  continue;
	}
	
	hasInfo++;
	if(opCode&BAM_CINS){
	  wpos--; //insertion/deletion is bound to the last position of the read
	  tmpNode = dn.nds[wpos]?dn.nds[wpos]:(dn.nds[wpos]=initNodeT(UPPILE_LEN)); 
	  if(tmpNode->l2 >= tmpNode->m2){//check if we need to make a new vector because we have another insertion
	    tmpNode->m2++;
	    kroundup32(tmpNode->m2);
	    tmpNode->insert =(tNode**) realloc(tmpNode->insert,tmpNode->m2*sizeof(tNode*));
	  }
	  tNode *ins = initNodeT(opLen);
	  tmpNode->insert[tmpNode->l2] = ins; 
	  for(int ii=0;ii<opLen;ii++){
	    char c = bam_nt16_rev_table[bam_seqi(seq, seq_pos)];
	    tmpNode->insert[tmpNode->l2]->seq[ii] = bam_is_rev(rd)? tolower(c) : toupper(c);
	    tmpNode->insert[tmpNode->l2]->posi[ii] = seq_pos + 1;
	    tmpNode->insert[tmpNode->l2]->isop[ii] =rd->core.l_qseq- seq_pos - 1;
	    tmpNode->insert[tmpNode->l2]->qs[ii] = quals[seq_pos];
	    if(trim5&& (!bam_is_rev(rd))&&tmpNode->insert[tmpNode->l2]->posi[ii]<trim5)
	      tmpNode->insert[tmpNode->l2]->seq[ii] ='N';
	    if(trim3&& (bam_is_rev(rd))&&tmpNode->insert[tmpNode->l2]->isop[ii]<trim3)
	      tmpNode->insert[tmpNode->l2]->seq[ii] ='n';

	    if( quals[seq_pos]<minQ || seq_pos + 1 < trim || rd->core.l_qseq- seq_pos - 1 < trim)  
	      tmpNode->insert[tmpNode->l2]->seq[ii] =  bam_is_rev(rd)? tolower('n') : toupper('N');
	    tmpNode->insert[tmpNode->l2]->mapQ[ii] = mapQ;
	    seq_pos++;// <- important, must be after macro
	  }
	  tmpNode->l2++;//incrementor!
	  wpos++;
	}else {//this is the deletion part
	  
	  for(int ii=0;ii<opLen;ii++){
	    tmpNode = dn.nds[wpos+ii]?dn.nds[wpos+ii]:(dn.nds[wpos+ii]=initNodeT(UPPILE_LEN));  
	    //	    tmpNode = dn.nds[wpos+ii];
	    tmpNode->refPos = wpos +ii + offs;
	    tmpNode->deletion ++;
	  }
	  wpos += opLen;
	}

      }else if(opCode==BAM_CSOFT_CLIP){
	//occurs only at the ends of the read
	if(seq_pos == 0){
	  //then we are at beginning of read and need to write mapQ
	  tmpNode = dn.nds[wpos]?dn.nds[wpos]:(dn.nds[wpos]=initNodeT(UPPILE_LEN)); 
	  //	  tmpNode = dn.nds[wpos];
	  seq_pos += opLen;
	}else//we are at the end of read, then break CIGAR loop
	  break;
      }else if(opCode==BAM_CMATCH||opCode==BAM_CEQUAL||opCode==BAM_CDIFF) {
	hasInfo++;
	for(int fix=wpos ;wpos<(fix+opLen) ;wpos++) {
	  tmpNode = dn.nds[wpos]?dn.nds[wpos]:(dn.nds[wpos]=initNodeT(UPPILE_LEN)); 
	  // tmpNode =  dn.nds[wpos];
	  tmpNode->refPos=wpos+offs;
	  if(tmpNode->l>=tmpNode->m){//shouldnt need to realloc each member in struct
	    tmpNode->m = tmpNode->m*2;
	    tmpNode->seq =(char *) realloc(tmpNode->seq,tmpNode->m);
	    tmpNode->qs =(char *) realloc(tmpNode->qs,tmpNode->m);
	    tmpNode->posi =(suint *) realloc(tmpNode->posi,sizeof(suint)*tmpNode->m);
	    tmpNode->isop =(suint *) realloc(tmpNode->isop,sizeof(suint)*tmpNode->m);
	    tmpNode->mapQ = (unsigned char *) realloc(tmpNode->mapQ,tmpNode->m);
	  }
	  

	  char c = bam_nt16_rev_table[bam_seqi(seq, seq_pos)];
	  tmpNode->seq[tmpNode->l] = bam_is_rev(rd)? tolower(c) : toupper(c);
	  tmpNode->qs[tmpNode->l] =  quals[seq_pos];
	  // fprintf(stderr,"tmpNode->l:%d tmpNode->m:%d seq_pos:%d\n",tmpNode->l,tmpNode->m,seq_pos);
	  tmpNode->posi[tmpNode->l] = seq_pos;
	  tmpNode->isop[tmpNode->l] =rd->core.l_qseq- seq_pos-1;
	  tmpNode->mapQ[tmpNode->l] = mapQ;
	  if(trim5&& (!bam_is_rev(rd))&&tmpNode->posi[tmpNode->l]<trim5)
	    tmpNode->seq[tmpNode->l] ='N';
	  if(trim3&& (bam_is_rev(rd))&&tmpNode->isop[tmpNode->l]<trim3)
	    tmpNode->seq[tmpNode->l] ='n';
	  
	  if( quals[seq_pos]<minQ || seq_pos  < trim || rd->core.l_qseq - seq_pos - 1 < trim )
	    tmpNode->seq[tmpNode->l] = bam_is_rev(rd)? tolower('n') : toupper('N');
	  
	  tmpNode->l++;
	  seq_pos++;
	}
      }else if(opCode==BAM_CREF_SKIP) {
	for(int fix=wpos;wpos<fix+opLen;wpos++){
	  tmpNode = dn.nds[wpos]?dn.nds[wpos]:(dn.nds[wpos]=initNodeT(UPPILE_LEN)); 
	  //	    tmpNode = dn.nds[wpos];
	  tmpNode->refPos=wpos+offs;
	}
      }else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP) {
	//dont care
      }else{
	fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
      }
      //      exit(0);
    }
    //after end of read/parsing CIGAR always put the endline char
  }

  //plug the reads back up //FIXME maybe do list type instead

  int miss= sgl->l-r;
  for(int i=0;i<(sgl->l-r);i++){
    sgl->first[i] =sgl->first[i+r];
    sgl->last[i] =sgl->last[i+r];
    std::swap(sgl->reads[i],sgl->reads[i+r]);
  }
  sgl->l = miss;
  //copy the part not meant for printing in this round into the buffer



  int lastSecureIndex =regionLen;
  int tailPos = dn.nds[regionLen-1]->refPos;
  if(tailPos>lastSecure)
    lastSecureIndex = regionLen-tailPos+lastSecure;
  dn.l = lastSecureIndex;
  
  if(regionLen-lastSecureIndex+4>np->m){
    //delete [] np->nds;
    np->m=regionLen+4;
    kroundup32(np->m);
    np->nds =(tNode**) realloc(np->nds,np->m*sizeof(tNode*));
  }
  assert(regionLen-lastSecureIndex+4<=np->m);
  //  fprintf(stderr,"np->m:%d\n",np->m)
  np->l=0;
  for(int i=lastSecureIndex;i<regionLen;i++)
    np->nds[np->l++] = dn.nds[i];


  if(np->l!=0){
    np->first = np->nds[0]->refPos;
    np->last = np->nds[np->l-1]->refPos;
  }

  if(dn.l!=0){
    dn.first = dn.nds[0]->refPos;
    dn.last = dn.nds[dn.l-1]->refPos;
#if 0
    fprintf(stderr,"[%s] first=%d last=%d\n",__FUNCTION__,ret->first,ret->last);
#endif
  }
  return dn;
}



/*
  old function.
  from:= min(dn);
  to:=min(max{dn[1],dn[2],...})

*/


void get_span_all_samples(const nodePool *dn,int nFiles,int &from,int &to){
  for(int i=0;0&&i<nFiles;i++)
    fprintf(stderr,"[%s]i=%d (%d,%d)=%d\n",__FUNCTION__,i,dn[i].first,dn[i].last,dn[i].l);

  for(int i=0;i<nFiles;i++)
    if(dn[i].l!=0){
      from = dn[i].first;
      to = dn[i].last;
      break;
    }
  
  for(int i=0;i<nFiles;i++){//might skip the 'i' used above
    if(dn[i].l==0)
      continue;
    if(dn[i].last>to)
      to = dn[i].last;
    if(dn[i].first<from)
      from = dn[i].first;
  }
}


void get_span_all_samplesT(const nodePoolT *dn,int nFiles,int &from,int &to){
  for(int i=0;0&&i<nFiles;i++)
    fprintf(stderr,"[%s]i=%d (%d,%d)=%d\n",__FUNCTION__,i,dn[i].first,dn[i].last,dn[i].l);

  for(int i=0;i<nFiles;i++)
    if(dn[i].l!=0){
      from = dn[i].first;
      to = dn[i].last;
      break;
    }
  
  for(int i=0;i<nFiles;i++){//might skip the 'i' used above
    if(dn[i].l==0)
      continue;
    if(dn[i].last>to)
      to = dn[i].last;
    if(dn[i].first<from)
      from = dn[i].first;
  }
}


typedef std::map<int,tNode **> umapT; 

chunkyT *slow_mergeAllNodes_new(nodePoolT *dn,int nFiles){
  //  fprintf(stderr,"starting [%s]\n",__FUNCTION__);
  
  umapT ret;
  umapT::iterator it;
  for(int f=0;f<nFiles;f++) {
    nodePoolT sm = dn[f];
    for(int l=0;l<sm.l;l++) {
      tNode **perSite = NULL;
      int thepos =sm.nds[l]->refPos;
      it=ret.find(thepos);
      if(it!=ret.end())
	perSite = it->second;
      else{
	perSite = (tNode**)calloc(nFiles,sizeof(tNode*));
	ret.insert(std::make_pair(thepos,perSite));
      }
      perSite[f] = sm.nds[l];
    }

  }
  //now we perform the backrool
  int nnodes = ret.size();
  int *refPos = new int [ret.size()];
  chunkyT *chk =new chunkyT;  
  chk->nd = new tNode**[nnodes];
  int p=0;
  for(it = ret.begin();it!=ret.end();++it){
    chk->nd[p++] = it->second;
    refPos[p-1] = it->first;
  }

  chk->nSamples=nFiles;
  chk->nSites=nnodes;
  chk->refPos = refPos;
  return chk;  
}


chunkyT *mergeAllNodes_new(nodePoolT *dn,int nFiles) {

  int *depth = NULL;
  int *refPos2 = NULL;
  tNode ***super = NULL;
  
  if(dn==NULL){//cleanup called after end of looping through all files.
    return NULL;
  }
  int first,last;
  get_span_all_samplesT(dn,nFiles,first,last);

  int rlen = last-first+1;
  assert(rlen>=0);
  if(rlen>BUG_THRES)
    return slow_mergeAllNodes_new(dn,nFiles);

  super = new tNode**[rlen];
  depth = new int[rlen];
  refPos2 = new int[rlen];
  
  memset(depth,0,rlen*sizeof(int));
  int offs = first;
  //looping through the different samples
  for(int n=0;n<nFiles;n++) {
    nodePoolT sm = dn[n];
    int i;
    //looping through all nodes
    for( i=0;((i<sm.l)&&(sm.nds[i]->refPos <= std::min(sm.last,last) ));i++) {
      
      int posi = sm.nds[i]->refPos-offs;
      if(depth[posi]==0){
	//	fprintf(stderr,"posi=%d\n",posi);
	super[posi] =(tNode**) calloc(nFiles,sizeof(tNode*));
	refPos2[posi]=sm.nds[i]->refPos;
      }
      depth[posi]++;
      super[posi][n] = sm.nds[i];
      //      super[posi][n].len++;
    }
  }
  int nnodes =0;
  for(int i=0;i<rlen;i++)
    if(depth[i]!=0)
      nnodes++;
  
  //  fprintf(stderr,"YYYY nnnosed=%d\n",nnodes);
  chunkyT *chk =new chunkyT;  
  chk->nd = new tNode**[nnodes];
  int p=0;
  int *refPos = new int[nnodes];
  for(int i=0;i<nnodes;i++)
    refPos[i] = -9;
  for(int i=0;i<rlen;i++){
    if(depth[i]!=0){
      //   fprintf(stderr,"refpos2[%d]=%d\n",i,refPos2[i]);
      refPos[p] = refPos2[i];
      chk->nd[p++] = super[i];
    }
    //   delete [] super[i];
  }

  delete [] super;
  delete [] refPos2;
  delete [] depth;
  chk->nSamples=nFiles;
  chk->nSites=nnodes;
  chk->refPos = refPos;
  return chk;
}


typedef std::map<int,node *> umap; 

chunky *slow_mergeAllNodes_old(nodePool *dn,int nFiles){
  //  fprintf(stderr,"starting [%s]\n",__FUNCTION__);
  
  umap ret;
  umap::iterator it;
  for(int f=0;f<nFiles;f++){
    nodePool sm = dn[f];
    for(int l=0;l<sm.l;l++) {
      node *perSite = NULL;
      int thepos =sm.nds[l].refPos;
      it=ret.find(thepos);
      if(it!=ret.end())
	perSite = it->second;
      else{
	perSite = new node[nFiles];
	for(int ii=0;ii<nFiles;ii++)
	  perSite[ii] = initNode(UPPILE_LEN,thepos,thepos);
	
	ret.insert(std::make_pair(thepos,perSite));
		
      }
      perSite[f] = sm.nds[l];
    }
  }
  //now we perform the backrool
  int nnodes = ret.size();
  int *refPos = new int [ret.size()];
  //  fprintf(stderr,"nnodes=%d\n",nnodes);
  chunky *chk =new chunky;  
  chk->nd = new node*[nnodes];
  int p=0;
  for(it = ret.begin();it!=ret.end();++it){
    chk->nd[p++] = it->second;
    refPos[p-1] = chk->nd[p-1][0].refPos;
  }

  chk->nSamples=nFiles;
  chk->nSites=nnodes;
  chk->refPos = refPos;
  chk->first = refPos[0];
  chk->last = refPos[nnodes-1];
  return chk;  
}




chunky *mergeAllNodes_old(nodePool *dn,int nFiles) {
  int first,last;
  get_span_all_samples(dn,nFiles,first,last);
  int rlen = last-first+1; 

  if(rlen>BUG_THRES) 
    return slow_mergeAllNodes_old(dn,nFiles);
  
  int depth[rlen];
  
  node **super;
  try{
    //    fprintf(stderr,"rlen=%d\n",rlen);
    fflush(stderr);
    super = new node*[rlen];
  }
  catch(char *str){
    fprintf(stderr,"rlen=%d\n",rlen);
    fprintf(stderr,"problems allocating:%s\n",str);
  }
  int *refPos2 = new int[rlen];
  memset(depth,0,sizeof(int)*rlen);
  int offs = first;
  //looping through the different samples
  for(int n=0;n<nFiles;n++){
    nodePool sm = dn[n];
    int i;
    //looping through all nodes
    for( i=0;((i<sm.l)&&(sm.nds[i].refPos <= std::min(sm.last,last) ));i++) {
      
      int posi = sm.nds[i].refPos-offs;
      if(depth[posi]==0){
	super[posi] = new node[nFiles];
	for(int ii=0;ii<nFiles;ii++){
	  super[posi][ii] = initNode(UPPILE_LEN,posi+offs,posi);
	}
	
	refPos2[posi]=sm.nds[i].refPos;
      }
      depth[posi]++;
      super[posi][n] = sm.nds[i];
      super[posi][n].len++;
    }
    
  }
  int nnodes =0;
  for(int i=0;i<rlen;i++)
    if(depth[i]!=0)
      nnodes++;

  //  fprintf(stderr,"YYYY nnnosed=%d\n",nnodes);
  chunky *chk =new chunky;  
  chk->nd = new node*[nnodes];
  int p=0;
  int *refPos = new int[nnodes];
  for(int i=0;i<nnodes;i++)
    refPos[i] = -9;
  for(int i=0;i<rlen;i++){
    if(depth[i]!=0){
      //   fprintf(stderr,"refpos2[%d]=%d\n",i,refPos2[i]);
      refPos[p] = refPos2[i];
      chk->nd[p++] = super[i];
    }
    //   delete [] super[i];
  }
  //  fprintf(stderr,"nnodes=%d p=%d\n",nnodes,p);
  delete [] super;
  delete [] refPos2;

  chk->nSamples=nFiles;
  chk->nSites=nnodes;
  chk->refPos = refPos;
  chk->first = refPos[0];
  chk->last = refPos[nnodes-1];
  return chk;
}





int getSglStop5(readPool *sglp,int nFiles,int pickStop) {
#if 0
    fprintf(stderr,"[%s]\n",__FUNCTION__);
    for(int i=0;i<nFiles;i++){
      if(sglp[i].l!=0)
	fprintf(stderr,"ret.l=%d (%d,%d)\n",sglp[i].l,sglp[i].first[0],sglp[i].last[sglp[i].l-1]);
      else
	fprintf(stderr,"ret.l=%d (-1,-1)\n",sglp[i].l);
      for(int s=0;0&&s<sglp[i].l;s++)
	fprintf(stderr,"%d %d\n",sglp[i].first[s],sglp[i].last[s]);
    }

    fprintf(stderr,"\n");
#endif
    int lowestStart = pickStop ;
    if(lowestStart==-1){
      for(int i=0;i<nFiles;i++)
	if(sglp[i].l>0&&(getmax(sglp[i].first,sglp[i].l)<lowestStart))
	  lowestStart = getmax(sglp[i].first,sglp[i].l);
      assert(lowestStart!=pickStop);
    }

  lowestStart--;
  //fprintf(stderr,"loweststart=%d\n",lowestStart);
  /*
    now we looptrhough each readpool, and find 2 values.
    1) the readnumber of the last read to be included for uppiling,this
    2) the lastposition for which we know that the uppiling is complete
  */
  int tmp =0;
  for(int i=0;i<nFiles;i++){
    sglp[i].lowestStart = lowestStart;
    int j;
    for( j=0;j<sglp[i].l;j++)
      if(sglp[i].first[j]>lowestStart)
	break;
    sglp[i].readIDstop=j;
    tmp += j;
  }
  return tmp;
}


/*
  this function just returns the highest position along all buffered nodepools and readpools
 */

void getMaxMax2(readPool *sglp,int nFiles,nodePool *nps){
#if 0
    fprintf(stderr,"[%s].nFiles=%d\n",__FUNCTION__,nFiles);
    for(int i=0;i<nFiles;i++){
      fprintf(stderr,"%d=sglp.l=%d\n",i,sglp[i].l);
    }
#endif
    
    int last_bp_in_chr = -1;
    
    for(int i=0;i<nFiles;i++){
      readPool *tmp = &sglp[i];
      if(tmp->l>0 && (getmax(tmp->last,tmp->l)>last_bp_in_chr) )
	last_bp_in_chr = getmax(tmp->last,tmp->l);
      if(nps[i].last>last_bp_in_chr)
	last_bp_in_chr = nps[i].last;
    }
    //  assert(last_bp_in_chr!=-1); This assertion should not be active, since it might be true if we don't have data for a whole chromosomes
    for(int i=0;i<nFiles;i++){
      sglp[i].lowestStart = last_bp_in_chr;
      sglp[i].readIDstop=sglp[i].l;
    }
    
}

void getMaxMax2(readPool *sglp,int nFiles,nodePoolT *nps){
#if 0
  fprintf(stderr,"[%s].nFiles=%d\n",__FUNCTION__,nFiles);
  for(int i=0;i<nFiles;i++)
    fprintf(stderr,"%d=sglp.l=%d\n",i,sglp[i].l);

#endif

  int last_bp_in_chr = -1;

  for(int i=0;i<nFiles;i++){
    readPool *tmp = &sglp[i];
    if(tmp->l>0 && (getmax(tmp->last,tmp->l)>last_bp_in_chr) )
      last_bp_in_chr = getmax(tmp->last,tmp->l);
    if(nps[i].last>last_bp_in_chr)
      last_bp_in_chr = nps[i].last;
  }
  //  assert(last_bp_in_chr!=-1); This assertion should not be active, since it might be true if we don't have data for a whole chromosomes
  for(int i=0;i<nFiles;i++){
    sglp[i].lowestStart = last_bp_in_chr;
    sglp[i].readIDstop=sglp[i].l;
  }

}


void getOffsets(htsFile *fp,char *fn,const bam_hdr_t *hd,hts_idx_t **idx,hts_itr_t **itr,int ref,int start,int stop,bam_hdr_t *hdr){

  if(*idx==NULL)
    *idx = sam_index_load(fp,fn);

  if (*idx == NULL) { // index is unavailable
    fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files. (file: \'%s\' )\n",fn);
    exit(0);
  }
  char tmp[1024];
  snprintf(tmp,1024,"%s:%d-%d",hd->target_name[ref],start+1,stop);
  if(*itr)
    hts_itr_destroy(*itr);
  *itr = sam_itr_querys(*idx, hdr, tmp);
  if (*itr == NULL) { // reference name is not found
    fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n",tmp);
    exit(0);
  }

}



void *setIterator1(void *args){
  bufReader *rd = (bufReader *) args;

  int start,stop,ref;
  start = rd->regions.start;
  stop = rd->regions.stop;
  ref = rd->regions.refID;
  getOffsets(rd->fp,rd->fn,rd->hdr,&rd->idx,&rd->itr,ref,start,stop,rd->hdr);//leak
  return EXIT_SUCCESS;
}


 void setIterator1_thread(void *args,int index){
   //   fprintf(stderr,"[%s]\n",__FUNCTION__);
  bufReader *rds = (bufReader *) args;
  bufReader *rd = &rds[index];
  // free(rd->it.off);
  int start,stop,ref;
  start = rd->regions.start;
  stop = rd->regions.stop;
  ref = rd->regions.refID;
  getOffsets(rd->fp,rd->fn,rd->hdr,&rd->idx,&rd->itr,ref,start,stop,rd->hdr);
}

void setIterators(bufReader *rd,regs regions,int nFiles,int nThreads){
  for(int i=0;1&&i<nFiles;i++){
    rd[i].regions=regions;
    if(nThreads==1)
      setIterator1(&rd[i]);
  }
  if(nThreads>1){
    pthread_t myT[nThreads];
    int cnt=0;
    while(cnt<nFiles){
      int nTimes;
      if(nFiles-cnt-nThreads>=0)
	nTimes = nThreads;
      else
	nTimes = nFiles-cnt;
      for(int i=0;0&&i<nTimes;i++)
	fprintf(stderr,"cnt:%d i:%d\n",cnt+i,i);
      for(int i=0;i<nTimes;i++)
	pthread_create(&myT[i],NULL,setIterator1,&rd[cnt+i]);
      for(int i=0;i<nTimes;i++)
	pthread_join(myT[i], NULL);

      cnt+=nTimes;
    }
  }
}


void callBack_bambi(fcb *fff);//<-angsd.cpp
//type=1 -> samtool mpileup textoutput
//type=0 -> callback to angsd

void destroy_tnode_pool(){
  if(tnodes==NULL)
    return;
  flush_queue();
  //  fprintf(stderr,"\n\t-> npools:%zu unfreed tnodes before clean:%lu \n",tnodes->npools,currentnodes);
  for(int i=0;i<tnodes->npools;i++){
    //fprintf(stderr,"%d/%d\n",i,tnodes->npools);    
    tNode **nd =(tNode**)tnodes->pools[i].pool;
    
    int nitem = tnodes->pools[i].used/tnodes->dsize;

    for(int j=0;j<nitem;j++){
      tNode *tn =(tNode*) nd+j;
      if(tn->m){
	//	fprintf(stdout,"nd[%d] :%p tn.refPos:%d tn.l:%d tn.m:%d\n",j,tn,tn->refPos,tn->l,tn->m);
	free(tn->seq);
	free(tn->qs);
	free(tn->posi);
	free(tn->isop);
	free(tn->mapQ);
	tn=NULL;
      }
    }
  }



  //  delete [] ary;
  tpool_destroy(tnodes);
}


//f Typenames with T are the ones for the callback, in main angsd
//Most likely this can be written more beautifull by using templated types. But to lazy now.
int uppile(int show,int nThreads,bufReader *rd,int nLines,int nFiles,std::vector<regs> &regions,abcGetFasta *gf) {
  
  assert(nLines&&nFiles);
  fprintf(stderr,"\t-> Parsing %d number of samples \n",nFiles);
  fflush(stderr);
  if(show!=1)
    pthread_mutex_lock(&mUpPile_mutex);//just to make sure, its okey to clean up

  readPool *sglp= new readPool[nFiles];
  
  nodePool *nps =NULL;
  nodePoolT *npsT = NULL;
  if(show)
    nps = new nodePool[nFiles];// <- buffered nodes
  else
    npsT = new nodePoolT[nFiles];// <- buffered nodes

  for(int i=0;i<nFiles;i++){
    sglp[i] = makePoolb(nLines);
    if(show)
      nps[i] = allocNodePool(MAX_SEQ_LEN);
    else
      npsT[i] = allocNodePoolT(MAX_SEQ_LEN);
  }

  int itrPos = -1;
  int sumEof = 0;//sumof readerobject that has eof

  int theRef=-1;//<- this is current referenceId,
  while( SIG_COND) {

    int notDone = nFiles;
    sumEof =0;
    //reset the done region flag, while checking that any file still has data
    for(int i=0;i<nFiles;i++){
      sumEof += rd[i].isEOF;
      rd[i].regionDone =0;
    }
#if 0
    for(int i=0;i<nFiles;i++)
      fprintf(stderr,"[%s] i=%d eof=%d regionDone=%d\n",__FUNCTION__,i,rd[i].isEOF,rd[i].regionDone);
#endif

    //break loop if all files are done
    if(sumEof==nFiles&&regions.size()==0)
      break;

    if(regions.size()!=0) {
      if(itrPos+1==(int)regions.size()){
	break;
      }else {
	itrPos++;
	fprintf(stderr,"\t-> Region lookup %d/%lu\n",itrPos+1,regions.size());
	fflush(stderr);
	setIterators(rd,regions[itrPos],nFiles,nThreads);
	//fprintf(stderr,"done region lookup %d/%lu\n",itrPos+1,regions.size());
	//fflush(stderr);///BAME

	theRef = rd[0].itr->tid; 
	//	fprintf(stderr,"theRef:%d %p %p\n",theRef,rd[0].it.off,rd[0].it.hts_itr->off);
	//validate we have offsets;
	int gotData = 0;
	for(int i=0;i<nFiles;i++)
	  if((rd[i].itr &&(rd[i].itr->off!=NULL))){
	    gotData =1;
	    break;
	  }
	//reset eof and region donepointers.
	for(int i=0;i<nFiles;i++){
	  rd[i].isEOF = 0;
	  rd[i].regionDone =0;
	  sglp[i].l =0;
	  if(sglp[i].bufferedRead)
	    bam_destroy1(sglp[i].bufferedRead);
	  sglp[i].bufferedRead=NULL;
	}
      }
    }else{
      if(theRef==-1){//init
	theRef =0;

      }else if(theRef==rd[0].hdr->n_targets-1){
	break;//then we are done
      }else{
	int minRef = rd[0].hdr->n_targets;
	for(int i=0;i<nFiles;i++)
	  if(sglp[i].bufferedRead && sglp[i].bufferedRead->core.tid<minRef)
	    minRef = sglp[i].bufferedRead->core.tid;
	if(minRef==-2||minRef == rd[0].hdr->n_targets){
	  theRef++;
	}else
	  theRef = minRef;
      }
    }

    if(theRef==rd[0].hdr->n_targets)
      break;
    assert(theRef>=0 && theRef<rd[0].hdr->n_targets);//ref should be inrange [0,#nref in header]

    //load fasta for this ref if needed
    void waiter(int);
    waiter(theRef);//will wait for exisiting threads and then load stuff relating to the chromosome=theRef;
    if(gf->ref!=NULL && theRef!=gf->ref->curChr)
      gf->loadChr(gf->ref,rd[0].hdr->target_name[theRef],theRef);
    
     

    //now we have changed chromosome we should plug in buffered, if buffered is same as the new chr
    //this is not needed if we use the indexing

    if(regions.size()==0){
      for(int i=0;i<nFiles;i++){
	extern abcGetFasta *gf;
	if(sglp[i].bufferedRead ==NULL)//<- no buffered no nothing
	  continue;
	if(sglp[i].bufferedRead->core.tid==theRef) {//buffered is on correct chr
	  int doCpy =1;
	  if(gf->ref!=NULL){
	    //this will modify the quality and mapQ if -baq >0 or adjustMapQ= SAMtools -c
	    doCpy = restuff(sglp[i].bufferedRead);
	  }
	  if(doCpy){
	    std::swap(sglp[i].reads[sglp[i].l],sglp[i].bufferedRead);
	    sglp[i].first[sglp[i].l] = sglp[i].reads[sglp[i].l]->core.pos;
	    sglp[i].last[sglp[i].l] =   bam_endpos(sglp[i].reads[sglp[i].l]);
	    sglp[i].l++;
	  }
	  //if we haven't performed the copy its because adjustedMap<minMapQ, in either case we don't need the read anymore
	  if(sglp[i].bufferedRead)
	    bam_destroy1(sglp[i].bufferedRead);
	  sglp[i].bufferedRead = NULL;
	}else // buffered was not on correct chr say that the 'region' is done
	  rd[i].regionDone=1;
      }
    }

    //below loop will continue untill eoc/oef/eor <- funky abbrev =end of chr, file and region
    while(notDone) {

      //before collecting reads from the files, lets first check if we should pop reads from the buffered queue in each sgl
      for(int i=0;i<nFiles;i++){
	if(sglp[i].bufferedRead&&sglp[i].bufferedRead->core.tid==theRef){
	  std::swap(sglp[i].reads[sglp[i].l] , sglp[i].bufferedRead);
	  sglp[i].first[sglp[i].l] = sglp[i].reads[sglp[i].l]->core.pos;
	  sglp[i].last[sglp[i].l] = bam_endpos(sglp[i].reads[sglp[i].l]);
	  sglp[i].l++;
	  if(sglp[i].bufferedRead)
	    bam_destroy1(sglp[i].bufferedRead);
	  sglp[i].bufferedRead = NULL;
	}
	
      }
      
      int pickStop=MAX_SEQ_LEN;
      
#ifndef __NEW__
      int doFlush =  (collect_reads(rd,nFiles,notDone,sglp,nLines,theRef,pickStop)==nFiles)?1:0 ;
#else
      int doFlush = (collect_reads2(rd,nFiles,notDone,sglp,nLines,theRef,pickStop)==nFiles)?1:0 ;
#endif
#if 0
      for(int i=0;i<nFiles;i++)
	fprintf(stderr,"2) [%s] sgl[%d].l=%d (%d:%d) l.0refid:%d\n",__FUNCTION__,i,sglp[i].l,sglp[i].first[0],sglp[i].last[sglp[i].l-1],sglp[i].reads[0]->core.tid);
      
#endif
      
      //prepare uppiling for all data, if end of chromosome or all files are EOF
      if(doFlush||notDone==0){
	//fprintf(stderr,"[%s]Last chunk Will flush all remaining reads\n",__FUNCTION__);
	if(show)
	  getMaxMax2(sglp,nFiles,nps);
	else
	  getMaxMax2(sglp,nFiles,npsT);
      }else{
	//getSglStop returns the sum of reads to parse.
	int hasData = getSglStop5(sglp,nFiles,pickStop);
	//	fprintf(stderr,"hasData=%d\n",hasData);
	if(hasData==0&&notDone!=0){
	  continue;
	}
      }

      //make uppile nodes
      nodePool dn[nFiles];
      nodePoolT *dnT = NULL;
      int tmpSum = 0;

      if(show!=1){
	dnT = new nodePoolT[nFiles];
	for(int i=0;i<nFiles;i++){
	  dnT[i] = mkNodes_one_sampleTb(&sglp[i],&npsT[i],theRef);
	  tmpSum += dnT[i].l;
	}
      }else
	for(int i=0;i<nFiles;i++) {
	  dn[i] = mkNodes_one_sampleb(&sglp[i],&nps[i],gf);
	  tmpSum += dn[i].l;
	}
      //fprintf(stderr,"tmpSum:%d\n",tmpSum);
#if 0
	for(int i=0;i<nFiles;i++)
	  if(show!=1)
	    fprintf(stderr,"[%s.%s():%d] l=%d first=%d last=%d\n",__FILE__,__FUNCTION__,__LINE__,dn[i].l,dn[i].first,dn[i].last);
	  else
	    fprintf(stderr,"[%s.%s():%d] l=%d first=%d last=%d\n",__FILE__,__FUNCTION__,__LINE__,dnT[i].l,dnT[i].first,dnT[i].last);
#endif
            
      //simple sanity check for validating that we have indeed something to print.

      if(tmpSum==0){
	if(regions.size()!=0)
	  notDone=1;
	else{
	  fprintf(stderr,"No data for chromoId=%d chromoname=%s\n",theRef,rd[0].hdr->target_name[theRef]);
	  fprintf(stderr,"This could either indicate that there really is no data for this chromosome\n");
	  fprintf(stderr,"Or it could be problem with this program regSize=%lu notDone=%d\n",regions.size(),notDone);
	}
	delete [] dnT;
	break;
      }

      int regStart,regStop;
      if(itrPos!=-1){
	regStart = regions[itrPos].start;
	regStop = regions[itrPos].stop;
      }else{
	regStart = -1;
	regStop = rd[0].hdr->target_len[theRef];
      }
      if(show){
	//merge the perFile upnodes.FIXME for large regions (with gaps) this will allocate to much...
	chunky *chk =mergeAllNodes_old(dn,nFiles);
	assert(chk->refPos[0]<=regStop);
	
	chk->regStart = regStart;
	chk->regStop = regStop;
	chk->refName = rd[0].hdr->target_name[theRef];
	chk->refId = theRef;
	printChunky2(chk,stdout,chk->refName,gf);
	
	for(int i=0;1&&i<nFiles;i++){
	  if(dn[i].l!=0)
	    delete [] dn[i].nds;
	}
      }else{

	fcb *f = new fcb; //<-for call back
	f->dn=dnT; f->nFiles = nFiles; f->regStart = regStart; f->regStop = regStop; f->refId = theRef;
	callBack_bambi(f);
      }
      if(SIG_COND==0){
	for(int i=0;i<nFiles;i++){
	  rd[i].isEOF = 1;
	  rd[i].regionDone =1;
	}
      }
      
      //if we flush, its due to end of chr/region or end of files
      if(doFlush)
	break;
    }

  }
  //below can be written nice but is just a copy paste, time is 10pm to lazy now
  if(show==0) {
    callBack_bambi(NULL);//cleanup signal
    pthread_mutex_lock(&mUpPile_mutex);//just to make sure, its okey to clean up
    for(int i=0;1&&i<nFiles;i++){
      free(npsT[i].nds);
      dalloc(&sglp[i]);
    }
    delete [] npsT;
    delete [] sglp;
    pthread_mutex_unlock(&mUpPile_mutex);//just to make sure, its okey to clean up
  }else{
    delete [] nps;
    delete [] sglp;
  }
  void waiter(int);
  waiter(theRef);//will wait for exisiting threads and then load stuff relatin

  return 0;
}

