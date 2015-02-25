#include <vector>
#include <cstring>
#include <cstdlib>
#include <ctype.h>
#include <pthread.h>
#include <cassert>
#include "bams.h"
#include <htslib/kstring.h>
#include "mUpPile.h"
#include "abcGetFasta.h"
#include "analysisFunction.h"
#include "bams.h"
#include <htslib/hts.h>
extern int SIG_COND;
extern int minQ;
extern int trim;
static const char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";


/*
  When merging different nodes, we are using a greedy but fast approach,
  however if we have a large genomic region with no data we are are allocating huge chunks.
  So the cutoff below will make a fallback to a slower but memsecure version
 */
#define BUG_THRES 1000000 //<- we allow 1mio sites to be allocated,

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


void dalloc_node(tNode &n){
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);
  if(n.l!=0||n.m!=0){
    free(n.seq);
    free(n.qs);
    free(n.posi);
    free(n.mapQ);
    free(n.isop);
  }
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
  np.nds = new tNode[np.m];
  return np;
}


void dalloc_node(node &n){
  free(n.seq.s);
  free(n.qs.s);
  free(n.pos.s);
}


void dalloc_nodePool (nodePool& np){
  for(int i=0;i<np.l;i++)
    dalloc_node(np.nds[i]);
  
}

void dalloc_nodePoolT (nodePoolT& np){
  for(int i=0;i<np.l;i++)
    dalloc_node(np.nds[i]);
  
}

 

void cleanUpChunkyT(chunkyT *chk){
  for(int s=0;s<chk->nSites;s++) {
    for(int i=0;i<chk->nSamples;i++) {
      if(chk->nd[s][i].l2!=0)
	for(int j=0;j<chk->nd[s][i].l2;j++){
	  dalloc_node(chk->nd[s][i].insert[j]);
	}
      free(chk->nd[s][i].insert);
      dalloc_node(chk->nd[s][i]);
    }
    delete [] chk->nd[s];
  }
  delete [] chk->nd;
  delete [] chk->refPos;
  delete chk;
  
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


tNode initNodeT(int l){
  tNode d;
  d.l = d.l2 = d.m2 = 0;
  d.m = l;
  kroundup32(d.m);
  if(d.m==0){
    d.seq=d.qs=NULL;
    d.posi=d.isop=d.mapQ=NULL;
  }else{
    d.seq=(char *)malloc(d.m);
    d.qs=(char *)malloc(d.m);
    d.posi=(unsigned char *)malloc(d.m);
    d.isop=(unsigned char *)malloc(d.m);
    d.mapQ=(unsigned char *)malloc(d.m);
  }
  d.refPos= -999;
  d.insert = NULL;
  return d;
}


void realloc(char **c,int l,int m){
  char *tmp =(char *) malloc(m);
  memcpy(tmp,*c,l);
  free(*c);
  *c= tmp;
}

void realloc(unsigned char **c,int l,int m){
  unsigned char *tmp =(unsigned char *) malloc(m);
  memcpy(tmp,*c,l);
  free(*c);
  *c= tmp;
}
void realloc(tNode *d,int newsize){
  kroundup32(newsize);
  if(newsize<=d->m)
    fprintf(stderr,"[%s] problems newsize should be bigger than oldesize\n",__FUNCTION__);
  d->m=newsize;
  realloc(&(d->seq),d->l,d->m);
  realloc(&(d->qs),d->l,d->m);
  realloc(&(d->posi),d->l,d->m);
  realloc(&(d->isop),d->l,d->m);
  realloc(&(d->mapQ),d->l,d->m);
}


void realloc(sglPool &ret,int l){
  ret.m =l;
  kroundup32(ret.m);
  aRead *tmp = new aRead[ret.m];
  memcpy(tmp,ret.reads,sizeof(aRead)*ret.l);
  delete [] ret.reads;
  ret.reads = tmp;
  
  int *tmpI = new int[ret.m];
  memcpy(tmpI,ret.first,sizeof(int)*ret.l);
  delete [] ret.first;
  ret.first = tmpI;
  
    
  tmpI = new int[ret.m];
  memcpy(tmpI,ret.last,sizeof(int)*ret.l);
  delete [] ret.last;
  ret.last = tmpI;
  
}



void read_reads_usingStop(htsFile *fp,int nReads,int &isEof,sglPool &ret,int refToRead,iter_t *it,int stop,int &rdObjEof,int &rdObjRegionDone,bam_hdr_t *hdr) {
#if 0
  fprintf(stderr,"[%s]\n",__FUNCTION__);
#endif 

  //if should never be in this function if we shouldnt read from the file.
  assert(rdObjRegionDone!=1 &&rdObjEof!=1 );
  
  if((nReads+ret.l)>ret.m)
    realloc(ret,nReads+ret.l);


  //this is the awkward case this could cause an error in some very unlikely scenario.
  //whith the current buffer empty and the first read being a new chromosome. 
  //this is fixed good job monkey boy
  while(ret.l==0){
    aRead &b = ret.reads[0];
    int tmp;
    b.vDat = new uint8_t[RLEN];
    if((tmp=bam_iter_read1(fp,it,b,hdr))<0){//FIXME
      if(tmp==-1){
	rdObjEof =1;
	//	ret.isEOF =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }
    if(b.nCig!=0){
      ret.first[ret.l] = b.pos;
      ret.last[ret.l] = bam_calend(b,getCig(&b));
      ret.l++;
      nReads--;//breaking automaticly due to ret.l++
    }
  }

   /*
    now do the general loop,
    1) read an alignemnt calulate start and end pos
    check
    a) if same chromosome that a new read has a geq FINE
    b) if same but < then UNSORTED
    
  */ 
  int i=0;
  int tmp;
  int buffed=0;
  while(1) {
    if(i+ret.l>=(ret.m-1)){
      ret.l += i;
      i=0;
      realloc(ret,ret.m+1);//will double buffer
    }
    
    aRead &b = ret.reads[i+ret.l];
    b.vDat = new uint8_t[RLEN];
    if((tmp=bam_iter_read1(fp,it,b,hdr))<0){
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;

      delete [] b.vDat;
      break;
    }if(b.nCig==0){
      delete [] b.vDat;
      continue;
    }
    //general fine case
    
    if((refToRead==b.refID)&&( b.pos >= ret.first[ret.l+i-1])&&(b.pos<stop)){
      ret.first[ret.l+i] = b.pos;
      ret.last[ret.l+i] = bam_calend(b,getCig(&b));
    }else if((refToRead==b.refID)&&b.pos>=stop){
      buffed=1;
      ret.bufferedRead = b;
      break;
    }
    else if(b.refID>refToRead){
      buffed=1;
      ret.bufferedRead = b;
      rdObjRegionDone =1;
      break;
    }else if(ret.first[ret.l+i-1]>b.pos){
      fprintf(stderr,"unsorted file detected will exit\n");
      fflush(stderr);
      exit(0);
    }
    i++;
  }
  ret.nReads =ret.l+i;
  ret.l = ret.nReads;
}


//nothing with buffered here
void read_reads_noStop(htsFile *fp,int nReads,int &isEof,sglPool &ret,int refToRead,iter_t *it,int &rdObjEof,int &rdObjRegionDone,bam_hdr_t *hdr) {
#if 0
  fprintf(stderr,"\t->[%s] buffRefid=%d\trefToRead=%d\n",__FUNCTION__,ret.bufferedRead.refID,refToRead);
#endif
  assert(rdObjEof==0 && ret.bufferedRead.refID==-2);
  
  if((nReads+ret.l)>ret.m)
    realloc(ret,nReads+ret.l);
 
  //start by allocating the memeory we need
  for(int i=0;i<nReads;i++){
    aRead b;
    b.vDat=new uint8_t[RLEN];
    ret.reads[ret.l+i] =b;
  }

  //this is the awkward case this could cause an error in some very unlikely scenario.
  //whith the current buffer empty and the first read being a new chromosome. 
  //this is fixed good job monkey boy
  while(ret.l==0){
    aRead &b = ret.reads[0];
    int tmp;
    if((tmp=bam_iter_read1(fp,it,b,hdr))<0){//FIXME
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }
    //check that the read is == the refToread
    if(b.refID!=refToRead){
       ret.bufferedRead = b;
       rdObjRegionDone =1;
       return;
    }

    if(b.nCig!=0){
      ret.first[ret.l] = b.pos;
      ret.last[ret.l] = bam_calend(b,getCig(&b));
      ret.l++;
      nReads--;
    }
  }


  /*
    now do the general loop,
    1) read an alignemnt calulate start and end pos
    check
    a) if same chromosome that a new read has a geq FINE
    b) if same but < then UNSORTED
    
  */ 
  int i;
  int tmp;
  int buffed=0;
  for( i=0;i<nReads;i++){
    aRead &b = ret.reads[i+ret.l];

    if((tmp=bam_iter_read1(fp,it,b,hdr))<0){
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }if(b.nCig==0){//only get reads with CIGAR
      i--;
      continue;
    }
    //general fine case
    if((refToRead==b.refID)&&( b.pos >= ret.first[ret.l+i-1])){
      ret.first[ret.l+i] = b.pos;
      ret.last[ret.l+i] = bam_calend(b,getCig(&b));
    }else if(b.refID>refToRead){
      buffed=1;
      ret.bufferedRead = b;
      //      fprintf(stderr,"[%s] new chromosome detected will temporarliy stop reading at read # :%d\n",__FUNCTION__,i);
      rdObjRegionDone =1;
      break;
    }else if(ret.first[ret.l+i-1]>b.pos){
      // fprintf(stderr,"unsorted file detected will exit new=%d new(-1)=%d\n",ret.first[ret.l+i-1],b.pos);
      fprintf(stderr,"[%s] unsorted file detected will exit new(-1)=%d new=%d,ret.l=%d i=%d\n",__FUNCTION__,ret.first[ret.l+i-1],b.pos,ret.l,i);
      exit(0);
    }
    
  }


  if(i!=nReads){
    for(int s=ret.l+i+buffed;s<ret.l+nReads;s++)
      delete [] ret.reads[s].vDat;
  }

  ret.nReads =ret.l+i;
  ret.l = ret.nReads;
}


//function will read data from all bamfiles, return value is the number of 'done' files
int collect_reads(bufReader *rd,int nFiles,int &notDone,sglPool *ret,int &readNlines,int ref,int &pickStop) {
#if 0
  fprintf(stderr,"\t[%s] Reading from referenceID=%d\n",__FUNCTION__,ref);
#endif

  int usedPicker=-1;
  for(int ii=0;ii<nFiles;ii++) {
    extern int *bamSortedIds;
    int i=bamSortedIds[ii];
    if(rd[i].regionDone||rd[i].isEOF){//if file is DONE with region go to next
      rd[i].regionDone = 1;
      continue;
    }
    
    int pre=ret[i].l;//number of elements before
    //function reads readNlines reads, from rd[i].fp and modifies isEOF and regionDone
    read_reads_noStop(rd[i].fp,readNlines,notDone,ret[i],ref,&rd[i].it,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
    
    //first check if reading caused and end of region event to occur
    if(rd[i].regionDone||rd[i].isEOF) 
      rd[i].regionDone = 1;
    
    if(ret[i].l>pre){//we could read data
      usedPicker = i;
      pickStop = ret[i].first[ret[i].l-1];
      break;
    }
  }
#if 0
  fprintf(stderr,"usedPicker:%d\n",usedPicker);
  exit(0);
#endif
 
  if(usedPicker==-1)  //<-this means we are done with the current chr/region
    return nFiles;

  //at this point we should have picked a picker. Thats a position=pickstop from a fileID=usedPicker that will be used as 'stopping point'
  for(int i=0;i<nFiles;i++) {
    if(rd[i].isEOF || rd[i].regionDone||i==usedPicker)
      continue;

#if 0
    fprintf(stderr,"i=%d regdone=%d\tiseof=%d\n",i,rd[i].regionDone,rd[i].isEOF);
    fprintf(stderr," getpool on i=%d\n",i);
#endif
    if(ret[i].l>0&&ret[i].first[ret[i].l-1]>pickStop)
      continue;
    read_reads_usingStop(rd[i].fp,readNlines,notDone,ret[i],ref,&rd[i].it,pickStop,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
  } 
  int nDone =0;
  for(int i=0;i<nFiles;i++)
    if(rd[i].regionDone)
      nDone++;
  
  return nDone;

}


/*
  What does this do?
  returns the number of basepairs covered by np and sgl
*/
int coverage_in_bp(nodePool *np, sglPool *sgl){
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

int coverage_in_bpT(nodePoolT *np, sglPool *sgl){
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
 
nodePool mkNodes_one_sample(sglPool *sgl,nodePool *np,abcGetFasta *gf) {
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
    aRead rd = sgl->reads[r];
    //    fprintf(stderr,"r=%d\tpos=%d\n",r,rd.pos);
    if(sgl->first[r] > last){
      int diffs = (rd.pos-last);
      offs = offs + diffs;
      last = sgl->last[r];
    }else
      last = std::max(sgl->last[r],last);

    char *seq =(char *) getSeq(&rd);
    char *quals =(char *) getQuals(&rd);
    int nCig = rd.nCig;

    uint32_t *cigs = getCig(&rd);
    int seq_pos =0; //position within sequence
    int wpos = rd.pos-offs;//this value is the current position assocatied with the positions at seq_pos
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
	  if(rd.mapQ!=255)
	    kputc(rd.mapQ+33, &tmpNode->seq);
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
	    char c = bam_nt16_rev_table[bam1_seqi(seq, seq_pos)];
	    kputc(bam1_strand(&rd)? tolower(c) : toupper(c), &tmpNode->seq);
	    kputw(seq_pos+1,&tmpNode->pos);
	    kputc(',',&tmpNode->pos);
	    seq_pos++;
	  }
	  wpos++;
	}else {//this is the deletion part
	  if(i!=0){
	    if(gf->ref==NULL)
	      for(int ii=0;ii<opLen;ii++)
		kputc(bam1_strand(&rd)? tolower('N') : toupper('N'),&tmpNode->seq);
	    else
	      for(int ii=0;ii<opLen;ii++)
		kputc(bam1_strand(&rd)? tolower(gf->ref->seqs[offs+wpos+ii+1]) : toupper(gf->ref->seqs[offs+wpos+ii+1]),&tmpNode->seq);
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
	  if(rd.mapQ!=255)
	    kputc(rd.mapQ+33, &tmpNode->seq);
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
	    if(rd.mapQ!=255)
	      kputc(rd.mapQ+33, &tmpNode->seq);
	    else
	      kputc('~', &tmpNode->seq);
	  }
	  char c = bam_nt16_rev_table[bam1_seqi(seq, seq_pos)];

	  if(gf->ref==NULL ||gf->ref->chrLen<wpos+offs)//prints the oberved allele
	    kputc(bam1_strand(&rd)? tolower(c) : toupper(c), &tmpNode->seq);
	  else{
	    if(refToInt[c]==refToInt[gf->ref->seqs[wpos+offs]])
	      kputc(bam1_strand(&rd)? ',' : '.', &tmpNode->seq);
	    else
	      kputc(bam1_strand(&rd)? tolower(c) : toupper(c), &tmpNode->seq);
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
	    bam1_strand(&rd)?kputc('<',&tmpNode->seq):kputc('>',&tmpNode->seq);
	    kputw(seq_pos+1, &tmpNode->pos);
	    kputc(quals[seq_pos]+33, &tmpNode->qs);
	  }
      }else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP) {
	//dont care
      }else{
	fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
	printErr();//unknown CIGAR
      }
    }
    //after end of read/parsing CIGAR always put the endline char
    //  fprintf(stderr,"printing endpileup for pos=%d\n",rd->pos);
    if(nCig!=0)
      kputc('$', &tmpNode->seq);
    delete [] rd.vDat;
  }

  //plug the reads back up //FIXME maybe do list type instead

  int miss= sgl->l-r;
  for(int i=0;i<(sgl->l-r);i++){
    sgl->reads[i] =sgl->reads[i+r];
    sgl->first[i] =sgl->first[i+r];
    sgl->last[i] =sgl->last[i+r];
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


nodePoolT mkNodes_one_sampleT(sglPool *sgl,nodePoolT *np) {
  int regionLen = coverage_in_bpT(np,sgl);//true covered regions
  nodePoolT dn;
  dn.l =0;
  if(regionLen==0)
    return dn;
  int lastSecure = sgl->lowestStart;
  
  dn = allocNodePoolT(regionLen);
  nodePoolT *ret = &dn;
  tNode *nds = ret->nds;
  int offs = np->first;//first position from buffered
  int last = np->last+1;//because we are comparing with calc_end

  //plug in the old buffered nodes
  for(int i=0;i<np->l;i++)
    nds[np->nds[i].refPos-offs] = np->nds[i];

  //initialize new nodes for the rest of the region
  for(int i=np->l;i<regionLen;i++)
    nds[i] = initNodeT(UPPILE_LEN);

  if(np->l==0){//if we didn't have have anybuffered then
    offs = sgl->first[0];
    last = sgl->last[0];
  }

  //parse all reads
  int r;
  
  for( r=0;r<sgl->readIDstop;r++) {

    aRead rd = sgl->reads[r];
    int mapQ = rd.mapQ;
    if(mapQ>=255) mapQ = 20;

    if(sgl->first[r] > last){
      int diffs = (rd.pos-last);
      offs = offs + diffs;
      last = sgl->last[r];
    }else
      last = std::max(sgl->last[r],last);

    char *seq =(char *) getSeq(&rd);
    char *quals =(char *) getQuals(&rd);
    int nCig = rd.nCig;

    uint32_t *cigs = getCig(&rd);
    int seq_pos =0; //position within sequence
    int wpos = rd.pos-offs;//this value is the current position assocatied with the positions at seq_pos
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
	  tmpNode = &nds[wpos];  
	  if(tmpNode->l2 >= tmpNode->m2){
	    tmpNode->m2++;
	    kroundup32(tmpNode->m2);
	    tNode *dddd =(tNode *) malloc(tmpNode->m2*sizeof(tNode));
	    for(int ddddd=0;ddddd<tmpNode->l2;ddddd++)
	      dddd[ddddd] = tmpNode->insert[ddddd];
	    free(tmpNode->insert);
	    tmpNode->insert = dddd;
	  }
	  
	  tmpNode->insert[tmpNode->l2] = initNodeT(opLen);
	  for(int ii=0;ii<opLen;ii++){
	    char c = bam_nt16_rev_table[bam1_seqi(seq, seq_pos)];
	    tmpNode->insert[tmpNode->l2].seq[ii] = bam1_strand(&rd)? tolower(c) : toupper(c);
	    tmpNode->insert[tmpNode->l2].posi[ii] = seq_pos + 1;
	    tmpNode->insert[tmpNode->l2].isop[ii] =rd.l_seq- seq_pos - 1;
	    tmpNode->insert[tmpNode->l2].qs[ii] = quals[seq_pos];
	    if( quals[seq_pos]<minQ || seq_pos + 1 < trim || rd.l_seq- seq_pos - 1 < trim)  
	      tmpNode->insert[tmpNode->l2].seq[ii] =  bam1_strand(&rd)? tolower('n') : toupper('N');
	    tmpNode->insert[tmpNode->l2].mapQ[ii] = mapQ;
	    seq_pos++;// <- important, must be after macro
	  }
	  tmpNode->l2++;//incrementor!
	  wpos++;
	}else {//this is the deletion part
	  
	  for(int ii=0;ii<opLen;ii++){
	    tmpNode = &nds[wpos+ii];
	    tmpNode->refPos = wpos +ii + offs;
	    tmpNode->deletion ++;
	  }
	  wpos += opLen;
	}

      }else if(opCode==BAM_CSOFT_CLIP){
	//occurs only at the ends of the read
	if(seq_pos == 0){
	  //then we are at beginning of read and need to write mapQ
	  tmpNode = &nds[wpos];
	  seq_pos += opLen;
	}else//we are at the end of read, then break CIGAR loop
	  break;
      }else if(opCode==BAM_CMATCH||opCode==BAM_CEQUAL||opCode==BAM_CDIFF) {
	hasInfo++;
	for(int fix=wpos ;wpos<(fix+opLen) ;wpos++) {
	  tmpNode =  &nds[wpos];
	  tmpNode->refPos=wpos+offs;
	  if(tmpNode->l>=tmpNode->m){
	    tmpNode->m = tmpNode->m*2;
	    tmpNode->seq =(char *) realloc(tmpNode->seq,tmpNode->m);
	    tmpNode->qs =(char *) realloc(tmpNode->qs,tmpNode->m);
	    tmpNode->posi =(unsigned char *) realloc(tmpNode->posi,tmpNode->m);
	    tmpNode->isop =(unsigned char *) realloc(tmpNode->isop,tmpNode->m);
	    tmpNode->mapQ = (unsigned char *) realloc(tmpNode->mapQ,tmpNode->m);
	  }
	  

	  char c = bam_nt16_rev_table[bam1_seqi(seq, seq_pos)];
	  tmpNode->seq[tmpNode->l] = bam1_strand(&rd)? tolower(c) : toupper(c);
	  tmpNode->qs[tmpNode->l] =  quals[seq_pos];
	
	  tmpNode->posi[tmpNode->l] = seq_pos;
	  tmpNode->isop[tmpNode->l] =rd.l_seq- seq_pos-1;
	  tmpNode->mapQ[tmpNode->l] = mapQ;
	  if( quals[seq_pos]<minQ || seq_pos  < trim || rd.l_seq - seq_pos - 1 < trim )
	    tmpNode->seq[tmpNode->l] = bam1_strand(&rd)? tolower('n') : toupper('N');
	  
	  tmpNode->l++;
	  seq_pos++;
	}
      }else if(opCode==BAM_CREF_SKIP) {
	  for(int fix=wpos;wpos<fix+opLen;wpos++){
	    tmpNode = &nds[wpos];
	    tmpNode->refPos=wpos+offs;
	  }
      }else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP) {
	//dont care
      }else{
	fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
	printErr();//unknown CIGAR
      }
      //      exit(0);
    }
    //after end of read/parsing CIGAR always put the endline char
    delete [] rd.vDat;
  }

  //plug the reads back up //FIXME maybe do list type instead

  int miss= sgl->l-r;
  for(int i=0;i<(sgl->l-r);i++){
    sgl->reads[i] =sgl->reads[i+r];
    sgl->first[i] =sgl->first[i+r];
    sgl->last[i] =sgl->last[i+r];
  }
  sgl->l = miss;
  //copy the part not meant for printing in this round into the buffer



  int lastSecureIndex =regionLen;
  //fprintf(stderr,"lastSecureIndex=%d\tregionLen=%d\t ret->nds.refpos=%d\n",lastSecureIndex,regionLen,ret->nds[regionLen-1].refPos);
  for(int i=0;0&&i<regionLen;i++)
    fprintf(stderr,"TYTYYTTYYTYT: refpos[%d]=%d seq.l=%d\n",i,ret->nds[i].refPos,ret->nds[i].l);
  int tailPos = ret->nds[regionLen-1].refPos;
  if(tailPos>lastSecure)
    lastSecureIndex = regionLen-tailPos+lastSecure;
  ret->l = lastSecureIndex;
  
  if(regionLen-lastSecureIndex+4>np->m){
    delete [] np->nds;
    np->m=regionLen+4;
    kroundup32(np->m);
    np->nds = new tNode[np->m];
  }
  assert(regionLen-lastSecureIndex+4<=np->m);
  //  fprintf(stderr,"np->m:%d\n",np->m)
  np->l=0;
  for(int i=lastSecureIndex;i<regionLen;i++)
    np->nds[np->l++] = nds[i];


  if(np->l!=0){
    np->first = np->nds[0].refPos;
    np->last = np->nds[np->l-1].refPos;
  }

  if(ret->l!=0){
    ret->first = ret->nds[0].refPos;
    ret->last = ret->nds[ret->l-1].refPos;
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


typedef std::map<int,tNode *> umapT; 

chunkyT *slow_mergeAllNodes_new(nodePoolT *dn,int nFiles){
  //  fprintf(stderr,"starting [%s]\n",__FUNCTION__);
  
  umapT ret;
  umapT::iterator it;
  for(int f=0;f<nFiles;f++) {
    nodePoolT sm = dn[f];
    for(int l=0;l<sm.l;l++) {
      tNode *perSite = NULL;
      int thepos =sm.nds[l].refPos;
      it=ret.find(thepos);
      if(it!=ret.end())
	perSite = it->second;
      else{
	perSite = new tNode[nFiles];
	for(int ii=0;ii<nFiles;ii++){
	  perSite[ii].l=perSite[ii].l2=perSite[ii].m=0;
	  perSite[ii].insert =NULL;
	  perSite[ii].refPos = thepos;
	}

	
	ret.insert(std::make_pair(thepos,perSite));
		
      }
      perSite[f] = sm.nds[l];
    }

  }
  //now we perform the backrool
  int nnodes = ret.size();
  int *refPos = new int [ret.size()];
  //  fprintf(stderr,"nnodes=%d\n",nnodes);
  chunkyT *chk =new chunkyT;  
  chk->nd = new tNode*[nnodes];
  int p=0;
  for(it = ret.begin();it!=ret.end();++it){
    chk->nd[p++] = it->second;
    refPos[p-1] = chk->nd[p-1][0].refPos;
  }

  chk->nSamples=nFiles;
  chk->nSites=nnodes;
  chk->refPos = refPos;
  return chk;  
}


chunkyT *mergeAllNodes_new(nodePoolT *dn,int nFiles) {

  int *depth = NULL;
  int *refPos2 = NULL;
  tNode **super = NULL;
  
  if(dn==NULL){//cleanup called after end of looping through all files.
    return NULL;
  }
  int first,last;
  get_span_all_samplesT(dn,nFiles,first,last);

  int rlen = last-first+1;
  assert(rlen>=0);
  if(rlen>BUG_THRES)
    return slow_mergeAllNodes_new(dn,nFiles);

  super = new tNode*[rlen];
  depth = new int[rlen];
  refPos2 = new int[rlen];
  
  memset(depth,0,rlen*sizeof(int));
  int offs = first;
  //looping through the different samples
  for(int n=0;n<nFiles;n++) {
    nodePoolT sm = dn[n];
    int i;
    //looping through all nodes
    for( i=0;((i<sm.l)&&(sm.nds[i].refPos <= std::min(sm.last,last) ));i++) {
      
      int posi = sm.nds[i].refPos-offs;
      if(depth[posi]==0){
	//	fprintf(stderr,"posi=%d\n",posi);
	super[posi] = new tNode[nFiles];
	for(int ii=0;ii<nFiles;ii++){
	  super[posi][ii].l=super[posi][ii].l2=super[posi][ii].m=0;
	  super[posi][ii].insert =NULL;
	  super[posi][ii].refPos = posi+offs;
	  //super[posi][ii] = initNodeT(UPPILE_LEN,posi+offs,posi);
	}
	
	refPos2[posi]=sm.nds[i].refPos;
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
  chk->nd = new tNode*[nnodes];
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

pthread_mutex_t mUpPile_mutex = PTHREAD_MUTEX_INITIALIZER;


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



sglPool makePool(int l){
  sglPool ret;
  ret.l=0;
  ret.m=l;
  kroundup32(ret.m);
  if(0){
    ret.reads=(aRead *)malloc(ret.m*sizeof(aRead));
    ret.first=(int *)malloc(ret.m*sizeof(int));
    ret.last=(int *)malloc(ret.m*sizeof(int));
  }else{
    ret.reads= new aRead[ret.m];
    ret.first=new int[ret.m];
    ret.last=new int[ret.m];

  }

  ret.bufferedRead.refID = -2;
  //  fprintf(stderr,"allocing pool with l=%d and m=%d\n",l,ret.m);
  return ret;
}



int getSglStop5(sglPool *sglp,int nFiles,int pickStop) {
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
void getMaxMax2(sglPool *sglp,int nFiles,nodePool *nps){
#if 0
    fprintf(stderr,"[%s].nFiles=%d\n",__FUNCTION__,nFiles);
    for(int i=0;i<nFiles;i++){
      fprintf(stderr,"%d=sglp.l=%d\n",i,sglp[i].l);
    }
#endif
    
    int last_bp_in_chr = -1;
    
    for(int i=0;i<nFiles;i++){
      sglPool *tmp = &sglp[i];
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

void getMaxMax2T(sglPool *sglp,int nFiles,nodePoolT *nps){
#if 0
  fprintf(stderr,"[%s].nFiles=%d\n",__FUNCTION__,nFiles);
  for(int i=0;i<nFiles;i++)
    fprintf(stderr,"%d=sglp.l=%d\n",i,sglp[i].l);

#endif

  int last_bp_in_chr = -1;

  for(int i=0;i<nFiles;i++){
    sglPool *tmp = &sglp[i];
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


void getOffsets(htsFile *fp,char *fn,const bam_hdr_t *hd,iter_t &iter,int ref,int start,int stop,bam_hdr_t *hdr){
  if(iter.hts_idx==NULL)
    iter.hts_idx = sam_index_load(fp,fn);
  if (iter.hts_idx == 0) { // index is unavailable
    fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
    exit(0);
  }
  char tmp[1024];
  snprintf(tmp,1024,"%s:%d-%d",hd->target_name[ref],start+1,stop);
  if(iter.hts_itr)
    hts_itr_destroy(iter.hts_itr);
  iter.hts_itr = sam_itr_querys(iter.hts_idx, hdr, tmp);
  if (iter.hts_itr == NULL) { // reference name is not found
    fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n",tmp);
    exit(0);
  }
  //  fprintf(stderr,"HEJSA MED IDG%p\n",iter.hts_itr);
  hts_itr_t *idx=iter.hts_itr;
}



void *setIterator1(void *args){
  bufReader *rd = (bufReader *) args;
  // free(rd->it.off);
  //  fprintf(stderr,"[%s]\n",__FUNCTION__,rd->regions[rd->itrPos]);
  int start,stop,ref;
  start = rd->regions.start;
  stop = rd->regions.stop;
  ref = rd->regions.refID;
  getOffsets(rd->fp,rd->fn,rd->hdr,rd->it,ref,start,stop,rd->hdr);//leak
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
  getOffsets(rd->fp,rd->fn,rd->hdr,rd->it,ref,start,stop,rd->hdr);
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
    
#if 0
    threadpool *tp = init_threadpool(nThreads,setIterator1_thread,nFiles,1);
    threadpool_iter(tp,rd,nFiles);
    wait_kill(tp);
    tp = NULL;
#endif
  }
}


void printChunky2(const chunky* chk,FILE *fp,char *refStr); //<-bammer_main
void callBack_bambi(fcb *fff);//<-angsd.cpp
//type=1 -> samtool mpileup textoutput
//type=0 -> callback to angsd

//function below is a merge between the original TESTOUTPUT and the angsdcallback. Typenames with T are the ones for the callback
//Most likely this can be written more beautifull by using templated types. But to lazy now.
int uppile(int show,int nThreads,bufReader *rd,int nLines,int nFiles,std::vector<regs> regions) {
  assert(nLines&&nFiles);
  fprintf(stderr,"\t-> Parsing %d number of samples \n",nFiles);
 fflush(stderr);
  if(show!=1)
    pthread_mutex_lock(&mUpPile_mutex);//just to make sure, its okey to clean up
  extern abcGetFasta *gf;

  sglPool *sglp= new sglPool[nFiles];
  
  nodePool *nps =NULL;
  nodePoolT *npsT = NULL;
  if(show)
    nps = new nodePool[nFiles];// <- buffered nodes
  else
    npsT = new nodePoolT[nFiles];// <- buffered nodes

  for(int i=0;i<nFiles;i++){
    sglp[i] = makePool(nLines);
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
	//fprintf(stderr,"region lookup %d/%lu\n",itrPos+1,regions.size());
	//fflush(stderr);
	setIterators(rd,regions[itrPos],nFiles,nThreads);
	//fprintf(stderr,"done region lookup %d/%lu\n",itrPos+1,regions.size());
	//fflush(stderr);///BAME

	theRef = rd[0].it.hts_itr->tid; 
	//	fprintf(stderr,"theRef:%d %p %p\n",theRef,rd[0].it.off,rd[0].it.hts_itr->off);
	//validate we have offsets;
	int gotData = 0;
	for(int i=0;i<nFiles;i++)
	  if((rd[i].it.hts_itr &&(rd[i].it.hts_itr->off!=NULL))){
	    gotData =1;
	    break;
	  }
	//reset eof and region donepointers.
	for(int i=0;i<nFiles;i++){
	  rd[i].isEOF = 0;
	  rd[i].regionDone =0;
	  sglp[i].l =0;
	  sglp[i].bufferedRead.refID=-2;
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
	  if(sglp[i].bufferedRead.refID!=-2 && sglp[i].bufferedRead.refID<minRef)
	    minRef = sglp[i].bufferedRead.refID;
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
	if(sglp[i].bufferedRead.refID==-2)//<- no buffered no nothing
	  continue;
	if(sglp[i].bufferedRead.refID==theRef) {//buffered is on correct chr
	  int doCpy =1;
	  if(gf->ref!=NULL){
	    //this will modify the quality and mapQ if -baq >0 or adjustMapQ= SAMtools -c
	    doCpy = restuff(sglp[i].bufferedRead);
	  }
	  if(doCpy){
	    sglp[i].reads[sglp[i].l] = sglp[i].bufferedRead;
	    sglp[i].first[sglp[i].l] = sglp[i].reads[sglp[i].l].pos;
	    sglp[i].last[sglp[i].l] =   bam_calend(sglp[i].reads[sglp[i].l],getCig(&sglp[i].reads[sglp[i].l]));
	    sglp[i].l++;
	  }
	  //if we haven't performed the copy its because adjustedMap<minMapQ, in either case we don't need the read anymore
	  sglp[i].bufferedRead.refID = -2;
	}else // buffered was not on correct chr say that the 'region' is done
	  rd[i].regionDone=1;
      }
    }

    //below loop will continue untill eoc/oef/eor <- funky abbrev =end of chr, file and region
    while(notDone&& SIG_COND) {

      //before collecting reads from the files, lets first check if we should pop reads from the buffered queue in each sgl
      for(int i=0;i<nFiles;i++){
	if(sglp[i].bufferedRead.refID==theRef){
	  sglp[i].reads[sglp[i].l] = sglp[i].bufferedRead;
	  sglp[i].first[sglp[i].l] = sglp[i].reads[sglp[i].l].pos;
	  sglp[i].last[sglp[i].l] = bam_calend(sglp[i].reads[sglp[i].l],getCig(&sglp[i].reads[sglp[i].l]));
	  sglp[i].l++;
	  sglp[i].bufferedRead.refID = -2;
	}
	
      }
      
      int pickStop=-1;
      
      int doFlush = (collect_reads(rd,nFiles,notDone,sglp,nLines,theRef,pickStop)==nFiles)?1:0 ;
#if 0
      for(int i=0;i<nFiles;i++)
	fprintf(stderr,"[%s] sgl[%d].l=%d (%d,%d)\n",__FUNCTION__,i,sglp[i].l,sglp[i].first[0],sglp[i].last[sglp[i].l-1]);
      
#endif
      
      //prepare uppiling for all data, if end of chromosome or all files are EOF
      if(doFlush||notDone==0){
	//fprintf(stderr,"[%s]Last chunk Will flush all remaining reads\n",__FUNCTION__);
	if(show)
	  getMaxMax2(sglp,nFiles,nps);
	else
	  getMaxMax2T(sglp,nFiles,npsT);
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
	dnT = new nodePoolT[nFiles];//<- this can leak now
	for(int i=0;i<nFiles;i++){
	  dnT[i] = mkNodes_one_sampleT(&sglp[i],&npsT[i]);
	  tmpSum += dnT[i].l;
	}
      }else
	for(int i=0;i<nFiles;i++) {
	  dn[i] = mkNodes_one_sample(&sglp[i],&nps[i],gf);
	  tmpSum += dn[i].l;
	}

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
	printChunky2(chk,stdout,chk->refName);
	
	for(int i=0;1&&i<nFiles;i++){
	  if(dn[i].l!=0)
	    delete [] dn[i].nds;
	}
      }else{

	fcb *f = new fcb; //<-for call back
	f->dn=dnT; f->nFiles = nFiles; f->regStart = regStart; f->regStop = regStop; f->refId = theRef;
	callBack_bambi(f);
      }

      //if we flush, its due to end of chr/region or end of files
      if(doFlush)
	break;
    }

  }
  //below can be written nice but is just a copy paste, time is 10pm to lazy now
  if(show==0){
    callBack_bambi(NULL);//cleanup signal
    pthread_mutex_lock(&mUpPile_mutex);//just to make sure, its okey to clean up
    for(int i=0;1&&i<nFiles;i++){
      dalloc_nodePoolT(npsT[i]);
      delete [] npsT[i].nds;
      dalloc(sglp[i]);
    }
    delete [] npsT;
    delete [] sglp;
    pthread_mutex_unlock(&mUpPile_mutex);//just to make sure, its okey to clean up
  }else{
    //clean up
    for(int i=0;1&&i<nFiles;i++){
      dalloc_nodePool(nps[i]);
      delete [] nps[i].nds;
      dalloc(sglp[i]);
    }

    delete [] nps;
    delete [] sglp;
  }
  return 0;
}

