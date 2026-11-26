#include "makeReadPool.h"
#include "aio.h"
#define bam_dup(b) bam_copy1(bam_init1(), (b))
readPool makePoolb(int l){
  readPool ret;
  ret.l=0;
  ret.m=l;
  kroundup32(ret.m);
  
  ret.reads=(bam1_t**)malloc(ret.m*sizeof(bam1_t*));// new bam1_t*[ret.m];
  for(int i=0;i<ret.m;i++)
    ret.reads[i] = bam_init1();
  ret.first=(int *)malloc(ret.m*sizeof(int));//new int[ret.m];
  ret.last=(int *)malloc(ret.m*sizeof(int));//new int[ret.m];
  ret.bufferedRead=NULL;
  ret.lastTid = ret.lastPos = -1;

  return ret;
}

//pop1_read, and without -r/-rf check that the reads are sorted by position. Unsorted reads were
//put on the wrong chromosome or dropped, or angsd stopped with 'ReferenceID ... is not' or crashed (#377)
//With -r/-rf the file is indexed, so sorted, and the regions may go back on purpose
static int read1(htsFile *fp,hts_itr_t *itr,bam1_t *b,bam_hdr_t *hdr,readPool &ret){
  int r = pop1_read(fp,itr,b,hdr);
  if(r<0||itr!=NULL)
    return r;
  if(b->core.tid<ret.lastTid||(b->core.tid==ret.lastTid&&b->core.pos<ret.lastPos)){
    fprintf(stderr,"\t-> \'%s\' is not sorted by position: read \'%s\' at %s:%lld comes after a read at %s:%d\n",fp->fn,bam_get_qname(b),hdr->target_name[b->core.tid],(long long)b->core.pos+1,hdr->target_name[ret.lastTid],ret.lastPos+1);
    fprintf(stderr,"\t-> Sort it with \'samtools sort\'\n");
    exit(1);
  }
  ret.lastTid = b->core.tid;
  ret.lastPos = b->core.pos;
  return r;
}

void dalloc (readPool *ret){
  for(int i=0;i<ret->m;i++)
    bam_destroy1(ret->reads[i]);
  free(ret->reads);
  free(ret->first);
  free(ret->last);
  //not needed: bufferedRead is freed when changing region and is NULL at the end of the files,
  //so it should always be NULL here, we could assert that instead
  if(ret->bufferedRead)
    bam_destroy1(ret->bufferedRead);
  ret->bufferedRead = NULL;
}

void realloc(readPool *ret,int l){
  int old=ret->m;
  ret->m =l;
  kroundup32(ret->m);
  ret->reads =(bam1_t**) realloc(ret->reads,sizeof(bam1_t*)*ret->m);
  for(int i=old;i<ret->m;i++)
    ret->reads[i]= bam_init1();
  ret->first =(int*) realloc(ret->first,sizeof(int)*ret->m);
  ret->last =(int*) realloc(ret->last,sizeof(int)*ret->m);
}

void read_reads_usingStop(htsFile *fp,int nReads,int &isEof,readPool &ret,int refToRead,hts_itr_t *itr,int stop,int &rdObjEof,int &rdObjRegionDone,bam_hdr_t *hdr) {

  //if should never be in this function if we shouldnt read from the file.
  aio::doAssert(rdObjRegionDone!=1 &&rdObjEof!=1 ,1,AT,"");
  
  if((nReads+ret.l)>ret.m)
    realloc(&ret,nReads+ret.l);


  //this is the awkward case this could cause an error in some very unlikely scenario.
  //whith the current buffer empty and the first read being a new chromosome. 
  //this is fixed good job monkey boy

  while(ret.l==0){
    bam1_t *b = ret.reads[0];
    int tmp;

    if((tmp=read1(fp,itr,b,hdr,ret))<0){
      if(tmp==-1){
	rdObjEof =1;
	//	ret.isEOF =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      //return, the loop below would read again and count the end of the file twice. Then the chromosome/region
      //was stopped while the other files still had reads, when a file had no reads (in the region) or all were filtered
      return;
    }
    if(b->core.tid>refToRead){
      //buffed=1;
      ret.bufferedRead = bam_dup(b);
      rdObjRegionDone =1;
      return;
    }
    
    ret.first[ret.l] = b->core.pos;
    ret.last[ret.l] = bam_endpos(b);
    ret.l++;
    nReads--;//breaking automaticly due to ret.l++
    
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
  //  int buffed=0;
  while(1) {
    if(i+ret.l>=(ret.m-1)){
      ret.l += i;
      i=0;
      realloc(&ret,ret.m+1);//will double buffer
    }
    
    bam1_t *b = ret.reads[i+ret.l];
    if((tmp=read1(fp,itr,b,hdr,ret))<0) {
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }
 
    //general fine case
    
    if((refToRead==b->core.tid)&&( b->core.pos >= ret.first[ret.l+i-1])&&(b->core.pos<stop)){
      ret.first[ret.l+i] = b->core.pos;
      ret.last[ret.l+i] = bam_endpos(b);
    }else if((refToRead==b->core.tid)&&b->core.pos>=stop){
      //buffed=1;
      ret.bufferedRead = bam_dup(b);
      break;
    }
    else if(b->core.tid>refToRead){
      //buffed=1;
      ret.bufferedRead = bam_dup(b);
      rdObjRegionDone =1;
      break;
    }else if(ret.first[ret.l+i-1]>b->core.pos){
      fprintf(stderr,"\t-> \'%s\' is not sorted by position, sort it with \'samtools sort\'\n",fp->fn);
      exit(1);
    }
    i++;
  }
  ret.l += i;
}


void read_reads(htsFile *fp,int nReads,int &isEof,readPool &ret,int refToRead,hts_itr_t *itr,int stop,int &rdObjEof,int &rdObjRegionDone,bam_hdr_t *hdr) {
  // fprintf(stderr,"stop:%d nreads:%d reftoread:%d ret.m:%d ret.l:%d\n",stop,nReads,refToRead,ret.m,ret.l);
  //if should never be in this function if we shouldnt read from the file.
  aio::doAssert(rdObjRegionDone!=1 &&rdObjEof!=1 ,1,AT,"");
  

  /*
      1) read an alignemnt calulate start and end pos
    check
    a) if same chromosome that a new read has a geq FINE
    b) if same but < then UNSORTED
    
  */ 
  int i=0;
  int tmp;

  while(1) {
  
    if(i+ret.l>=(ret.m-4)){//minus four will ofcourse work.
      realloc(&ret,ret.m+ret.l+i);
    }
    //fprintf(stderr,"ret.l:%d,i:%d\n",ret.l,i);
    bam1_t *b = ret.reads[i+ret.l];
    if((tmp=read1(fp,itr,b,hdr,ret))<0) {
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }

    //general fine case
    
    if(refToRead==b->core.tid){
      ret.first[ret.l+i] = b->core.pos;
      ret.last[ret.l+i] = bam_endpos(b);
      //      fprintf(stderr,"(%d,%d)\n",ret.first[ret.l+i],ret.last[ret.l+i]);
    }else if(b->core.tid!=refToRead){
      //buffed=1;
      ret.bufferedRead = bam_dup(b);
      rdObjRegionDone =1;
      break;
    }

    i++;
    //fprintf(stderr,"i:%d first:%d\n",i,ret.first[ret.l+i-1]);
    if(i>nReads&&ret.first[ret.l+i-1]>stop)
      break;
  }
  ret.l += i;
}


//nothing with buffered here
void read_reads_noStop(htsFile *fp,int nReads,int &isEof,readPool &ret,int refToRead,hts_itr_t *itr,int &rdObjEof,int &rdObjRegionDone,bam_hdr_t *hdr) {
#if 0
  fprintf(stderr,"\t->[%s] buffRefid=%d\trefToRead=%d\n",__FUNCTION__,ret.bufferedRead.refID,refToRead);
#endif
  aio::doAssert(rdObjEof==0 && ret.bufferedRead ==NULL,1,AT,"");
 
  if((nReads+ret.l)>ret.m)
    realloc(&ret,nReads+ret.l);
 
  //this is the awkward case this could cause an error in some very unlikely scenario.
  //whith the current buffer empty and the first read being a new chromosome. 
  //this is fixed good job monkey boy
  while(ret.l==0){
    bam1_t *b = ret.reads[0];
    int tmp;
    if((tmp=read1(fp,itr,b,hdr,ret))<0){//FIXME
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }
   
    //check that the read is == the refToread
    if(b->core.tid!=refToRead){
      ret.bufferedRead = bam_dup(b);
      rdObjRegionDone =1;
      return;
    }
    ret.first[ret.l] = b->core.pos;
    ret.last[ret.l] = bam_endpos(b);
    ret.l++;
    nReads--;
   
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
    bam1_t *b = ret.reads[i+ret.l];

    if((tmp=read1(fp,itr,b,hdr,ret))<0){
      if(tmp==-1){
	rdObjEof =1;
	isEof--;
      }else
	rdObjRegionDone =1;
      break;
    }

    //general fine case
    if((refToRead==b->core.tid)&&( b->core.pos >= ret.first[ret.l+i-1])){
      ret.first[ret.l+i] = b->core.pos;
      ret.last[ret.l+i] = bam_endpos(b);
    }else if(b->core.tid>refToRead){
      buffed=1;
      ret.bufferedRead = bam_dup(b);
      //      fprintf(stderr,"[%s] new chromosome detected will temporarliy stop reading at read # :%d\n",__FUNCTION__,i);
      rdObjRegionDone =1;
      break;
    }else if(ret.first[ret.l+i-1]>b->core.pos){
      fprintf(stderr,"\t-> \'%s\' is not sorted by position, sort it with \'samtools sort\'\n",fp->fn);
      exit(1);
    }
    
  }

  ret.l += i;

}



//function will read data from all bamfiles, return value is the number of 'done' files
int collect_reads(bufReader *rd,int nFiles,int &notDone,readPool *ret,int &readNlines,int ref,int &pickStop) {
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
    read_reads_noStop(rd[i].fp,readNlines,notDone,ret[i],ref,rd[i].itr,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
    //first check if reading caused and end of region event to occur
    if(rd[i].regionDone||rd[i].isEOF) 
      rd[i].regionDone = 1;
    
    if(ret[i].l>pre){//we could read data
      usedPicker = i;
      pickStop = ret[i].first[ret[i].l-1];
      break;
    }
  }
 
  if(usedPicker==-1)  //<-this means we are done with the current chr/region
    return nFiles;

  //at this point we should have picked a picker. Thats a position=pickstop from a fileID=usedPicker that will be used as 'stopping point'
  for(int i=0;i<nFiles;i++) {
    if(rd[i].isEOF || rd[i].regionDone||i==usedPicker)
      continue;
    
    if(ret[i].l>0&&ret[i].first[ret[i].l-1]>pickStop)
      continue;
    read_reads_usingStop(rd[i].fp,readNlines,notDone,ret[i],ref,rd[i].itr,pickStop,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
  } 
  int nDone =0;
  for(int i=0;i<nFiles;i++)
    if(rd[i].regionDone)
      nDone++;
  
  return nDone;

}

int collect_reads2(bufReader *rd,int nFiles,int &notDone,readPool *ret,int &readNlines,int ref,int &pickStop) {
  int usedPicker=-1;
  const int stopOffset = pickStop;//each picker candidate reads until stopOffset past its last read
  for(int ii=0;ii<nFiles;ii++) {
    //    fprintf(stderr,"read reads[%d]\n",ii);
    extern int *bamSortedIds;
    int i=bamSortedIds[ii];
    if(rd[i].regionDone||rd[i].isEOF){//if file is DONE with region go to next
      rd[i].regionDone = 1;
      continue;
    }
    
    int pre=ret[i].l;//number of elements before
    //function reads readNlines reads, from rd[i].fp and modifies isEOF and regionDone
    //from the offset, it was added to the stop of the files before that had no reads, so a file could read far ahead
    pickStop = stopOffset;
    if(ret[i].l>0)
      pickStop += ret[i].first[ret[i].l-1];
    // fprintf(stderr,"reading until:%d\n",pickStop);
    read_reads(rd[i].fp,readNlines,notDone,ret[i],ref,rd[i].itr,pickStop,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
    //first check if reading caused and end of region event to occur
    if(rd[i].regionDone||rd[i].isEOF) 
      rd[i].regionDone = 1;
    
    if(ret[i].l>pre){//we could read data
      usedPicker = i;
      pickStop = ret[i].first[ret[i].l-1];
      break;
    }
  }
  //  fprintf(stderr,"pickstop:%d\n",pickStop);
  if(usedPicker==-1)  //<-this means we are done with the current chr/region
    return nFiles;

  //at this point we should have picked a picker. Thats a position=pickstop from a fileID=usedPicker that will be used as 'stopping point'
  for(int i=0;i<nFiles;i++) {
    //    fprintf(stderr,"using stop[%d] \n",i);
    if(rd[i].isEOF || rd[i].regionDone||i==usedPicker)
      continue;
    
    if(ret[i].l>0&&ret[i].first[ret[i].l-1]>pickStop)
      continue;
    read_reads_usingStop(rd[i].fp,readNlines,notDone,ret[i],ref,rd[i].itr,pickStop,rd[i].isEOF,rd[i].regionDone,rd[i].hdr);
  } 
  int nDone =0;
  for(int i=0;i<nFiles;i++)
    if(rd[i].regionDone)
      nDone++;
  
  return nDone;

}

