//based on SAMtools code 1.17

#include <cstdio>
#include <stdint.h>
#include <map>
#include <cstring>
#include <cstdlib>

#include "bams.h"
#include <algorithm>
#include "argStruct.h"
//linear index per chromosome/reference

#define BAM_MAX_BIN 37450 // =(8^6-1)/7+1
#define BAM_LIDX_SHIFT    14


typedef struct{
  int32_t n_intv;
  uint64_t *ioffset;//length is n_intv
}lindex_t;

typedef struct{
  int32_t n_chunk;
  pair64_t *chunky;
}aBin;


//binned index per chromosome/reference
typedef struct{
  uint32_t n_bin;//nbin per chromosome
  std::map<uint32_t,aBin> myMap;//size is n_bin
}bindex_t;


struct __tindex{
  int32_t loadedRef;
  lindex_t lindex;//loadedRef long
  bindex_t bindex;//loadedRef long
};



void dalloc(tindex idx){
  if(idx==NULL)
    return;
  delete [] idx->lindex.ioffset;
  std::map<uint32_t,aBin> myMap = idx->bindex.myMap;
  for(std::map<uint32_t,aBin>::iterator it=myMap.begin();it!=myMap.end();++it)
    delete [] it->second.chunky;
  delete idx;
}



int parse_region(char *extra,const aHead *hd,int &ref,int &start,int &stop,const aMap *revMap) {
  aMap::const_iterator it;
   if(strrchr(extra,':')==NULL){//only chromosomename
     if((it = revMap->find(extra))==revMap->end()){
       fprintf(stderr,"[%s.%s():%d] Problems finding chromo:%s\n",__FILE__,__FUNCTION__,__LINE__,extra);
       fflush(stderr);
       exit(0);
       return -1;
     }
     ref = it->second;
     start =0;
     stop = hd->l_ref[ref];
     return 1;
   }

   char *tok = strtok(extra,":");
   if((it =revMap->find(tok))==revMap->end()){
       fprintf(stderr,"[%s.%s():%d] (-r) Problems finding chromo:%s\n",__FILE__,__FUNCTION__,__LINE__,extra);
       fflush(stderr);
       exit(0);
       return -1;
   }
   ref = it->second;

   start =0;
   stop = hd->l_ref[ref];
   tok = extra+strlen(tok)+1;//tok now contains the rest of the string

   if(strlen(tok)==0)//not start and/or stop ex: chr21:
     return 1;


   if(tok[0]=='-'){//only contains stop ex: chr21:-stop
     tok =strtok(tok,"-");

     stop = atoi(tok);
   }else{
     //catch single point
     int isProper =0;
     for(size_t i=0;i<strlen(tok);i++)
       if(tok[i]=='-'){
	 isProper=1;
	 break;
       }
     //fprintf(stderr,"isProper=%d\n",isProper);
     if(isProper){
       tok =strtok(tok,"-");
       start = atoi(tok)-1;//this is important for the zero offset
       tok = strtok(NULL,"-");
       if(tok!=NULL)
	 stop = atoi(tok);
     }else{
       //single point
       stop = atoi(tok);
       start =stop -1;

     }
     
   }
   if(stop<start){
     fprintf(stderr,"endpoint:%d is larger than startpoint:%d\n",start,stop);
     exit(0);
     
   }
   if(0){
     fprintf(stderr,"[%s] ref=%d,start=%d,stop=%d\n",__FUNCTION__,ref,start,stop);
     exit(0);
   }
   return 1;
 }




lindex_t getLindex(FILE *fp,int keep){
   lindex_t id;
   if(sizeof(int32_t)!=fread(&id.n_intv,1,sizeof(int32_t),fp))
     printErr();
   //   fprintf(stderr,"[%s] nlin index: %d\n",__FUNCTION__,id.n_intv);
   id.ioffset = NULL;
   if(id.n_intv!=0){
     id.ioffset = new uint64_t[id.n_intv];
     if(sizeof(uint64_t)!=fread(id.ioffset,id.n_intv,sizeof(uint64_t),fp))
       printErr();
   }
   if(!keep)
     delete [] id.ioffset;
   for(int ii=0;0&&ii<id.n_intv;ii++)
     fprintf(stderr,"iooffset[%d]=%lu\n",ii,(size_t)id.ioffset[ii]);
   return id;
 }

bindex_t getBindex(FILE *fp,int keep){
   bindex_t id;
   if(1!=fread(&id.n_bin,sizeof(int32_t),1,fp))//read number of bins
     printErr();
   // fprintf(stderr,"n_bin=%d\n",id.n_bin);

   for(uint32_t j=0;  j<id.n_bin;j++) {
     uint32_t bin_key;
     aBin abin;
     if(sizeof(uint32_t)!=fread(&bin_key,1,sizeof(uint32_t),fp))//read bin_key
       printErr();
     if(sizeof(int32_t)!=fread(&abin.n_chunk,1,sizeof(int32_t),fp))//read number of chunks associated with bin_key
       printErr();
     //fprintf(stderr,"bin[%d] of n_bin[%u] =%d,key=%u \n",j,id.n_bin,abin.n_chunk,bin_key);
     abin.chunky = new pair64_t[abin.n_chunk];

     for(int s=0; s < abin.n_chunk;s++) {
       if(sizeof(uint64_t)!=fread(&abin.chunky[s].chunk_beg,1,sizeof(uint64_t),fp))
	 printErr();
       if(sizeof(uint64_t)!=fread(&abin.chunky[s].chunk_end,1,sizeof(uint64_t),fp))
	 printErr();

       //fprintf(stderr," start:%lu-%lu\n",abin.chunky[s].chunk_beg,abin.chunky[s].chunk_end);
     }
     if(keep)
       id.myMap.insert(std::pair<uint32_t,aBin>(bin_key,abin));
     else
       delete [] abin.chunky;
   }
   return id;
 }

tindex getIdx(const char *fname,int chooseRef) {
  //  fprintf(stderr,"bainame is: = %s\n",fname);
   FILE *fp = NULL;
   if(NULL==(fp=fopen(fname,"r"))){
     fprintf(stderr,"\t-> Problems opening indexfile: %s\n",fname);
     exit(0);
   }
   char magic[4];
   if(4!=fread(magic,1,4,fp))
     printErr();
   //fprintf(stderr,"%s\n",magic);
   if(strncmp(magic,"BAI\1",4)!=0){
     fprintf(stderr,"bai file: %s looks corrupted will exit\n",fname);
     exit(0);
   }   
   tindex id = new __tindex;
   if(1!=fread(&id->loadedRef,sizeof(int32_t),1,fp))
     printErr();
   //   fprintf(stderr,"n_ref=%d\n",id->n_ref);
   // id->lindex = new lindex_t[id->n_ref];
   //id->bindex = new bindex_t[id->n_ref];
   for(int i=0;i<id->loadedRef;i++)  {//loop through each chr
     //read binned index
     id->bindex = getBindex(fp,chooseRef==i);
     //now read linear index
     id->lindex = getLindex(fp,chooseRef==i);
     if(i==chooseRef){
       id->loadedRef=chooseRef;
       break;
     }
   }
   if(fp) fclose(fp);
   return id;
 }




 static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[BAM_MAX_BIN]){
   int i = 0, k;
   if (beg >= end) return 0;
   if (end >= 1u<<29) end = 1u<<29;
   --end;
   list[i++] = 0;
   for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
   for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
   for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
   for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
   for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
   return i;
 }

 inline bool operator<(const pair64_t& a, const pair64_t& b)
 {
   return a.chunk_beg < b.chunk_beg;
 }

 


void buildOffsets(tindex id,int ref,int beg, int end,iter_t &iter) {
  //fprintf(stderr,"[%s] start: tid=%d, beg=%d, end=%d\n",__FUNCTION__,ref,beg,end);
   if(id==NULL)
     fprintf(stderr,"write msg in fun =%s\n",__FUNCTION__);


   uint16_t *bins;
   int i, n_bins, n_off;

   uint64_t min_off;
   //  iter_t iter;// = 0;

   iter.tid = ref, iter.beg = beg, iter.end = end; iter.i = -1;
   iter.from_first =0;iter.finished =0;iter.curr_off=0;iter.off=NULL;

   bindex_t bindex = id->bindex;
   lindex_t lindex = id->lindex;

   // initialize iter
   //iter = calloc(1, sizeof(struct __bam_iter_t));
   
   //

   //do linear index;
   if (lindex.n_intv > 0) {
     min_off = (beg>>BAM_LIDX_SHIFT >= lindex.n_intv)? lindex.ioffset[lindex.n_intv-1]
       : lindex.ioffset[beg>>BAM_LIDX_SHIFT];
     if (min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
       int n = beg>>BAM_LIDX_SHIFT;
       if (n > lindex.n_intv) 
	 n = lindex.n_intv;
       for (i = n - 1; i >= 0; --i)
	 if (lindex.ioffset[i] != 0) 
	   break;
       if (i >= 0) min_off = lindex.ioffset[i];
     }
   } else min_off = 0; // tabix 0.1.2 may produce such index files

   //do binned index
   bins = (uint16_t*)calloc(BAM_MAX_BIN, 2);
   n_bins = reg2bins(beg, end, bins);

   std::map<uint32_t,aBin>::iterator it;
   for (i = n_off = 0; i < n_bins; ++i) {
     it = bindex.myMap.find( bins[i]);
     if (it != bindex.myMap.end())
       n_off += it->second.n_chunk;
   }
   if (n_off == 0) {
     free(bins); return ;
   }

   pair64_t *off = (pair64_t*)calloc(n_off, 16);

   for (i = n_off = 0; i < n_bins; ++i) {
     it = bindex.myMap.find( bins[i]);
     if (it != bindex.myMap.end()) {
       for (int j = 0; j < it->second.n_chunk; ++j)
	 if (it->second.chunky[j].chunk_end > min_off) 
	   off[n_off++] = it->second.chunky[j];
     }
   }
   free(bins);
   if (n_off == 0) {
     free(off); return;
   }
   {

     int l;
     std::sort(off,off+n_off);

     // resolve completely contained adjacent blocks
     for (i = 1, l = 0; i < n_off; ++i)
       if (off[l].chunk_end < off[i].chunk_end)
	 off[++l] = off[i];
     n_off = l + 1;
     //          fprintf(stderr,"LLL: n_off=%d\n",n_off);
     
     // resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
     for (i = 1; i < n_off; ++i) //exists
       if (off[i-1].chunk_end >= off[i].chunk_beg) 
	 off[i-1].chunk_end = off[i].chunk_beg;//exits
     { // merge adjacent blocks
       for (i = 1, l = 0; i < n_off; ++i) {
	 //	 fprintf(stderr,"%d %d\n",i,l);
	 if (off[l].chunk_end>>16 == off[i].chunk_beg>>16) 
	   off[l].chunk_end = off[i].chunk_end;
	 else{
	   off[++l] = off[i];
	 }
       }
       n_off = l + 1;

     }

   }
  
   iter.n_off = n_off; iter.off = off;
#if 0
   fprintf(stderr,"printing offsets:\n");
   fprintf(stderr,"n_off=%d\n",iter.n_off);
   for(int i=0;i<iter.n_off;i++)
     fprintf(stderr,"%lu-%lu\n",iter.off[i].chunk_beg,iter.off[i].chunk_end);  
#endif
   //return ;
 }


void getOffsets(char *bainame,const aHead *hd,iter_t &ITER,int ref,int start,int stop){

  if(ITER.dasIndex==NULL||ITER.dasIndex->loadedRef !=ref){
    //    fprintf(stderr,"reding bai");
    dalloc(ITER.dasIndex);
    ITER.dasIndex = getIdx(bainame,ref);//index file
  }
    
  buildOffsets(ITER.dasIndex,ref,start,stop,ITER);
}
