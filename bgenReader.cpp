#include <stdio.h>
#include <vector>
#include <cmath>
//#include <string>
#include <errno.h>
#include <cassert>
#include "analysisFunction.h"
#include "bgenReader.h"


#include <cstdio>//fprintf etc
#include <cstdlib>// calloc etc
#include <zstd.h>// zstandard
#include <zlib.h>// zlib
#include <cstring> //memcpy



//dumb little function to convert char allele to int allele (0->A, 1->C, 2->G, 3->T, else 4)
int refToInt2(char c){
  
  int allele;
   switch (c) {
   case 'A':
     allele=0;
     break;
   case 'C':
     allele=1;
     break;
   case 'G':
     allele=2;
     break;
   case 'T':
     allele=3;
     break;
     //if cannot match the char
   default:     
     allele=4;
   }

   return(allele);
   
}

//returns 0 if not indel higher values indicates indel
int isindel(bgenLine *bgen){
  
  int isindel = 0;
  
  //we should stop if it is not a SNP
  if(bgen->la_l1 > 1 | bgen->la_l2 > 1){
    isindel = 1;
  }
  
  return isindel;
  
}

//returns chr
int getChr(bgenLine *bgen){

  //see if string is empty - meaning there is no line to read in
  if(bgen->Lchr[0]=='\0'){
    return -9;
  }
  
  //convering it to chr if not valid number throws error
  int chr = atoi(bgen->Lchr);    
  if(chr==0){
    fprintf(stderr,"\t-> Chromosome is not a valid number! Chromosome is: \'%s\'\n",bgen->Lchr);
    exit(0);
  }

  return chr;
  
}




header *bgenReader::parseheader(FILE *fp){

  header *hd = new header;
  //header block
  fread(&hd->LH,sizeof(unsigned),1,fp);
  fread(&hd->M,sizeof(unsigned),1,fp);
  fread(&hd->N,sizeof(unsigned),1,fp);
  //in order to have space for terminating char
  char magic[5]={'\0','\0','\0','\0','\0'};
  fread(magic,sizeof(char),4,fp);
  //assert compare magic value
  assert(strcmp(magic,"bgen") or strcmp(magic,"0000"));

  unsigned int fda_l = hd->LH-20;
  //always use calloc
  char *fda =(char*) calloc(fda_l,1);//free data area
  fread(fda,sizeof(char),fda_l,fp);
  unsigned flags;
  fread(&flags,sizeof(unsigned),1,fp);
  //for getting the first 2 bits - and with 3, for getting 2 first bits (max value of first 2 bits)
  hd->compressed = flags & 3;
  assert(hd->compressed>=0&&hd->compressed<3);
  //take the five first bits and move 2 to the right (2-5 bit)
  hd->layout = (flags>>2) & 15;
  if(hd->layout==0){
    fprintf(stderr,"file has unsupported layout\n");
    assert(0!=1);
  } else if(hd->layout!=1 & hd->layout!=2){
     fprintf(stderr,"file has unsupported layout\n");
     exit(0);
  }
    
  //move them 31 places to the right (leftmost bit) - could also have used &
  hd->si = flags>>31;
  if(hd->si==1){
    unsigned LSI;
    fread(&LSI,sizeof(unsigned),1,fp);
    unsigned N2;
    fread(&N2,sizeof(unsigned),1,fp);
    assert(hd->N==N2);
    hd->sampleids= new char*[hd->N];
    for(uint i=0;i<hd->N;i++){
      //short is 2 bytes
      unsigned short LSI_l;
      fread(&LSI_l,sizeof(unsigned short),1,fp);
      hd->sampleids[i] =(char*) calloc(LSI_l+1,sizeof(char));
      fread(hd->sampleids[i],sizeof(char),LSI_l,fp);
    }
  }

  return hd;
}

void cleanBgen(bgenLine *bgen){

  free(bgen->Lchr);
  free(bgen->Lrsid);
  free(bgen->Lid);
  free(bgen->la1);
  free(bgen->la2);
  
}

//implement this for .bgen simply parse a line and store it
// reads stuff into funkyPars *r
bgenLine *bgenReader::parseline(FILE *fp,header *hd){

  bgen->indis = hd->N;

  unsigned short Lid_l;
  char *Lid;
  //read in how many variant ID is
  fread(&Lid_l,sizeof(unsigned short),1,fp);
  bgen->Lid =(char*) calloc(Lid_l+1,sizeof(char));
  fread(bgen->Lid,sizeof(char),Lid_l,fp);
  
  unsigned short Lrsid_l;
  fread(&Lrsid_l,sizeof(unsigned short),1,fp);
  bgen->Lrsid =(char*) calloc(Lrsid_l+1,sizeof(char));
  fread(bgen->Lrsid,sizeof(char),Lrsid_l,fp);

  unsigned short Lchr_l;
  fread(&Lchr_l,sizeof(unsigned short),1,fp);
  bgen->Lchr =(char*) calloc(Lchr_l+1,sizeof(char));
  fread(bgen->Lchr,sizeof(char),Lchr_l,fp);

  fread(&bgen->vpos,sizeof(unsigned),1,fp);    
  fread(&bgen->nal,sizeof(unsigned short),1,fp);
  char *la1;
  char *la2;
  fread(&bgen->la_l1,sizeof(unsigned),1,fp);

  bgen->la1 =(char*) calloc(bgen->la_l1+1,sizeof(char));
  fread(bgen->la1,sizeof(char),bgen->la_l1,fp);
  fread(&bgen->la_l2,sizeof(unsigned),1,fp);
  
  bgen->la2 =(char*) calloc(bgen->la_l2+1,sizeof(char));
  fread(bgen->la2,sizeof(char),bgen->la_l2,fp);
        
  //done reading variant information;
  // uncompressed size of genotype probs
  unsigned C;
  fread(&C,sizeof(unsigned),1,fp);
  // compressed size of genotype probs
  unsigned D;
  if(hd->compressed !=0)
    fread(&D,sizeof(unsigned),1,fp);
  else
    D=C;
  
  //now do magic decompression
  // probability data storage,
  // we round up to nearst 32 bit
  // division removes everyhting after "." - same as ceiling function
  // we allocate 10 extra bytes, so that when we read in 8 bytes, we will not go into unallocated memory
  unsigned char *pds=(unsigned char*) calloc(4*(D/4+10),sizeof(unsigned char));
  
  unsigned char *upds=(unsigned char*) calloc(C-4,sizeof(unsigned char));
  
  if(hd->compressed==0){
    fread(pds,sizeof(unsigned char),C,fp);
  } else if(hd->compressed==1){
    //added this as otherwise, does not read from current position in file
    fread(upds,sizeof(unsigned char),C-4,fp);  
    //zlib https://gist.github.com/arq5x/5315739
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;

    // setup "b" as the input and "c" as the compressed output    
    infstream.avail_in = C-4;
    infstream.next_in = (Bytef *)upds; // input char array
    infstream.avail_out = D; // size of output
    infstream.next_out = (Bytef *)pds; // output char array
    
    // the actual DE-compression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
  } else if(hd->compressed==2) {
    //zstd
    fread(upds,sizeof(unsigned char),C-4,fp);  
    ZSTD_decompress(pds,D,upds,C-4);
  }

  //data is now in pds;//<- probability datastorage
  unsigned int N2;memcpy(&N2,pds,4);

  assert(hd->N==N2); 

  unsigned short nal2;memcpy(&nal2,pds+4,2);
  assert(bgen->nal==nal2);
  char pmin=pds[6];
  char pmax=pds[7];
  char pmis[hd->N];//missingness and ploidy
  memcpy(pmis,pds+8,hd->N);
  for(uint ii=0;ii<hd->N;ii++){
    //extract MSB - if this indi has missing
    int missing = pmis[ii]>>7;
    // extract first 6 bits ploidy of this individual
    int ploidy = pmis[ii] & 63;
  }
  
  //jump 8+N bytes in pds array, and copy to phased variable
  unsigned char phased;memcpy(&phased,pds+8+hd->N,1);
  //number of bits used for each genotype prob - copied like phased variable
  // can be from 1 to 32 bits
  unsigned char B;memcpy(&B,pds+8+hd->N+1,1);
  
  // where genotype probs start
  unsigned char *X=pds+10+hd->N;  
  double scal = pow(2,B)-1;  
  uint to_vals[2];
  int at=0;
  unsigned int mask = pow(2,B)-1;

  for(int ii=0;ii<2*hd->N;ii++){
    
    // which bits are we reading in, for which indi (ii goes through 2*indis)
    // point dingdong points to closest byte (before the bit)
    size_t *dingdong =(size_t *) (X +ii*B/8);//<- parantheses are important
    size_t dingdong_val = *dingdong; //copy value not pointer.
    
    //now we work with dingdong_val
    int shift = ii*B-8*(ii*B/8);//shift is how much we need to shift such that our data stats at bit0. NB paranthesis is important
    //fprintf(stderr,"shift_how_much[%d]:%d\n",ii,shift);
    dingdong_val = dingdong_val >> shift ; //<- notice that shift right >>means shifting towards the least signifcant bit
    
    uint tmptmp = 0;
    memcpy(&tmptmp,&dingdong_val,4);//only copy 32bits
    tmptmp = tmptmp &mask ;
    to_vals[ii %2 ] = tmptmp;
    
    //below should be the same
    if( ii %2 ){

      uint p11 = to_vals[0];
      // genotype prob HE allele 1 - genotype is 4 bits or 0.5 bytes
      uint p12 = to_vals[1];

      bgen->probs[at] = (1.0*p11)/scal;
      bgen->probs[at+1] = (1.0*p12)/scal;
      bgen->probs[at+2] = 1-(1.0*p11)/scal-(1.0*p12)/scal;
      
      if(bgen->probs[at+2]<0) //<-needed?
	bgen->probs[at+2] = 0;
      
      double sum = bgen->probs[at]+bgen->probs[at+1]+bgen->probs[at+2];
      
      bgen->probs[at] = bgen->probs[at] / sum;
      bgen->probs[at+1] = bgen->probs[at+1] / sum;
      bgen->probs[at+2] = bgen->probs[at+2] / sum;

      at+=3;
      
    }
        
  }

  free(pds);
  free(upds);

  return bgen;
  
}


void bgenReader::funkyCopy(bgenLine *bgen, funkyPars *r, int &balcon){    
  
  //to translate chr from fai indexing to normal indexing 1,2,3,4,5,...
  int wacko[26] = {0,11,15,16,17,18,19,20,21,1,2,3,4,5,6,7,8,9,10,12,13,14,22,23,24};
  //because 0 indexed
  r->refId = wacko[atoi(bgen->Lchr)-1];

  //because it is assumed to be 0-indexed in ANGSD
  r->posi[balcon] = bgen->vpos-1;
  //we do not know which one is major and minor, but probaly like beagle files
  //TO DO check which one is major and minor
  r->major[balcon] = refToInt2(bgen->la1[0]);
  r->minor[balcon] = refToInt2(bgen->la2[0]);
  //so that site is analysed
  r->keepSites[balcon] = 1;

  for(int at=0;at<3*bgen->indis;at+=3){
    //jumps 3 at a time
    r->post[balcon][at] = bgen->probs[at];
    r->post[balcon][at+1] = bgen->probs[at+1];
    r->post[balcon][at+2] = bgen->probs[at+2];
  }

  balcon++;
  
}


funkyPars *bgenReader::fetch(int chunksize){

  funkyPars *r = funkyPars_init();
  
  //using refId for chr
  //TO DO add chr field in funkyPars struct
  //r->refId=new int[chunkSize];
  r->posi=new int[chunksize];
  
  r->post=new double*[chunksize];
  r->major = new char[chunksize];
  r->minor = new char[chunksize];
  r->keepSites = new int[chunksize];
  r->refId = 0;

  r->post = new double*[chunksize];
  r->nInd = nInd;
  
  for(int s=0;s<chunksize;s++)
    r->post[s] = new double[nInd*3];
  
  //keeps track of which site we are at (we might skip some non diallelic sites)
  int balcon=0;

  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  
  //go to statement
 never_ever: //haha
  
  //if acpy not null - read read that has survived
  //member variable of vcfReader class 
  if(readAgain){

    fprintf(stderr,"READ AGAIN!\n");
    
    //read stored line of file    
    //bgen = parseline(bgenFile,hd);
    //curChr=getChr(bgen);
        
    funkyCopy(bgen,r,balcon);

    cleanBgen(bgen);
    
    readAgain = 0;

    //after read line we have update prevChr - current chr will be updated at the beginning of the loop
    prevChr=curChr;
        
  }

  //TO DO what if fewer sites than one chunck??
  //balcon how many ok sites
  while(balcon<chunksize && balcon<sites && sitesRead<sites) {

    bgen = parseline(bgenFile,hd);

    // to get chr
    curChr=getChr(bgen);

    //keep track of how many sites read in total
    sitesRead++;    
    
    //meaning there is no more data to be read no chr to read..
    if(curChr==-9){
      break;
    }
    
    if(prevChr==-1){
      prevChr=curChr;
    }
      
    n++;
        
    //skip nonsnips - gives back 0 if indel
    //reads first X bytes of line to see if a SNP, then jumps back X bytes
    if(isindel(bgen)){
      prevChr=curChr;
      if(onlyPrint>0){
	fprintf(stderr,"\t Skipping due to non snp pos:%d (this message will be silenced after 10 sites)\n",r->posi+1);
	onlyPrint--;
      }      
      continue;
    }
        
    nsnp++;

    //if we are changing chromosomes
    if(prevChr!=curChr){

      readAgain=1;
  
      //I do not think I have to take a copy I can just
      if(balcon==0){
	goto never_ever;
      }
      //acpy=bgen_dup(rec);
      break;
    }

    //copy data from line into funkyPars
    funkyCopy(bgen,r,balcon);

    cleanBgen(bgen);
    
    prevChr=curChr;

  }

  //to delete sites, if we do not have a whole chunck, but we have allocated a whole chunck
  if(balcon<chunksize){    
    for(int s=balcon;s<chunksize;s++)
      delete [] r->post[s];
  }
  
  // to catch when no more data for this chr and chunk is over
  //if(balcon==0&&bcf_retval==0)
  //  goto never_ever;

  //  fprintf(stderr, "\t-> [file=\'%s\'][chr=\'%s\'] Read %i records %i of which were SNPs number of sites with data:%lu\n",fname,seek, n, nsnp,mygl.size()); 
  
  r->numSites=balcon;
  if(r->numSites==0){    
    funkyPars_destroy(r);
    r=NULL;
  }

  return r;
 
}


