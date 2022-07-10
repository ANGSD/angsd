#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <errno.h>
#include <cassert>
#include <inttypes.h>
#include "analysisFunction.h"
#include "vcfReader.h"

bam_hdr_t *bcf_hdr_2_bam_hdr_t (htsstuff *hs){
  bam_hdr_t *ret = bam_hdr_init();
  ret->l_text = 0;
  ret->text =NULL;

  int nseq=0;

  for (int i=0; i<hs->hdr->nhrec; i++){
    bcf_hrec_t *hrec=hs->hdr->hrec[i];
    if(strcmp(hrec->key,"contig")==0)
      nseq++;
  }
  
  ret->n_targets = nseq;
  ret->target_len = (uint32_t*) malloc(sizeof(uint32_t)*nseq);
  ret->target_name = (char**) malloc(sizeof(char*)*nseq);
  int at=0;
  for (int i=0; i<hs->hdr->nhrec; i++){
    bcf_hrec_t *hrec=hs->hdr->hrec[i];
    if(strcmp(hrec->key,"contig")==0){
      //      fprintf(stderr,"%d) hrec->value:%s key:%s\n",i,hrec->value,hrec->key);
      int newlen =-1;
      char *chrnam=NULL;
      for(int j=0;j<hrec->nkeys;j++){
	if(strcmp("ID",hrec->keys[j])==0)
	  chrnam = strdup(hrec->vals[j]);
	if(strcmp("length",hrec->keys[j])==0)
	  newlen = atoi(hrec->vals[j]);
	//fprintf(stderr,"i:%d j:%d keys:%s vals:%s\n",i,j,hrec->keys[j],hrec->vals[j]);
      }
      //fprintf(stderr,"at: %d ID:%s len:%d\n",at,chrnam,newlen);
      ret->target_len[at] = newlen;
      ret->target_name[at] = chrnam;
      at++;
    }
  }

  return ret;
}

void htsstuff_seek(htsstuff *hs,char *seek){
  assert(seek);
  if(seek){
    fprintf(stderr,"\t-> Setting iterator to: %s\n",seek);fflush(stderr);
    if(hs->idx==NULL)
      hs->idx=bcf_index_load(hs->fname);
    if(hs->idx==NULL){
      fprintf(stderr,"\t-> Problem opening index file of file: \'%s\'\n",hs->fname);
      exit(0);
    }
    if(hs->hts_file->format.format!=bcf){
      fprintf(stderr,"\t-> File are required to be vcf for using region specification\n");
      exit(0);
      
    }
    if(hs->iter)
      hts_itr_destroy(hs->iter);
    hs->iter=bcf_itr_querys(hs->idx,hs->hdr,seek);
    assert(hs->iter);
  }
}

htsstuff *htsstuff_init(char *fname,char *seek){
  //  fprintf(stderr,"[%s] seek:%s\n",__FUNCTION__,seek);
  htsstuff *hs = new htsstuff;
  hs->fname=NULL;hs->fname=strdup(fname);
  hs->seek=NULL;
  if(seek) hs->seek=strdup(seek);
  hs->hts_file=NULL;
  hs->hdr=NULL;
  hs->idx=NULL;
  hs->iter=NULL;
  hs->nsamples = 0;
  hs->hts_file=hts_open(hs->fname,"r");
  assert(hs->hts_file);
  
  hs->hdr = bcf_hdr_read(hs->hts_file);
  assert(hs->hdr);

  hs->nsamples = bcf_hdr_nsamples(hs->hdr);
  if(hs->seek)
    htsstuff_seek(hs,hs->seek);
  return hs;
}

void htsstuff_destroy(htsstuff *hs){
  free(hs->fname);
  free(hs->seek);
  
  if(hs->hdr)
    bcf_hdr_destroy(hs->hdr);
  if(hs->hts_file)
    bcf_close(hs->hts_file);

  if(hs->iter)
    hts_itr_destroy(hs->iter);
  if(hs->idx)
    hts_idx_destroy(hs->idx);
  

  delete hs;
}


// https://github.com/samtools/htslib/blob/57fe419344cb03e2ea46315443abd242489c32f2/vcf.c#L53
// uint32_t bcf_float_missing    = 0x7F800001;
// uint32_t bcf_float_vector_end = 0x7F800002;

// double pl2ln[256];

const char *wildcar= "<*>";

//mindfuck2000
void buildreorder(int swap[10],char **alleles,int len){
  char bases[] = {'A','C','G','T'};
  for(int i=0;0&&i<len;i++)
    fprintf(stderr,"%d):%s\n",i,alleles[i]);
  int at = 0;//<- yes sir
  char tmp[4]={'N','N','N','N'};
  for(int b=0;b<4;b++) {
    // fprintf(stderr,"b:%d\n",b);
    int i;
    for(i=0;i<len;i++){
      // fprintf(stderr,"cmp: %c %c\n",alleles[i][0],bases[b]);
      if(alleles[i][0]==bases[b])
	break;
    }
    //  fprintf(stderr,"i:%d\n",i);
    if(i==len)
      tmp[at++]=bases[b];
  }
  for(int i=0;0&&i<4;i++)
    fprintf(stderr,"tmp[%d]:%c\n",i,tmp[i]);
  for(int i=0;i<10;i++)
    swap[i]= -1;
 
  int adder =0;
  for(int en=0;en<len;en++){
    for(int to=0;to<=en;to++){
      // fprintf(stderr,"RAW[%d]: en:%s to:%s \n",adder,alleles[en],alleles[to]);
      if(strcmp(alleles[en],wildcar)!=0&&strcmp(alleles[to],wildcar)!=0){
	//	fprintf(stderr,"en:%s to:%s\n",alleles[en],alleles[to]);
	swap[angsd::majorminor[refToInt[alleles[en][0]]][refToInt[alleles[to][0]]]] = adder ;
      }
      else if(strcmp(alleles[en],wildcar)==0&&strcmp(alleles[to],wildcar)!=0){
	for(int i=0;tmp[i]!='N';i++){
	  //	  fprintf(stderr,"Een=<*> en:%s to:%s tre:%c\n",alleles[en],alleles[to],tmp[i]);
	  swap[angsd::majorminor[refToInt[alleles[to][0]]][refToInt[tmp[i]]]] = adder;
	}
      }
      else if(strcmp(alleles[en],wildcar)!=0&&strcmp(alleles[to],wildcar)==0){
	for(int i=0;tmp[i]!='N';i++){
	  //	  fprintf(stderr,"Eto=<*> en:%s to:%s tre:%c\n",alleles[en],alleles[to],tmp[i]);
	  swap[angsd::majorminor[refToInt[alleles[en][0]]][refToInt[tmp[i]]]] = adder;
	}
      }
      else if(strcmp(alleles[en],wildcar)==0&&strcmp(alleles[to],wildcar)==0){
	for(int i=0;tmp[i]!='N';i++)
	  for(int j=i;tmp[j]!='N';j++){
	    //	    fprintf(stderr,"ASDF en:%s to:%s tre:%c fire:%c\n",alleles[en],alleles[to],tmp[i],tmp[j]);
	    swap[angsd::majorminor[refToInt[tmp[i]]][refToInt[tmp[j]]]] = adder;
      	    
	  }
      }
      adder++;
    }
  }
#if 0
    fprintf(stderr,"AA swap[%d]:%d\n",0,swap[0]);
    fprintf(stderr,"AC swap[%d]:%d\n",1,swap[1]);
    fprintf(stderr,"AG swap[%d]:%d\n",2,swap[2]);
    fprintf(stderr,"AT swap[%d]:%d\n",3,swap[3]);
    fprintf(stderr,"CC swap[%d]:%d\n",4,swap[4]);
    fprintf(stderr,"CG swap[%d]:%d\n",5,swap[5]);
    fprintf(stderr,"CT swap[%d]:%d\n",6,swap[6]);
    fprintf(stderr,"GG swap[%d]:%d\n",7,swap[7]);
    fprintf(stderr,"GT swap[%d]:%d\n",8,swap[8]);
    fprintf(stderr,"TT swap[%d]:%d\n",9,swap[9]);
#endif
}

//returns which info field is the INDEL
int whichisindel(const bcf_hdr_t *h){
  int hit =-1;
  //fprintf(stderr,"ninfo:%d\n",h->n[BCF_DT_ID]);
  for(int i=0;i<h->n[BCF_DT_ID];i++){
    //    fprintf(stderr,"i:%d key:%s\n",i,h->id[BCF_DT_ID][i].key);
    if(h->id[BCF_DT_ID][i].key&&strcmp(h->id[BCF_DT_ID][i].key,"INDEL")==0)
      return i;
  }
  if(hit==-1)
    fprintf(stderr,"\t-> No indel tag in vcf/bcf file, will therefore not be able to filter out indels\n");
  return hit;
}

//returns 0 if not indel higher values indicates indel
int isindel(const bcf_hdr_t *h, const bcf1_t *v){
  bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
  static int hit = whichisindel(h);//only calculated once
  //  assert(hit!=-1);//we assume INDEL exists in INFO field
  if(hit==-1)
    return 0;
  int returnvalue = 0;
  if (v->n_info) {//loop over all INFO fields
    for (int i = 0; i < v->n_info; ++i) {
      bcf_info_t *z = &v->d.info[i];
      if ( !z->vptr ) continue;
      if (z->key >= h->n[BCF_DT_ID]) {
	fprintf(stderr,"Invalid BCF, the INFO index is too large");
	errno = EINVAL;
	return -1;
      }
      //fprintf(stderr,"z-key:%d  %s\n",z->key,h->id[BCF_DT_ID][z->key].key);
      if(z->key==hit)//if key hits the INDEL 
	returnvalue++;
    }
  }
  return returnvalue;
}


//type=0 -> PL
//type=1 -> GL
//type=2 -> GP

int dumpcounterverbose[2] = {0,0};
int vcfReader::parseline(bcf1_t *rec,htsstuff *hs,funkyPars *r,int &balcon,int type){
  assert(type>=0&&type<=2);
  int n;
  if(isindel(hs->hdr,rec)!=0)
    return 0;

  r->major[balcon] = refToChar[rec->d.allele[0][0]];
  r->minor[balcon] = refToChar[rec->d.allele[1][0]];
  if(r->major[balcon]==4)//<- means N
    return 0;
  
  // pl data for each call

  int myreorder[10];
  //funky below took ridiculus long to make
  buildreorder(myreorder,rec->d.allele,rec->n_allele);
  
  
  double *dupergl = NULL;
  double *dupergp = NULL;
  if(type==0||type==1){
    dupergl= new double[10*hs->nsamples];
    for(int i=0;i<10*hs->nsamples;i++)
      dupergl[i] = log(0);
  }else if(type==2){
    dupergp= new double[3*hs->nsamples];
    for(int i=0;i<3*hs->nsamples;i++)
      dupergp[i] = 0;
  }
  if(type==2) {
    n = bcf_get_format_float(hs->hdr, rec, "GP", &farr, &mfarr);
    if(n>=ln_gl_m){
      ln_gl_m  = n;
      kroundup32(ln_gl_m);
      ln_gl = (float *) realloc(ln_gl,sizeof(float)*ln_gl_m);
    }
    if(n<0){
      // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
      // no PL tag is available
      fprintf(stderr, "BAD SITE %s:%" PRId64 " return code:%d" "while fetching GP tag rec->rid:%d\n", bcf_seqname(hs->hdr,rec), rec->pos, n,rec->rid);
      return 0;
    }{
      // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151

      for(int nsample = 0 ;nsample < r->nInd; nsample++){
	if (bcf_float_is_missing(farr[nsample*3])){
	  if(dumpcounterverbose[0]++ <10)
	    fprintf(stderr,"[bcf_float_is_missing error], this might be an error or missing data, this msg is printed %d times more \n",10-dumpcounterverbose[0]);
	}else if ( farr[nsample*3]==bcf_float_vector_end ){
	  if(dumpcounterverbose[1]++ <10)
	    fprintf(stderr,"[bcf_float_vector_end], this might be an error or missing data, this msg is printed %d times more\n",10-dumpcounterverbose[1]);
	  exit(0);
	} else{
	  dupergp[nsample*3] = farr[nsample*3];
	  dupergp[nsample*3+1] = farr[nsample*3+1];
	  dupergp[nsample*3+2] = farr[nsample*3+2];
	}
      }
    }
  } else  if(type==0) {
    //parse PL
    n = bcf_get_format_int32(hs->hdr, rec, "PL", &iarr, &miarr);
    if(n>=ln_gl_m){
      ln_gl_m  = n;
      kroundup32(ln_gl_m);
      ln_gl = (float *) realloc(ln_gl,sizeof(float)*ln_gl_m);
    }
    if(n<0){
      // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
      // no PL tag is available
      fprintf(stderr, "BAD SITE %s:%" PRId64 ". return code:%d while fetching PL tag rec->rid:%d\n", bcf_seqname(hs->hdr,rec), rec->pos, n,rec->rid);
      delete [] dupergl;
      return 0;
    }else{
      // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151
      for (int i=0; i<n; i++){
	if ( iarr[i]==bcf_int32_missing ){
	  bcf_float_set_missing(ln_gl[i]);
	} else if ( iarr[i]==bcf_int32_vector_end ){
	  bcf_float_set_vector_end(ln_gl[i]);
	} else{
	  ln_gl[i] = pl2ln_f(iarr[i]);
	}
	//         fprintf(stderr, "%d %f\n", pl[i], ln_gl[i]);
      }
    }
   
  }else if(type==1) {
    n = bcf_get_format_float(hs->hdr, rec, "GL", &farr, &mfarr);
    if(n>=ln_gl_m){
      ln_gl_m  = n;
      kroundup32(ln_gl_m);
      ln_gl = (float *) realloc(ln_gl,sizeof(float)*ln_gl_m);
    }
    
    if(n<0){
      // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
      // no PL tag is available
      fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching GL tag rec->rid:%d\n", bcf_seqname(hs->hdr,rec),(int) rec->pos, n,rec->rid);
      return 0;
    }{
      // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151
      for (int i=0; i<n; i++){
	if (bcf_float_is_missing(farr[i]) ){
	  bcf_float_set_missing(ln_gl[i]);
	} else if ( farr[i]==bcf_float_vector_end ){
	  bcf_float_set_vector_end(ln_gl[i]);
	} else{
	  ln_gl[i] = farr[i]/M_LOG10E;
	  //	  if(ln_gl[i]==0.0) 	    ln_gl[i] = -0.0;
	}
	//	fprintf(stderr, "%f %f\n", farr[i], ln_gl[i]);
      }
    }
  } else {
    fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, PL or GL\n");
    return 0;
  }
  /*
    We assume we have data and set all to -inf
    Then we check if something is nan, if one of them are nan we set all gls to missing -0.0
    Secondly we check that the gls we have plugged in contains information (difference between them), is not set to miggin
   */
  if(type==0||type==1){
    if(n>=0){//case where we have data
     // fprintf(stderr,"data\n");
      assert((n % hs->nsamples)==0 );
      int myofs=n/hs->nsamples;
      for(int ind=0;ind<hs->nsamples;ind++){
	int rollback =0;
	int at =0;
	double hastaken[10];
	for(int o=0;o<10;o++){
	  if(myreorder[o]!=-1){
	    dupergl[ind*10+o] = ln_gl[ind*myofs+myreorder[o]];
	    hastaken[at++] = ln_gl[ind*myofs+myreorder[o]];
	  }
	  if(std::isnan(dupergl[ind*10+o])){
	    rollback = 1;
	    break;
	  }
	}



	for(int o=0;rollback&&o<10;o++)
	    dupergl[ind*10+o] = 0.0;

	rollback =0;

	for(int i=1;at>1&&i<at;i++)
	  if(hastaken[0]!=hastaken[i])
	    rollback =1;

	if(rollback==0)
	  for(int o=0;o<10;o++)
	    dupergl[ind*10+o] = 0.0;

      }
    }
  }
  if(0&&type==2){//FIX
    for(int ind=0;ind<hs->nsamples;ind++){
      double tsum = dupergp[3*ind]+dupergp[3*ind+1]+dupergp[3*ind+2];
      dupergp[3*ind] /= tsum;
      dupergp[3*ind+1] /= tsum;
      dupergp[3*ind+2] /= tsum;
    }
  }
  if(type==0||type==1){
    r->likes[balcon] = dupergl;
#if 0
    for(int i=0;i<hs->nsamples;i++){
      for(int j=0;j<10;j++)
	fprintf(stderr,"%f\t",r->likes[balcon][i*10+j]);
      fprintf(stderr,"\n");
    }
#endif
    
  } else if(type==2){
    r->post[balcon] = dupergp;
  }
  
  r->refId=rec->rid;
  r->posi[balcon] = rec->pos;
  r->keepSites[balcon] = hs->nsamples;
  balcon++;
  return 1;//<- what should this function return?
}

//stupid little function to wrap the iter part. This could be done much more cleverly
int vcfReader::vcfReaderwrap_reader(htsstuff *hts,bcf1_t *rec){
  int ret;
  static int whichRegion =-1;
  static int hasSeeked =0;
  if(regions->size()==0){
    ret  = bcf_read(hts->hts_file,hts->hdr,rec);
    return ret;
  }else{
    if(hasSeeked==0){
      whichRegion++;
      if(whichRegion>=regions->size())
	return -1;
      itrname.l =0;
      int start=regions->at(whichRegion).start;
      int stop=regions->at(whichRegion).stop;
      int ref=regions->at(whichRegion).refID;
      ksprintf(&itrname,"%s:%d-%d",bamhdr->target_name[ref],start+1,stop);
      //fprintf(stderr,"ksprintf :%s\n",itrname.s);
      seek(itrname.s);
      hasSeeked =1;
      return vcfReaderwrap_reader(hts,rec);//calling again now that things has been seeked.
    }
    ret=bcf_itr_next(hts->hts_file,hts->iter, rec);
    if(ret<0){
      hasSeeked = 0;
      return vcfReaderwrap_reader(hts,rec);
    }
    return ret;
  }
}
//

int onlyprint = 10;
funkyPars *vcfReader::fetch(int chunkSize) {
  funkyPars *r = funkyPars_init();
  r->nInd = hs->nsamples;
  if(pl_gl_gp==0||pl_gl_gp==1)
    r->likes=new double*[chunkSize];
  else if(pl_gl_gp==2)
    r->post=new double*[chunkSize];
  r->keepSites = new int[chunkSize];
  r->major = new char[chunkSize];
  r->minor = new char[chunkSize];  
  for(int i=0;i<chunkSize;i++){
    if(pl_gl_gp==0||pl_gl_gp==1){
      memset(r->likes,0,chunkSize*sizeof(double*));
    } else if(pl_gl_gp==2){
      memset(r->post,0,chunkSize*sizeof(double*));
    }
  }
  
  r->posi=new int[chunkSize];
  
  bcf1_t *rec = bcf_init();

  //   http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  
  int balcon=0;
  int bcf_retval =0;
 never: //haha

  if(acpy){
    curChr=acpy->rid;
    parseline(acpy,hs,r,balcon,pl_gl_gp);
    bcf_destroy(acpy);
    acpy=NULL;
  }
  while(balcon<chunkSize) {
    //either parse a read with region or nextread
    if(((bcf_retval=vcfReaderwrap_reader(hs,rec)))!=0)
      break;
    
    n++;
    //skip nonsnips
    if(isindel(hs->hdr,rec)){
      if(onlyprint>0){
	fprintf(stderr,"\t Skipping due to non snp pos:%d (this message will be silenced after 10 sites)\n",(int)rec->pos+1);
	onlyprint--;
      }
      continue;
    }

    nsnp++;

    //initialize
    if(curChr==-1)
      curChr=rec->rid;
    
    //if we are changing chromosomes
    if(rec->rid!=curChr){
      acpy=bcf_dup(rec);
      break;
    }

    //regular case
    parseline(rec,hs,r,balcon,pl_gl_gp);
  }
  if(balcon==0&&bcf_retval==0)
    goto never;

  //  fprintf(stderr, "\t-> [file=\'%s\'][chr=\'%s\'] Read %i records %i of which were SNPs number of sites with data:%lu\n",fname,seek, n, nsnp,mygl.size()); 
  bcf_destroy(rec);
  r->numSites=balcon;
  if(r->numSites==0){
    funkyPars_destroy(r);
    r=NULL;
  }
  return r;
}

bcf_hdr_t *vcfreader_hs_bcf_hdr  = NULL;//nasty hack global dragon
vcfReader::vcfReader(char *fname,char *seek,int pl_or_gl_a,std::vector<regs> *regions_a){
  itrname.s=NULL;itrname.l=itrname.m =0;
  regions = regions_a;
  farr=NULL;
  iarr=NULL;
  mfarr=0;
  miarr=0;
  ln_gl_m = 8;
  ln_gl =(float *) malloc(sizeof(float)*ln_gl_m);
  pl_gl_gp = pl_or_gl_a;
  hs=htsstuff_init(fname,seek);
  vcfreader_hs_bcf_hdr = hs->hdr;
  bamhdr = bcf_hdr_2_bam_hdr_t(hs);
  acpy=NULL;
  for(int i=0;i<PHREDMAX;i++)
    pl2ln[i] = log(pow(10.0,-0.1*i));
  curChr=-1;
  pl=NULL;
}
