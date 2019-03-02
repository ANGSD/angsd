#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <errno.h>
#include "analysisFunction.h"
#include "vcfReader.h"


void htsstuff_seek(htsstuff *hs,char *seek){
  assert(seek);
  if(seek){
    fprintf(stderr,"\t-> Setting iterator to: %s\n",seek);fflush(stderr);
    hs->idx=bcf_index_load(hs->fname);
    if(hs->idx==NULL){
      fprintf(stderr,"\t-> Problem opening index file of file: \'%s\'\n",hs->fname);
      exit(0);
    }
    if(hs->hts_file->format.format!=bcf){
      fprintf(stderr,"\t-> File are required to be vcf for using region specification\n");
      exit(0);
    }
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
#if 0
//copied from vcf.c in htslib to understand accessing of these variables
int vcf_format_tsk(const bcf_hdr_t *h, const bcf1_t *v)
{
  kstring_t *s=(kstring_t *)malloc(sizeof(kstring_t));
  int i;
  bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
  kputs(h->id[BCF_DT_CTG][v->rid].key, s); // CHROM
  kputc('\t', s); kputw(v->pos + 1, s); // POS
  kputc('\t', s); kputs(v->d.id ? v->d.id : ".", s); // ID
  kputc('\t', s); // REF
  if (v->n_allele > 0) kputs(v->d.allele[0], s);
  else kputc('.', s);
  kputc('\t', s); // ALT
  if (v->n_allele > 1) {
    for (i = 1; i < v->n_allele; ++i) {
      if (i > 1) kputc(',', s);
      kputs(v->d.allele[i], s);
    }
  } else kputc('.', s);
  kputc('\t', s); // QUAL
  if ( bcf_float_is_missing(v->qual) ) kputc('.', s); // QUAL
  else kputd(v->qual, s);
  kputc('\t', s); // FILTER
  if (v->d.n_flt) {
    for (i = 0; i < v->d.n_flt; ++i) {
      if (i) kputc(';', s);
      kputs(h->id[BCF_DT_ID][v->d.flt[i]].key, s);
    }
  } else kputc('.', s);
  kputc('\t', s); // INFO
  if (v->n_info) {
    int first = 1;
    for (i = 0; i < v->n_info; ++i) {
      bcf_info_t *z = &v->d.info[i];
      if ( !z->vptr ) continue;
      if ( !first ) kputc(';', s);
      first = 0;
            if (z->key >= h->n[BCF_DT_ID]) {
	      hts_log_error("Invalid BCF, the INFO index is too large");
	      errno = EINVAL;
	      return -1;
            }
            kputs(h->id[BCF_DT_ID][z->key].key, s);
	    fprintf(stderr,"%s: \n",h->id[BCF_DT_ID][z->key].key);
            if (z->len <= 0) continue;
            kputc('=', s);
            if (z->len == 1)
            {
	      switch (z->type)
                {
		case BCF_BT_INT8:  if ( z->v1.i==bcf_int8_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
		case BCF_BT_INT16: if ( z->v1.i==bcf_int16_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
		case BCF_BT_INT32: if ( z->v1.i==bcf_int32_missing ) kputc('.', s); else kputw(z->v1.i, s); break;
		case BCF_BT_FLOAT: if ( bcf_float_is_missing(z->v1.f) ) kputc('.', s); else kputd(z->v1.f, s); break;
		case BCF_BT_CHAR:  kputc(z->v1.i, s); break;
		default: hts_log_error("Unexpected type %d", z->type); exit(1); break;
                }
            }
            else bcf_fmt_array(s, z->len, z->type, z->vptr);
    }
    if ( first ) kputc('.', s);
  } else kputc('.', s);
  // FORMAT and individual information

    kputc('\n', s);
    fprintf(stderr,"kstirngts: %s\n",s->s);
    free(s->s);free(s);
      
    return 0;
}
#endif
//returns which info field is the INDEL
int whichisindel(const bcf_hdr_t *h){
  int hit =-1;
  //fprintf(stderr,"ninfo:%d\n",h->n[BCF_DT_ID]);
  for(int i=0;i<h->n[BCF_DT_ID];i++){
    //    fprintf(stderr,"i:%d key:%s\n",i,h->id[BCF_DT_ID][i].key);
    if(h->id[BCF_DT_ID][i].key&&strcmp(h->id[BCF_DT_ID][i].key,"INDEL")==0)
      return i;
  }
  return hit;
}

//returns 0 if not indel higher values indicates indel
int isindel(const bcf_hdr_t *h, const bcf1_t *v){
  bcf_unpack((bcf1_t*)v, BCF_UN_ALL);
  static int hit = whichisindel(h);//only calculated once
  assert(hit!=-1);//we assume INDEL exists in INFO field
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



int vcfReader::parseline(bcf1_t *rec,htsstuff *hs,funkyPars *r,int &balcon){
  if(isindel(hs->hdr,rec)!=0)
    return 0;

  // pl data for each call
  int npl_arr = 0;
  int npl     = 0;

  std::string vcf_format_field = "PL"; 
  int myreorder[10];
  //funky below took ridiculus long to make
  buildreorder(myreorder,rec->d.allele,rec->n_allele);
  
  
  double *dupergl = new double[10*hs->nsamples]; 
  for(int i=0;i<10*hs->nsamples;i++)
    dupergl[i] = log(0);
  //parse PL
  if(vcf_format_field == "PL") {
    npl = bcf_get_format_int32(hs->hdr, rec, "PL", &pl, &npl_arr);
    float ln_gl[npl];
    if(npl<0){
      // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
      // no PL tag is available
      fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching PL tag rec->rid:%d\n", bcf_seqname(hs->hdr,rec), rec->pos, npl,rec->rid);
      return 0;
    }
    // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151
    for (int i=0; i<npl; i++){
      if ( pl[i]==bcf_int32_missing ){
	bcf_float_set_missing(ln_gl[i]);
      } else if ( pl[i]==bcf_int32_vector_end ){
	bcf_float_set_vector_end(ln_gl[i]);
      } else{
	ln_gl[i] = pl2ln_f(pl[i]);
      }
      //         fprintf(stderr, "%d %f\n", pl[i], ln_gl[i]);
    }
    assert((npl % hs->nsamples)==0 );
    int myofs=npl/hs->nsamples;
    for(int ind=0;ind<hs->nsamples;ind++){
      for(int o=0;o<10;o++)
	if(myreorder[o]!=-1)
	  dupergl[ind*10+o] = ln_gl[ind*myofs+myreorder[o]]; 
    }
    r->likes[balcon] = dupergl;
  } else {
    fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, GT or PL\n");
    return 0;
   
  }
  r->refId=rec->rid;
  r->posi[balcon] = rec->pos;
  r->keepSites[balcon] = hs->nsamples;

  balcon++;
  //naf = bcf_get_info_float(hdr, rec, vcf_allele_field.c_str(), &af, &naf_arr);//maybe include this<-
}



funkyPars *vcfReader::fetch(int chunkSize){
  funkyPars *r = funkyPars_init();
  r->nInd = hs->nsamples;
  r->likes=new double*[chunkSize];
  r->post=new double*[chunkSize];
  r->keepSites = new int[chunkSize];
  for(int i=0;i<chunkSize;i++){
    memset(r->likes,0,chunkSize*sizeof(double*));
    memset(r->post,0,chunkSize*sizeof(double*));
  }
  
  r->posi=new int[chunkSize];
  
  bcf1_t *rec = NULL;rec=bcf_init();assert(rec);

  //   http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
  // counters
  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  int nseq = 0;  // number of sequences
  
  int balcon=0;
  
  if(acpy){
    parseline(acpy,hs,r,balcon);
    bcf_destroy(acpy);
    acpy=NULL;
  }
  
  while(balcon<chunkSize) {
    //either parse a read with region or nextread
    if(hs->iter==NULL){
      if(bcf_read(hs->hts_file,hs->hdr,rec)!=0)	
	break;
    }else{
      if(bcf_itr_next(hs->hts_file, hs->iter, rec)!=0)
	break;
    }

    n++;
    //skip nonsnips

    if(isindel(hs->hdr,rec)){
      //      fprintf(stderr,"skipping due to non snp\n");
      continue;
    }

    nsnp++;

    //initialize
    if(curChr==-1){
      curChr=rec->rid;
    }

    //if we are changing chromosomes
    if(rec->rid!=curChr){
      //remove old
      if(acpy){
	bcf_destroy(acpy);
	acpy=NULL;
      }
      //copy new
      acpy=bcf_dup(rec);
      curChr=rec->rid;
      break;
    }

    //regular case
    parseline(rec,hs,r,balcon);
  }
  //  fprintf(stderr, "\t-> [file=\'%s\'][chr=\'%s\'] Read %i records %i of which were SNPs number of sites with data:%lu\n",fname,seek, n, nsnp,mygl.size()); 

  bcf_destroy(rec);
  r->numSites=balcon;
  if(r->numSites==0){
    funkyPars_destroy(r);
    r=NULL;
  }
  return r;
 
}


