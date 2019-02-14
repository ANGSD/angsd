/*
  example modified from http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html

  compile and run with:
  g++ vcf.cpp -I../htslib/ ../htslib/libhts.a -lz -D__WITH_MAIN__

./a.out my.bcf
 */

#include <stdio.h>
#include <vector>
#include <htslib/vcf.h>
#include <cmath>
#include <string>
#include "analysisFunction.h"

typedef struct satan_t{
  char*fname;
  std::string vcf_format_field;
  std::string vcf_allele_field;
  char *seek;
  std::vector<double *> mygl;
  std::vector<double> freqs;
  int nind;
}satan;

std::vector<satan> jobs;


// https://github.com/samtools/htslib/blob/57fe419344cb03e2ea46315443abd242489c32f2/vcf.c#L53
// uint32_t bcf_float_missing    = 0x7F800001;
// uint32_t bcf_float_vector_end = 0x7F800002;

const int PHREDMAX=256;
float pl2ln[PHREDMAX];

// double pl2ln[256];
float pl2ln_f(int32_t & val){
  if(val>=PHREDMAX){
    return log(pow(10.0,-0.1*val));
  } else {
    return pl2ln[val];
  }
}

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

size_t getgls(char*fname,std::vector<double *> &mygl, std::string &vcf_format_field, std::string &vcf_allele_field,char *seek){
  for(int i=0;i<PHREDMAX;i++){    
    pl2ln[i] = log(pow(10.0,-0.1*i));
  }
  htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);
  bcf_hdr_t *hdr = bcf_hdr_read(inf);
  bcf1_t *rec = NULL;rec=bcf_init();assert(rec);
  hts_idx_t *idx=NULL;
  hts_itr_t *iter=NULL;

  if(seek){
    fprintf(stderr,"\t-> Setting iterator to: %s\n",seek);fflush(stderr);
    idx=bcf_index_load(fname);
    iter=bcf_itr_querys(idx,hdr,seek);
  }
  //   http://wresch.github.io/2014/11/18/process-vcf-file-with-htslib.html
  // counters
  int n    = 0;  // total number of records in file
  int nsnp = 0;  // number of SNP records in file
  int nseq = 0;  // number of sequences
  int nsamples = 0;

  // pl data for each call
  int npl_arr = 0;
  int npl     = 0;
  int32_t *pl = NULL;

   // gt data for each call
  int32_t ngt_arr = 0;
  int32_t ngt     = 0;
  int32_t *gt     = NULL;

  
  // af1/af data for each call
  int naf_arr = 0;
  int naf     = 0;
  float *af     = NULL;

  // read header
  nsamples = bcf_hdr_nsamples(hdr);
  //  fprintf(stderr, "\t-> File %s contains %i samples\n", fname, nsamples);
  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(hdr, &nseq); assert(seqnames);//bcf_hdr_id2name(hdr,i)

  char *chr;
  while(1){
    if(seek==NULL){
      if(bcf_read(inf,hdr,rec)!=0)	
	break;
    }else{
      if(bcf_itr_next(inf, iter, rec)!=0)
	break;
    }
    n++;
    if(!bcf_is_snp(rec))
      continue;
    nsnp++;
    fprintf(stderr,"[%d] rec->n_allele:%d ",rec->pos,rec->n_allele);

    for(int i=0;i<rec->n_allele;i++)
      fprintf(stderr,"(%d: %s)",i,rec->d.allele[i]);
    fprintf(stderr,"\n");
    int myreorder[10];
    //funky below took ridiculus long to make
    buildreorder(myreorder,rec->d.allele,rec->n_allele);
    

    double dupergl[10*nsamples]; 

    if(vcf_format_field == "PL") {
      npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
      float ln_gl[npl];
      if(npl<0){
        // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
        // no PL tag is available
        fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching PL tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
        continue;
      }
      // https://github.com/samtools/bcftools/blob/e9c08eb38d1dcb2b2d95a8241933daa1dd3204e5/plugins/tag2tag.c#L151
      fprintf(stderr,"npl:%d\n",npl);
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

      for(int ind=0;ind<nsamples;ind++){
	for(int o=0;o<10;o++)
	  dupergl[ind*10+o] = ln_gl[myreorder[o]]; 
      }
    } else {
      fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, GT or PL\n");
    }
 
    //    fprintf(stderr,"keepind:%d\n",keepInd);

    
    naf = bcf_get_info_float(hdr, rec, vcf_allele_field.c_str(), &af, &naf_arr);
    fprintf(stderr,"rec->pos:%d npl:%d naf:%d rec->n_allele:%d\n",rec->pos,npl,naf,rec->n_allele);    
    //if multiple alt alleles then n_allele>3. We only care about diallelic ref/alt alleless
    //		if(rec->n_allele==4) fprintf(stdout,"\n%s\n",rec->d.allele[2]);
    //ok this is a bit messed up. apparantly sometime the allele is <*> sometimes not.
    // just use the first two alleles now and discard the rest of the alleles.

    //should matter, program should never run on such low freqs, this is just for validation between formats
   
    //filtering
#ifdef __WITH_MAIN__
      fprintf(stdout,"%s\t%i\t%s\t%s\tqual:%f n_info:%d n_allele:%d n_fmt:%d n_sample:%d\n",
	      seqnames[rec->rid],
	      rec->pos+1,
	      rec->d.allele[0],
	      rec->d.allele[1],
	      rec->qual,
	      rec->n_info,
	      rec->n_allele,
	      rec->n_fmt,
	      rec->n_sample);
      
#endif
      
  }
  fprintf(stderr, "\t-> [file=\'%s\'][chr=\'%s\'] Read %i records %i of which were SNPs number of sites with data:%lu\n",fname,seek, n, nsnp,mygl.size()); 

  free(pl);
  free(gt);
  bcf_hdr_destroy(hdr);
  bcf_close(inf);
  bcf_destroy(rec);

  if(iter)
    hts_itr_destroy(iter);
  if(idx)
    hts_idx_destroy(idx);

  //for(int i=0;i<nseq;i++)
    free(seqnames);

  return nsamples;
 
}

void *wrap(void *ptr){
  satan *god = (satan*) ptr;
  god->nind=getgls(god->fname, god->mygl,god->vcf_format_field,god->vcf_allele_field,god->seek);
  //  pthread_exit(NULL);//this is sometimes called without thread
  return NULL;
}

double ** readbcfvcf(char*fname,int &nind, std::string vcf_format_field, std::string vcf_allele_field,char *seek){
  fprintf(stderr,"\t-> readbcfvcf seek:%s nind:%d\n",seek,nind);
  htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);  
  bcf_hdr_t *hdr = NULL;hdr=bcf_hdr_read(inf);assert(hdr);
  int isbcf=0;

  satan god;
  god.fname=fname;
  god.vcf_format_field=vcf_format_field;
  god.vcf_allele_field=vcf_allele_field;
  god.seek=seek;
  

  if(inf->format.format==bcf){
    isbcf=1;
  }if(seek&&isbcf==0){
    fprintf(stderr,"\t-> if choosing region then input file has to be bcf\n");
    exit(0);
  }
 
  double **gls=NULL;
  wrap(&god);
  nind=god.nind;
  gls=new double *[god.mygl.size()];
  for(int i=0;i<god.mygl.size();i++)
    gls[i] = god.mygl[i];

  fprintf(stderr,"Done reading everything we have nsites:%lu for samples:%d\n",god.mygl.size(),nind);
  //merge results
  
  if(hdr) bcf_hdr_destroy(hdr);
  hts_close(inf);
  return gls;
}


#ifdef __WITH_MAIN__
int main(int argc, char **argv) {
  if (argc == 1)
    return 1;

  double **gls=NULL;
  std::string pl=std::string("PL");
  std::string fr=std::string("AFngsrelate");
  char *reg = NULL;
  if(argc==3)
    reg=strdup(argv[2]);
  fprintf(stderr,"reg:%s\n",reg);
  int nind;
  gls = readbcfvcf(argv[1],nind,pl,fr,reg);
  return 0;
}

#endif
