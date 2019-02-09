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
#include <limits>
#include <string>
#include <pthread.h>

#define diskio_threads 24
int std_queue = 0;


pthread_mutex_t mymut = PTHREAD_MUTEX_INITIALIZER;
int mycounter = 0; //my semaphore

typedef struct satan_t{
  char*fname;
  int minind;
  double minfreq;
  std::string vcf_format_field;
  std::string vcf_allele_field;
  char *seek;
  std::vector<double *> mygl;
  std::vector<double> freqs;
  int nind;
}satan;

std::vector<satan> jobs;

//populates a vector with the names of which we have data
std::vector<char *> hasdata(char *fname){
  htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);
  bcf_hdr_t *hdr = NULL;hdr=bcf_hdr_read(inf);assert(hdr);
  bcf1_t *rec = NULL;rec=bcf_init();assert(rec);
  hts_idx_t *idx=NULL;idx=bcf_index_load(fname);assert(idx);
  hts_itr_t *iter=NULL;
  int nseq = 0;  // number of sequences
  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(hdr, &nseq); assert(seqnames);//bcf_hdr_id2name(hdr,i)  
  std::vector<char *> ret;
  for(int i=0;i<nseq;i++){
    char buf[strlen(seqnames[i])+100];
    snprintf(buf,strlen(seqnames[i])+100,"%s:1-1000000000",seqnames[i]);
    iter=bcf_itr_querys(idx,hdr,buf);
    if(bcf_itr_next(inf, iter, rec)==0)
      ret.push_back(strdup(seqnames[i]));
    hts_itr_destroy(iter);
  }
  free(seqnames);
  bcf_destroy1(rec);
  bcf_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  hts_close(inf);
  fprintf(stderr,"\t-> Done with preliminary parsing of file: we have data for %lu out of %d reference sequences\n",ret.size(),nseq);
  return ret;
}


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

template <class T>
bool same_val_vcf(T a, T b) {
  return std::fabs(a - b) < std::numeric_limits<T>::epsilon();  
}


// https://en.cppreference.com/w/c/numeric/math/isnan
bool is_nan_vcf(double x) { return x != x; }

//from angsd
double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
     
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;

  
  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      if(std::isnan(sum))
	fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f) W(%f,%f,%f) sum=%f\n",i,loglike[i*3],loglike[i*3+1],loglike[i*3+2],W0,W1,W2,sum);
    }
    
    p=sum/keepInd;
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }
  
  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);
    
    for(int ii=0;1&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1){
	//	fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
	fprintf(stderr,"\n");
      }
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=exp(loglike[i*3+0])*pow(1-p,2);
      W1=exp(loglike[i*3+1])*2*p*(1-p);
      W2=exp(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"[%s.%s():%d] p=%f W %f\t%f\t%f sum=%f loglike: %f\n",__FILE__,__FUNCTION__,__LINE__,p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
      break;
    }
    p=-999;
    assert(p!=999);
    return p;
  }

  return(p);
}


size_t getgls(char*fname,std::vector<double *> &mygl, std::vector<double> &freqs,int minind,double minfreq, std::string &vcf_format_field, std::string &vcf_allele_field,char *seek){
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
    continue;
    //if(rec->n_allele>=3||rec->n_allele==1)//last case shouldt happen
	// continue;
    
    float ln_gl[3*nsamples];    

    if(vcf_format_field == "PL") {
      npl = bcf_get_format_int32(hdr, rec, "PL", &pl, &npl_arr);
      if(npl<0){
        // return codes: https://github.com/samtools/htslib/blob/bcf9bff178f81c9c1cf3a052aeb6cbe32fe5fdcc/htslib/vcf.h#L667
        // no PL tag is available
        fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching PL tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
        continue;
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
    } else if(vcf_format_field == "GT"){
       int ngts = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
       if ( ngts<0 ){
         fprintf(stderr, "BAD SITE %s:%d. return code:%d while fetching GT tag\n", bcf_seqname(hdr,rec), rec->pos, npl);
         continue;
       }
       for(int ns=0; ns<nsamples;ns++){
         int32_t *curr_ptr = gt + ns*2;
         float *ln_gl_ptr = ln_gl + ns*3;
         if ( bcf_gt_is_missing(curr_ptr[0]) ||
              bcf_gt_is_missing(curr_ptr[1]) ){ // obs genotype is missing
           // missing
           ln_gl_ptr[0] = 0;
           ln_gl_ptr[1] = 0;
           ln_gl_ptr[2] = 0;
           
         } else if (bcf_gt_allele(curr_ptr[0])!=bcf_gt_allele(curr_ptr[1])){ // this is obs genotype
           // het
	   ln_gl_ptr[0] = -INFINITY;
	   ln_gl_ptr[1] = 0;
	   ln_gl_ptr[2] = -INFINITY;
	 } else if(bcf_gt_allele(curr_ptr[0])==1){ // this is obs genotype
           // hom alt
	   ln_gl_ptr[0] = -INFINITY;
	   ln_gl_ptr[1] = -INFINITY;
	   ln_gl_ptr[2] = 0;
	 } else{ // this is obs genotype
           // hom ref
	   ln_gl_ptr[0] = 0;
	   ln_gl_ptr[1] = -INFINITY;
	   ln_gl_ptr[2] = -INFINITY;
	 }
       }
    } else {
      fprintf(stderr, "\t\t-> BIG TROUBLE. Can only take one of two tags, GT or PL\n");
    }
    
    int keepInd=0;
    char keep[nsamples];
    double *tmp = new double[3*nsamples];    
    for(int ns=0;ns<nsamples;ns++){
      float *ary= ln_gl+ns*3;
      if ((is_nan_vcf(ary[0]) || is_nan_vcf(ary[1]) || is_nan_vcf(ary[2])) ||(same_val_vcf(ary[0], ary[1]) && same_val_vcf(ary[0], ary[2]))){
        keep[ns]=0;
      }else{
	keep[ns]=1;
	keepInd++;
      }
      tmp[ns*3] = ary[0];
      tmp[ns*3+1] = ary[1];
      tmp[ns*3+2] = ary[2];
      // fprintf(stderr, "TMP: %d %d: %f %f %f\n", ns+1, rec->pos+1, tmp[ns*3], tmp[ns*3+1], tmp[ns*3+2]);
    }
    //    fprintf(stderr,"keepind:%d\n",keepInd);

    
    naf = bcf_get_info_float(hdr, rec, vcf_allele_field.c_str(), &af, &naf_arr);
    fprintf(stderr,"rec->pos:%d npl:%d naf:%d rec->n_allele:%d\n",rec->pos,npl,naf,rec->n_allele);    
    //if multiple alt alleles then n_allele>3. We only care about diallelic ref/alt alleless
    //		if(rec->n_allele==4) fprintf(stdout,"\n%s\n",rec->d.allele[2]);
    //ok this is a bit messed up. apparantly sometime the allele is <*> sometimes not.
    // just use the first two alleles now and discard the rest of the alleles.

    double freq;
    if(naf==1){
      freq = af[0];
    }else{
      freq = emFrequency(tmp,nsamples,50,0.05,keep,keepInd);
    }
    //should matter, program should never run on such low freqs, this is just for validation between formats
    if(freq>0.999)
      freq=1;
    if(freq<0.001)
      freq=0;
    //fprintf(stderr,"freq:%f minfreq:%f keepInd:%d minind:%d\n",freq,minfreq,keepInd,minind);
    //filtering
    if(keepInd>minind&&freq>=minfreq && freq<= (1-minfreq)) {

#ifdef __WITH_MAIN__
      fprintf(stdout,"%s\t%i\t%s\t%s\tqual:%f n_info:%d n_allele:%d n_fmt:%d n_sample:%d n_samplws_with_data:%d freq:%f",
	      seqnames[rec->rid],
	      rec->pos+1,
	      rec->d.allele[0],
	      rec->d.allele[1],
	      rec->qual,
	      rec->n_info,
	      rec->n_allele,
	      rec->n_fmt,
	      rec->n_sample,
	      keepInd,
	      freq
	      );
      for(int i=0;i<3*nsamples;i++)
	fprintf(stdout," %f",ln_gl[i]);
      fprintf(stdout,"\n");
#endif
      for(int ns=0;ns<3*nsamples;ns++)
	tmp[ns]=exp(tmp[ns]);
      mygl.push_back(tmp);
      freqs.push_back(freq);
      //populate debug names

    } else {
      delete [] tmp;
    }
    // fprintf(stderr,"rec->pos:%d npl:%d naf:%d rec->n_allele:%d af[0]:%f\n",rec->pos,npl,naf,rec->n_allele,freq);
    // exit(0);
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
  god->nind=getgls(god->fname, god->mygl,god->freqs, god->minind, god->minfreq,god->vcf_format_field,god->vcf_allele_field,god->seek);
  //  pthread_exit(NULL);//this is sometimes called without thread
}


void *wrap2(void *){
  while(1){
    pthread_mutex_lock(&mymut);
    int myvar = mycounter;
    mycounter++;
    pthread_mutex_unlock(&mymut);
    if(myvar>=jobs.size())
      pthread_exit(NULL);
    satan *god = (satan*) &jobs[myvar];
    god->nind=getgls(god->fname, god->mygl,god->freqs, god->minind, god->minfreq,god->vcf_format_field,god->vcf_allele_field,god->seek);
  }
}


double ** readbcfvcf(char*fname,int &nind, std::vector<double> &freqs,int minind,double minfreq, std::string vcf_format_field, std::string vcf_allele_field,char *seek){
  fprintf(stderr,"\t-> readbcfvcf seek:%s nind:%d\n",seek,nind);
  htsFile * inf = NULL;inf=hts_open(fname, "r");assert(inf);  
  bcf_hdr_t *hdr = NULL;hdr=bcf_hdr_read(inf);assert(hdr);
  int isbcf=0;
  std::vector<char *> hd;

  satan god;
  god.fname=fname;
  god.minind=minind;
  god.minfreq=minfreq;
  god.vcf_format_field=vcf_format_field;
  god.vcf_allele_field=vcf_allele_field;
  god.seek=seek;
  

  if(inf->format.format==bcf){
    isbcf=1;
    hd = hasdata(fname);
  }if(seek&&isbcf==0){
    fprintf(stderr,"\t-> if choosing region then input file has to be bcf\n");
    exit(0);
  }

  int nseq = 0;  // number of sequences
  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(hdr, &nseq); 
  assert(seqnames);
  
  double **gls=NULL;
  if(seek!=NULL||isbcf==0){//single run
    wrap(&god);
    nind=god.nind;
   
    gls=new double *[god.mygl.size()];
    for(int i=0;i<god.mygl.size();i++){
      gls[i] = god.mygl[i];
    }
    freqs=god.freqs;
  }else{
      for(int i=0;i<hd.size();i++){
	jobs.push_back(god);
	jobs[i].seek=hd[i];
      }
    
    if(diskio_threads==1||isbcf==0){
      for(int i=0;i<jobs.size();i++){
	wrap(&jobs[i]);
      }
    }else{
      int at =0;

      while(at<hd.size()){
	//	fprintf(stderr,"at:%d hdsize:%lu\n",at,hd.size());
	int howmany=std::min(diskio_threads,(int)hd.size()-at);
	pthread_t mythd[howmany];
	for(int i=0;i<howmany;i++){
	  if(std_queue)
	    assert (pthread_create(&mythd[i],NULL,wrap,&jobs[i+at])==0);
	  else
	    assert (pthread_create(&mythd[i],NULL,wrap2,NULL)==0);
	}
	for(int i=0;i<howmany;i++)
	  assert(pthread_join(mythd[i],NULL)==0);
	at+=howmany;
      }
    }
    int nsites =0;
    for(int i=0;i<jobs.size();i++)
      nsites += jobs[i].mygl.size();
    nind=jobs[0].nind;
    //    fprintf(stderr,"Done reading everything we have nsites:%d for samples:%d\n",nsites,nind);
    //merge results
    gls=new double *[nsites];
    freqs.reserve(nsites);
    int at =0;
    for(int i=0;i<jobs.size();i++){
      for(int j=0;j<jobs[i].mygl.size();j++)
	gls[at++] = jobs[i].mygl[j];
      freqs.insert(freqs.end(),jobs[i].freqs.begin(),jobs[i].freqs.end());
    }
    for(int i=0;i<hd.size();i++)
      free(hd[i]);
  }
  free(seqnames);
  if(hdr) bcf_hdr_destroy(hdr);
  hts_close(inf);
  return gls;
}


#ifdef __WITH_MAIN__

int main(int argc, char **argv) {
  if (argc == 1)
    return 1;

  double **gls=NULL;
  int nind;
  std::vector<double> freqs;
  std::string pl=std::string("PL");
  std::string fr=std::string("AFngsrelate");
  char *reg = NULL;
  if(argc==3)
    reg=strdup(argv[2]);
  fprintf(stderr,"reg:%s\n",reg);
  gls = readbcfvcf(argv[1],nind,freqs,2,0.04,pl,fr,reg);
  return 0;
}

#endif
