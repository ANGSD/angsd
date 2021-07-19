/*
  This is a class that dumps BCF output. Its mainly a wrapper around other classes here in angsd. Maybe this should be streamlined abit
*/

#include <assert.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "analysisFunction.h"
#include "shared.h"
#include "abcFreq.h"
#include "abcCallGenotypes.h"
#include "abcMajorMinor.h"
#include "abcWriteBcf.h"
#include "aio.h"
#include "abcMcall.h"
#include "abcCounts.h"

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

void abcWriteBcf::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doBcf\t%d\n",doBcf);
  fprintf(argFile,"\t1:  (still beta, not really working)\n");
  fprintf(argFile,"\n\tNB This is a wrapper around -gl -domajorminor and -dopost -dogeno\n");
}

void abcWriteBcf::run(funkyPars *pars){
  if(doBcf==0)
    return ;
}

//last is from abc.h
void print_bcf_header(htsFile *fp,bcf_hdr_t *hdr,argStruct *args,kstring_t &buf,const bam_hdr_t *bhdr,int domcall){
  if(args->ref){
    buf.l = 0;
    ksprintf(&buf, "##reference=file://%s\n", args->ref);
    bcf_hdr_append(hdr,buf.s);
  }
  if (args->anc){
    buf.l = 0;
    ksprintf(&buf, "##ancestral=file://%s\n", args->anc);
    bcf_hdr_append(hdr,buf.s);
  }

  // Translate BAM @SQ tags to BCF ##contig tags
  // todo: use/write new BAM header manipulation routines, fill also UR, M5
  
  for (int i=0; i<bhdr->n_targets; i++){
    buf.l = 0; 
    ksprintf(&buf, "##contig=<ID=%s,length=%d>", bhdr->target_name[i], bhdr->target_len[i]);
    bcf_hdr_append(hdr, buf.s);
  }
  buf.l = 0;
  bcf_hdr_append(hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
  bcf_hdr_append(hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
  bcf_hdr_append(hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of raw reads supporting an indel\">");
  bcf_hdr_append(hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of raw reads supporting an indel\">");
  bcf_hdr_append(hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
  bcf_hdr_append(hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");
  bcf_hdr_append(hdr,"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
  bcf_hdr_append(hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
  bcf_hdr_append(hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
  if(domcall){
    bcf_hdr_append(hdr,"##INFO=<ID=QS_angsd,Number=R,Type=Float,Description=\"Sum of quality scores for A,C,G,T,N \">");
    bcf_hdr_append(hdr,"##INFO=<ID=Quals_angsd ,Number=R,Type=Float,Description=\"Quality score for alternative and null hypothesis\">");
  }
  bcf_hdr_append(hdr,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Minor Allele Frequency\">");
  bcf_hdr_append(hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
  bcf_hdr_append(hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");
  bcf_hdr_append(hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">");
  bcf_hdr_append(hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"scaled Genotype Likelihoods (loglikeratios to the most likely (in log10))\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=EBD,Number=4,Type=Float,Description=\"The effective basedepth\">");
  //fprintf(stderr,"samples from sm:%d\n",args->sm->n);
  buf.l =0;
  assert(args);
  ksprintf(&buf, "##angsdVersion=%s",args->version);
  bcf_hdr_append(hdr, buf.s);
  
  buf.l = 0;
  ksprintf(&buf, "##angsdCommand=%s;Date=%s",args->cmdline,args->datetime);
  bcf_hdr_append(hdr, buf.s);
  buf.l=0;

  for(int i=0;i<args->sm->n;i++){
    bcf_hdr_add_sample(hdr, args->sm->smpl[i]);
  }
  bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
  
  if ( bcf_hdr_write(fp, hdr)!=0 )
    fprintf(stderr,"Failed to write bcf\n");
  buf.l=0;
}

void abcWriteBcf::clean(funkyPars *pars){
  if(doBcf==0)
    return;
}

void abcWriteBcf::print(funkyPars *pars){
  if(doBcf==0)
    return;
 
  //  lh3struct *lh3 = (lh3struct*) pars->extras[5];
  freqStruct *freq = (freqStruct *) pars->extras[7];
  genoCalls *geno = (genoCalls *) pars->extras[11];
  angsd_mcall *mcall = NULL;
  if(domcall)
    mcall = (angsd_mcall *) pars->extras[5];
  for(int s=0;s<pars->numSites;s++) {
    if(pars->keepSites[s]==0)
      continue;
    rec->rid = bcf_hdr_name2id(hdr,header->target_name[pars->refId]);
    rec->pos = pars->posi[s];//<- maybe one index?
    //    bcf_update_id(hdr, rec, "rs6054257");
    if(domcall!=NULL){
      // 
      char alleles[16];
      int goa =0;
      alleles[goa++] = intToRef[mcall->als[4*s]];
      for(int i=1;i<4;i++)
	if(mcall->als[4*s+i]==4)
	  break;
	else{
	  alleles[goa++] = ',';
	  alleles[goa++] = intToRef[mcall->als[4*s+i]];
	}
      alleles[goa] = '\0';
      //    fprintf(stderr,"ALS: %s\n",alleles);
      //exit(0);
      //    ={intToRef[pars->major[s]],',',intToRef[pars->minor[s]],'\0'};
      bcf_update_alleles_str(hdr, rec, alleles);
    }else{
      char majmin[4]={intToRef[pars->major[s]],',',intToRef[pars->minor[s]],'\0'};
      bcf_update_alleles_str(hdr, rec, majmin);
    }
    if(freq==NULL&&mcall==NULL)
      rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = pars->keepSites[s];
    assert(bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1)==0);
    if(mcall){

      float *tmp = mcall->QS+5*s;
      //      for(int i=0;i<5;i++)	fprintf(stderr,"qs[%d]: %f\n",i,tmp[i]);
      assert(bcf_update_info_float(hdr, rec, "QS_angsd", tmp, 5)==0);
      tmp = mcall->quals+2*s;
      assert(bcf_update_info_float(hdr, rec, "Quals_angsd", tmp, 2)==0);
      if(mcall->isvar[s]>0)
	rec->qual = tmp[0];
      else
	rec->qual = tmp[1];
      
    }
    
    if(pars->counts){
      int depth = 0;
      for(int i=0; i<4*pars->nInd; i++)
	depth += pars->counts[s][i];
      tmpi = depth;
      assert(bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1)==0);
    }
    if(freq){
      float tmpf = freq->freq_EM[s];
      assert(bcf_update_info_float(hdr, rec, "AF", &tmpf, 1)==0);
    }
    
    // .. FORMAT
    // assert(geno);
    //    fprintf(stderr,"bcf_hdr_nsamples(hdr): %d\n",bcf_hdr_nsamples(hdr));
    if(geno&&mcall==NULL){
      int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int32_t));
      for(int i=0; i<pars->nInd;i++){
	if(geno->dat[s][i]==-1){
	  tmpia[2*i+0] = bcf_gt_missing;
	  tmpia[2*i+1] = bcf_gt_missing;
	}else if(geno->dat[s][i]==0){
	  tmpia[2*i+0] = bcf_gt_unphased(0);
	  tmpia[2*i+1] = bcf_gt_unphased(0);
	}else if(geno->dat[s][i]==1){
	  tmpia[2*i+0] = bcf_gt_unphased(0);
	  tmpia[2*i+1] = bcf_gt_unphased(1);
	}  else if(geno->dat[s][i]==2){
	  tmpia[2*i+0] = bcf_gt_unphased(1);
	  tmpia[2*i+1] = bcf_gt_unphased(1);
	}
      }
      bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2); 
      free(tmpia);
    }
    else if(mcall!=NULL){
      int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int32_t));
      for(int i=0; i<pars->nInd;i++){
	if(mcall->gcdat[s][2*i]==-1)
	  assert(mcall->gcdat[s][2*i]==mcall->gcdat[s][2*i+1]);

	if(mcall->gcdat[s][2*i]==-1){
	  tmpia[2*i+0] = bcf_gt_missing;
	  tmpia[2*i+1] = bcf_gt_missing;
	}else {
	  tmpia[2*i+0] = bcf_gt_unphased(mcall->gcdat[s][2*i+0]);
	  tmpia[2*i+1] = bcf_gt_unphased(mcall->gcdat[s][2*i+1]);
	}
      }
      bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2); 
      free(tmpia);
    }
    if(pars->counts){
      int32_t *tmpfa = (int32_t*)malloc(sizeof(int32_t)*bcf_hdr_nsamples(hdr));
      suint *ary=pars->counts[s];
      for(int i=0;i<bcf_hdr_nsamples(hdr);i++)
	tmpfa[i] = ary[i*4+0]+ary[i*4+1]+ary[i*4+2]+ary[i*4+3];
      bcf_update_format_int32(hdr, rec, "DP", tmpfa,bcf_hdr_nsamples(hdr) );
      free(tmpfa);
    }
    if(pars->likes[s]){
      float *tmpfa  =   (float*)malloc(3*bcf_hdr_nsamples(hdr)*sizeof(float  ));
      int32_t *tmpi = (int32_t*)malloc(3*bcf_hdr_nsamples(hdr)*sizeof(int32_t));

      int major = pars->major[s];
      int minor = pars->minor[s];
      assert(major!=4&&minor!=4);
      for(int i=0;i<bcf_hdr_nsamples(hdr);i++){
	double val[3];
	val[0] = pars->likes[s][i*10+angsd::majorminor[major][major]];
	val[1] = pars->likes[s][i*10+angsd::majorminor[major][minor]];
	val[2] = pars->likes[s][i*10+angsd::majorminor[minor][minor]];
	
	//fixed awkward case where all gls are -Inf, should only happen with -gl 6
	if(std::isinf(val[0])&&std::isinf(val[1])&&std::isinf(val[2]))
	  val[0]=val[1]=val[2]=0;
	//angsd::logrescale(val,3);	
	for(int j=0;j<3;j++){
	  tmpfa[i*3+j] = val[j]/M_LN10;
	  double tmptmp = -log10(exp(val[j]))*10.0;
	  tmpi[i*3+j] = round(tmptmp);
	}
      }
      bcf_update_format_float(hdr, rec, "GL", tmpfa,3*bcf_hdr_nsamples(hdr) );
      bcf_update_format_int32(hdr, rec, "PL", tmpi,3*bcf_hdr_nsamples(hdr) );
      free(tmpfa);
      free(tmpi);
    }
    if(pars->extras[2]){
      counts *cnts = (counts*) pars->extras[2];
      float *ebd = cnts->ebd[s];
      bcf_update_format_float(hdr, rec, "EBD",ebd ,4*bcf_hdr_nsamples(hdr) );
    }
    if(pars->post&&pars->post[s]){
      float *tmpfa  =   (float*)malloc(3*bcf_hdr_nsamples(hdr)*sizeof(float  ));

      for(int i=0;i<bcf_hdr_nsamples(hdr);i++){
	double *val = pars->post[s] +i*3;
	//	fprintf(stderr,"va: %f %f %f\n",val[0],val[1],val[2]);
	//fixed awkward case where all gls are -Inf, should only happen with -gl 6
	if(std::isinf(val[0])&&std::isinf(val[1])&&std::isinf(val[2]))
	  val[0]=val[1]=val[2]=0;

	//angsd::logrescale(val,3);	
	for(int j=0;j<3;j++){
	  tmpfa[i*3+j] = -log10(val[j])*10.0; 
	}
	double mmax = std::min(tmpfa[i*3],std::min(tmpfa[3*i+1],tmpfa[3*i+2]));
	//	fprintf(stderr,"min:%f\n",mmax);
	for(int j=0;j<3;j++)
	  tmpfa[i*3+j] -= mmax;
	//angsd::logrescale(tmpfa+i*3,3);	
      }
      bcf_update_format_float(hdr, rec, "GP", tmpfa,3*bcf_hdr_nsamples(hdr) );

      free(tmpfa);

    }

    if ( bcf_write1(fp, hdr, rec)!=0 ){
      fprintf(stderr,"Failed to write %s:%d\n",__FILE__,__LINE__);
      exit(0);
    }
    //    fprintf(stderr,"------\n");
    bcf_clear1(rec);
  }
}

void abcWriteBcf::getOptions(argStruct *arguments){

  doBcf=angsd::getArg("-doBcf",doBcf,arguments);
  
  if(doBcf==0)
    return;
  else{
    fprintf(stderr,"\n[bcfoutput] \tPlease add the following parameters \n\t\t '-gl 1 -dopost 1 -domajorminor 1 -domaf 1 -dobcf 1 --ignore-RG 0 -dogeno 1 -docounts 1'\n\n");

  }
  refName = ancName = NULL;
  gl = doMajorMinor = doCounts = doMaf = doPost = doGeno = 0;
  ancName = angsd::getArg("-anc",ancName,arguments);
  refName = angsd::getArg("-ref",refName,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  gl=angsd::getArg("-gl",gl,arguments);
  doMaf = angsd::getArg("-doMaf",doMaf,arguments);
  doCounts = angsd::getArg("-doCounts",doCounts,arguments);
  doGeno = angsd::getArg("-doGeno",doGeno,arguments);
  domcall = angsd::getArg("-doMcall",domcall,arguments);

  if(doPost==0||doMajorMinor==0||gl==0){
    fprintf(stderr,"\nPotential problem. -doBcf is a wrapper around -doMajorMinor -doPost and -gl. These values must therefore be set\n\n");
    //    exit(0);
  }
  
  
}
extern bcf_hdr_t *vcfreader_hs_bcf_hdr;
//constructor
abcWriteBcf::abcWriteBcf(const char *outfiles_a,argStruct *arguments,int inputtype) {
  rec=NULL;
  hdr=NULL;
  fp=NULL;
  doBcf =0;
  args=arguments;
  outfiles=outfiles_a;
  domcall = 0;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doBcf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);

  if(doBcf==0){
    shouldRun[index] =0;
    return;
  }else
    shouldRun[index] =1;
  printArg(arguments->argumentFile);
  
  kstring_t buf;
  buf.s=NULL;
  buf.l=buf.m=0;
  fp=aio::openFileHtsBcf(outfiles,".bcf");
  
  rec    = bcf_init1();
  if(arguments->inputtype==INPUT_BAM){
    hdr = bcf_hdr_init("w");
    print_bcf_header(fp,hdr,args,buf,header,domcall);
  } else if(arguments->inputtype==INPUT_VCF_GP||arguments->inputtype==INPUT_VCF_GL){
    fprintf(stderr,"\t-> Building bcf header\n");
    hdr = bcf_hdr_dup(vcfreader_hs_bcf_hdr);
    bcf_hdr_append(hdr,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"scaled Genotype Likelihoods (loglikeratios to the most likely (in log10))\">");
    bcf_hdr_append(hdr,"##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">");
    bcf_hdr_append(hdr,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    if(doMaf)
      bcf_hdr_append(hdr,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Minor Allele Frequency\">");
    ksprintf(&buf, "##angsdVersion=%s",args->version);
    bcf_hdr_append(hdr, buf.s);
    buf.l = 0;
    ksprintf(&buf, "##angsdCommand=%s;Date=%s",args->cmdline,args->datetime);
    bcf_hdr_append(hdr, buf.s);
    buf.l=0;
    if ( bcf_hdr_write(fp, hdr)!=0 )
      fprintf(stderr,"Failed to write bcf %s:%d\n",__FILE__,__LINE__);
  }
  if(buf.s)
    free(buf.s);
}


abcWriteBcf::~abcWriteBcf(){
  if(rec)
    bcf_destroy(rec);
  if(hdr)
    bcf_hdr_destroy(hdr);
  if(fp!=NULL) hts_close(fp);
}
