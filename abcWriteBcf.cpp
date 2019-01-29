/*
  This is a class that dumps plink file output
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
  fprintf(argFile,"\n\tNB This is a wrapper around -gl -domajorminor and -dopost\n");
}

void abcWriteBcf::run(funkyPars *pars){
  if(doBcf==0)
    return ;
 
}

//last is from abc.h
void print_bcf_header(htsFile *fp,bcf_hdr_t *hdr,argStruct *args,kstring_t &buf,const bam_hdr_t *bhdr){
  assert(args);
  ksprintf(&buf, "##angsdVersion=%s+htslib-%s\n",angsd_version(),hts_version());
  bcf_hdr_append(hdr, buf.s);
  
  buf.l = 0;
  ksprintf(&buf, "##angsdCommand=");
  for (int i=1; i<args->argc; i++)
    ksprintf(&buf, " %s", args->argv[i]);
  kputc('\n', &buf);
  bcf_hdr_append(hdr, buf.s);
  buf.l=0;

  if (args->ref){
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

  bcf_hdr_append(hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");

  bcf_hdr_append(hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
  bcf_hdr_append(hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
  bcf_hdr_append(hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
  bcf_hdr_append(hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
  bcf_hdr_append(hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
  bcf_hdr_append(hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
  bcf_hdr_append(hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths\">");
  bcf_hdr_append(hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand\">");
  bcf_hdr_append(hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand\">");
  //fprintf(stderr,"samples from sm:%d\n",args->sm->n);
  for(int i=0;i<args->sm->n;i++){
    bcf_hdr_add_sample(hdr, args->sm->smpl[i]);
  }
  bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
  
  if ( bcf_hdr_write(fp, hdr)!=0 )
    fprintf(stderr,"Failed to write bcf\n");
}



void abcWriteBcf::clean(funkyPars *pars){
  if(doBcf==0)
    return;
    

}

void abcWriteBcf::print(funkyPars *pars){
  if(doBcf==0)
    return;
  if(fp==NULL){
    fp=aio::openFileHts(outfiles,".bcf");
    hdr = bcf_hdr_init("w");
    rec    = bcf_init1();
    print_bcf_header(fp,hdr,args,buf,header);
  }
  lh3struct *lh3 = (lh3struct*) pars->extras[5];
  freqStruct *freq = (freqStruct *) pars->extras[6];
  genoCalls *geno = (genoCalls *) pars->extras[10];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
  }
#if 0
    rec->rid = bcf_hdr_name2id(hdr,header->target_name[pars->refId]);
    rec->pos = pars->posi[s]+1;//<- maybe one index?
    //    bcf_update_id(hdr, rec, "rs6054257");
    char majmin[4]={intToRef[pars->major[s]],',',intToRef[pars->minor[s]],'\0'};
    bcf_update_alleles_str(hdr, rec, majmin);
    rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = 3;
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 14;
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    tmpi = -127;
    bcf_update_info_int32(hdr, rec, "NEG", &tmpi, 1);
    float tmpf = 0.5;
    bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
    
        // .. FORMAT
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_unphased(0);
    tmpia[1] = bcf_gt_unphased(0);
    tmpia[2] = bcf_gt_unphased(1);
    tmpia[3] = bcf_gt_unphased(0);
    tmpia[4] = bcf_gt_unphased(1);
    tmpia[5] = bcf_gt_unphased(1);
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
    char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3);
    if ( bcf_write1(fp, hdr, rec)!=0 )
      error("Failed to write to \n");
  
#endif

   rec->rid = bcf_hdr_name2id(hdr, "20");
    // .. POS
    rec->pos = 14369;
    // .. ID
    bcf_update_id(hdr, rec, "rs6054257");
    // .. REF and ALT
    bcf_update_alleles_str(hdr, rec, "G,A");
    // .. QUAL
    rec->qual = 29;
    // .. FILTER
    int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
    bcf_update_filter(hdr, rec, &tmpi, 1);
    // .. INFO
    tmpi = 3;
    bcf_update_info_int32(hdr, rec, "NS", &tmpi, 1);
    tmpi = 14;
    bcf_update_info_int32(hdr, rec, "DP", &tmpi, 1);
    tmpi = -127;
    bcf_update_info_int32(hdr, rec, "NEG", &tmpi, 1);
    float tmpf = 0.5;
    bcf_update_info_float(hdr, rec, "AF", &tmpf, 1);
    bcf_update_info_flag(hdr, rec, "DB", NULL, 1);
    bcf_update_info_flag(hdr, rec, "H2", NULL, 1);
    // .. FORMAT
    int32_t *tmpia = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
    tmpia[0] = bcf_gt_phased(0);
    tmpia[1] = bcf_gt_phased(0);
    tmpia[2] = bcf_gt_phased(1);
    tmpia[3] = bcf_gt_phased(0);
    tmpia[4] = bcf_gt_unphased(1);
    tmpia[5] = bcf_gt_unphased(1);
    bcf_update_genotypes(hdr, rec, tmpia, bcf_hdr_nsamples(hdr)*2);
    tmpia[0] = 48;
    tmpia[1] = 48;
    tmpia[2] = 43;
    bcf_update_format_int32(hdr, rec, "GQ", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 1;
    tmpia[1] = 8;
    tmpia[2] = 5;
    bcf_update_format_int32(hdr, rec, "DP", tmpia, bcf_hdr_nsamples(hdr));
    tmpia[0] = 51;
    tmpia[1] = 51;
    tmpia[2] = 51;
    tmpia[3] = 51;
    tmpia[4] = bcf_int32_missing;
    tmpia[5] = bcf_int32_missing;
    bcf_update_format_int32(hdr, rec, "HQ", tmpia, bcf_hdr_nsamples(hdr)*2);
    char *tmp_str[] = {"String1","SomeOtherString2","YetAnotherString3"};
    bcf_update_format_string(hdr, rec, "TS", (const char**)tmp_str, 3);
    if ( bcf_write1(fp, hdr, rec)!=0 )
      error("Failed to write to %s");
#if 0
  for(int s=0;s<pars->numSites;s++){
    //chr pos id
       kputc(intToRef[pars->minor[s]],kstr);kputc('\t',kstr);
    ksprintf(kstr,".\tPASS\tNS=%d", pars->keepSites[s]);
    // Total per site depth
    if(doCounts != 0){
      int depth = 0;
      for(int i=0; i<4*pars->nInd; i++)
	depth += pars->counts[s][i];
      ksprintf(kstr,";DP=%d", depth);
    }
    if(refName != NULL)
      ksprintf(kstr,";RA=%c", intToRef[pars->ref[s]]);
    if(ancName != NULL)
      ksprintf(kstr,";AA=%c", intToRef[pars->anc[s]]);

    // MAF
    if(doMaf != 0)
      ksprintf(kstr,";AF=%f", freq->freq_EM[s]);
    // GP and GL
    kputc('\t',kstr);
    if(doGeno != 0)
      ksprintf(kstr,"GT:");
    if(doCounts != 0)
      ksprintf(kstr, "DP:AD:");
    ksprintf(kstr,"GP:GL");
    // Per-indiv data
    for(int i=0; i<pars->nInd;i++){
      kputc('\t',kstr);
      if(doGeno != 0){
	int g = geno->dat[s][i];
	int gg[2] = {'0','1'};
	if(g == 0)
	  gg[0] = gg[1] = '0';
	else if(g == 1);
	else if(g == 2)
	  gg[0] = gg[1] = '1';
	else
	  gg[0] = gg[1] = '.';
	ksprintf(kstr, "%c/%c:", gg[0], gg[1]);
      }
      if(doCounts != 0)
        ksprintf(kstr,"%d:%d,%d:", pars->counts[s][i*4+pars->major[s]]+pars->counts[s][i*4+pars->minor[s]], pars->counts[s][i*4+pars->major[s]], pars->counts[s][i*4+pars->minor[s]]);
      ksprintf(kstr,"%f,%f,%f:",pars->post[s][i*3+0],pars->post[s][i*3+1],pars->post[s][i*3+2]);
      ksprintf(kstr,"%f,%f,%f",lh3->lh3[s][i*3+0]/M_LN10,lh3->lh3[s][i*3+1]/M_LN10,lh3->lh3[s][i*3+2]/M_LN10);
    }
    ksprintf(kstr,"\n");
  }
#endif
  // aio::bgzf_write(fp,kstr->s,kstr->l);kstr->l=0;
}

void abcWriteBcf::getOptions(argStruct *arguments){

  doBcf=angsd::getArg("-doBcf",doBcf,arguments);
  
  if(doBcf==0)
    return;
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

  if(doPost==0||doMajorMinor==0||gl==0){
    fprintf(stderr,"\nPotential problem. -doBcf is a wrapper around -doMajorMinor -doPost and -gl. These values must therefore be set\n\n");
    exit(0);
  }
  
  
}
//constructor
abcWriteBcf::abcWriteBcf(const char *outfiles_a,argStruct *arguments,int inputtype){
  fp=NULL;
  doBcf =0;
  args=arguments;
  outfiles=outfiles_a;
  buf.s=NULL;buf.l=buf.m=0;
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
  
  
  
}


abcWriteBcf::~abcWriteBcf(){
  if(fp!=NULL) hts_close(fp);
}
