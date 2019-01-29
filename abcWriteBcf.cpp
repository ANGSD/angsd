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
  fprintf(stderr,"-dobcf2000:%d\n",doBcf);
  if(doBcf==0)
    return ;
   
}

void print_bcf_header(htsFile *fp,bcf_hdr_t *hdr,argStruct *args){

  bcf_hdr_append(hdr, "##fileDate=20090805");
  bcf_hdr_append(hdr, "##FORMAT=<ID=UF,Number=1,Type=Integer,Description=\"Unused FORMAT\">");
  bcf_hdr_append(hdr, "##INFO=<ID=UI,Number=1,Type=Integer,Description=\"Unused INFO\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=Flt,Description=\"Unused FILTER\">");
  bcf_hdr_append(hdr, "##unused=<XX=AA,Description=\"Unused generic\">");
  bcf_hdr_append(hdr, "##unused=unformatted text 1");
  bcf_hdr_append(hdr, "##unused=unformatted text 2");
  bcf_hdr_append(hdr, "##contig=<ID=Unused,length=62435964>");
  bcf_hdr_append(hdr, "##source=myImputationProgramV3.1");
  bcf_hdr_append(hdr, "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta");
  bcf_hdr_append(hdr, "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>");
  bcf_hdr_append(hdr, "##phasing=partial");
  bcf_hdr_append(hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
  bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
  bcf_hdr_append(hdr, "##INFO=<ID=NEG,Number=.,Type=Integer,Description=\"Test Negative Numbers\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">");
  bcf_hdr_append(hdr, "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">");
  bcf_hdr_append(hdr, "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=q10,Description=\"Quality below 10\">");
  bcf_hdr_append(hdr, "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=TS,Number=1,Type=String,Description=\"Test String\">");

  for(int i=0;i<args->sm->n;i++){
    fprintf(stderr,"arguments->sm->sampl[9]:%s\n",args->sm->smpl[i]);
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
  fprintf(stderr,"-dobcf:%d\n",doBcf);
  if(doBcf==0)
    return;
  if(fp==NULL){
    fp=aio::openFileHts(outfiles,".bcf");
    hdr = bcf_hdr_init("w");
    rec    = bcf_init1();
    print_bcf_header(fp,hdr,args);
  }
  lh3struct *lh3 = (lh3struct*) pars->extras[5];
  freqStruct *freq = (freqStruct *) pars->extras[6];
  genoCalls *geno = (genoCalls *) pars->extras[10];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    //chr pos id
    ksprintf(kstr,"%s\t%d\t.\t",header->target_name[pars->refId],pars->posi[s]+1);
    kputc(intToRef[pars->major[s]],kstr);kputc('\t',kstr);
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

abcWriteBcf::abcWriteBcf(const char *outfiles_a,argStruct *arguments,int inputtype){
  fp=NULL;
  doBcf =0;
  kstr=NULL;
  outfiles=outfiles_a;
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
  if(kstr && kstr->s)
    free(kstr->s);
}
