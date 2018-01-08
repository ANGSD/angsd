/*
  This is a class that dumps plink file output

 */
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abcFreq.h"
#include "abcCallGenotypes.h"
#include "abcMajorMinor.h"
#include "abcWriteVcf.h"
void abcWriteVcf::printArg(FILE *argFile){
  fprintf(argFile,"------------------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doVcf\t%d\n",doVcf);
  fprintf(argFile,"\t1:  (still beta, not really working)\n");
  fprintf(argFile,"\n\tNB This is a wrapper around -gl -domajorminor and -dopost\n");
}

void abcWriteVcf::run(funkyPars *pars){
  if(doVcf==0)
    return ;
   
}

void abcWriteVcf::clean(funkyPars *pars){
  if(doVcf==0)
    return;
    

}

void abcWriteVcf::print(funkyPars *pars){
  if(doVcf==0)
    return;

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

  aio::bgzf_write(fp,kstr->s,kstr->l);kstr->l=0;
}

void abcWriteVcf::getOptions(argStruct *arguments){

  doVcf=angsd::getArg("-doVcf",doVcf,arguments);
  if(doVcf==0)
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
    fprintf(stderr,"\nPotential problem. -doVcf is a wrapper around -doMajorMinor -doPost and -gl. These values must therefore be set\n\n");
    exit(0);
  }
  
  
}

abcWriteVcf::abcWriteVcf(const char *outfiles,argStruct *arguments,int inputtype){
  fp=NULL;
  doVcf =0;
  kstr=NULL;
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doVcf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);

  if(doVcf==0){
    shouldRun[index] =0;
    return;
  }
  printArg(arguments->argumentFile);

  kstr =(kstring_t*) calloc(1,sizeof(kstring_t));
  //format is taken from: http://faculty.washington.edu/browning/beagle/intro-to-vcf.html
  const char *hdstring= "##fileformat=VCFv4.2(angsd version)\n##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n##INFO=<ID=RA,Number=1,Type=String,Description=\"Reference Allele (included since ANGSD places the MAJOR allele under REF)\">\n##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Minor Allele Frequency\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the major and minor alleles in the order listed\">\n##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">\n##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"Phred-scaled Genotype Likelihoods\">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"scaled Genotype Likelihoods (loglikeratios to the most likely (in log10))\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  ksprintf(kstr,"%s",hdstring);
  for(int i=0;i<arguments->nInd;i++)
    ksprintf(kstr,"\tind%d",i);
  ksprintf(kstr,"\n");

  fp=aio::openFileBG(outfiles,".vcf.gz");
  aio::bgzf_write(fp,kstr->s,kstr->l);kstr->l=0;
}


abcWriteVcf::~abcWriteVcf(){
  if(fp!=NULL) bgzf_close(fp);
  if(kstr && kstr->s)
    free(kstr->s);
}
