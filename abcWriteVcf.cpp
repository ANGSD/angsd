/*
  This is a class that dumps plink file output

 */
#include <assert.h>

#include "analysisFunction.h"
#include "shared.h"
#include <htslib/kstring.h>
#include "abcCallGenotypes.h"
#include "abcMajorMinor.h"
#include "abcWriteVcf.h"
#include <zlib.h>
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

  lh3struct *lh3 =(lh3struct*) pars->extras[5];

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0)
      continue;
    //chr pos id
    ksprintf(&kstr,"%s\t%d\t.\t",header->target_name[pars->refId],pars->posi[s]+1);
    kputc(intToRef[pars->major[s]],&kstr);kputc('\t',&kstr);
    kputc(intToRef[pars->minor[s]],&kstr);kputc('\t',&kstr);
    ksprintf(&kstr,".\tPASS\t.\tGL:GP\t");
    for(int i=0;i<pars->nInd;i++){
      ksprintf(&kstr,"%f,%f%f:",pars->post[s][i*3+0],pars->post[s][i*3+1],pars->post[s][i*3+2]);
      ksprintf(&kstr,"%f,%f%f",lh3->lh3[s][i*3+0],lh3->lh3[s][i*3+1],lh3->lh3[s][i*3+2]);
      if(i<pars->nInd-1)
	ksprintf(&kstr,"\t");
    }
    ksprintf(&kstr,"\n");
  }
  
  gzwrite(fp,kstr.s,kstr.l);kstr.l=0;
}

void abcWriteVcf::getOptions(argStruct *arguments){

  doVcf=angsd::getArg("-doVcf",doVcf,arguments);
  if(doVcf==0)
    return;
  int doPost = 0;
  int doMajorMinor =0;
  int gl =0;
  doPost=angsd::getArg("-doPost",doPost,arguments);
  doMajorMinor=angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  gl=angsd::getArg("-gl",gl,arguments);

  if(doPost==0||doMajorMinor==0||gl==0){
    fprintf(stderr,"\nPotential problem. -doVcf is a wrapper around -gl -doPost and -gl. These values must therefore be set\n\n");
    exit(0);
  }
  
  
}

abcWriteVcf::abcWriteVcf(const char *outfiles,argStruct *arguments,int inputtype){
  fp=Z_NULL;
  doVcf =0;
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-doVcf")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doVcf==0){
    shouldRun[index] =0;
    return;
  }
  kstr.l=kstr.m=0;kstr.s=NULL;
  //format is taken from: http://faculty.washington.edu/browning/beagle/intro-to-vcf.html
  const char *hdstring= "##fileformat=VCFv4.2(angsd version)\n##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">\n##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"Phred-scaled Genotype Likelihoods\">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"scaled Genotype Likelihoods (these are really llh eventhough they sum to one)\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  ksprintf(&kstr,"%s",hdstring);
  for(int i=0;i<arguments->nInd;i++)
    ksprintf(&kstr,"\tind%d",i);
  ksprintf(&kstr,"\n");

  fp=aio::openFileGz(outfiles,".vcf.gz",GZOPT);
  gzwrite(fp,kstr.s,kstr.l);kstr.l=0;
}


abcWriteVcf::~abcWriteVcf(){
  if(fp!=Z_NULL) gzclose(fp);
}


