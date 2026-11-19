/*
  -GL 5 (phys): genotype likelihoods with a per sample model of sequencing
  errors and ancient DNA damage, read from <bam>.phys for every bam file.
  See phys_genolike_calc.h.

  DEPRECATED: will be removed in a future version.
 */

#include <cstdio>
#include <cstdlib>
#include "analysisFunction.h"
#include "phys_likes.h"
#include "phys_genolike_calc.h"

static phys_genolike_calc *glc = NULL;

void phys_init(std::vector<char *> bamnames){
  fprintf(stderr,"\t-> -GL 5 (phys) is deprecated and will be removed in a future version\n");
  glc = new phys_genolike_calc();
  char phys_refpath[4096];

  for(size_t i=0;i<bamnames.size();i++){
    snprintf(phys_refpath,4096,"%s.phys",bamnames[i]);
    if(angsd::fexists(bamnames[i])==0){
      fprintf(stderr,"\t-> bamfile: %s doesnt exist\n",bamnames[i]);
      exit(1);
    }
    if(angsd::fexists(phys_refpath)==0){
      fprintf(stderr,"\t-> coefficient file: %s doesnt exist\n",phys_refpath);
      exit(1);
    }
    glc->read_coef(phys_refpath);
  }
}

void phys_destroy(){
  delete glc;
  glc = NULL;
}

void call_phys(chunkyT *chk,double **lk,int trim){
  for(int s=0;s<chk->nSites;s++){
    lk[s] = new double[10*chk->nSamples];
    for(int ii=0;ii<10*chk->nSamples;ii++)
      lk[s][ii] = -0.0;//set default values
    for(int i=0;i<chk->nSamples;i++)
      glc->get_genolikes( s, i, chk, trim, lk[s]+10*i );
  }
}
