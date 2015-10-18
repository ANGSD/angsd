/*
  The GLs directly from the qscore, assuming diploid.

  This is the one used in the original GATK paper. And is most likely very different from the newest gatk

 */


#include <cmath>
#include <cstdio>
#include "bambi_interface.h"
#include "analysisFunction.h"
#include "phys_likes.h"
#include "phys_genolike_calc.h"

extern int refToInt[256];
double **phys_probs = NULL;
//char *phys_refpath;
phys_genolike_calc *glc = NULL;
/*
  allocate the prob matrix with
  prob[1][qs] = homohit
  prob[2][qs] = hethit
  prob[3][qs] = homoNonhit

  I'm perfecly aware that these are not genotypes, but what should these be called
 */

double **genLikes_phys(int len){
  double **ret = new double*[3];
  for(int i=0;i<3;i++)
    ret[i] = new double[256];
  
  for(int i=0;i<len;i++){
    double p = 1-pow(10,-i/10.0);
    double p3 = (1-p)/3.0;
    ret[0][i] = log(p*0.5+p*0.5);
    ret[1][i] = log(p*0.5+p3*0.5);
    ret[2][i] = log(0.5*p3+0.5*p3);
  }
  return ret;
}

int offsetss[4][10]={
  {0,1,2,3,  4,5,6,7,8,9},//AA,AC,AG,AT,therest
  {4,1,5,6,  0,2,3,7,8,9},//CC,AC,CG,CT,therest
  {7,2,5,8,  0,1,3,4,6,9},//GG,AG,CG,GT,therest
  {9,3,6,8,  0,1,2,4,5,7},//TT,AT,CT,GT,therest
};


void phys_init(std::vector<char *> bamnames){
  glc = new phys_genolike_calc();
  char phys_refpath[4096];

  for(size_t i=0;i<bamnames.size();i++){
    snprintf(phys_refpath,4096,"%s.phys",bamnames[i]);
    fprintf(stderr,"Reading pars bam:\n\t\t%s\n\t Assuming pars name is:\n\t\t%s\n",bamnames[i],phys_refpath);
    if(angsd::fexists(bamnames[i])==0){
      fprintf(stderr,"bamfile:%s doesnt exists",bamnames[i]);
      exit(0);
    }
    if(angsd::fexists(phys_refpath)==0){
      fprintf(stderr,"coefficient filefile:%s doesnt exists\n",phys_refpath);
      exit(0);
    }
    glc->read_coef(phys_refpath);
  }
  phys_probs = genLikes_phys(256);

  //  

  

}

void phys_destroy(){

  for(int i=0;i<3;i++)
    delete [] phys_probs[i];
  delete [] phys_probs;
  phys_probs=NULL;
  //  delete [] phys_refpath;
}

void call_phys(chunkyT *chk,double **lk,int trim){

  //  new phys_genolike_calc( phys_refpath );
  //  glc->set_debug( true );

  for(int s=0;s<chk->nSites;s++){
    for(int i=0;i<chk->nSamples;i++){

      //allocate for all samples, for first sample
      if(i==0){
        lk[s] = new double[10*chk->nSamples];
        for(int ii=0;ii<10*chk->nSamples;ii++)
          lk[s][ii] = -0.0;//set default values
      }

      // Get pointer to current 'sample' genotypes
      double *likes1 = lk[s]+10*i;

      // Fill geno_probs array with the 10 values found for site s and smaple i
      glc->get_genolikes( s, i, chk, likes1 );

      
    }
  }

  // delete glc;

}
