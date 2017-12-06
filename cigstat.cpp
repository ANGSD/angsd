#include <cassert>
#include <cstdio>
#include <cstring>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include "analysisFunction.h"
#include "cigstat.h"
#define max_rlen 1024

size_t **cig_counts = NULL;
BGZF *bg = NULL;


char **popnams(){
char **BAM_NAM = NULL;
BAM_NAM=new char*[10];
int at=0;
BAM_NAM[at++] = strdup("BAM_CMATCH");
BAM_NAM[at++] = strdup("BAM_CINS");
BAM_NAM[at++] = strdup("BAM_CDEL");     
BAM_NAM[at++] = strdup("BAM_CREF_SKIP");
BAM_NAM[at++] = strdup("BAM_CSOFT_CLIP");
BAM_NAM[at++] = strdup("BAM_CHARD_CLIP");
BAM_NAM[at++] = strdup("BAM_CPAD");
BAM_NAM[at++] = strdup("BAM_CEQUAL");
BAM_NAM[at++] = strdup("BAM_CDIFF");
BAM_NAM[at++] = strdup("BAM_CBACK");
return BAM_NAM;
}


int cigstat_init(const char *fname){
  fprintf(stderr,"\t-> %s.%s():%d\n",__FILE__,__FUNCTION__,__LINE__);
bg = aio::openFileBG(fname,".cigstat.gz");
fprintf(stderr,"\t-> Will do basic cigstat\n");
cig_counts = new size_t*[10];
for(int i=0;i<10;i++){
  cig_counts[i] = new size_t[max_rlen];
memset(cig_counts[i],0,sizeof(size_t)*max_rlen);
}
  return 0;
}
int cigstat_calc(bam1_t *rd){
  //    fprintf(stderr,"\t-> %s.%s():%d\n",__FILE__,__FUNCTION__,__LINE__);
 int nCig = rd->core.n_cigar;

    uint32_t *cigs = bam_get_cigar(rd);

int cum =0;
    for(int i=0;i<nCig;i++) {
      //      if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
      int opCode = cigs[i]&BAM_CIGAR_MASK; //what to do
      int opLen = cigs[i]>>BAM_CIGAR_SHIFT; //length of what to do
for(int s=0;s<opLen;s++){
  assert(s+cum<1024);
cig_counts[opCode][s+cum]++;
}
cum += opLen;


}

return 0;
}
int cigstat_close(){
  fprintf(stderr,"\t-> %s.%s():%d\n",__FILE__,__FUNCTION__,__LINE__);
  kstring_t *ks=(kstring_t*)calloc(1,sizeof(kstring_t));
  ksprintf(ks,"CIGAR");
  for(int i=0;i<max_rlen;i++)
    ksprintf(ks,"\tPos%d",i);
  //ksprintf(ks,"\n");
  char **bam_nam=popnams();
  for(int ncig=0;ncig<10;ncig++){
    ksprintf(ks,"\n%s",bam_nam[ncig]);
    for(int p=0;p<1024;p++)
      ksprintf(ks,"\t%lu",cig_counts[ncig][p]);
  }
  ksprintf(ks,"\n");
  assert(ks->l==bgzf_write(bg,ks->s,ks->l));
  bgzf_close(bg);
  return 0;
}
