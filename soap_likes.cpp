/*
  This file implements soapsnp GL estimators.
  
  This is a 2 step precedure
  1) estimate a calibration matrix, based on a mismatch matrix (qs,position,refe,observed)
  2) from number one estimate GL, based on observed data
  
  it will dump the countsmatrix and the calibration matrix

 */
#include <map>
#include <cmath>
#include <fstream>
#include <libgen.h>//basename
#include <sys/stat.h>//mkdir
#include <sys/types.h>//mkdir
#include <vector>
#include <cassert>
#include <limits> //<- for setting nan
#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <pthread.h>

#include "bambi_interface.h"
#include "soap_likes.h"
#include "analysisFunction.h"

//calc_offset into an array of dim 256 256 4 4
int co(int a, int b, int c, int d){
  return (a<<12 | b<<4 | c<<2|d );
}

#define LENSS 100000

int conv[4] = {0,1,3,2};

extern int refToInt[256];

typedef struct{
  char qs;
  char pos;
  int isUpper;
}tmpStruct;

//traverse qs(hight->low),posi(low->high),BIGLETTER->smallletters

struct cmp_tmpStruct {
  bool operator()(const tmpStruct& first,const  tmpStruct& second) const {
    if(first.qs!=second.qs)
      return first.qs>second.qs;
    else if (first.pos!=second.pos)
      return first.pos<second.pos;
    else
      return first.isUpper >second.isUpper;
  }
};



typedef std::multimap<tmpStruct,int ,cmp_tmpStruct> soapMap;

int mySum(int ary[4]){
  int tsum =0;
  for(int i=0;i<4;i++)
    tsum += ary[i];
  return tsum;
}



//this function does the GL estimation, must supply a calibration matrix
int calc_gl(chunkyT *chk,double *p_matrix,int which_sample,double **lk,int trim) {
  int read_len = 256;
  //  fprintf(stderr,"[%s]\n",__FUNCTION__);


  double pcr_dependency = log10(0.5);
  double global_dependency = log10(0.9);

  for(int s=0;s<chk->nSites;s++) {
    if(which_sample==0){
      //      fprintf(stderr,"doubleling=%d\n",s);
      lk[s] = new double[10*chk->nSamples];
      for(int i=0;i<10*chk->nSamples;i++)
	lk[s][i] = -0.0;
    }
    int global_dep_count[4] = {-1,-1,-1,-1};
    int pcr_dep_count[4][2*read_len];
    for(int ii=0;ii<4;ii++) 
      memset(pcr_dep_count[ii],0,sizeof(int)*2*read_len);
    tNode *nd = chk->nd[s][which_sample];
    double type_likely[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

    soapMap aMap;

    //roll values into the map
for(int j=0;nd&&j<nd->l;j++) {
      //soapsnp only allows 255s
      if(aMap.size()>=0xFF)
	break;
      tmpStruct t;
      if(nd->qs[j]<0){
	fprintf(stderr,"Meaningless low quality\n");
	continue;
      }
      if(nd->posi[j]<trim||nd->isop[j]<trim|| refToInt[nd->seq[j]]==4 ){
	continue;
      }


      t.qs=nd->qs[j];
      if(1&&isupper(nd->seq[j]))
	t.pos= nd->posi[j];
      else
	t.pos= nd->isop[j];
      t.isUpper = isupper(nd->seq[j]);
      aMap.insert(std::pair<tmpStruct,int>(t,j));
    }


    //    for(int j=0;j<nd.l;j++) {
    for(soapMap::iterator it=aMap.begin();it!=aMap.end();++it) {
      int j=it->second;
      int workId = isupper(nd->seq[j])?nd->posi[j]:read_len+nd->isop[j];
      int alleleID = conv[refToInt[nd->seq[j]]];//{0,1,2,3}
      
      if(pcr_dep_count[alleleID][workId]==0)
	global_dep_count[alleleID] += 1;

      pcr_dep_count[alleleID][workId]++;

      double led1= (pcr_dep_count[alleleID][workId]-1)*pcr_dependency;
      double led2 =  global_dep_count[alleleID]*global_dependency;
      double led3 =nd->qs[j]*pow(10,led1+led2)+0.5;//+0.5000001;//this small eps addition is to make it work like soapsnp.

      int q_adjusted = int(led3);
      
      if(q_adjusted < 1)
	q_adjusted = 1;
      //      fprintf(stderr,"qadj=%d\n",q_adjusted);
      int coord;
      if(isupper(nd->seq[j]))
	coord = nd->posi[j];
      else
	coord = nd->isop[j];
      //      fprintf(stderr,"coord=%d\n",coord);
      for(int allele1 = 0;allele1 != 4;allele1++ ) {
	for(int allele2 = allele1; allele2 != 4; allele2++) {
	  //	  fprintf(stderr,"co1=%d co2=%d aqdj=%d coord=%d al1=%d al2= alID=%d\n",co(q_adjusted,coord,allele1,alleleID),co(q_adjusted,coord,allele2,alleleID),q_adjusted,coord,allele1,allele2,alleleID);
	  type_likely[allele1][allele2] += log10(0.5*p_matrix[co(q_adjusted,coord,allele1,alleleID)] +0.5*p_matrix[co(q_adjusted,coord,allele2,alleleID)]);
	  
	}
      }
    }
    if(mySum(global_dep_count)==-4)
      continue;

    //plugin likelihoods according, alloc for all samples if working with firstsample

    //soap snp order printf("#AA,CC,GG,TT,AC,AG,AT,CG,CT,GT\n")
    double *tmpLk = lk[s]+10*which_sample;
    tmpLk[0] = type_likely[0][0];//AA
    tmpLk[1] = type_likely[0][1];//AC
    tmpLk[2] = type_likely[0][3];//AG
    tmpLk[3] = type_likely[0][2];//AT
    tmpLk[4] = type_likely[1][1];//CC
    tmpLk[5] = type_likely[1][3];//CG
    tmpLk[6] = type_likely[1][2];//CT
    tmpLk[7] = type_likely[3][3];//GG
    tmpLk[8] = type_likely[2][3];//GT
    tmpLk[9] = type_likely[2][2];//TT

    for(int i=0;i<10;i++)//swap to normal log
      tmpLk[i] *= log(10.0);

    double mmax = tmpLk[0];//rescale to most likely

  }
  return 1;
}

void count_qual(size_t *count_matrix,const char*fname) {
  double *p_matrix = new double [256*256*4*4];
  for(int i=0;i<256*256*4*4;i++)
    p_matrix[i] = std::numeric_limits<double>::quiet_NaN();

  
  //  size_t coord;
  size_t o_base/*o_based base*/;
  size_t t_base/*theorecical(supposed) base*/;
  size_t sum[4];
  size_t  same_qual_count_by_type[4][4];// = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
  size_t  same_qual_count_by_t_base[4];// = {0,0,0,0};
  size_t  same_qual_count_total;
  size_t  same_qual_count_mismatch;
  
  const size_t sta_pow=10; // minimum number to say statistically powerful

  for(int q_char=0; q_char<255 ;q_char++)  {
    memset(same_qual_count_by_t_base,0,sizeof(size_t)*4);
    for(int i=0;i<4;i++)
      memset(same_qual_count_by_type[i],0,sizeof(size_t)*4);

    
    same_qual_count_total = 0;
    same_qual_count_mismatch = 0;
    for(int coord=0; coord < 255 ; coord++) {
      for(int t1=0;t1<4;t1++)
        for(int t2=0;t2<4;t2++){
          // If the sample is small, then we will not consider the effect of read cycle.
	  int funkyOffset = co(q_char,coord,t1,t2);
          same_qual_count_by_type[t1][t2] += count_matrix[funkyOffset];
          same_qual_count_by_t_base[t1] += count_matrix[funkyOffset];
          same_qual_count_total += count_matrix[funkyOffset];
          if(t1!= t2) 
	    // Mismatches
	    same_qual_count_mismatch += count_matrix[funkyOffset ];
        
        }
    }


    for(int coord=0; coord <255 ; coord++) {
      memset(sum, (size_t)0, sizeof(size_t)*4);
      // Count of all ref base at certain coord and quality
      for(int t1=0;t1<4;t1++)
        for(int t2=0;t2<4;t2++)
          sum[t1] += count_matrix[co(q_char,coord,t1,t2)]; // (type>>2)&3: the ref base

      for(t_base=0; t_base!=4; t_base++) {
        for(o_base=0; o_base!=4; o_base++) {
	  int fo = co( q_char,coord,t_base, o_base);
	  //	  fprintf(stderr,"fo=%d,qchar=%d coord=%d t_base=%lu obase_t=%lu\n",fo,q_char,coord,t_base,o_base);
	  if (count_matrix[fo] > sta_pow) {
	    // Statistically powerful
            p_matrix[fo] = ((double)count_matrix[fo]) /((double) sum[t_base]);
          }
          else if (same_qual_count_by_type[t_base][o_base] > sta_pow) {
	    // Smaller sample, given up effect from read cycle
            p_matrix [fo] =  ((double)same_qual_count_by_type[t_base][o_base]) /((double) same_qual_count_by_t_base[t_base]);
          }
          else if (same_qual_count_total > 0){
	    // Too small sample, given up effect of mismatch types
            if (o_base == t_base) {
	      p_matrix [fo] = ((double)(same_qual_count_total-same_qual_count_mismatch))/((double)same_qual_count_total);
            }
            else {
	      p_matrix [fo] = ((double)same_qual_count_mismatch)/same_qual_count_total;
            }
          }
          else {
	    //DRAGON FIX, deviation from std soapsnp
	    p_matrix[fo]  = std::numeric_limits<double>::quiet_NaN()
            ;
          }

          // For these cases like:
          // Ref: G o_base: G x10 Ax5. When calculate the probability of this allele to be A,
          // If there's no A in reference gives observation of G, then the probability will be zero,
          // And therefore exclude the possibility of this pos to have an A
          // These cases should be avoid when the dataset is large enough
          // If no base with certain quality is o_based, it also doesn't matter
	  
	  //tsk isnan check is a deviation from soapsnp, to avoid valgrind problems
          if(std::isnan(p_matrix[fo])|| (p_matrix [fo]==0) || (p_matrix [fo] ==1)) {
            if (o_base == t_base) {
              p_matrix [fo] = (1-pow(10, -((q_char)/10.0)));
	      if(p_matrix [fo]<0.25) {
                p_matrix [fo] = 0.25;
              }
            }
            else {
              p_matrix [fo] = (pow(10, -((q_char)/10.0))/3);
	      if(p_matrix [fo]>0.25) {
                p_matrix [fo] = 0.25;
              }
            }
          }
        }
      }

    }
  }
  if(1) {
    std::fstream fp;
    //    fprintf(stderr,"[%s] dumping file:%s\n",__FUNCTION__,fname);
    fp.open(fname,std::fstream::out);
    for(int qc = 0;qc<=255;qc++) {
      for(int i=0;i<255;i++) {
	fp<< qc << '\t' << i;
	for(int a1=0;a1<4;a1++)
	  for(int a2=0;a2<4;a2++)
	    fp << '\t' << p_matrix[co(qc,i,conv[a1],conv[a2])];
      fp << std::endl;
      }
    }
    fp.close();
  }
  delete [] p_matrix;
}


int actuallyrun =0;
void soap_likes::gen_counts(chunkyT *chk,size_t *count_matrix,int whichSample,char *refs) {
  //qs,seqpos,ref,nuc
  int printOrder[4] = {0,1,3,2};


  actuallyrun += chk->nSites;

  for(int s=0;s<chk->nSites;s++) {
    //    fprintf(stderr,"%c\n",refs[s]);
    if(refs[s]==4)
      continue;
    
    for(int l=0;chk->nd[s][whichSample]&&l<chk->nd[s][whichSample]->l;l++){
      char nuc = chk->nd[s][whichSample]->seq[l];
      if(nuc =='n' || nuc=='N')
        continue;
      /*
        increment:
        [qs][position in read][reference][nucleotide]
      */
      int funkyOffset =0;
      if(isupper(nuc))//qs not +33 anymore
	funkyOffset = co(chk->nd[s][whichSample]->qs[l], chk->nd[s][whichSample]->posi[l],refs[s],refToInt[chk->nd[s][whichSample]->seq[l]]);
      else
	funkyOffset = co(chk->nd[s][whichSample]->qs[l], chk->nd[s][whichSample]->isop[l],refs[s],refToInt[chk->nd[s][whichSample]->seq[l]]);
      
      count_matrix[funkyOffset]++;//increment the correct offset.
	
    }
  }
}






void dump_count_mat(const char *fname,size_t *count_mat){
  //  fprintf(stderr,"[%s] dumping file:%s\n",__FUNCTION__,fname);
  int printOrder[4] = {0,1,3,2};
  FILE *fpp =NULL;
  fpp = fopen(fname,"w");
    for(int i=0;i<255;i++) {
      for(int n=0;n<255;n++){
	fprintf(fpp,"%d\t%d\t",i, n);
	for(int a1=0;a1<4;a1++)
	  for(int a2=0;a2<4;a2++)
	    fprintf(fpp,"%lu\t",count_mat[co(i,n,printOrder[a1],printOrder[a2])]);//+33 to conform with soapsnp version, just debug..
	fprintf(fpp,"\n");
      }
    }
    if(fpp) fclose(fpp);
    
}


void soap_likes::run(chunkyT *chk,double **lk,char *refs,int trim) {
  assert(myMuts!=NULL);
  assert(chk!=NULL);

  if(doRecal==0){
    for(int i=0;i<chk->nSamples  ;i++)
      calc_gl(chk,p_matrix[i],i,lk,trim);
  }else{
    if(refs==NULL){
      fprintf(stderr,"\t-> Must supply reference for SOAPsnp calibration\n");
      exit(0);
    }
    for(int i=0;i<chk->nSamples;i++){
      pthread_mutex_lock(myMuts+i);//do persample lock maybe not nescearry due to atomic operations
      gen_counts(chk,count_mat[i],i,refs);
      pthread_mutex_unlock(myMuts+i);
    }
  }
}


void soap_likes::setCaliNames(int nInd){
  //  fprintf(stderr,"setcalinames\n");
  for(size_t i=0;i<nInd;i++){
    char *ret1 = new char[128];
    char *ret2 = new char[128];
    snprintf(ret1,128,"%s/%zu.qual%zu",tmpdir,i,i);
    snprintf(ret2,128,"%s/%zu.counts%zu",tmpdir,i,i);
    if(strlen(ret1)>100||strlen(ret1)>100){
      fprintf(stderr,"Might be problems with to long filenames\n");
      fprintf(stderr,"Please check/modify this function: %s in file:%s\n",__FUNCTION__,__FILE__);
    }
    //    fprintf(stderr,"ret1=%s ret2=%s\n",ret1,ret2);
    p_mat_names.push_back(ret1);
    count_mat_names.push_back(ret2);
  }
  
}

int soap_likes::p_mat_exists(){
  
  for(size_t i=0;i<p_mat_names.size();i++){
    //fprintf(stderr,"checking [%s]\n",p_mat_names[i]);
    if(aio::fexists(p_mat_names[i])==0)
      return 0;
  }
  return 1;
}

//simple filereader of a precalculated calibration matrix
double *get_p_matrix(const char*fname){
  //  fprintf(stderr,"[%s] reading file: %s\n",__FUNCTION__,fname);

  double *p_matrix = new double [256*256*4*4];
  for(size_t i=0;i<256*256*4*4;i++)
    p_matrix[i] = 0.0;//maybe set to nan?

  //now plug in values
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  char buf[LENSS];

  
  while(fgets(buf,LENSS,fp)){
    //size to to avoid typecasting when bitshitting... maybe not nescarry
    size_t qc = atoi(strtok(buf,"\t\n "));
    size_t rlen = atoi(strtok(NULL,"\t\n "));
    //    fprintf(stderr,"qc=%lu rlen=%lu shifted=%lu\n",qc,rlen,(qc<<12) | (rlen<<4) |16);
    for(int al2=0;al2<16;al2++)
      p_matrix[(qc<<12) | (rlen<<4) |al2   ] = atof(strtok(NULL,"\t \n"));

    //	p_matrix[i][j][al2/4][al2 %4] = atof(strtok(NULL,"\t \n"));

  }
  if(fp) fclose(fp);
  return p_matrix;
}


void soap_likes::init(int nInd_a,const char*wdir){
  tmpdir = strdup(wdir);
  nInd = nInd_a;
  mkdir(tmpdir,0777);

  setCaliNames(nInd);
  if(p_mat_exists()==0){
    //    fprintf(stderr,"Calibration matrix doesn't exist now allocating important stuff\n");
    doRecal =1;
    count_mat = new size_t*[nInd];
    for(int i=0;i<nInd;i++){
      count_mat[i] = new size_t[256*256*4*4];
      memset(count_mat[i],0,sizeof(size_t)*256*256*4*4);
    }
    
  }else{
    p_matrix = new double*[nInd];
    for(int i=0;i<nInd;i++)
      p_matrix[i] = get_p_matrix(p_mat_names[i]); 
    doRecal =0;
  }

  myMuts = new pthread_mutex_t[nInd];
  for(int i=0;i<nInd;i++)
    if(pthread_mutex_init(myMuts+i,NULL))
      fprintf(stderr,"problems initializing mutex\n");
  
  //  fprintf(stderr,"doRecal=%d\n",doRecal);
}

soap_likes::~soap_likes(){
  //  fprintf(stderr,"soaplikes destructor\n");
  
  if(doRecal){
    for(int i=0;i<nInd;i++){//dumpcount and qual
      dump_count_mat(count_mat_names[i],count_mat[i]);
      count_qual(count_mat[i],p_mat_names[i]);
    }
    fprintf(stderr,"Dumping persample recalibrations matrices in dir:%s\n",tmpdir);
  }
  if(p_matrix!=NULL){
    for(int i=0;i<nInd;i++)
      delete [] p_matrix[i];
    delete [] p_matrix;
  }
  for(size_t i=0;i<count_mat_names.size();i++){
    delete [] count_mat_names[i];
    delete [] p_mat_names[i];
  }
  free(tmpdir);
  delete [] myMuts;
}
