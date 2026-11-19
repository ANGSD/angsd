/*
  Genotype likelihoods for -GL 5 (phys). DEPRECATED: will be removed in
  a future version. See phys_genolike_calc.h for the model and the
  format of the coefficient file.
*/
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "phys_genolike_calc.h"

using namespace std;

extern int refToInt[256];

// Probability that a read base is drawn from a genotype, [base][genotype]
// with genotypes AA,AC,AG,AT,CC,CG,CT,GG,GT,TT. Every genotype sums to
// 2/5, which is a constant factor and does not change the likelihood
// ratios.
static const float m_base_geno[4][10] = {
  {2.0/5,1.0/5,1.0/5,1.0/5,0    ,0    ,0    ,0    ,0    ,0    },
  {0    ,1.0/5,0    ,0    ,2.0/5,1.0/5,1.0/5,0    ,0    ,0    },
  {0    ,0    ,1.0/5,0    ,0    ,1.0/5,0    ,2.0/5,1.0/5,0    },
  {0    ,0    ,0    ,1.0/5,0    ,0    ,1.0/5,0    ,1.0/5,2.0/5}
};

// A genotype probability below this is set to this value, so a read that
// does not fit a genotype lowers its likelihood instead of being skipped
static const float MIN_GT_PROB = 0.0000001;

void phys_genolike_calc::read_coef(const char *parspath){
  ifstream infile(parspath);
  if(!infile.is_open()){
    fprintf(stderr,"\t-> Problem opening coefficient file: %s\n",parspath);
    exit(1);
  }

  physPars p;
  int nQs = 0;
  int nPars = 0;
  string readline;
  while(getline(infile,readline)){
    istringstream ss(readline);
    string name;
    if(!(ss>>name))
      continue;
    vector<double> vals;
    double v;
    while(ss>>v)
      vals.push_back(v);
    if(!ss.eof()){
      fprintf(stderr,"\t-> Problem parsing line starting with '%s' in coefficient file: %s\n",name.c_str(),parspath);
      exit(1);
    }

    if(name.find("QscoreCorrections")!=string::npos){
      if(vals.size()<1||vals.size()>61){
	fprintf(stderr,"\t-> QscoreCorrections must have between 1 and 61 values, %s has %zu\n",parspath,vals.size());
	exit(1);
      }
      nQs = vals.size();
      for(int i=0;i<61;i++)//quality scores above the last value use the last value
	p.qscore_corr[i] = vals[i<nQs?i:nQs-1];
    }else if(name.find("ModelParameters")!=string::npos){
      if(vals.size()!=18){
	fprintf(stderr,"\t-> ModelParameters must have 18 values, %s has %zu\n",parspath,vals.size());
	exit(1);
      }
      nPars = vals.size();
      for(int i=0;i<18;i++)
	p.parlist[i] = vals[i];
    }
  }
  if(nQs==0||nPars==0){
    fprintf(stderr,"\t-> Coefficient file %s must have both a QscoreCorrections and a ModelParameters line\n",parspath);
    exit(1);
  }
  fprintf(stderr,"\t-> Read %d quality score corrections and %d model parameters from %s\n",nQs,nPars,parspath);

  // forward strand, reverse strand is the complemented matrix
  float *pl = p.parlist;
  float fw[4][4] = {{ 1.0  , pl[0] , pl[1] , pl[2] },
		    { pl[3], 1.0   , pl[4] , pl[5] },
		    { pl[6], pl[7] , 1.0   , pl[8] },
		    { pl[9], pl[10], pl[11], 1.0   }};
  for(int a=0;a<4;a++)
    for(int b=0;b<4;b++){
      p.pararray[0][a][b] = fw[a][b];
      p.pararray[1][a][b] = fw[3-a][3-b];
    }
  pars.push_back(p);
}

// Probabilities that the true base is A,C,G,T given the read base at depth
void phys_genolike_calc::update_pbase(int depth,const physPars &p,tNode *nd,float *base_prob){
  const float *parlist = p.parlist;
  int qscore = nd->qs[depth];
  if(qscore>60)
    qscore = 60;
  int posi   = nd->posi[depth];
  int isop   = nd->isop[depth];

  int strand = isupper( nd->seq[depth] ) == 0;

  int allele = refToInt[nd->seq[depth]];

  bool use_posi = (posi <= isop);
  int dir = use_posi ? 0 : 1;
  int abs_position = use_posi ? posi : isop; // Absolute position from nearest end

  float p_uncorr = pow( 10, -0.1 * qscore ) / 3.0 ;

  base_prob[0] = 1.0; base_prob[1] = 1.0; base_prob[2] = 1.0; base_prob[3] = 1.0;

  for( int ia=0; ia<4; ia++ ){
    if( ia == allele ) continue;

    float p_transversion = p_uncorr * p.qscore_corr[qscore] * p.pararray[strand][ia][allele];
    float p_transition = 0.0;

    if ((dir == 0 && ia == 1 && allele == 3) || (dir == 1 && ia == 2 && allele == 0)) {
      p_transition = parlist[12]*exp(-parlist[13]*abs_position) + parlist[14];
    }
    else if ((dir == 1 && ia == 1 && allele == 3) || (dir == 0 && ia == 2 && allele == 0)) {
      p_transition = parlist[15]*exp(-parlist[16]*abs_position) + parlist[17];
    }

    base_prob[ia] = p_transversion + p_transition - p_transversion * p_transition;

    if( base_prob[ia] < 0.0 ) base_prob[ia] = p_transition;
    if( base_prob[ia] > 1.0 ) base_prob[ia] = 0.999999;
  }

  // probability to conserve base defined as 1 - prob( transition U transversion )
  for( int ia=0; ia<4; ia++ ){
    if( ia == allele ) continue;
    base_prob[allele] -= base_prob[ia];
  }
}

// Adds the log likelihoods of the 10 genotypes to return_likes
void phys_genolike_calc::get_genolikes(int site,int sample,chunkyT *chk,int trim,double *return_likes){
  tNode *nd = chk->nd[site][sample];
  if( !nd ) return;

  float base_prob[4];
  for( int idepth=0; idepth<nd->l; idepth++ ){
    if( refToInt[nd->seq[idepth]] >= 4 ) continue;
    if( nd->posi[idepth]<trim || nd->isop[idepth]<trim ) continue;

    update_pbase( idepth, pars[sample], nd, base_prob );

    for( int igt=0; igt<10; igt++ ){
      float gt_prob_temp = 0.0;
      for( int ibase=0; ibase<4; ibase++ )
	gt_prob_temp += base_prob[ibase] * m_base_geno[ibase][igt];

      if( gt_prob_temp < MIN_GT_PROB ) gt_prob_temp = MIN_GT_PROB;
      return_likes[igt] += log( gt_prob_temp );
    }
  }
}
