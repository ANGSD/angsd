// Dear emacs, this is -*- c++ -*-
#ifndef PHYS_GENOLIKE_CALC_H
#define PHYS_GENOLIKE_CALC_H

//Angsd includes
#include "bambi_interface.h"

// C++ includes
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <math.h>

using namespace std;

extern int refToInt[256];
class phys_genolike_calc {

 private : 
  
  // Information stored in chunkyT
  
  

  bool debug;

  // Use thie to store probabilties for the resulting genotypes
  // implemented with double precision in angsd
  // double geno_likes[10];
  
  // Probability matrix that a base is from a given geno-type
  //float m_base_geno[4][10];
  std::vector<float**> m_base_geno_vec;

  // Contains q score corrections
  //float qscore_corr[61];
  std::vector<float*> qscore_corr_vec;  
  // model parameters
  //float parlist[18];
  std::vector<float*> parlist_vec;
  //  float pararray[2][2][4][4];

  void init_p_base();

 public : 

  // Use this to store the probabilities for the individual bases
  float base_prob[4];
 
  // Default constructor
  //phys_genolike_calc( char *parspath );

  // Default destructor
  ~phys_genolike_calc();

  // 
  void update_pbase( int depth,float *qscore_cor,float *parlist, tNode *nd );

  void get_genolikes( int site, int sample, chunkyT *chk, double *return_likes );

  // Helper functions to get results
  void get_base_prob_str( char *results_str );
  void get_genolikes_str( char *results_str );

  void set_debug( bool db );
  void read_coef(char *);

};

#endif
