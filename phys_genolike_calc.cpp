#include "phys_genolike_calc.h"

// Corrections to q-score determined from Yana sample
//const float phys_genolike_calc::qscore_corr[61] = {1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
//						   0.623,1.082,0.930,0.979,3.552,4.249,5.041,1.837,2.817,2.272,2.460,1.854,1.263,
//						   1.841,1.466,1.994,0.972,7.822,1.150,2.499,1.468,1.681,1.962,2.306,2.372,3.564,.501,3.056,3.471,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
//						   1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000};





// Fitted model parameters, also not to be hardcoded
//const float phys_genolike_calc::parlist[18] = { 1.40519, 1.03243, 0.85508, 1.15778, 0.21355, -0.15465, 0.33590, 0.50438, 0.47351 , 
//						0.99191, 0.72090, 0.50917, 0.09018, 0.40386, 0.01166, -0.01656, 0.02668, 0.01919 };


//const float phys_genolike_calc::pararray[2][2][4][4] = {{{{ 1.0    , parlist[0] , parlist[1] , parlist[2] },{ parlist[3], 1.0     , parlist[4], parlist[5] },
//							  { parlist[6], parlist[7] , 1.0     , parlist[8] },{ parlist[9], parlist[10],parlist[11], 1.0     }},   // 4X4 Matrix
//							 {{ 1.0    , parlist[11], parlist[10], parlist[9] },{ parlist[8], 1.0     , parlist[7], parlist[6] },
//							  { parlist[5], parlist[4] , 1.0     , parlist[3] },{ parlist[2], parlist[1] , parlist[0], 1.0     }}},  // 4X4 Matrix
//							{{{ 1.0    , parlist[0] , parlist[1] , parlist[2] },{ parlist[3], 1.0     , parlist[4], parlist[5] },
//							  { parlist[6], parlist[7] , 1.0     , parlist[8] },{ parlist[9], parlist[10],parlist[11], 1.0     }},   // 4X4 Matrix
//							 {{ 1.0    , parlist[11], parlist[10], parlist[9] },{ parlist[8], 1.0     , parlist[7], parlist[6] },
//							  { parlist[5], parlist[4] , 1.0     , parlist[3] },{ parlist[2], parlist[1] , parlist[0], 1.0     }}}}; // 4x4 Matrix

vector<string> split_string( string instring, char split_by );
double str_to_double( string textstring );

// --------------------------------------------------------------------------------
// 
void phys_genolike_calc::read_coef( char *parspath ) {
    
  // Don't spam output to stderr unless specified
  debug     = false;

  // Use this to store the probabilities for the individual bases
  for( int ibase=0; ibase<4; ibase++ ){
    base_prob[ibase] = 0.0;
  }
  
  // Inititialise array containing probabilities for individual geno types
#if 0
  for( int igt=0; igt<10; igt++ ){
    geno_likes[igt] = 0.0;
  }
#endif
  // Inititialise 'matrix' containing probablilities that a base come
  // from a specific genotype
  
  // 0 : AA
  // 1 : AC
  // 2 : AG
  // 3 : AT
  //
  // 4 : CC
  // 5 : CG
  // 6 : CT
  // 
  // 7 : GG
  // 8 : GT
  // 
  // 9 : TT
  float **m_base_geno = new float*[4];
  for( int ibase=0; ibase<4; ibase++ ){
    m_base_geno[ibase] = new float[10];
    for( int igt=0; igt<10; igt++ ){
      m_base_geno[ibase][igt] = 0.0;
    }
  }
  
  // A from XY 
  m_base_geno[0][0] = 2.0/5;
  m_base_geno[0][1] = 1.0/5;
  m_base_geno[0][2] = 1.0/5;
  m_base_geno[0][3] = 1.0/5;
  
  // C from XY 
  m_base_geno[1][1] = 1.0/5;
  m_base_geno[1][4] = 2.0/5;
  m_base_geno[1][5] = 1.0/5;
  m_base_geno[1][6] = 1.0/5;
  
  // G from XY 
  m_base_geno[2][2] = 1.0/5;
  m_base_geno[2][5] = 1.0/5;
  m_base_geno[2][7] = 2.0/5;
  m_base_geno[2][8] = 1.0/5;
  
  // T from XY 
  m_base_geno[3][3] = 1.0/5;
  m_base_geno[3][6] = 1.0/5;
  m_base_geno[3][8] = 1.0/5;
  m_base_geno[3][9] = 2.0/5;

  // Print out 
  if( debug ){
    fprintf( stderr, " Prob matrix for base from geno-type \n" );
    fprintf( stderr, "         \tAA       AC       AG       AT       CC       CG       CT       GG       GT       TT \n" );
    for(int ibase=0; ibase<4; ibase++){
      fprintf( stderr, "base(%d) \t", ibase );

      for( int igt=0; igt<10; igt++ ){
	fprintf( stderr, "%6.4f   ", m_base_geno[ibase][igt] );
      }

      // fprintf( stderr, "\n", ibase );
    }
  }


  ifstream infile;
  string readline;
  fprintf(stderr,"parsparth:%s\n",parspath);
  infile.open( parspath );
  float *qscore_corr = new float[61];
  float *parlist = new float[18];
  while( !infile.eof() && infile.good() ){

    getline(infile, readline);
    fprintf( stderr, "Reading Parameter line %s \n ", readline.c_str() );

    vector<double> qscores;
    vector<double> params;

    //Start by splitting string by tabs
    vector<string> format_line = split_string( readline, ' ' );
    if( format_line.size() < 2 ) continue;

    //Remove empty strings
    for( unsigned int i=0; i<format_line.size(); i++ ){
      format_line[i].erase( remove(format_line[i].begin(), format_line[i].end(), ' '), format_line[i].end() );
    }


    if( format_line[0].find( "QscoreCorrections" ) != string::npos ){
      fprintf( stderr, "Reading Qscores \n" );
      int counter = 0;
      for( unsigned int i=1; i<format_line.size(); i++ ){
	if( format_line[i].empty() ) continue;

	qscore_corr[counter] = str_to_double( format_line[i] );

	fprintf( stderr, "\t %d %5.3f \n",i, qscore_corr[counter] );

	counter++;
      }
    }
 
    if( format_line[0].find( "ModelParameters" ) != string::npos ){
      fprintf( stderr, "Reading Model Parameters \n" );
      int counter = 0;
      for( unsigned int i=1; i<format_line.size(); i++ ){
	if( format_line[i].empty() ) continue;
	
	parlist[counter] = str_to_double( format_line[i] );

	fprintf( stderr, "\t %d %5.3f \n",i, parlist[counter] );

	counter++;
      }
    }

  }

  infile.close();

  // Initialise probability model
  // It is likely that once we stop hardcoding the model that 
  // the arrays will be initialised by the constructor instead?
  // init_p_base();
  m_base_geno_vec.push_back(m_base_geno);
  qscore_corr_vec.push_back(qscore_corr);
  parlist_vec.push_back(parlist);
} //End constructor

// --------------------------------------------------------------------------------
// Destructor
phys_genolike_calc::~phys_genolike_calc(){

  //  delete chk;
  //  delete nd;

} //End destructor

// // --------------------------------------------------------------------------------
// // Method to initialise pbase histogram
// void phys_genolike_calc::init_p_base(){
// 
// 
//   if( debug ){
//     fprintf( stderr, "Initialising models \n" );
//   }
// 
// // Needs to be not hardcoded at some point all this junk
// 
// 
// 
// 
// 
//   if( debug ){
//     fprintf( stderr, "Q-score corrections ( qscore, corrected p_val, p_val, correction factor \n" );
//     for( int iqs=0; iqs<42; iqs++ ){
//       float p_uncorr = pow( 10, -0.1 * iqs );
//       float p_transversion = p_uncorr * qscore_corr[iqs];
// 
//       fprintf( stderr, "%2d, %7.5f, %7.5f, %5.3f \n", iqs, p_transversion, p_uncorr, qscore_corr[iqs] );
//     }
//   }
// 
// } // End init_p_base

// --------------------------------------------------------------------------------
// Update probabilities (4) for each base, for observed base at depth ''depth''
void phys_genolike_calc::update_pbase( int depth, float *qscore_corr, float *parlist, tNode *nd ){

  //junk MOVE TO CONSTRUCTOR!!!!
  float pararray[2][2][4][4] = {{{{ 1.0    , parlist[0] , parlist[1] , parlist[2] },{ parlist[3], 1.0     , parlist[4], parlist[5] },
				  { parlist[6], parlist[7] , 1.0     , parlist[8] },{ parlist[9], parlist[10],parlist[11], 1.0     }},   // 4X4 Matrix
				 {{ 1.0    , parlist[11], parlist[10], parlist[9] },{ parlist[8], 1.0     , parlist[7], parlist[6] },
				  { parlist[5], parlist[4] , 1.0     , parlist[3] },{ parlist[2], parlist[1] , parlist[0], 1.0     }}},  // 4X4 Matrix
				{{{ 1.0    , parlist[0] , parlist[1] , parlist[2] },{ parlist[3], 1.0     , parlist[4], parlist[5] },
				  { parlist[6], parlist[7] , 1.0     , parlist[8] },{ parlist[9], parlist[10],parlist[11], 1.0     }},   // 4X4 Matrix
				 {{ 1.0    , parlist[11], parlist[10], parlist[9] },{ parlist[8], 1.0     , parlist[7], parlist[6] },
				  { parlist[5], parlist[4] , 1.0     , parlist[3] },{ parlist[2], parlist[1] , parlist[0], 1.0     }}}}; // 4x4 Matrix


  // Retrieve information from node
  int qscore = nd->qs[depth];
  int posi   = nd->posi[depth];
  int isop   = nd->isop[depth];

  int strand = isupper( nd->seq[depth] ) == 0;

  int allele = refToInt[nd->seq[depth]];

  bool use_posi = (posi <= isop);
  int dir = use_posi ? 0 : 1;
  int abs_position = use_posi ? posi : isop; // Absolute position from nearest end

  // Start by correcting qscore and convert to a probability
  //    float corr_qscore = qscore * qscore_corr[qscore];
  float p_uncorr = pow( 10, -0.1 * qscore ) / 3.0 ;
  
  // Restart individual base probabilities
  base_prob[0] = 1.0; base_prob[1] = 1.0; base_prob[2] = 1.0; base_prob[3] = 1.0; 


  
  //Loop over different bases
  for( int ia=0; ia<4; ia++ ){
    
    // For now ignore the probability of conserving
    if( ia == allele ) continue;
    
    float p_transversion = p_uncorr * qscore_corr[qscore] * pararray[dir][strand][ia][allele];
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

    /*
    if( 0.0 > base_prob[ia] || 1.0 < base_prob[ia] ){
      fprintf( stderr, "WARNING !!! P Not defined : %10.7f = %10.7f + %10.7f - %10.7f \n", base_prob[ia], p_transversion, p_transition, p_transition * p_transversion );

      fprintf( stderr, "Calc p(base) (%2d) / allele (%d -> %d) : q-score( %d ), posi (%2d), isop(%2d), dir(%1d)", 
	       depth, ia, allele, qscore, posi, isop, dir );

      fprintf( stderr, " :  P(A)=%5.3f, P(C)=%5.3f, P(G)=%5.3f, P(T)=%5.3f \n", 
	       base_prob[0], base_prob[1], base_prob[2], base_prob[3] );
    }
    */

  }// done looping over non observed allele

  // probability to conserve base defined as 1 - prob( transition U transversion )
  for( int ia=0; ia<4; ia++ ){

    if( ia == allele ) continue;
    base_prob[allele] -= base_prob[ia];

  }

  
  if( debug ){
    fprintf( stderr, "Calc p(base) (%2d) / allele (%d) : q-score( %d ), posi (%2d), isop(%2d)", 
	     depth, allele, qscore, posi, isop );

    fprintf( stderr, " :  P(A)=%5.3f, P(C)=%5.3f, P(G)=%5.3f, P(T)=%5.3f \n", 
	     base_prob[0], base_prob[1], base_prob[2], base_prob[3] );
  }

}// end update p-base

// --------------------------------------------------------------------------------
// Function to return probabilities for each of the 10 individual geno-types
// Takes array with 10 floating points numbers
void phys_genolike_calc::get_genolikes( int site, int sample, chunkyT *chk, double *return_likes ){
  float **m_base_geno = m_base_geno_vec[sample];
  // Get reference to current tNode
  tNode *nd = chk->nd[site][sample];
  if( !nd ) return;

  //Reset geno_likes values TSK THESE ARE RESET AT ALLOCATION
#if 0
  for( int igt=0; igt<10; igt++ ){
    geno_likes[igt] = 0.0;
  }
#endif
  // Loop over node depth
  if( debug ) fprintf( stderr, "GenoTs   \tAA       AC       AG       AT       CC       CG       CT       GG       GT       TT \n" );

  for( int idepth=0; idepth<nd->l; idepth++ ){
    if( refToInt[nd->seq[idepth]] >= 4 ) continue;

    update_pbase( idepth,qscore_corr_vec[sample],parlist_vec[sample], nd );

    if( debug ) fprintf( stderr, "Dep (%2d) \t", idepth );
    for( int igt=0; igt<10; igt++ ){

      float gt_prob_temp = 0.0;

      for( int ibase=0; ibase<4; ibase++ ){
	gt_prob_temp += base_prob[ibase] * m_base_geno[ibase][igt];
      }//base loop

      if( gt_prob_temp < 0.0000001 ) continue;
      return_likes[igt] += log( gt_prob_temp );

    } //geno-type loop
    if( debug ) fprintf( stderr, "\n");

  } // Depth loop
#if 0
  for( int igt=0; igt<10; igt++ ){
    return_likes[igt] = geno_likes[igt];
  }
#endif
  
}// get_genoprobs

// --------------------------------------------------------------------------------
// Put base probabilities into the results_str
void phys_genolike_calc::get_base_prob_str( char *results_str ){

  int lbuffer=0;
  for( int ibase=0; ibase<4; ibase++ ){
    lbuffer += sprintf( results_str+lbuffer, "%7.5f   ", base_prob[ibase] );
  }

} // end get_base_prob_str
#if 0
// --------------------------------------------------------------------------------
// Put genotype probabilities into the results_str
void phys_genolike_calc::get_genolikes_str( char *results_str ){

  int lbuffer=0;
  for( int igt=0; igt<10; igt++ ){
    lbuffer += sprintf( results_str+lbuffer, "%6.4f   ", geno_likes[igt] );
  }


} // end get_geno_likes_str
#endif
// --------------------------------------------------------------------------------
// Use this to change to debug mode
void phys_genolike_calc::set_debug( bool db ){
  debug = db;
}


//---------------------------------------------------------------------------
// format input reduce_string to match RooAbsData and RooSimultaneous
vector<string> split_string( string instring, char split_by ){
  istringstream ss( instring );
  string token;

  vector<string> return_vec;
  if( instring.size() == 0 ) return return_vec;

  while( getline( ss, token, split_by ) ) {

    return_vec.push_back( token );

  } // 
  return return_vec;

}//end format_cut

//---------------------------------------------------------------------------
// Convert string to double
double str_to_double( string textstring ){
  stringstream ss; //create a stringstream
  double value;
  ss << textstring;    //add number to the stream
  ss >> value;
  return value;
}// end str to double
