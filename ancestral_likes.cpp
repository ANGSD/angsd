#include <cassert>
#include <cmath>
#include <cstdlib>   /* atoi */
#include <cstring>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <libgen.h>//basename
#include <limits> //<- for setting nan
#include <pthread.h>
#include <sstream>
#include <stdlib.h>     
#include <sys/stat.h>//mkdir
#include <sys/types.h>//mkdir
#include <vector>

#include "bambi_interface.h"
#include "ancestral_likes.h"
#include "analysisFunction.h"


extern int refToInt[256];
const int PRIMES = 3;
const int NUCLEOTIDES = 4;
const int STRANDS = 2;

void anc_likes::checkfilehandle(std::ifstream &fh, const std::string &filename) {
  if (!fh.is_open()) {
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}
void anc_likes::checkfilehandle(std::ofstream &fh, const std::string &filename) {
  if (!fh.is_open()) {
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}

void anc_likes::load_error_mat(const int & nInd, const int & readpos, const int & quals ){
  // sample, readpos, prime(5p,3p,c), strand ,qual, ancbase, errorbase: earrorrates
  // sample, readpos, prime(5p,3p,c), strand ,qual, allele1, allele2, obs: log(genotype_likelihood)
  errorrates.resize(nInd);
  genotyperates.resize(nInd);
  for (int sample=0; sample<nInd; sample++){
    errorrates[sample].resize(readpos + 1); // + 1 as we have 'c' as well
    genotyperates[sample].resize(readpos + 1); // + 1 as we have 'c' as well

    for (int p = 0; p < (readpos + 1); p++) {
      errorrates[sample][p].resize(PRIMES);
      genotyperates[sample][p].resize(PRIMES);

      for (int pr = 0; pr < PRIMES; pr++) {
        errorrates[sample][p][pr].resize(STRANDS);
        genotyperates[sample][p][pr].resize(STRANDS);

        for (int strand = 0; strand < STRANDS; strand++) {
          errorrates[sample][p][pr][strand].resize(quals);
          genotyperates[sample][p][pr][strand].resize(quals);

          for (int qu = 0; qu < quals; qu++) {
            errorrates[sample][p][pr][strand][qu].resize(NUCLEOTIDES);
            genotyperates[sample][p][pr][strand][qu].resize(NUCLEOTIDES);

            for (int allele1 = 0; allele1 < NUCLEOTIDES; allele1++) {
              errorrates[sample][p][pr][strand][qu][allele1].resize(NUCLEOTIDES);
              genotyperates[sample][p][pr][strand][qu][allele1].resize(NUCLEOTIDES);

              for (int allele2 = allele1; allele2 < NUCLEOTIDES; allele2++) {
                genotyperates[sample][p][pr][strand][qu][allele1][allele2].resize(NUCLEOTIDES);
                
              }
            }
          }
        }
      }
    }
  }

  std::ifstream error_fh(error_filename.c_str(), std::ios::in);
  checkfilehandle(error_fh, error_filename);

  std::string line;
  int sample, pos, prime, strand, qual, ancbase, errorbase;
  double errorrate;
  while ( getline(error_fh, line) ){
      std::istringstream ss (line);
      ss >> sample >> pos >> prime >> strand >>
        qual >> ancbase >> errorbase >> errorrate;
      
      errorrates[sample][pos][prime]
            [strand][qual]
            [ancbase][errorbase] = errorrate;
  }
  error_fh.close();
  
  double error1, error2;
  for (int sample=0; sample<nInd; sample++){
    for (int p = 0; p < (readpos + 1); p++) {
      for (int pr = 0; pr < PRIMES; pr++) {
        for (int strand = 0; strand < STRANDS; strand++) {
          for (int qu = 0; qu < quals; qu++) {
            for (int allele1 = 0; allele1 < NUCLEOTIDES; allele1++) {
              for (int allele2 = allele1; allele2 < NUCLEOTIDES; allele2++) {
                for (int obsallele = 0; obsallele < NUCLEOTIDES; obsallele++) {
                  error1 = errorrates[sample][p][pr][strand][qu]
                                      [allele1][obsallele];
                  error2 = errorrates[sample][p][pr][strand][qu]
                                      [allele2][obsallele];

                  if(error1 > 1 || error1<0){
                    fprintf(stderr,
                            "\nERROR: Sample: %d; readposition: %d; prime: %d; strand: "
                            "%d; qualitybin: %d; allele: %d; obsallele: %d "
                            "\n\terror2: %f; Values <0 and >1 are not possible. Check the error estimates",
                            sample, p, pr, strand, qu, allele1,
                            obsallele, error1);
                    exit(0);
                  }
                  if(error2 > 1 || error2<0){
                    fprintf(stderr,
                            "\nERROR: Sample: %d; readposition: %d; prime: %d; strand: "
                            "%d; qualitybin: %d; allele: %d; obsallele: %d "
                            "\n\terror2: %f; Values <0 and >1 are not possible. Check the error estimates",
                            sample, p, pr, strand, qu, allele2,
                            obsallele, error2);
                    exit(0);
                  }
                  if(error1<1e-9){
                    // fprintf(stderr, "these are all the double zeros\n");
                    error1=1e-9;
                  }
                  if(error2<1e-9){
                    // fprintf(stderr, "these are all the double zeros\n");
                    error2=1e-9;
                  }
                  
                  genotyperates[sample][p][pr][strand][qu][allele1]
                    [allele2][obsallele] = std::log( (error1 * 0.5) + (error2 * 0.5) );
                }
              }
            }
          }
        }
      }
    }
  }
}


bool anc_likes::check_filename_exists(const std::string & filename){
  // https://stackoverflow.com/a/19841704
  std::ifstream infile(filename.c_str());
  return infile.good();
}

std::vector<int> anc_likes::parse_qual_bins(const std::string &row,
                                            const std::string &delim) {
  std::vector<int> results;
  int token;
  size_t last = 0, next = 0;
  while ((next = row.find_first_of(delim, last)) != std::string::npos) {

    // token = std::stoi(row.substr(last, next - last));
    token = std::atoi(row.substr(last, next - last).c_str());
    results.push_back(token);
    last = next + 1;
  }
  // token = std::stoi(row.substr(last)); // to the end
  token = std::atoi(row.substr(last).c_str()); // atoi
  results.push_back(token);
  return results;
}

void anc_likes::setreadpos(const int & val){
  readpos = val;
  readpos_0based = val - 1;
}

void anc_likes::setfilenames(){
  mat_filename = tmpdir+ "/" + "anc_likes.errors";
  error_filename = mat_filename+".est";
}

std::vector<int> anc_likes::get_qual_converter(const std::vector<int> & quals_vector){

  std::vector<int> qual_converter;
  size_t maxqual=100, qmax;

  qual_converter.reserve(maxqual);
  for(size_t i=0; i<quals_vector.size(); i++){
      if(i==0){
        qmax=quals_vector[i];
      } else {
        qmax = quals_vector[i]-quals_vector[i-1];
      }

    for(size_t q=0; q < qmax; q++){
      qual_converter.push_back(i);
    }
  }

  while (qual_converter.size()<maxqual){
    qual_converter.push_back(quals_vector.size());
  }
  
  return qual_converter;
}

void anc_likes::init_count_mat(const int & nInd, const int & readpos, const int & quals){
  // int readpos = 20;
  // readpos, prime(5p,3p,c), strand ,qual, outgroup, perfect, sample
  count_mat *count_mat_sample;  // a pointer to a single count mat
  count_mat_samples.resize(nInd);
  for (int ind = 0; ind < nInd; ind++) {
    count_mat_sample = &count_mat_samples[ind];
    count_mat_sample->resize(readpos + 1); // + 1 as we have 'c' as well
    for (int p = 0; p < (readpos + 1); p++) {
      count_mat_sample->at(p).resize(PRIMES);

      for (int pr = 0; pr < PRIMES; pr++) {
        count_mat_sample->at(p)[pr].resize(STRANDS);

        for (int strand = 0; strand < STRANDS; strand++) {
          count_mat_sample->at(p)[pr][strand].resize(quals);

          for (int qu = 0; qu < quals; qu++) {
            count_mat_sample->at(p)[pr][strand][qu].resize(NUCLEOTIDES);

            for (int out = 0; out < NUCLEOTIDES; out++) {
              count_mat_sample->at(p)[pr][strand][qu][out].resize(NUCLEOTIDES);

              for (int per = 0; per < NUCLEOTIDES; per++) {
                count_mat_sample->at(p)[pr][strand][qu][out][per].resize(
                    NUCLEOTIDES);
              }
            }
          }
        }
      }
    }
  }
}

void anc_likes::init(int nInd_a, const char *wdir){
  tmpdir = strdup(wdir);
  mkdir(tmpdir.c_str(), 0777);
  setfilenames(); 
  nInd = nInd_a;

  setreadpos(20);  // e.g. 0-19, should NOT be hardcoded
  quals = "20,30"; // e.g. 0:0-19,1:20-29,2:30-100; should NOT be hardcoded

  quals_vector = parse_qual_bins(quals, ",");
  if(quals_vector[0]==0 || quals_vector[quals_vector.size()-1] > 100){
    std::cerr << "First quality base has to be >0 and <=100" << '\n';
    exit(EXIT_FAILURE);
  }

  qual_converter = get_qual_converter(quals_vector);

  if(! check_filename_exists(error_filename) ){
    if ( check_filename_exists(mat_filename)){
      fprintf(stderr,
              "\t %s is already created. Run 'python "
              "${ANGSD}/misc/est_ancestral_errors.py  %s ${cores}'. It will "
              "generate %s, now rerun this ANGSD command. This "
              "ANGSD process should be killed!!!",
              mat_filename.c_str(), mat_filename.c_str(),
              error_filename.c_str());
      fprintf(stderr, "[%s] EXITING", __FILE__);
      exit(0);
    }
    
    doRecal =1;
    init_count_mat(nInd, readpos, quals_vector.size()+1);

    myMuts = new pthread_mutex_t[nInd];
    for(int i=0;i<nInd;i++){
      if(pthread_mutex_init(myMuts+i,NULL)){
        fprintf(stderr,"problems initializing mutex\n");
      }
    }
  } else {
    // the error filename has been identified
    load_error_mat(nInd, readpos, quals_vector.size()+1);
  }
}

void anc_likes::gen_counts(const chunkyT *chk, count_mat *count_mat_sample,
                           const int &whichsample, char *refs, char *ancs, int *keepSites,
                           const int &trim) {

  int qual, strand, readpos, termini, dist_5p, dist_3p;
  for(int s=0;s<chk->nSites;s++) {
    
    if(refs[s]==4 || ancs[s]==4 || keepSites[s]==0){
      continue;
    }
    // we only handle one sample at the time.
    tNode *nd = chk->nd[s][whichsample];  // this is tNode at sites s for a single individual i

    if(nd==NULL)
      continue;

    for(int j=0;j<nd->l;j++){  // nd->l is depth
      
      int allele = refToInt[nd->seq[j]];
      int qs = nd->qs[j];
      //filter qscore, mapQ,trimming, and always skip n/N
      if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
        continue;
      }
      
      qual = qual_converter[qs];  // 0:0-19, 1:20-29, 2:30-99
      strand = isupper(nd->seq[j])==0; // true (1) if negative strand else false (0)
      if( strand ){
        // negative strand
        dist_5p = nd->isop[j];
        dist_3p = nd->posi[j];
      } else {
        // positive strand
        dist_5p = nd->posi[j];
        dist_3p = nd->isop[j];
      }
      
      if (( dist_3p > readpos_0based ) && ( dist_5p > readpos_0based )){
          readpos = readpos_0based+1;
          termini = 2; // c
      } else if (dist_3p < dist_5p) {
        readpos = dist_3p;
        termini = 1;  // 3p
      } else {
        readpos = dist_5p;
        termini = 0; // 5p
      }
      count_mat_sample->at(readpos)[termini][strand][qual][ancs[s]][refs[s]][allele]++;
    }
  }
}

void anc_likes::run(chunkyT *chk, double **lk, char *refs, char *ancs, int *keepSites, int trim){
  assert(chk!=NULL);
  if(doRecal==1){
    assert(myMuts!=NULL);

    if(refs==NULL){
      fprintf(stderr,"\t-> Must supply -ref (reference) for ancestral matrix generation\n");
      exit(0);
    }
    
    if(ancs==NULL){
      fprintf(stderr,"\t-> Must supply -anc (ancestral) for ancestral matrix generation\n");
      exit(0);
    }
    for(int i=0;i<chk->nSamples;i++){
      pthread_mutex_lock(myMuts+i);//do persample lock maybe not nescearry due to atomic operations
      gen_counts(chk, &count_mat_samples[i], i, refs, ancs, keepSites, trim);
      pthread_mutex_unlock(myMuts+i);
    }
  } else
    geno_ancestral(chk,lk,trim);
}

anc_likes::~anc_likes() {
    if (doRecal) {
      std::ofstream outfile(mat_filename.c_str(), std::fstream::out);
      checkfilehandle(outfile, mat_filename);
      for (int ind = 0; ind < nInd; ind++) {
        for (int p = 0; p < (readpos + 1); p++) {
          for (int pr = 0; pr < PRIMES; pr++) {
            for (int strand = 0; strand < STRANDS; strand++) {
              for (int qu = 0; qu < (quals_vector.size()+1); qu++) {
                for (int out = 0; out < NUCLEOTIDES; out++) {
                  for (int per = 0; per < NUCLEOTIDES; per++) {
                    for (int sample = 0; sample < NUCLEOTIDES; sample++) {
                      if(count_mat_samples[ind][p][pr][strand][qu][out][per][sample]>0){
                        outfile << ind << " " << p << " " << pr << " " << strand
                                << " " << qu << " " << out << " " << per << " "
                                << sample << " "
                                << count_mat_samples[ind][p][pr][strand][qu][out]
                                                    [per][sample]
                                << '\n';
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      outfile.close();

      fprintf(stderr,
              "->\tDumping persample recalibrations matrices in "
              "dir:%s\n->\tNow run 'python "
              "${ANGSD}/misc/est_ancestral_errors.py  %s ${cores}' to estimate "
              "the error rates (%s). Finally, rerun angsd.\n",
              tmpdir.c_str(), mat_filename.c_str(), error_filename.c_str());
    } // do recall ends

    if (myMuts != NULL) {
      delete[] myMuts;
    }
    // in cpp you dont have to collect when quitting if you dont use new
    errorrates.clear(); 
}

void anc_likes::geno_ancestral(chunkyT *chk,double **lk,int trim){
  // std::vector< std::vector < std::vector < double > > >  * t;
  int qual, sample, strand, readpos, termini, dist_5p, dist_3p;
  for(int s=0;s<chk->nSites;s++){  // this is the number of sites in the chunk
    for(int i=0;i<chk->nSamples;i++){  // this is also number of individuals
      //allocate for all samples, for first sample
      if(i==0){
	lk[s] = new double[10*chk->nSamples];  // lk[s] is per site and  10*numberofsamples
	for(int ii=0;ii<10*chk->nSamples;ii++)
	  lk[s][ii] = -0.0;//set default values
      }
      
      tNode *nd = chk->nd[s][i];  // this is tNode at sites s for a single individual i
      if(nd==NULL)
	continue;
      //calc like persample
      double *likes1 = lk[s]+10*i;  // 
      for(int j=0;j<nd->l;j++){  // nd->l is depth
        
	int allele = refToInt[nd->seq[j]];
	int qs = nd->qs[j];
	//filter qscore, mapQ,trimming, and always skip n/N
	if(nd->posi[j]<trim||nd->isop[j]<trim||allele==4){
	  continue;
	}

        qual = qual_converter[qs];  // 0, 1, 2 
        strand = isupper(nd->seq[j])==0; // true (1) if negative strand else false (0)
        sample = i; // see major loop
        if( strand ){
          // negative strand
          dist_5p = nd->isop[j];
          dist_3p = nd->posi[j];
        } else {
          // positive strand
          dist_5p = nd->posi[j];
          dist_3p = nd->isop[j];
        }

        if (( dist_3p > readpos_0based ) && ( dist_5p > readpos_0based )){
          readpos = readpos_0based+1;
          termini = 2; // c
        } else if (dist_3p < dist_5p) {
          readpos = dist_3p;
          termini = 1;  // 3p
        } else {
          readpos = dist_5p;
          termini = 0; // 5p
        }
        int index_counter = 0;
        // t = &genotyperates[sample][readpos][termini][strand][qual];
        for (int geno1 = 0; geno1 < NUCLEOTIDES; geno1++) {
          for (int geno2 = geno1; geno2 < NUCLEOTIDES; geno2++) {
            likes1[index_counter] += genotyperates[sample][readpos][termini][strand][qual][geno1][geno2][allele];
            // likes1[index_counter] += t->at(geno1)[geno2][allele];
            index_counter++;
          }
        }
      }
    }
  }
}
