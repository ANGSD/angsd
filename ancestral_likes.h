#include <string>
#include <vector>
#include <pthread.h>



// samplenumber, readpos, prime(5p,3p,c), strand ,qual, ancbase, errobaser

typedef std::vector<  // sample
    std::vector<                     // readpos
        std::vector<                 // prime
            std::vector<             // strand
                std::vector<         // qual
                    std::vector<     // ancbase
                        std::vector< // errorbase
                            double > > > > > > > error_mat;

typedef std::vector<       // readpos
    std::vector<                     // prime
        std::vector<                 // strand
            std::vector<             // qual
                std::vector<         // outgroup
                    std::vector<     // perfect
                        std::vector< // sample
                            unsigned long long int > > > > > > > count_mat;

typedef std::vector <                // sample
  std::vector<                       // readpos
    std::vector<                     // prime
        std::vector<                 // strand
            std::vector<             // qual
                std::vector<         // allele1
                    std::vector<     // allele2
                        std::vector< // observed allele
                          double > > > > > > > >  genotypes_mat;


class anc_likes {
private:
  // set in init
  std::string tmpdir;
  int nInd;
  int readpos, readpos_0based; // 0-19
  error_mat errorrates;
  genotypes_mat genotyperates;
  // quals vector is based on "20,30" -> (20, 30)
  std::vector<int> quals_vector;
  // this is 100 long and based on 20,30 from quals_vector
  std::vector<int> qual_converter;

  std::vector<int> parse_qual_bins(const std::string &row,
                                   const std::string &delim);
  void init_error_mat(const int & readpos, const int qualbins);

  std::vector<int> get_qual_converter(const std::vector<int> & quals_vector);

  // check in tmp dir for a hardcoded file -> anc_likes.errors.est
  // if doesnt exists but anc_likes.errors exists, remind the user to
  // run the python script for estimating errors
  // if negative set doRecal to 1 and let the fun begin

  std::string mat_filename, error_filename;
  void setfilenames();
  bool check_filename_exists(const std::string & filename);
  std::vector<count_mat> count_mat_samples;
  pthread_mutex_t *myMuts;
  void init_count_mat(const int & nInd, const int & readpos, const int & quals);
  void load_error_mat(const int & nInd, const int & readpos, const int & quals);
  void checkfilehandle(std::ifstream &fh, const std::string &filename);
  void checkfilehandle(std::ofstream &fh, const std::string &filename);
  void gen_counts(const chunkyT *chk, count_mat *count_mat_sample,
                  const int &whichsample, char *refs, char *ancs, int *keepSites,
                  const int &trim);
  void geno_ancestral(chunkyT *chk,double **lk,int trim);
public:

  std::string quals;
  int doRecal;
  void setreadpos(const int & val);
  anc_likes()
      : tmpdir("angsd_tempdir"), nInd(), readpos(), readpos_0based(),
        errorrates(), genotyperates(), quals_vector(), qual_converter(),
        mat_filename(), error_filename(), count_mat_samples(), myMuts(),
        quals(), doRecal(0){};

  // find out what analyses to do
  void init(int nInd, const char *wdir);
  // run calls gen_counts if dorecall==1 else calc_gl(already implemented)
  void run(chunkyT *chk, double **lk, char *refs, char *ancs, int *keepSites,int trim);
  // anc_likes print count_mat if dorecall else do noting
  ~anc_likes();
};
