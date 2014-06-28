

class soap_likes{
  char* tmpdir;
  std::vector<char *> p_mat_names;
  std::vector<char *> count_mat_names;
  void setCaliNames(int nInd);
  int p_mat_exists();
  int nInd;
  double **p_matrix;//the calibration matrices p_matrix;; each matrix is flat,
  size_t **count_mat;//the counts matrix;; each matrix is flat
  pthread_mutex_t *myMuts;
  void gen_counts(chunkyT *chk,size_t *count_matrix,int whichSample,char *refs);
public:
  int doRecal;
  soap_likes() {  p_matrix = NULL; count_mat = NULL; doRecal = 0;tmpdir=NULL;myMuts=NULL;};
  void init(int nInd,const char *wdir);
  void run(chunkyT *chk,double **lk,char *refs,int trim);
  ~soap_likes();
};
