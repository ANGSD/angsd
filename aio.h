#include <vector>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/hts.h>
#include <zlib.h>


/*
 * Macro:[AT]
 * Injects the file and line info as string
 */
#define STRINGIFY(x) #x
#define ASSTR(x) STRINGIFY(x)
#define AT __FILE__ ":" ASSTR(__LINE__)

/*
 * Macro:[ASSERT]
 * shortcut to evaluate an expression, works the same way as the C-macro assert
 */
#define ASSERT(expr) if (!(expr)) {fprintf(stderr,"\n\n*******\n[ERROR](%s:%d) %s\n*******\n",__FILE__,__LINE__,#expr);exit(1);}


//angsd io
namespace aio{
  size_t fsize(const char* fname);
  int fexists(const char* str);//{///@param str Filename given as a string.
  FILE *openFile(const char* a,const char* b);
  FILE *getFILE(const char*fname,const char* mode);
  BGZF *openFileBG(const char* a,const char* b);
  htsFile *openFileHts(const char * a, const char*b);
  htsFile *openFileHtsBcf(const char * a, const char*b);
  int isNewer(const char *newer,const char *older);
  ssize_t bgzf_write(BGZF *fp, const void *data, size_t length);
  int tgets(gzFile gz,char**buf,int *l);



/* [doAssert] evaluate a statement, works the same way as C macro assert
 *
 * @param exp_eval		expression to evaluate, exit if evaluates to false
 * @param exit_code 	code to exit with if expression evaluates to false
 * @error_location		specifies where the error is encountered at
 * 						macro AT (defined in aio.h) can be used to print
 * 						(which_file.cpp:which_line)
 * 
 * @exit_text			exit text to print before exit, use as in fprintf
 * 						followed by optional arguments
 *
 * @use
 * aio::doAssert(thisStatementShallBeTrue,1,AT,"Fancy error msg with a val %d",my_var);
 */

  void doAssert(int exp_eval, int exit_code, const char* error_location, const char* exit_text,...);
  void doAssert(int exp_eval, const char* error_location, const char* exit_text,...);
  void doAssert(int exp_eval, const char* error_location);
  void doAssert(int exp_eval);

}
