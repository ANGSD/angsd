#include "abc.h"

class abcCovar:public abc{
public:
  //none optional stuff
  FILE *outfile;
  int doCovar;
  
  void run(funkyPars *pars);
  void clean(funkyPars *pars);  
  void print(funkyPars *pars);  
  void addDefault(funkyPars *pars);
  void openfile(const char *outfiles);
  abcCovar(const char *outfiles,argStruct *arguments);
  void printArg(FILE *argfile);
  void getOptions(argStruct *arguments);
  ~abcCovar();
  //other stuff
  
};
