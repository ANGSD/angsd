#include <ctype.h>
#include "abc.h"
#include "bambi_interface.h"
#include "shared.h"

class fsb : abc{

  void run(funkyPars *p);
  void print(funkyPars *p);
  void clean();

};

void fsb::clean(){

}

void fsb::print(funkyPars *p){
 
}

char *calcStuff(funkyPars *p){
  char *ret = new char [p->numSites];
  for(int i=0;i<p->nInd;i++){
    for(int s=0;s<p->numSites;s++){
      int tmp[2]={0,0};
      tNode *tsk =p->chk->nd[i][s];
      if(tsk==NULL)
	continue;
      for (int c=0;c<tsk->l;c++){
	if(isupper(tsk->seq[c]))
	  tmp[0]++;
	else
	  tmp[1]++;
	
      }
      if(tmp[0]!=0 && tmp[1]!=0)
	ret[s] = 0;

    }

      
}
  return 0;
}

void fsb::run(funkyPars *p){

  char *strIsBad = calcStuff(p);


  for(int i=0;i<p->numSites;p++){
    if(strIsBad[i])
      p->keepSites[i] =0;
    else
      p->keepSites[i] =1;
  }


}
