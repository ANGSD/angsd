/*
  small class to filter positions and assign major minor.
  This class can be improved very much.

  Thorfinn 7march 2013

  On the most basic level then this class will remove those sites with an effective sample size<minInd

  It can also be used for filtering away those sites not included in the -sites file.keep

 */

#include <cassert>
#include <sys/stat.h>
#include "shared.h"
#include "abc.h"
#include "analysisFunction.h"
#include "abcFilter.h"

abcFilter::~abcFilter(){
  if(fl!=NULL){
    if(fl->keeps)
      free(fl->keeps);
    if(fl->major)
      free(fl->major);
    if(fl->minor)
      free(fl->minor);
    if(fl->fp) fclose(fl->fp);
    bgzf_close(fl->bg);
    delete fl;
  }
  free(fname);
}


abcFilter::abcFilter(argStruct *arguments){
  isBed = 0; //<- indicates that we don't assume bed
  //below if shared for all analysis classes
  shouldRun[index] = 1;
  header = arguments->hd;
  revMap = arguments->revMap;
  fl = NULL;
  //his is used by this class
  keepsChr = NULL;
  curChr = -1;
  fp = NULL;
  minInd = 0;
  fname = NULL;
  doMajorMinor =0;
  
  if(arguments->argc==2){
    if(!strcasecmp(arguments->argv[1],"-sites")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  //get options and print them
  getOptions(arguments);
  printArg(arguments->argumentFile);

}


void abcFilter::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr pos major minor))\n",fname);
  fprintf(argFile,"\t-sites\t\t%s\t(File containing sites to keep (chr posStart posStop mapScore))\n",fname);
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if atleast minInd of samples has data\n",minInd);
  fprintf(argFile,"\t-isBed\t\t%d\t-sites file is a bed file\n",isBed);
  fprintf(argFile,"\t1) You can force major/minor by -doMajorMinor 3\n\tAnd make sure file contains 4 columns (chr tab pos tab major tab minor)\n");
  fprintf(argFile,"\t2) You can also supply a bed file containing regions to use last column is a mappability score.\n\tSet this to a value differet from zero will then assume the input file is a bedfile. \n");
}

void abcFilter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-sites",fname,arguments);
  isBed=angsd::getArg("-isBed",isBed,arguments);
  if(fname!=NULL)  
    fl = filt_init(fname,revMap,header,isBed);
  if(fl!=NULL)
    fprintf(stderr,"\t-> [%s] -sites is still beta, use at own risk...\n",__FILE__);


  //1=bim 2=keep
  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && fl!=NULL&& fl->hasMajMin!=1){
    fprintf(stderr,"\t-> Must supply -sites with a file containing major and minor if -doMajorMinor 3\n");
  }
  if(doMajorMinor!=3 && fl!=NULL&& fl->hasMajMin==1){
    fprintf(stderr,"\t-> Filter file contains major/minor information to use these in analysis supper \'-doMajorMinor 3\'\n");
  }
  
  minInd = angsd::getArg("-minInd",minInd,arguments);
}


void abcFilter::run(funkyPars *p){
  //  fprintf(stderr,"nsites=%d\n",p->numSites);
  p->keepSites=new int[p->numSites];
  
  for(int s=0;s<p->numSites;s++){
    p->keepSites[s]=p->nInd;
    //    p->results->freq->keepInd[s]=nInd;  
  }


  if(fl!=NULL && fl->hasMajMin==1){
    //    fprintf(stderr,"aloocating for major and minor\n");
    p->major = new char [p->numSites];
    p->minor = new char [p->numSites];
    for(int i=0;i<p->numSites;i++){
      p->major[i] = 4;
      p->minor[i] = 4;
    }
  }

  if(fl!=NULL) {
    for(int s=0;s<p->numSites;s++){

      if(fl->keeps[p->posi[s]]==0){
	//	fprintf(stderr,"Plugging inf vals std\n");
	p->keepSites[s] =0;
      }
      if(p->keepSites[s] && fl->hasMajMin==1){
	//fprintf(stderr,"Plugging inf vals std majorminor\n");
	p->major[s] = fl->major[p->posi[s]];
	p->minor[s] = fl->minor[p->posi[s]];
      }
    }
  }
  //how set the keepsites according the effective sample size persite
  //if(0!=minInd){
    if(p->chk!=NULL){
      //loop over sites;
      for(int s=0;s<p->numSites;s++){
	if(p->keepSites[s]==0)
	  continue;
	int nInfo =0;
	tNode *tn = p->chk->nd[s];
	//loop over samples;
	for(int i=0;i<p->nInd;i++){
	  if(tn[i].l!=0)
	    nInfo++;
	}
	p->keepSites[s] =nInfo;
	if(0!=minInd){
	  if(minInd>nInfo)
	    p->keepSites[s] =0;
	}
    
    }
  }
}
void abcFilter::print(funkyPars *p){
}

void abcFilter::clean(funkyPars *p){
  
}


void abcFilter::readSites(int refId) {
  //  fprintf(stderr,"[%s].%s():%d -> refId:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);
  if(fl==NULL)
    return;
  std::map<int,asdf_dats> ::iterator it = fl->offs.find(refId);
  if(it==fl->offs.end()){
    fprintf(stderr,"\n\t-> Potential problem: The filereading has reached a chromsome: \'%s\', which is not included in your \'-sites\' file.\n\t-> Please consider limiting your analysis to the chromsomes of interest \n",header->name[refId]);
    fprintf(stderr,"\t-> see \'http://www.popgen.dk/angsd/index.php/Sites\' for more information\n");
    fprintf(stderr,"\t-> Program will continue reading this chromosome... \n");
    //exit(0);
    if(header->l_ref[refId]> fl->curLen){
      fl->keeps=(char*) realloc(fl->keeps,header->l_ref[refId]);
      fl->curLen = header->l_ref[refId];
      if(fl->hasMajMin==1){
	fl->major=(char*) realloc(fl->major,header->l_ref[refId]);
	fl->minor=(char*) realloc(fl->minor,header->l_ref[refId]);
      }
    }
    memset(fl->keeps,0,fl->curLen);
    return;
  }
  bgzf_seek(fl->bg,it->second.offs,SEEK_SET);
  if(it->second.len>fl->curLen) 
    fl->keeps=(char*) realloc(fl->keeps,it->second.len);
  bgzf_read(fl->bg,fl->keeps,it->second.len);

  if(fl->hasMajMin==1){
    if(it->second.len>fl->curLen) {
      fl->major = (char*) realloc(fl->major,it->second.len);
      fl->minor = (char*) realloc(fl->minor,it->second.len);
    }
    bgzf_read(fl->bg,fl->major,it->second.len);
    bgzf_read(fl->bg,fl->minor,it->second.len);
  }
 
}
