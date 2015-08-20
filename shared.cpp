/*
  remember to cleanup the threaded cleaning.
  it is not faster to thread the cleaning.
  

  idea for improvement.
  1. make blockwise threading of the indexed iteration? tsk I have no idea what this means

  
 */

#include <cstdio>
#include <vector>
#include <cstring>
#include <sys/stat.h>
#include <pthread.h>
#include <cassert>
#include <errno.h>
#include "bambi_interface.h"
#include "shared.h"
#include "abc.h" //<-new way, to add additional analysis to angsd
#include "analysisFunction.h" //<-included fancy parsing of argument

//small hack for the keepsites, this should be included allready in the bam parsing
#include "abcFilter.h"
#include "abcSmartCounts.h"
#include "abcWriteFasta.h"
#include "abcSaf.h"
#include "mUpPile.h"

pthread_attr_t attr;
int andersSux =0;
abc **allMethods = NULL;
char *shouldRun = NULL;
extern pthread_mutex_t mUpPile_mutex;

int maxThreads=1;//this value is given from the commandline -nThreads
pthread_t *thread = NULL;
funkyPars **fdata = NULL;
pthread_barrier_t   barrier; // barrier synchronization object
int currentChr=-1;
int chunkNumber =0;
int SIG_COND2 = 1;
// print nicely into files.
int isAtty =1;


const bam_hdr_t *header = NULL;

// howoften should we printout to screen default 100
int howOften = 100;


void collapse(funkyPars *p){
  fcb *f = p->for_callback;
  chunkyT *chk = mergeAllNodes_new(f->dn,f->nFiles);
  chk->regStart = f->regStart;
  chk->regStop = f->regStop;
  chk->refId = f->refId;

  p->refId = chk->refId;
  p->numSites=chk->nSites;
  p->nInd = chk->nSamples;  

  //now chk contains the merged data
  p->chk = chk;
  p->posi = new int[chk->nSites];
  for(int i=0;i<chk->nSites;i++)
    for(int j=0;j<chk->nSamples;j++)
      if(chk->nd[i][j]!=NULL){
	p->posi[i] = chk->nd[i][j]->refPos;
	break;
      }


  p->keepSites = new int [p->numSites];
  for(int i=0;i<p->numSites;i++)
    p->keepSites[i] = p->nInd;

  /*
    we are in running different threads now, so dont clean up nodes
    but modify keeplist
  */
  if((p->posi[0]<chk->regStart)||(p->posi[p->numSites-1]>chk->regStop)){
    for(int i=0;i<p->numSites;i++){
      if(p->posi[i]<chk->regStart)
	p->keepSites[i] = 0;
      if(p->posi[i]>=chk->regStop)
	p->keepSites[i] = 0;
    }
  }
    

}

int main_analysis(funkyPars *p) {
  assert(p);
  //  fprintf(stderr,"p->fhunknurmber:%d\n",p->chunkNumber);
  //first step is to make a chunk of data from the sample "uppiles"
  if(p->for_callback!=NULL) 
    collapse(p);
  
  //run all methods (ORDER is defined in general.cpp)
  if(p->numSites==0)
    return 0;
  for(int i=0;i<andersSux;i++){
    if(shouldRun[i])
      allMethods[i]->run(p);
  }
  return 1;
}

// this will run eternaly untill data[threadid] is NULL
void *BusyWork(void *t){
   int i;
   long tid;
   tid = (long)t;
   pthread_barrier_wait (&barrier);//we are locking all threads, untill we are certain we have data
   while(SIG_COND2){
     pthread_barrier_wait (&barrier);
     if(fdata[tid]!=NULL)
       main_analysis(fdata[tid]);
     pthread_barrier_wait (&barrier);
   }
   pthread_exit((void*) t);
}


void init(argStruct *arguments){
  if(!isatty(fileno(stderr))){
    isAtty = 0;
  }
  maxThreads=angsd::getArg("-nThreads",maxThreads,arguments);
  maxThreads=angsd::getArg("-P",maxThreads,arguments);
  howOften = angsd::getArg("-howOften",howOften,arguments);
  if(!isatty(fileno(arguments->argumentFile)))
    fprintf(arguments->argumentFile,"--------------------\n[%s:%s()]\n\t-nThreads\t%d\tNumber of threads to use\n\\n\t-howOften\t%d\tHow often should the program show progress\n",__FILE__,__FUNCTION__,maxThreads,howOften);

  //these are  analysis that might be performed
  allMethods = extra(andersSux,arguments->outfiles,arguments->inputtype,arguments);
  shouldRun = allMethods[0]->shouldRun;
  header = arguments->hd;
  pthread_attr_init(&attr);
  pthread_barrier_init(&barrier, NULL, maxThreads+1);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  thread = new pthread_t [maxThreads];
  fdata =(funkyPars**) calloc(maxThreads,sizeof(funkyPars*));
  for(size_t t=0; t<maxThreads; t++) {
    int rc = pthread_create(&thread[t], &attr, BusyWork, (void *)t); 
    if (rc) {
      fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }
  pthread_barrier_wait (&barrier);//we are locking all threads, untill we are certain we have data
}

void closethreads(){
  fprintf(stderr,"\n");//after last positions print add a neu line mutafuka
  pthread_barrier_wait (&barrier);
  //  fprintf(stderr,"Will close threads\n");fflush(stderr);
  SIG_COND2=0;
  pthread_barrier_wait (&barrier);
  void *status = NULL;
  for(int t=0; t<maxThreads; t++) {
    int rc = pthread_join(thread[t], &status);
    if (rc) {
      fprintf(stderr,"ERROR; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
    fprintf(stderr,"Main: completed join with thread %d having a status   of %ld\n",t,(long)status);
  }

}



void destroy(){
  fprintf(stderr,"\t-> Calling destroy\n");
  closethreads();

  

  fprintf(stderr,"\t-> Done waiting for threads\n");
  for(int i=0;i<andersSux;i++)
    delete allMethods[i];
  delete [] abc::shouldRun;
  delete [] allMethods;

}
void tnode_destroy(tNode*);
void cleanUptNodeArray(tNode **row,int nSamples){
  //  fprintf(stderr,"nodearray\n");
    for(int i=0;i< nSamples;i++) {
      if(row[i]==NULL)
	continue;
      if(row[i]->l2!=0){
	for(int j=0;j<row[i]->l2;j++)
	  tnode_destroy(row[i]->insert[j]);
	free(row[i]->insert);
      }
      
      tnode_destroy(row[i]);
    }
    free(row);
}



//only one instance at a time is running this function
void printFunky(funkyPars *p){
  if((p->chunkNumber%howOften)==0){
    if(isAtty){
      //	fprintf(stderr,"posi:%d\n",p->posi[0]);
      fprintf(stderr,"\r\t-> Printing at chr: %s pos:%d chunknumber %d ",header->target_name[p->refId],p->posi[0]+1,p->chunkNumber);
    }else
      fprintf(stderr,"\t-> Printing at chr: %s pos:%d chunknumber %d\n",header->target_name[p->refId],p->posi[0]+1,p->chunkNumber);
  }if(p->numSites!=0){
    for(int i=0;i<andersSux;i++)
      if(shouldRun[i])
	allMethods[i]->print(p);
  }
   
  deallocFunkyPars(p);
  
}





void interface(funkyPars *p){
  main_analysis(p);
  printFunky(p);
}


void changeChr(int refId){
  //fprintf(stderr,"[%s.%s():%d] refid:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);fflush(stderr);
  funkyPars p;p.chunkNumber=-1;
  int runner(funkyPars *as);
  runner(&p);
  ((abcFilter *)allMethods[0])->readSites(refId);
  ((abcWriteFasta *)allMethods[19])->changeChr(refId);//used when changing chr;
  ((abcSmartCounts *)allMethods[20])->changeChr(refId);//used when changing chr;
  ((abcSaf *)allMethods[11])->changeChr(refId);//used when changing chr;

}


void waiter(int refId){
  //  fprintf(stderr,"_%s_\n",__FUNCTION__);fflush(stderr);

  if(currentChr==-1||refId!=currentChr){
    currentChr=refId;
    changeChr(refId);
  }
  
}


int runner(funkyPars *as){
  assert(as!=NULL);
  static int batch=-1;//<- this is the inarray position id
  int i;

  //plugin data in array if we have data;
  if(as->chunkNumber!=-1){
    batch = as->chunkNumber %maxThreads;
    fdata[batch] = as;

  }
  //  fprintf(stderr,"[%s] chunknr:%d batch:%d\n",__FUNCTION__,as->chunkNumber,batch);fflush(stderr); 
  //case where we launch all analysis threads
  if(batch==maxThreads-1){
    //fprintf(stderr,"FULL Will launch all analysis:\n");fflush(stderr);
    pthread_barrier_wait (&barrier);//RUN ALL THREADS
    pthread_barrier_wait (&barrier);//MAKE SURE THEY ARE FINISHED
    //    fprintf(stderr,"Will print results\n");fflush(stderr);
    for(i=0;i<maxThreads;i++)
      if(fdata[i])
	printFunky(fdata[i]);
    batch=-1;
    memset(fdata,0,sizeof(funkyPars*)*maxThreads);
    //fprintf(stderr,"Done print results\n");fflush(stderr);
    //fprintf(stderr,"RESULTS-> %d) chunknr:%d\n",i,fdata[i]->chunkNumber);
    //fflush(stderr);
  }else if(as->chunkNumber==-1){
    //fprintf(stderr,"SUBSUB Will launch all analysis but only to bach:%d:\n",batch);fflush(stderr);
    for(i=batch+1;i<maxThreads;i++){
      //fprintf(stderr,"data[%d] is set to null\n",i);fflush(stderr);
      fdata[i] = NULL;
    }
    pthread_barrier_wait (&barrier);//RUN ALL THREADS
    pthread_barrier_wait (&barrier);//MAKE SURE THEY ARE FINISHED
    for(i=0;i<=batch;i++)
      if(fdata[i]!=NULL)
	printFunky(fdata[i]);
      //      fprintf(stderr,"RESULTS-> %d) chunknr:%d\n",i,fdata[i]->chunkNumber);
    //fflush(stderr);
    batch=-1;
    memset(fdata,0,sizeof(funkyPars*)*maxThreads);
  }
}




/*
  This is the function that determines whether or not to start threads
*/
void selector(funkyPars *p){
 
  if(p!=NULL){
    p->chunkNumber = chunkNumber++;
    
    if(maxThreads==1)
      interface(p);
    else
      runner(p);
  }else{
    fprintf(stderr,"selector is NULL\n");fflush(stderr);
    if(maxThreads>1){
      funkyPars pp;pp.chunkNumber=-1;
      runner(&pp);
    }
    fprintf(stderr,"selector is NULL2\n");fflush(stderr);
    pthread_mutex_unlock(&mUpPile_mutex);
  }
}


/*
  initialize all pointers to zero
*/

funkyPars *allocFunkyPars(){

  funkyPars *r = new funkyPars;
  r->numSites =0;
  r->extras = new void*[andersSux];//funky;

  r->counts = NULL;
  r->likes = NULL;
  r->post = NULL;

  r->major = NULL;
  r->minor = NULL;

  r->ref = NULL;
  r->anc= NULL;
  r->keepSites =NULL;
  r->chk = NULL;
  r->for_callback = NULL;
  r->posi = NULL;
  r->refId = -1;
  return r;
}





void deallocFunkyPars(funkyPars *p) {
  //clean up data related to each analysis class;
  if(p->numSites!=0){
    for(int i=0;i<andersSux;i++)
      if(shouldRun[i])
	allMethods[i]->clean(p);
  }
  
  //cleanup begin  

  delete[] p->keepSites;
  
  delete []  p->anc;
  delete [] p->ref;
  

  //cleanup extra stuff
  if(p->chk!=NULL){
    //    fprintf(stderr,"cleanuing up chunkyT\n");
    cleanUpChunkyT(p->chk);
    p->chk = NULL;
  }

  if(p->for_callback!=NULL){
    fcb *f = p->for_callback;
    for(int i=0;i<f->nFiles;i++)
      if(f->dn[i].l!=0)//this is nescesarry...
	free(f->dn[i].nds);
    delete [] f->dn;
    delete f;
    f=NULL;
  }
  delete [] p->extras;

  delete [] p->posi;//<-new way
  
  //below for beagle input
  if(p->major){
    delete [] p->major;p->major=NULL;
  }if(p->minor){
    delete [] p->minor;p->minor=NULL;
  }
  if(p->post){
    for(int i=0;i<p->numSites;i++)
      delete [] p->post[i];
    delete [] p->post;
    p->post=NULL;
  }
  if(p->likes){
    for(int i=0;i<p->numSites;i++)
      delete [] p->likes[i];
    delete [] p->likes;
    p->post=NULL;
  }
  delete p;

}

//plus one, plus 33
void printChunkyT(chunkyT *chk,double **liks,char *refs,FILE *fp){
  for(int s=0;s<chk->nSites;s++){
    fprintf(fp,"%d\t%d\t",chk->refId,chk->nd[s][0]->refPos+1);
    if(refs!=NULL)
      fprintf(fp,"%c\t",intToRef[refs[s]]);
    for(int n=0;n<chk->nSamples;n++){
      tNode *nd = chk->nd[s][n];
      if(nd==NULL)
	continue;
      fprintf(fp,"%d\t",nd->l);
      for(int i=0;i<nd->l;i++)
	fprintf(fp,"%c",nd->seq[i]);
      fprintf(fp,"\t");
      for(int i=0;i<nd->l;i++)
	fprintf(fp,"%c",nd->qs[i]+33);
      fprintf(fp,"\t");
      for(int i=0;i<10;i++)
	fprintf(fp,"%f ",liks[s][n*10+i]);
    }
    fprintf(fp,"\n");
  }
}
#ifdef __WITH_POOL__
extern int currentnodes;
#endif
 
