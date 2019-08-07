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
#include <htslib/khash_str2int.h>
#define THRESHOLD_FOR_NICEOUTPUT 0//200 //we normally print whenever we change a chr, but don't bother if we have more the 200, looks ugly for contig data

//small hack for the keepsites, this should be included allready in the bam parsing
#include "abcFilter.h"
#include "abcSmartCounts.h"
#include "abcWriteFasta.h"
#include "abcGL.h"
#include "abcSaf.h"
//class to keep track of chunk order when dumping results
#include "printRes.h"
#include "mUpPile.h"
#include "pooled_alloc.h"
#include "abcPSMC.h"
#include "cigstat.h"
extern tpool_alloc_t *tnodes;

pthread_attr_t attr;
int andersSux =0;
abc **allMethods = NULL;
char *shouldRun = NULL;
extern pthread_mutex_t mUpPile_mutex;

int maxThreads=1;//this value is given from the commandline -nThreads
pthread_t thread1;
pthread_cond_t cvMaxThread;
pthread_mutex_t counterMut;

int currentChr=-1;
int curRunning =0;
int chunkNumber =1;
printRes printer;

// print nicely into files.
int isAtty =1;


const bam_hdr_t *header = NULL;

//Max number of "finished" threads waiting to be printed
// -1 indicates no limit
int nQueueSize = -1;
// queueStop =1 shouldnt stop new threads
// queueStop =0 should spawn newthreads;
int queueStop =0;

// howoften should we printout to screen default 100
int howOften = 100;


void init(argStruct *arguments){
  if(!isatty(fileno(stderr))){
    isAtty = 0;
  }
  maxThreads=angsd::getArg("-nThreads",maxThreads,arguments);
  maxThreads=angsd::getArg("-P",maxThreads,arguments);
  nQueueSize = angsd::getArg("-nQueueSize",nQueueSize,arguments);
  howOften = angsd::getArg("-howOften",howOften,arguments);
  if(!isatty(fileno(arguments->argumentFile)))
    fprintf(arguments->argumentFile,"--------------------\n[%s:%s()]\n\t-nThreads\t%d\tNumber of threads to use\n\t-nQueueSize\t%d\tMaximum number of queud elements\n\t-howOften\t%d\tHow often should the program show progress\n",__FILE__,__FUNCTION__,maxThreads,nQueueSize,howOften);
  //  fprintf(stderr,"nQueue=%d howOften=%d\n",nQueueSize,howOften);

  //these are  analysis that might be performed
  allMethods = extra(andersSux,arguments->outfiles,arguments->inputtype,arguments);
  shouldRun = allMethods[0]->shouldRun;
  pthread_mutex_init(&counterMut, NULL);
  pthread_cond_init (&cvMaxThread, NULL);


  header = arguments->hd;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  tnodes = tpool_create(sizeof(tNode));
}


void destroy_shared(){
  //  fprintf(stderr,"\t-> Calling destroy\n");
  while(1){
    pthread_mutex_lock(&counterMut);
    if(curRunning==0&&printer.contains()==0)
      break;
    fprintf(stderr,"\t-> Waiting for nThreads:%d\n",curRunning);
    pthread_mutex_unlock(&counterMut);
    sleep(1);
  }
  pthread_mutex_unlock(&counterMut);
  fprintf(stderr,"\t-> Done waiting for threads\n");
  for(int i=0;i<andersSux;i++)
    delete allMethods[i];
  delete [] abc::shouldRun;
  delete [] allMethods;
  void destroy_tnode_pool();
  destroy_tnode_pool();
  extern int cigstat;
  if(cigstat)
    cigstat_close();
  extern void *rghash;
  if(rghash)
    khash_str2int_destroy_free(rghash);
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



void interface(funkyPars *p){

  if(p->killSig==0)
    main_analysis(p);//dont do analysis on the killsig chunk

  printer.prop(p);
}

void *slave(void *ptr){
  funkyPars *p = (funkyPars *) ptr;
  interface(p);

  pthread_mutex_lock( &counterMut );
  curRunning--;

  
  if(nQueueSize==-1){//no limit on queuesize, always signal
    pthread_cond_signal(&cvMaxThread);
  }else{
    // this is for limiting the number of elements in queue. 
    // Only a problem with large >1000 samples.    
    if(printer.contains()<(unsigned int)nQueueSize){
      pthread_cond_signal(&cvMaxThread);
      queueStop =0; //can
    }else
      queueStop =1; //should not
  }
  //  fprintf(stderr,"queueStop =%d contains=%lu\n",queueStop,printer.contains());
  pthread_mutex_unlock( &counterMut );
  pthread_exit(NULL);
}


void *slave2(void *ptr){
  funkyPars *p = (funkyPars *) ptr;
  if(p->killSig==0)
    main_analysis(p);
  printFunky(p);
  pthread_mutex_unlock( &counterMut );
  pthread_exit(NULL);
}

void master2(funkyPars *p){
  //  fprintf(stderr,"[%s] Number of threads running=%d\n",__FUNCTION__,curRunning);
  pthread_mutex_lock( &counterMut );
  if(pthread_create( &thread1, NULL, slave2, (void*) p)){
    fprintf(stderr,"[%s] Problem spawning thread\n%s at chunknumber:%d\n",__FUNCTION__,strerror(errno),p->chunkNumber);
    exit(0);
  }
  
  pthread_detach(thread1);
}



void master(funkyPars *p){
  //  fprintf(stderr,"[%s] Number of threads running=%d inqueue=%lu\n",__FUNCTION__,curRunning,printer.contains());
  pthread_mutex_lock( &counterMut );

  curRunning++;
  if(curRunning==maxThreads||queueStop==1){
    //fprintf(stderr,"We are running max number of threads will wait for finishing threads: %d\n",curRunning);
    pthread_cond_wait(&cvMaxThread, &counterMut);
    //fprintf(stderr,"Done waiting for finishing threads:%d\n",curRunning);
  }
  pthread_mutex_unlock( &counterMut );
  if(pthread_create( &thread1, &attr, slave, (void*) p)){
    fprintf(stderr,"[%s] Problem spawning thread\n%s at chunknumber:%d\n",__FUNCTION__,strerror(errno),p->chunkNumber);
    while(1){
      fprintf(stderr,"[%s] Problem spawning thread\n%s at chunknumber:%d\n",__FUNCTION__,strerror(errno),p->chunkNumber);
      
    }
    exit(0);
  }
  // pthread_detach(thread1);
}


void changeChr(int refId){
  fprintf(stderr,"[%s.%s():%d] refid:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);
  ((abcFilter *)allMethods[0])->readSites(refId);
  ((abcGL *)allMethods[4])->changeChr(refId);//used when changing chr;
  ((abcWriteFasta *)allMethods[19])->changeChr(refId);//used when changing chr;
  ((abcSmartCounts *)allMethods[20])->changeChr(refId);//used when changing chr;
  ((abcSaf *)allMethods[11])->changeChr(refId);//used when changing chr;
  ((abcPSMC *)allMethods[27])->changeChr(refId);//used when changing chr;

  void flush_queue();
  flush_queue();
}


void waiter(int refId){
  //fprintf(stderr,"_%s_\n",__FUNCTION__);fflush(stderr);

  //change of chr detected wait untill all threads are done
  if(allMethods[0]->header->n_targets<THRESHOLD_FOR_NICEOUTPUT)
    fprintf(stderr,"Change of chromo detected Waiting for nThreads:%d\n",curRunning);
  while(1){
    pthread_mutex_lock(&counterMut);
    if(curRunning==0&&printer.contains()==0){
      pthread_mutex_unlock(&counterMut);//need to do this because rest of while is not used
      break;
    }
    if(allMethods[0]->header->n_targets<THRESHOLD_FOR_NICEOUTPUT){
      fprintf(stderr,"Change of chromo detected Waiting for nThreads:%d printer.contains=%lu\n",curRunning,printer.contains());
      fflush(stderr);
    }
    pthread_mutex_unlock(&counterMut);
    sleep(1);
  }
  if(currentChr==-1||refId!=currentChr){
    currentChr=refId;
    changeChr(refId);
  }
  
}




/*
  This is the function that determines whether or not to start a new thread
*/
void selector(funkyPars *p){
  // fprintf(stderr,"funkypars\n");
  if(p==NULL){
    p = funkyPars_init();
    p->killSig = 1;
  }else{
    p->extras = new void*[andersSux];//funky
    memset(p->extras,0,sizeof(void*)*andersSux);
  }
  p->chunkNumber = chunkNumber++;

  if(maxThreads==1)
    interface(p);
  else if(maxThreads==2)//simple serialization, allow one chunk to be computated
    master2(p);
  else
    master(p);
}


/*
  initialize all pointers to zero
*/

funkyPars *funkyPars_init(){

  funkyPars *r = new funkyPars;
  r->numSites =0;
  r->extras = NULL;

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
  r->killSig =0;
  r->posi = NULL;
  r->refId = -1;
  return r;
}





void funkyPars_destroy(funkyPars *p) {
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
    delete [] p->major;
    p->major=NULL;
  }if(p->minor){
    delete [] p->minor;
    p->minor=NULL;
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

//only one instance at a time is running this function
size_t total_number_of_sites_unfiltred =0;
size_t total_number_of_sites_filtered =0;
void printFunky(funkyPars *p){
  //  fprintf(stderr,"printFunky killsig=%d nsites=%d refid:%d\n",p->killSig,p->numSites,p->refId);
  if(p->killSig==0) {//don't print the empty killSig chunk
    if((p->chunkNumber%howOften)==0){
      if(isAtty)
	fprintf(stderr,"\r\t-> Printing at chr: %s pos:%d chunknumber %d contains %d sites",header->target_name[p->refId],p->posi[0]+1,p->chunkNumber,p->numSites);
      else
	fprintf(stderr,"\t-> Printing at chr: %s pos:%d chunknumber %d contains %d sites\n",header->target_name[p->refId],p->posi[0]+1,p->chunkNumber,p->numSites);
    }if(p->numSites!=0){
      total_number_of_sites_unfiltred += p->numSites;
      for(int i=0;i<p->numSites;i++)
	if(p->keepSites[i])
	  total_number_of_sites_filtered++;
      for(int i=0;i<andersSux;i++)
	if(shouldRun[i])
	  allMethods[i]->print(p);
    }
   
    funkyPars_destroy(p);
    
  }else{
    funkyPars_destroy(p);
    pthread_mutex_unlock(&mUpPile_mutex);
    fprintf(stderr,"\n");//after last positions print add a neu line mutafuka
  }

}



const char *angsd_version(){
  return ANGSD_VERSION;
}
