#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define NUM_THREADS	3
pthread_t thread[NUM_THREADS];
pthread_attr_t attr;

int SIG_COND2 =1;//<- this could cause race condition, but lets assume printing takes some time

pthread_barrier_t   barrier; // barrier synchronization object

typedef struct{
  int chunknr;//<-chunknr -1 indicates no data
}astruct; //<- NULL indicates shutdown thread

astruct **data=NULL;

void doanal(int tid,astruct *as){
  if(as!=NULL){
    int sleepval=lrand48() % 2+1;
    sleep(sleepval);
    as->chunknr += 100;
  }else{
    fprintf(stderr,"donal is dryrun at:%d\n",tid);
    fflush(stderr);
  }
}

// this will run eternaly untill data[threadid] is NULL
void *BusyWork(void *t){
   int i;
   long tid;
   tid = (long)t;
   pthread_barrier_wait (&barrier);//we are locking all threads, untill we are certain we have data
   while(SIG_COND2){
     pthread_barrier_wait (&barrier);
     doanal(tid,data[tid]);
     pthread_barrier_wait (&barrier);//we are waiting here to make sure all threads are finished with analysis
   }
   fprintf(stderr,"Will shutdown threadid: %ld\n",tid);fflush(stderr);
   pthread_exit((void*) t);
}

void closethreads(){
  pthread_barrier_wait (&barrier);
  fprintf(stderr,"Will close threads\n");fflush(stderr);
  SIG_COND2=0;
  pthread_barrier_wait (&barrier);
}


int runner(astruct *as){
  assert(as!=NULL);
  static int batch=-1;//<- this is the inarray position id
  int i;

  //plugin data in array if we have data;
  if(as->chunknr!=-1){
    batch = as->chunknr %NUM_THREADS;
    data[batch] = as;

  }
  fprintf(stderr,"[%s] chunknr:%d batch:%d\n",__FUNCTION__,as->chunknr,batch);fflush(stderr); 
  //case where we launch all analysis threads
  if(batch==NUM_THREADS-1){
    //    fprintf(stderr,"FULL Will launch all analysis:\n");fflush(stderr);
    pthread_barrier_wait (&barrier);//RUN ALL THREADS
    pthread_barrier_wait (&barrier);//MAKE SURE THEY ARE FINISHED
    for(i=0;i<NUM_THREADS;i++){
      fprintf(stderr,"RESULTS-> %d) chunknr:%d\n",i,data[i]->chunknr);
      free(data[i]);
}    
fflush(stderr);
  memset(data,0,sizeof(astruct*)*NUM_THREADS);
batch=-1;
  }else if(as->chunknr==-1){
free(as);
    //fprintf(stderr,"SUBSUB Will launch all analysis but only to bach:%d:\n",batch);fflush(stderr);
    for(i=batch+1;i<NUM_THREADS;i++){
      fprintf(stderr,"data[%d] is set to null\n",i);fflush(stderr);
      data[i] = NULL;
    }
    pthread_barrier_wait (&barrier);//RUN ALL THREADS
    pthread_barrier_wait (&barrier);//MAKE SURE THEY ARE FINISHED
    for(i=0;i<=batch;i++){
      fprintf(stderr,"RESULTS-> %d) chunknr:%d\n",i,data[i]->chunknr);
free(data[i]);
}batch=-1;
    fflush(stderr);
  memset(data,0,sizeof(astruct*)*NUM_THREADS);
  }
}

int selector(int a){
  static int counter=0;
  astruct *as =calloc(1,sizeof(astruct)); 
  if(a!=-1)
    as->chunknr =counter++;
  else
    as->chunknr=-1;
  runner(as);
}



#ifdef __WITH_MAIN__

int main (int argc, char *argv[])
{
  

   int rc;
   long t;
   int i;
   void *status;

   /* Initialize and set thread detached attribute */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   data = calloc(NUM_THREADS,sizeof(astruct*));
   pthread_barrier_init (&barrier, NULL, NUM_THREADS+1);
   
   for(t=0; t<NUM_THREADS; t++) {
     rc = pthread_create(&thread[t], &attr, BusyWork, (void *)t); 
     if (rc) {
       fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }
   }
   pthread_barrier_wait (&barrier);//we are waiting here to make sure all threads are finished with analysis
   //sleep(100);
   for(i=0;i<7;i++){
     selector(i);
     int sleep_val = drand48()*1e6+1;
     usleep(sleep_val);
   }
   selector(-1);
   closethreads();
   
   /* Free attribute and wait for the other threads */
   pthread_attr_destroy(&attr);
   for(t=0; t<NUM_THREADS; t++) {
      rc = pthread_join(thread[t], &status);
      if (rc) {
	fprintf(stderr,"ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
         }
      fprintf(stderr,"Main: completed join with thread %ld having a status   of %ld\n",t,(long)status);
      }
 
   fprintf(stderr,"Main: program completed. Exiting.\n");
free(data);
   pthread_exit(NULL);
}
#endif
