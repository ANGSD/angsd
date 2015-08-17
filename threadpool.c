#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NUM_THREADS	4

pthread_barrier_t   barrier; // barrier synchronization object

int data[NUM_THREADS];

void doanal(int a){
  int sleepval=lrand48() % 10+1;
  fprintf(stderr,"INLOOP Thread %ld is taken:%lu is sleeping:%d\n",a,data[a], sleepval);   fflush(stderr);   
  sleep(sleepval);
}

void *BusyWork(void *t){
   int i;
   long tid;
   tid = (long)t;
   pthread_barrier_wait (&barrier);
   while(1){
     fprintf(stderr,"is wiaigin\n");fflush(stderr);
     pthread_barrier_wait (&barrier);
     doanal(tid);
   }
   pthread_exit((void*) t);
}

int runner(int chunknr){
  int batch = chunknr %NUM_THREADS;
  data[batch] = chunknr;
  fprintf(stderr,"[%s] chunknr:%d batch:%d\n",__FUNCTION__,chunknr,batch);
  if(batch==NUM_THREADS-1){
    fprintf(stderr,"THREAD[%d] will run\n",batch);fflush(stderr);
    pthread_barrier_wait (&barrier);
    doanal(batch);
    fprintf(stderr,"THREAD[%d] Will print also\n",batch);fflush(stderr);
  }
 }

#ifdef __WITH_MAIN__

int main (int argc, char *argv[])
{
   pthread_t thread[NUM_THREADS];
   pthread_attr_t attr;
   int rc;
   long t;
   int i;
   void *status;
   pthread_barrier_init (&barrier, NULL, NUM_THREADS);
   /* Initialize and set thread detached attribute */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   for(t=0; t<NUM_THREADS-1; t++) {
     data[t]=-1;
     rc = pthread_create(&thread[t], &attr, BusyWork, (void *)t); 
     if (rc) {
       fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
       exit(-1);
     }

   }
   pthread_barrier_wait (&barrier);
   for(i=0;1&&i<100;i++)
     runner(i);



   
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
pthread_exit(NULL);
}
#endif
