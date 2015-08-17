#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NUM_THREADS	10
pthread_barrier_t   barrier; // barrier synchronization object


size_t data[NUM_THREADS];


void *BusyWork(void *t)
{
   int i;
   long tid;
   tid = (long)t;
   pthread_barrier_wait (&barrier);
   while(1){
     int sleepval=lrand48() % 5+1;
     printf("INLOOP Thread %ld is taken:%lu is sleeping:%d\n",tid,data[tid], sleepval);
     sleep(sleepval);
     pthread_barrier_wait (&barrier);
     fprintf(stderr,"after barrier\n");
     fflush(stderr);   
   }
     
   pthread_exit((void*) t);
}

int runner(){
  static int i=0;
  

}



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

   for(t=0; t<NUM_THREADS; t++) {
      printf("Main: creating thread %ld\n", t);
      rc = pthread_create(&thread[t], &attr, BusyWork, (void *)t); 
      if (rc) {
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         exit(-1);
         }
      }
   for(i=0;i<100;i++)
     runner();
   /* Free attribute and wait for the other threads */
   pthread_attr_destroy(&attr);
   for(t=0; t<NUM_THREADS; t++) {
      rc = pthread_join(thread[t], &status);
      if (rc) {
         printf("ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
         }
      printf("Main: completed join with thread %ld having a status   of %ld\n",t,(long)status);
      }
 
printf("Main: program completed. Exiting.\n");
pthread_exit(NULL);
}
