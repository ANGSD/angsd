#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define nthreads 10

pthread_t thd[nthreads];

pthread_mutex_t mutex;
pthread_cond_t count_threshold_cv;

pthread_barrier_t barr;
void *inner(void *ptr){
  size_t threadid=(size_t)ptr;
  fprintf(stderr,"thread: %lu speaks\n",threadid);fflush(stderr);
  sleep(1);
  while(1){
    fprintf(stderr,"in lloop\n");fflush(stderr);
    pthread_mutex_lock(&mutex);
    pthread_cond_wait(&count_threshold_cv, &mutex);
    int sleepval = lrand48() % 5 +1;
    fprintf(stderr,"thread: %lu will wait:%d\n",threadid,sleepval);fflush(stderr);  
    sleep(sleepval);

    int rc = pthread_barrier_wait(&barr);
    fprintf(stderr,"rc:%d\n",rc);fflush(stderr);
    fprintf(stderr,"in lloop after barrier\n");fflush(stderr);
    pthread_mutex_unlock(&mutex);
  }
  fprintf(stderr,"never here\n");fflush(stderr);
}




int outer(size_t ntimes){
  pthread_t thread1;
  size_t i,n;
  for(i=0;i<nthreads;i++){
    if(pthread_create( &thread1, NULL, inner, (void*) i))
      fprintf(stderr,"Problems creating thread\n");
  }

  for(n=0;n<ntimes;n++){
    pthread_mutex_lock(&mutex);
    pthread_cond_broadcast(&count_threshold_cv);
    pthread_mutex_unlock(&mutex);
  }
  for (i=0; i<nthreads; i++) {
    pthread_join(thd[i], NULL);
  }
}
int main(){
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init (&count_threshold_cv, NULL);
  pthread_barrier_init (&barr, NULL, nthreads);
  outer(5);

}
