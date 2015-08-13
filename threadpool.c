#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define nthreads 10

pthread_t thd[nthreads];

pthread_mutex_t mutex;
pthread_cond_t count_threshold_cv;
//This is taken from here:
//http://blog.albertarmea.com/post/47089939939/using-pthread-barrier-on-mac-os-x
#ifdef __APPLE__

#ifndef PTHREAD_BARRIER_H_
#define PTHREAD_BARRIER_H_

#include <pthread.h>
#include <errno.h>

typedef int pthread_barrierattr_t;
typedef struct
{
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int tripCount;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count)
{
    if(count == 0)
    {
        errno = EINVAL;
        return -1;
    }

if(pthread_mutex_init(&barrier->mutex, 0) < 0)
    {
        return -1;
    }
    if(pthread_cond_init(&barrier->cond, 0) < 0)
    {
        pthread_mutex_destroy(&barrier->mutex);
        return -1;
    }

    barrier->tripCount = count;
    barrier->count = 0;

    return 0;
}
int pthread_barrier_destroy(pthread_barrier_t *barrier)
{
    pthread_cond_destroy(&barrier->cond);
    pthread_mutex_destroy(&barrier->mutex);
    return 0;
}

int pthread_barrier_wait(pthread_barrier_t *barrier)
{
    pthread_mutex_lock(&barrier->mutex);
    ++(barrier->count);
    if(barrier->count >= barrier->tripCount)
    {
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }
    else
    {
        pthread_cond_wait(&barrier->cond, &(barrier->mutex));
        pthread_mutex_unlock(&barrier->mutex);
        return 0;
    }
}

#endif // PTHREAD_BARRIER_H_
#endif // __APPLE__
pthread_barrier_t barr;
void *inner(void *ptr){
  size_t threadid=(size_t)ptr;
  fprintf(stderr,"thread: %lu speaks\n",threadid);fflush(stderr);
  while(1){
    fprintf(stderr,"in lloop\n");fflush(stderr);
    pthread_mutex_lock(&mutex);
    pthread_cond_wait(&count_threshold_cv, &mutex);
    int sleepval = lrand48() % 5 +1;
    fprintf(stderr,"thread: %lu will wait:%d\n",threadid,sleepval);fflush(stderr);  
    sleep(sleepval);
    pthread_mutex_unlock(&mutex);
    int rc = pthread_barrier_wait(&barr);
    fprintf(stderr,"in lloop after barrier\n");fflush(stderr);
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
    pthread_cond_broadcast(&count_threshold_cv);
    
  }
  
}
int main(){
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init (&count_threshold_cv, NULL);
  outer(5);

}
