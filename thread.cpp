#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>//used for setting number of decimals in posterior
#include <cstdlib>
#include <zlib.h>
#include <sstream>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#include <pthread.h>

int nThreads;

pthread_t *threads = NULL;
pthread_t *threads1 = NULL;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexPrint = PTHREAD_MUTEX_INITIALIZER;
int NumJobs;
int jobs;
int *printArray;
int cunt;
int *data;
int *result;
int dataRead=0;
int currentJob;
int maxPrintQue;
int lengthOfPrintArray;
int reader();
int doPrint=0;
int getNprintJobs();



int printer(){

  int cunt=0;
  pthread_mutex_lock(&mutexPrint);
  while(cunt<NumJobs){

     if(printArray[cunt % lengthOfPrintArray]==1){
       
      pthread_mutex_unlock(&mutexPrint);
    
      int sleepval=lrand48() % 240+1;
      unsigned int microseconds=5*sleepval;
      usleep(microseconds);

     fprintf(stdout,"%d ",result[cunt % lengthOfPrintArray]); 

 
      pthread_mutex_lock(&mutexPrint);
      printArray[cunt % lengthOfPrintArray]=0;
      //  pthread_mutex_unlock(&mutexPrint);
      cunt++;

      // fflush(stdout);
     
       
    }
    else{
      if(doPrint)
      fprintf(stderr," printer waiting for jobs: \n");
      pthread_mutex_unlock(&mutexPrint);    
       unsigned int microseconds=100000;
      usleep(microseconds);
   
      pthread_mutex_lock(&mutexPrint);
    
    }
  }
 pthread_mutex_unlock(&mutexPrint);    
    
  return 1;
}

void waitForPrinter(){
        
   while(1){
      pthread_mutex_lock(&mutexPrint);
      int nPrintJobs=getNprintJobs();
      pthread_mutex_unlock(&mutexPrint);
 
     if(nPrintJobs<maxPrintQue)
	break;
    if(doPrint)
        fprintf(stderr," waiting for print: %d\n",nPrintJobs);
      int sleepval=lrand48() % 4+1;
      unsigned int microseconds=500*sleepval;
      usleep(microseconds);
       
   
   }

}



void waitForPrintSpace(int job){
        
   while(1){
      pthread_mutex_lock(&mutexPrint);
      int isFree=printArray[job % lengthOfPrintArray]==0;
      pthread_mutex_unlock(&mutexPrint);
 
      if(isFree){
	//	fprintf(stderr,"isFree %d\n",job);
	break;

      }
      unsigned int microseconds=100000;
      usleep(microseconds);
    
      
        if(doPrint)
    fprintf(stderr,"waiting for print space %d\n",job);
   }

}

void waitForPrintSpaceRead(int job){
        
   while(1){
      pthread_mutex_lock(&mutexPrint);
      int isFree=printArray[job % lengthOfPrintArray]==0;
      if(!isFree){
    if(doPrint)
  	fprintf(stderr,"waiting for print space READ %d, %d\n",job % lengthOfPrintArray,job);
      }
      pthread_mutex_unlock(&mutexPrint);
 
      if(isFree){
	//	fprintf(stderr,"isFree %d\n",job);
	break;

      }

    unsigned int microseconds=100000;
    usleep(microseconds);
    

     
   }

}


int reader(){
 //read in data
 

  for (int j=0;j<NumJobs;j++){

    //   fprintf(stderr,"reader:(%d) ",j);
    waitForPrintSpaceRead(j);
    int sleepval=lrand48() % 50+1;
    unsigned int microseconds=50*sleepval;
    usleep(microseconds);
    
    data[j  % lengthOfPrintArray]=j;
    pthread_mutex_lock(&mutex1);
    dataRead++;
    pthread_mutex_unlock(&mutex1);

  
  }
  return 1;

}

int getNprintJobs(){
  
  int nJobs=0;
  for(int i=0;i++;i<lengthOfPrintArray)
    nJobs+=printArray[i];
  
  return nJobs;
}


void *functionC(void *a) //the a means nothing
{
 int running_job;

  pthread_mutex_lock(&mutex1);

  while (currentJob<NumJobs) {
    running_job = currentJob++;
    pthread_mutex_unlock(&mutex1);
    if(running_job==-2){
      // fprintf(stderr,"reader:(%d) ",currentJob);
       int isDone=reader(); 
    }
    else if(running_job==-1){
      //    fprintf(stderr,"print:(%d) ",currentJob);
      int isDonePrint=printer(); 
      }
    else{
      //   fprintf(stderr," Job starting: %d\n", running_job);
       int datRe;
      while(1){
	//see how far the data is
	pthread_mutex_lock(&mutex1);
	int datRe=dataRead;
	pthread_mutex_unlock(&mutex1);

	//wait untill there is data
	if(running_job>=datRe){
	    if(doPrint)
	      fprintf(stderr," analysis waiting for jobs\n");

	  unsigned int microseconds=1000;
	  usleep(microseconds);
	  continue;
	}
	else{ // Do the analysis
	  int theData = data[running_job  % lengthOfPrintArray]; 
	  int sleepval=lrand48() % 20+1;
	  unsigned int microseconds=2200*sleepval;
	  usleep(microseconds);
	  //analysis goes here
	  int TmpResults= theData;   
	  waitForPrintSpace(running_job);
 

	 
	  //add results to printer que (can be included in waitFor..)
	  pthread_mutex_lock(&mutexPrint);
	  result[running_job % lengthOfPrintArray]=TmpResults;  
	  printArray[running_job % lengthOfPrintArray]=1;
	  pthread_mutex_unlock(&mutexPrint);
	  break;
 
	}
      }
 
    }
    pthread_mutex_lock(&mutex1);
  }
  pthread_mutex_unlock(&mutex1);

  return NULL;

}



int main(int argc, char *argv[]){

  nThreads=20;
  NumJobs=10000;
  maxPrintQue=nThreads*10; // maximal print or data que
  lengthOfPrintArray=maxPrintQue*2;

  fprintf(stderr," lengthOfPrintArray %d\n", lengthOfPrintArray);

  currentJob=-2;
  pthread_t thread1[nThreads];
  cunt = 0;

  data = new int[NumJobs];
  printArray = new int[lengthOfPrintArray];
  result = new int[lengthOfPrintArray];
  for (int j=0;j<NumJobs;j++){
    data[j] = -1;
   
  }
  
  for (int j=0;j<lengthOfPrintArray;j++){
    printArray[j] = 0;
  }
  

  //create all Threads
  for (int i = 0; i < nThreads; i++)
    pthread_create(&thread1[i], NULL, &functionC, NULL);

  
  // Wait all threads to finish
  for (int i = 0; i < nThreads; i++)
    pthread_join(thread1[i], NULL);


  return(0);
}
