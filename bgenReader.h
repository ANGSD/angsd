#include <cstdio>
#include <cstring>

#include "argStruct.h"
#include "analysisFunction.h"


typedef struct{
  unsigned LH;
  int layout;
  unsigned M;
  unsigned N;
  int compressed;
  int si;//<- sample identifier. If this is set we have names below.
  char **sampleids;
}header;


class bgenReader{
  //fields of class
private:
  int readAgain;
  int curChr;
  int prevChr;
  int onlyPrint;
  void parseline(FILE *fp,funkyPars *r,int &balcon,header *hd);
public:
  //for reading in chunk of bgen file
  funkyPars *fetch(int chunkSize);
  FILE *bgenFile;
  int layout;
  int sites;
  int nInd;
  header *hd;
  header *parseheader(FILE *fp);

  //constructor of class
  bgenReader(char *fname,  int intName_a,int &nInd_a){

    FILE *fp=fopen(fname,"rb");
    assert(fp!=NULL);

    unsigned offset;
    //unsigned is 4 bytes
    fread(&offset,sizeof(unsigned),1,fp);
    
    //read header function
    hd = parseheader(fp);
    //jump offset bytes
    layout = hd->layout;
    sites = hd->M;
    nInd = hd->N;

    //emil jump offset + 4 bytes into file (bytes of offset), from start
    fseek(fp, offset+4, SEEK_SET);

    bgenFile = fp;
    //declare other helper stuff
    readAgain=0;
    curChr=-1;
    prevChr=-1;
    onlyPrint=0;
    
  }

  //destructor
  ~bgenReader(){
    
    fclose(bgenFile);

    if(hd->si==1){

      for(uint i=0;i<hd->N;i++){	
	free(hd->sampleids[i]);
      }
      delete [] hd->sampleids;
    }
        
    delete [] hd;
    //TO DO
    //delete header struct
  }
  
};
