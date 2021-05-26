#include <cstdio>
#include <cstring>
#include <cassert>

#include "argStruct.h"
#include "analysisFunction.h"
#include "shared.h"

typedef struct{
  unsigned LH;
  int layout;
  unsigned M;
  unsigned N;
  int compressed;
  int si;//<- sample identifier. If this is set we have names below.
  char **sampleids;
}header;


typedef struct{

  unsigned indis;
  char *Lid;
  char *Lrsid;
  char *Lchr;
  unsigned int vpos;
  unsigned short nal;
  unsigned la_l1;
  unsigned la_l2;
  char *la1;
  char *la2;
  double *probs;

}bgenLine;


class bgenReader{
  //fields of class
private:
  int readAgain;
  int curChr;
  int prevChr;
  int onlyPrint;
  const aMap *revMap;

public:
  //for reading in chunk of bgen file
  funkyPars *fetch(int chunkSize);
  FILE *bgenFile;
  int layout;
  int sites;
  int sitesRead;
  int nInd;
  header *hd;
  header *parseheader(FILE *fp);
  bgenLine *bgen;
  bgenLine *parseline(FILE *fp,header *hd);
  void funkyCopy(bgenLine *bgen, funkyPars *r, int &balcon);
  
  //constructor of class
  bgenReader(char *fname, const aMap *revMap_a,  int intName_a,int &nInd_a){

    FILE *fp=fopen(fname,"rb");
    assert(fp!=NULL);

    unsigned offset;
    //unsigned is 4 bytes
    assert(fread(&offset,sizeof(unsigned),1,fp)==1);
    
    //read header function
    hd = parseheader(fp);
    //jump offset bytes
    layout = hd->layout;
    sites = hd->M;
    nInd = hd->N;
    //keeps track of how many sites read
    sitesRead = 0;

    bgenLine *bgen2 = new bgenLine;
    bgen = bgen2;
    bgen->probs = new double[3*nInd];
    
    //emil jump offset + 4 bytes into file (bytes of offset), from start
    fseek(fp, offset+4, SEEK_SET);

    bgenFile = fp;
    //declare other helper stuff
    readAgain=0;
    curChr=-1;
    prevChr=-1;
    onlyPrint=0;

    revMap=revMap_a;
    
  }

  //destructor
  ~bgenReader(){
    
    fclose(bgenFile);

    delete [] bgen->probs;
    delete bgen;
    
    if(hd->si==1){

      for(uint i=0;i<hd->N;i++){	
	free(hd->sampleids[i]);
      }
      delete [] hd->sampleids;
    }
        
    delete hd;
    //TO DO
    //delete header struct
  }
  
};
