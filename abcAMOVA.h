#pragma once
#include "abc.h"

//typedef struct{
  //int **dat;
  //int *major;
//}haploCalls;

class abcAMOVA:public abc{
    private:
        int doAMOVA;
        char *refName;
        FILE *outfile;
        //non optional arguments
        int doHaploCall;
        int doCount;

        //optional arguments
        int maxMis;
        int minMinor;

        //out file
        BGZF* outfileZ;

        //functions
        void printHaplo(funkyPars *pars);
        void getHaplo(funkyPars *pars);

        //print buffer
        kstring_t bufstr;

    public:
        void run(funkyPars  *pars);
        void clean(funkyPars *pars);  
        void print(funkyPars *pars);  
        void openfile(const char *outfiles);
        void getOptions(argStruct *arguments);
        void printArg(FILE *argFile);
        abcAMOVA(const char *outfiles,argStruct *arguments,int inputtype);
        ~abcAMOVA();

};


