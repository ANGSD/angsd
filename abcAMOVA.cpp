/*

abcAMOVA: Analysis of Molecular Variance


notes:
-site should be diallelic

*/

#include <cassert>
#include <ctype.h>
#include <cmath> 
#include <stdlib.h>

#include "aio.h"
#include "shared.h" 
#include "analysisFunction.h" 
#include "abcCallGenotypes.h"

#include "abcAMOVA.h"
#include "abcHaploCall.h"



void abcAMOVA::printArg(FILE *argFile){
    fprintf(argFile,"------------------------\n%s:\n",__FILE__);
    fprintf(argFile,"-doAMOVA\t%d\n",doAMOVA);
    fprintf(argFile,"\t\t1: init\n");
}

void abcAMOVA::getOptions(argStruct *args){

    doAMOVA=angsd::getArg("-doAMOVA",doAMOVA,args);

    if(doAMOVA==0){
        shouldRun[index]=0;
        return;
    }

    if (doAMOVA == 1){


        int doMajorMinor = 0;
        doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,args);

        int doGeno=0;
        doGeno = angsd::getArg("-doGeno",doGeno,args);
        if(doGeno!=2){
            fprintf(stderr,"\t->\t-doAMOVA %d requires -doGeno 2\n",doAMOVA);
            exit(0);
        }


        //based on allele counts
        if(doMajorMinor!=2){
            fprintf(stderr,"\t->\t-doAMOVA %d requires -doMajorMinor 2\n",doAMOVA);
            exit(0);
        }

        int doCounts=0;
        doCounts=angsd::getArg("-doCounts",doCounts,args);

        if(doCounts!=1){
            fprintf(stderr,"\t->\t-doAMOVA %d requires -doCounts 1\n",doAMOVA);
            exit(0);
        }

    }else if (doAMOVA == 2){

        int doMajorMinor = 0;
        doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,args);

        if(doMajorMinor!=1){
            fprintf(stderr,"\t->\t-doAMOVA %d requires -doMajorMinor 1\n",doAMOVA);
            exit(0);
        }

        int GL = 0;
        GL = angsd::getArg("-GL",GL,args);

        if(GL!=1){
            fprintf(stderr,"\t->\t-doAMOVA %d requires -GL 1\n",doAMOVA);
            exit(0);
        }

    }else{
        fprintf(stderr,"\t->\tuse -doAMOVA 1 or 2\n");
        exit(0);
    }


}

//constructor
abcAMOVA::abcAMOVA(const char *outfiles,argStruct *arguments,int inputtype){
    doAMOVA = 0;
    refName = NULL;
    outfile = NULL;



    if (arguments->argc==2){
        if(!strcasecmp(arguments->argv[1],"-doAMOVA")){
            printArg(stdout);
            exit(0);
        } else
            return;
    }

    getOptions(arguments);

    if(doAMOVA==0){
        shouldRun[index] =0;
        return ;
    }

    printArg(arguments->argumentFile);
    outfile = aio::openFile(outfiles,".amovaout");

}


//destructor
abcAMOVA::~abcAMOVA(){
    if(doAMOVA==0)
        return;

    if(outfile!=NULL) 
        fclose(outfile);
}



void abcAMOVA::print(funkyPars *pars){

    /*
     * 0 = MM
     * 1 = Mm
     * 2 = mm
     *
     * M[i][j]
     * number of sites where ind1=i and ind2=j
     */
    int pairwiseGenotypeMatrix[3][3] = {{0}};

    if(doAMOVA==1){
        genoCalls *geno =(genoCalls *) pars->extras[11];

        for(int s=0;s<pars->numSites;s++){

            for(int i=0;i<pars->nInd;i++){

                //fprintf(stderr,"%d\t",geno->dat[s][i]);

                //print actual haplotypes
                //if(geno->dat[s][i]==0)
                //fprintf(stderr,"%c%c\n",intToRef[pars->major[s]],intToRef[pars->major[s]]);
                //else if(geno->dat[s][i]==1)
                //fprintf(stderr,"%c%c\n",intToRef[pars->major[s]],intToRef[pars->minor[s]]);
                //else if(geno->dat[s][i]==2)
                //fprintf(stderr,"%c%c\n",intToRef[pars->minor[s]],intToRef[pars->minor[s]]);
                //else if(geno->dat[s][i]==-1)
                //fprintf(stderr,"NN\n");

            }
            //fprintf(stderr,"%d %d\t",geno->dat[s][0],geno->dat[s][1]);

            //ind0 ind1
            pairwiseGenotypeMatrix[geno->dat[s][0]][geno->dat[s][1]]++;
        }

        //fprintf(stderr,"\n\n");
        for(int gi=0; gi<3;gi++){
            for (int gj=0; gj<3; gj++){
                //fprintf(stderr,"%d\t",pairwiseGenotypeMatrix[gi][gj]);
                fprintf(outfile,"%d\t",pairwiseGenotypeMatrix[gi][gj]);

            }
            //fprintf(stderr,"\n");
            fprintf(outfile,"\n");
        }
    }
}



void abcAMOVA::run(funkyPars *pars){

    if(doAMOVA==0)
        return;
}

void abcAMOVA::clean(funkyPars *pars){
}

