/*

abcAMOVA: Analysis of Molecular Variance

*/

#include <cassert>
#include <ctype.h>
#include <cmath> 
#include <stdlib.h>

#include "aio.h"
#include "shared.h" 
#include "analysisFunction.h" 

#include "abcAMOVA.h"



void abcAMOVA::printArg(FILE *argFile){
	fprintf(argFile,"------------------------\n%s:\n",__FILE__);
	fprintf(argFile,"-doAMOVA\t%d\n",doAMOVA);
	fprintf(argFile,"\t\t1: init\n");
}

void abcAMOVA::getOptions(argStruct *arguments){

	doAMOVA=angsd::getArg("-doAMOVA",doAMOVA,arguments);

	if(doAMOVA==0){
		shouldRun[index]=0;
		return;
	}

	if (doAMOVA == 1){
		int doMajorMinor = 0;
		doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);

		if(doMajorMinor!=1){
			fprintf(stderr,"\t->\t-doAMOVA requires -doMajorMinor 1\n");
			exit(0);
		}

		int GL = 0;
		GL = angsd::getArg("-GL",GL,arguments);

		if(GL!=1){
			fprintf(stderr,"\t->\t-doAMOVA requires -GL 1\n");
			exit(0);
		}
	}else{
		fprintf(stderr,"\t->\tuse -doAMOVA 1\n");
		exit(0);
	}


}

//constructor
abcAMOVA::abcAMOVA(const char *outfiles,argStruct *arguments,int inputtype){
	doAMOVA = 0;
	refName = NULL;
	outfile = NULL;

	int doMajorMinor = 0;
	int GL = 0;

	char *empty;
	//like_calc = new phys_genolike_calc( empty );

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
	fprintf(outfile,"chr\tpos\tmajor\tminor\tlikes");

}


//destructor
abcAMOVA::~abcAMOVA(){
	if(doAMOVA==0)
		return;

	if(outfile!=NULL) 
		fclose(outfile);
}



void abcAMOVA::print(funkyPars *pars){

	if(doAMOVA==1){

		chunkyT *chk = pars->chk;

		for(int s=0;s<pars->numSites;s++){

			for(int i=0;i<pars->nInd;i++){

				assert(pars->likes!=NULL);

				//print [chr pos major minor likes]
				fprintf(outfile,"\n%s %d %c %c %d\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->likes[s]);
				fprintf(stderr,"\n%s %d %c %c %d\n",header->target_name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->likes[s]);

				tNode *nd = chk->nd[s][i];

				//for ind_i i, loop through all alignments at pos site_i s
				for(int l=0;l<nd->l;l++){
					//base_i
					char c = nd->seq[l];
				}
			}
		}
	}
}

void abcAMOVA::run(funkyPars *pars){
	if(doAMOVA==0)
		return;
}

void abcAMOVA::clean(funkyPars *pars){
}

