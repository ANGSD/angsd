
#include <cassert>
#include "analysisFunction.h"
#include "beagleReader.h"
#include "aio.h"

beagle_reader::~beagle_reader(){
	free(buffer);
}


beagle_reader::beagle_reader(gzFile gz_a,const aMap *revMap_a,int intName_a,int &nInd_a){
	gz=gz_a;
	revMap=revMap_a;
	l=128;
	intName = intName_a;

	original=buffer =(char *) calloc(l,sizeof(char));
	const char *delims = "\t \n";

	int nCol=1;

	aio::tgets(gz,&buffer,&l);
	if(buffer!=original)
		original=buffer;
	strtok_r(buffer,delims,&buffer);
	while(strtok_r(NULL,delims,&buffer))
		nCol++;
	if(nCol % 3 ){
		fprintf(stderr,"\t-> Number of columns should be a multiple of 3, nCol=%d\n",nCol);
		exit(0);
	} 
	nInd_a=nCol/3-1;
	nInd = nInd_a;
}


void parsepost(char *buffer,double *post,int nInd,const char *delims){
	for(int i=0;i<nInd*3;i++){
		char *tsk = strtok_r(NULL,delims,&buffer);
		assert(tsk!=NULL);
		post[i] = atof(tsk);
	}
}


//parse beagle fields
//
//chr_A_B_C -> chr: chr_A_B ; pos: C
//bpost = to give parsepost
//bufpos = position info
//tok = chr info
//bufmaj = major
//bufmin = minor
void parse_beagle_fields(char *&buffer, char *&tok, char *&bufpos, char *&bufmaj, char *&bufmin, char *&bpost){

	static const char *delims = "\t \n";
	char *buf;
	char *bufchr;
	char *tok1;

	buf=strdup(buffer);
	// bpost = strchr(buffer,"\t ");
	// bpost = strchr(buffer,*delims);


	bpost=strchr (buffer, ' ');
	if(bpost == NULL){
		bpost=strchr (buffer,'\t');
	}


	
	bufmaj=&refToChar[strtok_r(NULL,delims,&bpost)[0]];
	bufmin=&refToChar[strtok_r(NULL,delims,&bpost)[0]];


	// tok = strtok_r(buffer,(const char *)" ",&buffer);
	tok = strtok_r(buffer,(const char *)delims,&buffer);
	tok1=strrchr(tok,'_');
	bufpos=strrchr(tok,'_');
	bufpos=bufpos+1;
	*tok1='\0';

	bufchr= tok;
	tok=bufchr;
	buffer=buf;
}


funkyPars *beagle_reader::fetch(int chunksize){

	static const char *delims = "\t \n";
	static const char *delims2 = "_\t \n";
	double **post = new double*[chunksize];

	char *bufpos;
	char *bufmaj;
	char *bufmin;
	char *bpost;

	funkyPars * myfunky =funkyPars_init();
	myfunky->posi = new int[chunksize];
	myfunky->major = new char[chunksize];
	myfunky->minor = new char[chunksize];


	for(int s=0;s<chunksize;s++){
		post[s] = new double[nInd*3];
	}

	int nSites=0;
	static int positions =0;//every site is a new position across different chunks
	static int lastRefId =-1;
	static int changed =0;
READAGAIN:
	if(changed){
			char *tok;
		//parse an entire site:
		// fprintf(stdout,"nSites %d\n",nSites);
		myfunky->refId = lastRefId;
		// fprintf(stdout,"BUFfer= %s\n",buffer);

		parse_beagle_fields(buffer,tok,bufpos,bufmaj,bufmin,bpost);

		myfunky->posi[nSites] = atoi(bufpos)-1;

		myfunky->major[nSites] = *bufmaj;
		myfunky->minor[nSites] = *bufmin;

		parsepost(bpost,post[nSites],nInd,delims);    

		nSites++;
		changed =0;
	}
	buffer=original;


	while(aio::tgets(gz,&buffer,&l)) {


		if(buffer!=original)
			original=buffer;


		if(intName){

			//identifier is chr_pos or snp name

			char *tok;
		parse_beagle_fields(buffer,tok,bufpos,bufmaj,bufmin,bpost);


			aMap ::const_iterator it = revMap->find(tok);

			if(it==revMap->end()){
				fprintf(stderr,"\t-> Problem finding chr:%s from faifile\n",tok);
				exit(0);
			}
			if(lastRefId==-1){
				lastRefId = it->second;
			}
			if(lastRefId!=it->second){
				changed =1;
				lastRefId = it->second;
				if(nSites==0){// if chromosome if finish then read next one
					goto READAGAIN; 
				}
				break;
			}
			lastRefId = it->second;
			myfunky->refId = lastRefId;
			myfunky->posi[nSites] = atoi(bufpos)-1;
		}

		else{

			char *tok;
		parse_beagle_fields(buffer,tok,bufpos,bufmaj,bufmin,bpost);
			myfunky->refId = 0;
			myfunky->posi[nSites] = positions++;
		}


		myfunky->major[nSites] = *bufmaj;
		myfunky->minor[nSites] = *bufmin;

		parsepost(bpost,post[nSites],nInd,delims);
		buffer=original;

		nSites++;
		if(nSites>=chunksize)
			break;
	}

	if(nSites<chunksize){
		for(int s=nSites;s<chunksize;s++)
			delete[] post[s];
	}
	myfunky->nInd=nInd;
	myfunky->post=post;
	myfunky->numSites = nSites;

	if(nSites==0 & changed==0){
		fprintf(stdout,"Done reading beagle\n");
		funkyPars_destroy(myfunky);
		return(NULL);

	}
	return(myfunky);
}
