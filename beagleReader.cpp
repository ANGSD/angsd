
#include "analysisFunction.h"
#include "beagleReader.h"
#include "aio.h"

beagle_reader::~beagle_reader(){
	free(buffer);
	free(buf);
}


beagle_reader::beagle_reader(gzFile gz_a,const aMap *revMap_a,int intName_a,int &nInd_a){
	gz=gz_a;
	revMap=revMap_a;
	l=128;
	intName = intName_a;

	original=buffer=(char *) calloc(l,sizeof(char));
	buf=(char *) calloc(l,sizeof(char));


	const char *delims = "\t \n";

	int nCol=1;

	int presize=l;
	aio::tgets(gz,&buffer,&l);
	if(buffer!=original)
		original=buffer;


	if(l>presize){
		buf=(char*)realloc(buf,strlen(buffer)+1);
		presize=l;
	}

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
		aio::doAssert(tsk!=NULL,1,AT,"");
		post[i] = atof(tsk);
	}
}


//parse beagle fields
//
//chr_A_B_C -> chr: chr_A_B ; pos: C
//bpost = to give parsepost
//bufpos = position info
//buf = contig id
//bufmaj = major
//bufmin = minor
void parse_beagle_fields(char *buffer, char **buf, char **bufpos, char **bufmaj, char **bufmin, char **bpost,const char *delims){

	
	// copy buffer to working buffer buf
	// we need to keep buffer as it is 
	// size of buf is reallocated together with buffer
	// buff is cleaned together with buffer
	strcpy(*buf,buffer);

	char *a;
	
	//line= chr_a_b_123 g c 0.05 0.05 0.05
	//point to the first delim
	//a=chr_a_b_123
	a = strtok_r(*buf,delims,bpost);

	//bpost is now after first delim
	//bpost=g c 0.05 0.05 0.05
	// fprintf(stderr,"\n\nbpost=%s\n",*bpost);
	if ((a=strrchr(a,'_'))!=NULL){
		//a= _123
		//a+1= 123

		*bufpos=a+1;
		//bufpos now stores the position info=123
	}

	//here buf=chr_a_b_123
	
	//terminate after contig id
	//write \0 to the start of position info
	*strrchr(*buf,'_')='\0';

	//here buf now contains contig id=chr_a_b
	//extract major and minor from bpost
	*bufmaj=&refToChar[strtok_r(*bpost,delims,bpost)[0]];
	*bufmin=&refToChar[strtok_r(NULL,delims,bpost)[0]];

}


funkyPars *beagle_reader::fetch(int chunksize){

	static const char *delims = "\t \n";
	double **post = new double*[chunksize];

	// double **postbuf = new double*[chunksize];

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
		//parse an entire site:
		// fprintf(stdout,"nSites %d\n",nSites);
		myfunky->refId = lastRefId;

		parse_beagle_fields(buffer,&buf,&bufpos,&bufmaj,&bufmin,&bpost,delims);

		myfunky->posi[nSites] = atoi(bufpos)-1;

		myfunky->major[nSites] = *bufmaj;
		myfunky->minor[nSites] = *bufmin;

		parsepost(bpost,post[nSites],nInd,delims);    

		nSites++;
		changed =0;
	}
	buffer=original;


	int presize=l;

	while(aio::tgets(gz,&buffer,&l)) {


		if(l>presize){
			buf=(char*)realloc(buf,strlen(buffer)+1);
			presize=l;
		}

		if(buffer!=original)
			original=buffer;


		if(intName){

			//identifier is chr_pos or snp name

			parse_beagle_fields(buffer,&buf,&bufpos,&bufmaj,&bufmin,&bpost,delims);

			//buf now contains contig id
			aMap ::const_iterator it = revMap->find(buf);

			if(it==revMap->end()){
				fprintf(stderr,"\t-> Problem finding chr:%s from faifile\n",buf);
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

		parse_beagle_fields(buffer,&buf,&bufpos,&bufmaj,&bufmin,&bpost,delims);
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
