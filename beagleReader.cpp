
#include <climits>
#include <cstdarg>
#include "analysisFunction.h"
#include "beagleReader.h"
#include "aio.h"

static const char *delims = "\t \r\n";

static void lineError(int lineNo,const char *fmt,...){
	fprintf(stderr,"\t-> Problem in beagle file at line %d: ",lineNo);
	va_list ap;
	va_start(ap,fmt);
	vfprintf(stderr,fmt,ap);
	va_end(ap);
	fprintf(stderr,". Will exit\n");
	exit(1);
}

beagle_reader::~beagle_reader(){
	free(buffer);
	free(buf);
}


beagle_reader::beagle_reader(gzFile gz_a,const aMap *revMap_a,int intName_a,int &nInd_a){
	gz=gz_a;
	revMap=revMap_a;
	intName = intName_a;
	lineNo=0;
	positions=0;
	changed=0;

	l=128;
	buffer=(char *) calloc(l,sizeof(char));
	if(aio::tgets(gz,&buffer,&l)==0){
		fprintf(stderr,"\t-> Beagle file is empty. Will exit\n");
		exit(1);
	}
	lineNo++;
	bufSize=l;
	buf=(char *) malloc(bufSize);
	strcpy(buf,buffer);

	//header: marker allele1 allele2 and three columns per individual
	int nCol=0;
	char *save=NULL;
	for(char *tok=strtok_r(buf,delims,&save);tok;tok=strtok_r(NULL,delims,&save))
		nCol++;
	if(nCol<6||nCol % 3){
		fprintf(stderr,"\t-> Beagle header should have marker, allele1, allele2 and three columns per individual, found %d columns. Will exit\n",nCol);
		exit(1);
	}
	nInd_a=nCol/3-1;
	nInd = nInd_a;
}


//parse the line in buffer. Returns 0 for an empty line.
//With intName the marker is chr_pos, where chr can contain '_': chr_A_B_123 -> chr: chr_A_B ; pos: 123
int beagle_reader::parseLine(int *refId,int *posi,char *major,char *minor,double *post){
	//tokenize a copy, buffer is kept for when the line starts a new chromosome
	if(bufSize<l){
		buf=(char*)realloc(buf,l);
		bufSize=l;
	}
	strcpy(buf,buffer);

	char *save=NULL;
	char *marker=strtok_r(buf,delims,&save);
	if(marker==NULL)
		return 0;
	char *allele1=strtok_r(NULL,delims,&save);
	char *allele2=strtok_r(NULL,delims,&save);
	if(allele2==NULL)
		lineError(lineNo,"expected marker, allele1, allele2 and %d genotype probabilities",3*nInd);

	if(intName){
		char *us=strrchr(marker,'_');
		if(us==NULL)
			lineError(lineNo,"marker \'%s\' is not of the form chr_pos (use -intName 0 for other marker names)",marker);
		*us='\0';
		char *end=NULL;
		long p=strtol(us+1,&end,10);
		if(end==us+1||*end!='\0'||p<1||p>INT_MAX)
			lineError(lineNo,"invalid position \'%s\' in marker",us+1);
		aMap::const_iterator it = revMap->find(marker);
		if(it==revMap->end())
			lineError(lineNo,"chromosome \'%s\' is not in the fai file",marker);
		*refId=it->second;
		*posi=p-1;
	}

	*major=refToChar[(unsigned char)allele1[0]];
	*minor=refToChar[(unsigned char)allele2[0]];

	for(int i=0;i<nInd*3;i++){
		char *tok=strtok_r(NULL,delims,&save);
		if(tok==NULL)
			lineError(lineNo,"found %d genotype probabilities, expected %d",i,3*nInd);
		char *end=NULL;
		post[i]=strtod(tok,&end);
		if(end==tok||*end!='\0')
			lineError(lineNo,"genotype probability \'%s\' is not a number",tok);
	}
	if(strtok_r(NULL,delims,&save)!=NULL)
		lineError(lineNo,"more than the %d genotype probabilities given by the header",3*nInd);
	return 1;
}


funkyPars *beagle_reader::fetch(int chunksize){
	funkyPars * myfunky =funkyPars_init();
	myfunky->posi = new int[chunksize];
	myfunky->major = new char[chunksize];
	myfunky->minor = new char[chunksize];
	double **post = new double*[chunksize];
	for(int s=0;s<chunksize;s++)
		post[s] = new double[nInd*3];

	int nSites=0;
	//buffer still holds a line from the previous call that started a new chromosome,
	//so parse it before reading more
	int haveLine = changed;
	changed=0;
	while(haveLine||aio::tgets(gz,&buffer,&l)){
		if(!haveLine)
			lineNo++;
		haveLine=0;

		int refId=0;
		int posi=0;
		if(!parseLine(&refId,&posi,&myfunky->major[nSites],&myfunky->minor[nSites],post[nSites]))
			continue;//empty line
		if(intName){
			if(nSites>0&&refId!=myfunky->refId){
				//new chromosome: keep the line in buffer for the next call
				changed=1;
				break;
			}
		}else
			posi=positions++;
		myfunky->refId=refId;
		myfunky->posi[nSites]=posi;
		nSites++;
		if(nSites>=chunksize)
			break;
	}

	for(int s=nSites;s<chunksize;s++)
		delete[] post[s];
	myfunky->nInd=nInd;
	myfunky->post=post;
	myfunky->numSites = nSites;

	if(nSites==0){
		fprintf(stderr,"\t-> Done reading beagle\n");
		funkyPars_destroy(myfunky);
		return(NULL);
	}
	return(myfunky);
}
