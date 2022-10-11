#include <cstring>
#include <stdarg.h>
#include <sys/stat.h>
#include "aio.h"

int aio::fexists(const char* str){///@param str Filename given as a string.
	struct stat buffer ;
	return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

size_t aio::fsize(const char* fname){
	struct stat st ;
	stat(fname,&st);
	return st.st_size;
}

std::vector <char *> dumpedFiles;//small hack for getting a nice vector of outputfiles
FILE *aio::openFile(const char* a,const char* b){
	if(0)
		fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
	char *c = new char[strlen(a)+strlen(b)+1];
	strcpy(c,a);
	strcat(c,b);
	//  fprintf(stderr,"\t-> Dumping file: %s\n",c);
	dumpedFiles.push_back(strdup(c));
	FILE *fp = NULL;
	fp = fopen(c,"w");
	if(fp==NULL){
		fprintf(stderr,"\t-> Problem opening file: \'%s\' check permissions\n",c);
		exit(0);
	}
	delete [] c;
	return fp;
}

BGZF *aio::openFileBG(const char* a,const char* b){

	char *c = new char[strlen(a)+strlen(b)+1];
	strcpy(c,a);
	strcat(c,b);
	dumpedFiles.push_back(strdup(c));
	BGZF *fp = bgzf_open(c,"w6h");
	delete [] c;
	return fp;
}

htsFile *aio::openFileHts(const char* a,const char* b){

	char *c = new char[strlen(a)+strlen(b)+1];
	strcpy(c,a);
	strcat(c,b);
	dumpedFiles.push_back(strdup(c));
	htsFile *fp = hts_open(c,"w");
	delete [] c;
	return fp;
}

htsFile *aio::openFileHtsBcf(const char* a,const char* b){

	char *c = new char[strlen(a)+strlen(b)+1];
	strcpy(c,a);
	strcat(c,b);
	dumpedFiles.push_back(strdup(c));
	htsFile *fp = hts_open(c,"wb");
	delete [] c;
	return fp;
}

FILE *aio::getFILE(const char*fname,const char* mode){
	int writeFile = 0;
	for(size_t i=0;i<strlen(mode);i++)
		if(mode[i]=='w')
			writeFile = 1;
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"\t-> Error opening FILE handle for file:%s exiting\n",fname);
		exit(0);
	}
	return fp;
}

//checks that newer is newer than older
int aio::isNewer(const char *newer,const char *older){
	if (strstr(older, "ftp://") == older || strstr(older, "http://") == older)
		return 0;
	//  fprintf(stderr,"newer:%s older:%s\n",newer,older);
	// return 0;
	struct stat one;
	struct stat two;
	stat(newer, &one );
	stat(older, &two );

	return one.st_mtime>=two.st_mtime;
}

ssize_t aio::bgzf_write(BGZF *fp, const void *data, size_t length){
	if(length>0)
		return ::bgzf_write(fp,data,length);
	return 0;
}


int aio::tgets(gzFile gz,char**buf,int *l){
	int rlen = 0;
neverUseGoto:
	char *tok = gzgets(gz,*buf+rlen,*l-rlen);
	if(!tok)
		return rlen;
	int tmp = tok?strlen(tok):0;
	if(tok[tmp-1]!='\n'){
		rlen += tmp;
		*l *= 2;
		*buf = (char*) realloc(*buf,*l);
		goto neverUseGoto;
	}
	rlen += tmp;
	return rlen;
}

//  Usage: aio::doAssert(something==NULL,1,AT,"");
void aio::doAssert(int EXIT, int EXIT_CODE, const char* error_location, const char* format,...){
	if(EXIT){
		va_list args;
		va_start (args, format);
		fprintf(stderr, "\n");
		fprintf(stderr, "*******\n");
		fprintf(stderr, "[ERROR](%s)\n",error_location);
		vfprintf (stderr, format, args);
		va_end (args);
		fprintf(stderr, "\n");
		fprintf(stderr, "*******\n");
		exit(EXIT_CODE);
	}else{
		return;
	}
}

