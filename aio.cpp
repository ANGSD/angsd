#include <cstring>
#include <cstdlib>
#include <stdarg.h>
#include <sys/stat.h>
#include "aio.h"

int aio::fexists(const char* str){///@param str Filename given as a string.
	struct stat buffer ;
	return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

size_t aio::fsize(const char* fname){
	struct stat st ;
	if(stat(fname,&st)!=0)
		return 0;
	return st.st_size;
}

//prefix+postfix as a new[]'ed string
static char *joinName(const char *a,const char *b){
	char *c = new char[strlen(a)+strlen(b)+1];
	strcpy(c,a);
	strcat(c,b);
	return c;
}

static void openFailed(const char *fname){
	fprintf(stderr,"\t-> Problem opening file: \'%s\' check permissions\n",fname);
	exit(1);
}

std::vector <char *> dumpedFiles;//small hack for getting a nice vector of outputfiles
FILE *aio::openFile(const char* a,const char* b){
	char *c = joinName(a,b);
	dumpedFiles.push_back(strdup(c));
	FILE *fp = fopen(c,"w");
	if(fp==NULL)
		openFailed(c);
	delete [] c;
	return fp;
}

BGZF *aio::openFileBG(const char* a,const char* b){
	char *c = joinName(a,b);
	dumpedFiles.push_back(strdup(c));
	BGZF *fp = bgzf_open(c,"w6h");
	if(fp==NULL)
		openFailed(c);
	delete [] c;
	return fp;
}

htsFile *aio::openFileHts(const char* a,const char* b){
	char *c = joinName(a,b);
	dumpedFiles.push_back(strdup(c));
	htsFile *fp = hts_open(c,"w");
	if(fp==NULL)
		openFailed(c);
	delete [] c;
	return fp;
}

htsFile *aio::openFileHtsBcf(const char* a,const char* b){
	char *c = joinName(a,b);
	dumpedFiles.push_back(strdup(c));
	htsFile *fp = hts_open(c,"wb");
	if(fp==NULL)
		openFailed(c);
	delete [] c;
	return fp;
}

FILE *aio::getFILE(const char*fname,const char* mode){
	FILE *fp;
	if(NULL==(fp=fopen(fname,mode))){
		fprintf(stderr,"\t-> Error opening FILE handle for file:%s exiting\n",fname);
		exit(1);
	}
	return fp;
}

//checks that newer was modified after older. Files with the same
//modification time (e.g. copied with cp -p, rsync -a or tar) are not
//newer. If either file can't be stat'ed, it is not newer either.
int aio::isNewer(const char *newer,const char *older){
	if (strstr(older, "ftp://") == older || strstr(older, "http://") == older)
		return 0;
	struct stat one;
	struct stat two;
	if(stat(newer, &one)!=0 || stat(older, &two)!=0)
		return 0;

	return one.st_mtime>two.st_mtime;
}

ssize_t aio::bgzf_write(BGZF *fp, const void *data, size_t length){
	if(length>0)
		return ::bgzf_write(fp,data,length);
	return 0;
}


//read a whole line, growing *buf (whose capacity is *l) until it fits
int aio::tgets(gzFile gz,char**buf,int *l){
	int rlen = 0;
	while(1){
		char *tok = gzgets(gz,*buf+rlen,*l-rlen);
		//a read error, e.g. a truncated gzip file, would otherwise look
		//like the end of the file or a last line without a newline
		int err;
		const char *msg = gzerror(gz,&err);
		if(err<0){
			fprintf(stderr,"\t-> Problem reading file: %s\n",msg);
			exit(1);
		}
		if(!tok)
			return rlen;//end of file
		int tmp = strlen(tok);
		rlen += tmp;
		//a complete line, or the last line of a file without a newline
		if(tmp==0 || tok[tmp-1]=='\n' || gzeof(gz))
			return rlen;
		*l *= 2;
		char *tmpbuf = (char*) realloc(*buf,*l);
		if(tmpbuf==NULL){
			fprintf(stderr,"\t-> Problem allocating %d bytes for a line\n",*l);
			exit(1);
		}
		*buf = tmpbuf;
	}
}


static void assertHeader(const char *error_location){
	fprintf(stderr, "\n");
	fprintf(stderr, "*******\n");
	if(error_location)
		fprintf(stderr, "[ERROR](%s)\n",error_location);
	else
		fprintf(stderr, "[ERROR]\n");
}

static void assertExit(int exit_code){
	fprintf(stderr, "*******\n");
	exit(exit_code);
}

void aio::doAssert(int exp_eval, int exit_code, const char* error_location, const char* exit_text,...){
	if(exp_eval)
		return;
	assertHeader(error_location);
	va_list args;
	va_start (args, exit_text);
	vfprintf (stderr, exit_text, args);
	va_end (args);
	fprintf(stderr, "\n");
	assertExit(exit_code);
}


// Function overload to avoid specifying exit code and exit text. If not specified exit 1
void aio::doAssert(int exp_eval, const char* error_location){
	if(exp_eval)
		return;
	assertHeader(error_location);
	assertExit(1);
}


// Function overload to avoid specifying anything, just evaluate and exit 1
void aio::doAssert(int exp_eval){
	if(exp_eval)
		return;
	assertHeader(NULL);
	assertExit(1);
}
