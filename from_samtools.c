#include "from_samtools.h"
#include <stdio.h>
#include <assert.h>
/*
  This file contains modified versions from SAMtools (1.5-1-g27b628e (using htslib 1.5-2-g2739558))
  thorfinn@binf.ku.dk 7aug 2017 cambridge
 */


void* add_read_group_single(char *name){
    char *d = strdup(name);
    void * retval = NULL;
    int ret = 0;

    if (d == NULL) goto err;

    retval = khash_str2int_init();
    if (retval == NULL) goto err;
    ret = khash_str2int_inc(retval,name);
    fprintf(stderr,"\t-> [READGROUP info] %s->%d\n",name,ret);
    if (ret == -1) goto err;
    //if (ret ==  0) free(d); /* Duplicate */ //doesnt work with str2int dragon
    return retval;

 err:
    fprintf(stderr, "Couldn't add \"%s\" to read group list: memory exhausted?", name);
    free(d);
    return NULL;
}

void* add_read_groups_file(char *fn)
{
    FILE *fp;
    char buf[1024];
    int ret = 0;


    void *retval = khash_str2int_init();

    fp = fopen(fn, "r");
    if (fp == NULL) {
        fprintf(stderr, "failed to open \"%s\" for reading", fn);
        return NULL;
    }

    while (ret != -1 && !feof(fp) && fscanf(fp, "%1023s", buf) > 0) {
      char *d = strdup(buf);
      if (d != NULL) 
	ret = khash_str2int_inc(retval,d);
      if(ret==-1)
	fprintf(stderr,"\t-> [READGROUP info] Problems adding readgroup: %s\n",d);
      else{
	//int tmp;
	//	assert(khash_str2int_get(retval,d,&tmp)==0);
	fprintf(stderr,"\t-> [READGROUP info] Added readgroup %s\n",d);
      }
    }
    if (ferror(fp)) ret = -1;
    if (ret == -1) {
        fprintf(stderr, "failed to read \"%s\"", fn);
    }
    fclose(fp);
    return retval;
}

