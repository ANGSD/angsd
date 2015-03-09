/*
Copyright (c) 2009 Genome Research Ltd.
Author: Rob Davies <rmd@sanger.ac.uk>

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

   3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
Institute nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "tpooled_alloc.h"

//#define TEST_MAIN

#define PSIZE 1024*1024

tpool_alloc_t *tpool_create(size_t dsize) {
    tpool_alloc_t *p;

    if (NULL == (p = (tpool_alloc_t *)malloc(sizeof(*p))))
	return NULL;

    /* Minimum size is a pointer, for free list */
    dsize = (dsize + sizeof(void *) - 1) & ~(sizeof(void *)-1);
    if (dsize < sizeof(void *))
	dsize = sizeof(void *);
    p->dsize = dsize;

    p->ntpools = 0;
    p->tpools = NULL;
    p->free  = NULL;

    return p;
}

static tpool_t *new_tpool(tpool_alloc_t *p) {
    size_t n = PSIZE / p->dsize;
    tpool_t *tpool;
    
    tpool = realloc(p->tpools, (p->ntpools + 1) * sizeof(*p->tpools));
    if (NULL == tpool) return NULL;
    p->tpools = tpool;
    tpool = &p->tpools[p->ntpools];

    tpool->tpool = malloc(n * p->dsize);
    if (NULL == tpool->tpool) return NULL;

    tpool->used = 0;

    p->ntpools++;

    return tpool;
}

void tpool_destroy(tpool_alloc_t *p) {
    size_t i;

    for (i = 0; i < p->ntpools; i++) {
        free(p->tpools[i].tpool);
    }
    free(p->tpools);
    free(p);
}

void *tpool_alloc(tpool_alloc_t *p) {
    tpool_t *tpool;
    void *ret;

    /* Look on free list */
    if (NULL != p->free) {
        ret = p->free;
	p->free = *((void **)p->free);
	return ret;
    }

    /* Look for space in the last tpool */
    if (p->ntpools) {
        tpool = &p->tpools[p->ntpools - 1];
        if (tpool->used + p->dsize < PSIZE) {
	    ret = ((char *) tpool->tpool) + tpool->used;
	    tpool->used += p->dsize;
	    return ret;
	}
    }

    /* Need a new tpool */
    tpool = new_tpool(p);
    if (NULL == tpool) return NULL;

    tpool->used = p->dsize;
    return tpool->tpool;
}

void tpool_free(tpool_alloc_t *p, void *ptr) {
    *(void **)ptr = p->free;
    p->free = ptr;
}

#ifdef TEST_MAIN
typedef struct {
    int x, y, z;
} xyz;

#define NP 10000
int main(void) {
    int i;
    xyz *item;
    xyz **items;
    tpool_alloc_t *p = tpool_create(sizeof(xyz));

    items = (xyz **)malloc(NP * sizeof(*items));

    for (i = 0; i < NP; i++) {
	item = tpool_alloc(p);
	item->x = i;
	item->y = i+1;
	item->z = i+2;
	items[i] = item;
    }

    for (i = 0; i < NP; i++) {
	item = items[i];
	if (i % 3)
	    tpool_free(p, item);
    }

    for (i = 0; i < NP; i++) {
	item = tpool_alloc(p);
	item->x = 1000000+i;
	item->y = 1000000+i+1;
	item->z = 1000000+i+2;
    }

    for (i = 0; i < NP; i++) {
	item = items[i];
	printf("%d\t%d\t%d\t%d\n", i, item->x, item->y, item->z);
	tpool_free(p, item);
    }

    return 0;
}
#endif
