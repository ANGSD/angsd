#pragma once
#include "bambi_interface.h"

//-GL 1, the SAMtools model from htslib's errmod
void bam_likes_init();
void call_bam(chunkyT *chk,double **likes,int trim,int *keepSites,int seed);
void bam_likes_destroy();
