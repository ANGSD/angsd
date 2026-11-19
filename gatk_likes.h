#pragma once
#include "bambi_interface.h"
void gatk_init();
void call_gatk(chunkyT *chk,double **lk,int trim,int *keepSites);
void gatk_destroy();
