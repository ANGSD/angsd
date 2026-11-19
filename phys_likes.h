#pragma once
#include <vector>
#include "bambi_interface.h"
//-GL 5 (phys). DEPRECATED: will be removed in a future version.
void phys_init(std::vector<char*> bamnames);
void call_phys(chunkyT *chk,double **lk,int trim);
void phys_destroy();
