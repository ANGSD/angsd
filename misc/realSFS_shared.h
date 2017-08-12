#pragma once
#include "realSFS_args.h"
#include "safreader.h"
int set_intersect_pos(std::vector<persaf *> &saf,char *chooseChr,int start,int stop,char **curChr,filt *fl);
size_t calc_nsites(std::vector<persaf *> &pp,args *ar);
size_t helper(persaf * pp,char *chr);
void delGloc(std::vector<persaf *> &saf,size_t nSites);
void setGloc(std::vector<persaf *> &saf,size_t nSites);
size_t parspace(std::vector<persaf *> &saf);
template <typename T>
size_t fsizes(std::vector<persaf *> &pp, size_t nSites);
void readSFS(const char*fname,size_t hint,double *ret);
