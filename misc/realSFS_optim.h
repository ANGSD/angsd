#include <vector>
#include "safreader.h"
size_t calc_nsites(std::vector<persaf *> &pp,args *ar);
size_t parspace(std::vector<persaf *> &saf);
void setGloc(std::vector<persaf *> &saf,size_t nSites);
void delGloc(std::vector<persaf *> &saf,size_t nSites);
template <typename T>
int main_opt(args *arg);

