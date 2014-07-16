#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include <ctype.h>

const char *hapfile=NULL,mapfile=NULL,icounts=NULL;

typedef unsigned char uchar;

uchar flip(uchar c){
  c = toupper(c);
  if(c=='A')
    return 'T';
  if(c=='T')
    return 'A';
  if(c=='G')
    return 'C';
  if(c=='C')
    return 'G';
}

int main(int argc,char**argv){
  hapfile="../RES/hapMapCeuXlift.map.gz";

}
