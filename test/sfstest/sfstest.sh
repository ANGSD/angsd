#!/bin/bash


WDIR=$1

GLFINPUT=input/sfstest.glf.gz

#absolute path, since the script cd's into input/ and output/ while logging
LOG="$(cd "$(dirname "$0")/.." && pwd)/testSFS.sh.log"

MD5=md5sum
if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5="$(cd "$(dirname "$0")/.." && pwd)/md5osx.sh"
fi

##rm old
rm -rf input/ output/ ${LOG}

##make workdirs
mkdir -p input
mkdir -p output

##cp checksums
cp md5/sfsinput2.md5sum input/
cp md5/*.md5sum output

##make and check inputfiles
$WDIR/misc/supersim -outfiles input/sfstest -npop 1 -nind 10 -nsites 50000 -depth 8 -pvar 1.0 >>${LOG} 2>&1 ||exit 1 
cd input
${MD5}  -c sfsinput2.md5sum >>${LOG} 2>&1 || exit 2
 

cd ..
PROGGY=$WDIR/angsd
$PROGGY -glf  $GLFINPUT -isSim 1 -nInd 10  -out output/sfstest -doSaf 1 -underFlowProtect 0 -fai ../hg19.fa.fai >>${LOG} 2>&1||exit 3


cd output/
${MD5}  -c results2.md5sum >>${LOG} 2>&1||exit 4

cd ..
$WDIR/misc/realSFS output/sfstest.saf.idx -P 4 -nSites 10000000 -m 1 -seed -1 >output/sfs.est  2>>${LOG} ||exit 5

cd output
${MD5}  -c mlRes3.md5sum >>${LOG} 2>&1 ||exit 7



#echo do GC test
cd ..
$PROGGY -glf  $GLFINPUT -nInd 10  -out output/sfstestGC -doPost 3 -underFlowProtect 0 -fai ../hg19.fa.fai -doMajorMinor 1 -pest output/sfs.est -domaf 1 -dogeno 8  >>${LOG} 2>&1 ||exit 8
cd output
${MD5}  -c sfsgc2.md5sum >>${LOG} 2>&1 ||exit 9
