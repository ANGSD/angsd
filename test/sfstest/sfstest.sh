#!/bin/bash


WDIR=$1

GLFINPUT=input/sfstest.glf.gz

LOG=../../testSFS.sh.log

##rm old
rm -rf input/ output/ ${LOG}

##make workdirs
mkdir -p input
mkdir -p output

##cp checksums
cp md5/sfsinput2.md5sum input/
cp md5/*.md5sum output

##make and check inputfiles
$WDIR/misc/supersim -outfiles input/sfstest -npop 1 -nind 10 -nsites 50000 -depth 8 -pvar 1.0 &>>${LOG} ||exit 1 
cd input
md5sum  -c sfsinput2.md5sum &>>${LOG} || exit 2
 

cd ..
PROGGY=$WDIR/angsd
$PROGGY -glf  $GLFINPUT -isSim 1 -nInd 10  -out output/sfstest -doSaf 1 -underFlowProtect 0 -fai ../hg19.fa.fai &>>${LOG}||exit 3


cd output/
md5sum  -c results2.md5sum &>>${LOG}||exit 4

cd ..
$WDIR/misc/realSFS output/sfstest.saf.idx -P 4 -nSites 10000000 -m 1 -seed -1 >output/sfs.est  2>>${LOG} ||exit 5

cd output
md5sum  -c mlRes3.md5sum &>>${LOG} ||exit 7



#echo do GC test
cd ..
$PROGGY -glf  $GLFINPUT -nInd 10  -out output/sfstestGC -doPost 3 -underFlowProtect 0 -fai ../hg19.fa.fai -doMajorMinor 1 -pest output/sfs.est -domaf 1 -dogeno 8  &>>${LOG} ||exit 8
cd output
md5sum  -c sfsgc2.md5sum &>>${LOG} ||exit 9
