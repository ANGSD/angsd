#!/bin/bash
if [ $# -eq 1 ] 
then
    WDIR=$1
fi

LOG=${0}.log
echo "Cleaning old output dir fst/output/" &>${LOG}
rm -rf fst/output ${LOG}
mkdir -p fst/output

MSMS=fst/input/msout.gz
gunzip -c ${MSMS} >fst/input/msout
MSMS=fst/input/msout
ODIR=fst/output/
##below is for simulating haplotypes
#mkdir -p tajima/input
#msms -ms 40 100 -t 900 -r 400 -SAA 1000 -SaA 500 -N 10000 -SF 0 -Sp .5 -oTPi -oAFS >${MSMS}
#grep -i afs ${MSMS} |grep -i sum|cut -f2 -d":"|tr " " "\n" >${MSMS}.sfs

echo "Generating genotype likelihood based on the haplotypes" >>${LOG} 2>&1
${WDIR}/misc/msToGlf -in ${MSMS} -out ${ODIR}/glout -err 0.005 -depth 8 -singleOut 1 -regLen 0 >>${LOG} 2>&1
echo "Splitting gl file into different populations"
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 1 6 >${ODIR}/pop1.glf.gz
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 7 13 >${ODIR}/pop2.glf.gz
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 14 22 >${ODIR}/pop3.glf.gz
echo 1 250000000 >${ODIR}/fai.fai
echo "Calculating perpop saf"
${WDIR}/angsd -glf ${ODIR}/pop1.glf.gz -nind 6 -doSaf 1 -out ${ODIR}/pop1 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop2.glf.gz -nind 7 -doSaf 1 -out ${ODIR}/pop2 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop3.glf.gz -nind 9 -doSaf 1 -out ${ODIR}/pop3 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}

echo "Calculating perpop sfs based on the perpop saf"
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx >${ODIR}/pop1.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop2.saf.idx >${ODIR}/pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop3.saf.idx >${ODIR}/pop3.saf.idx.ml 2>>${LOG}
echo "Calculating parwise 2d sfs"
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx >${ODIR}/pop1.pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx >${ODIR}/pop1.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx >${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
echo "Calculating fst index for 3 pairwise, and multi fst"
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx -fstout ${ODIR}/pop1.pop2 -sfs ${ODIR}/pop1.pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop3 -sfs ${ODIR}/pop1.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop2.pop3 -sfs ${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop2.pop3 -sfs ${ODIR}/pop1.pop2.saf.idx.ml -sfs ${ODIR}/pop1.pop3.saf.idx.ml -sfs ${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
echo "Calculating fst stats for 3 pairwise, and multi fst"
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.fst.idx >${ODIR}/pop1.pop2.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop3.fst.idx >${ODIR}/pop1.pop3.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop2.pop3.fst.idx >${ODIR}/pop2.pop3.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.pop3.fst.idx >${ODIR}/pop1.pop2.pop3.fst.idx.res 2>>${LOG}


##when generated the results:
##md5sum tajima/output/*pestPG tajima/output/*.ml tajima/output/*.saf tajima/output/*.gz  >tajima/md5/pestPG.md5sum
md5sum  -c fst/md5/fst.md5sum >>${LOG} 2>&1 || exit 1
