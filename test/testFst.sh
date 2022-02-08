#!/bin/bash
MD5=md5sum

if [ $# -eq 1 ] 
then
    WDIR=$1
fi


if [[ "$OSTYPE" == "darwin"* ]]; then
    MD5=./md5osx.sh
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
##ms 44 10 -t 930 -r 400 -I 3 12 14 18 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1  [3.2rc Build:162]


echo "Generating genotype likelihood based on the haplotypes" >>${LOG} 2>&1
${WDIR}/misc/msToGlf -in ${MSMS} -out ${ODIR}/glout -err 0.005 -depth 8 -singleOut 1 -regLen 0 -simpleRand 0  >>${LOG} 2>&1
echo "Splitting gl file into different populations" >>${LOG} 2>&1
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 1 6 >${ODIR}/pop1.glf.gz 2>>${LOG}
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 7 13 >${ODIR}/pop2.glf.gz 2>>${LOG}
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 22 14 22 >${ODIR}/pop3.glf.gz 2>>${LOG}
echo 1 250000000 >${ODIR}/fai.fai
echo "Calculating perpop saf" >>${LOG} 2>&1
${WDIR}/angsd -glf ${ODIR}/pop1.glf.gz -nind 6 -doSaf 1 -out ${ODIR}/pop1 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop2.glf.gz -nind 7 -doSaf 1 -out ${ODIR}/pop2 -fai ${ODIR}/fai.fai -issim 1 2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop3.glf.gz -nind 9 -doSaf 1 -out ${ODIR}/pop3 -fai ${ODIR}/fai.fai -issim 1 2>>${LOG}

echo "Calculating perpop sfs based on the perpop saf" >>${LOG} 2>&1
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx -seed -1 >${ODIR}/pop1.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop2.saf.idx -seed -1 >${ODIR}/pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop3.saf.idx -seed -1 >${ODIR}/pop3.saf.idx.ml 2>>${LOG}
echo "Calculating parwise 2d sfs" >>${LOG} 2>&1
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx -seed -1 >${ODIR}/pop1.pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx -seed -1 >${ODIR}/pop1.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -seed -1 >${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
echo "Calculating fst index for 3 pairwise, and multi fst" >>${LOG} 2>&1
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx -fstout ${ODIR}/pop1.pop2 -sfs ${ODIR}/pop1.pop2.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop3 -sfs ${ODIR}/pop1.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop2.pop3 -sfs ${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop2.pop3 -sfs ${ODIR}/pop1.pop2.saf.idx.ml -sfs ${ODIR}/pop1.pop3.saf.idx.ml -sfs ${ODIR}/pop2.pop3.saf.idx.ml 2>>${LOG}
echo "Calculating fst stats for 3 pairwise, and multi fst" >>${LOG} 2>&1
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.fst.idx >${ODIR}/pop1.pop2.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop3.fst.idx >${ODIR}/pop1.pop3.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop2.pop3.fst.idx >${ODIR}/pop2.pop3.fst.idx.res 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.pop3.fst.idx >${ODIR}/pop1.pop2.pop3.fst.idx.res 2>>${LOG}


##when generated the results:
##md5sum tajima/output/*pestPG tajima/output/*.ml tajima/output/*.saf tajima/output/*.gz  >tajima/md5/pestPG.md5sum
#ls fst/output/*|grep -v arg|xargs -n1 md5sum >fst/md5/fst3.md5sum 
${MD5}  -c fst/md5/fst3.md5sum >>${LOG} 2>&1 || exit 1
