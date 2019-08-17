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
echo "Cleaning old output dir fst_folded/output/" &>${LOG}
rm -rf fst_folded/output ${LOG}
mkdir -p fst_folded/output

##below is for simulating haplotypes
#MS=/willerslev/software/msms/bin/msms
#$MS -ms 90 10 -t 930 -r 400 -I 3 26 30 34 -n 1 1.682020 -n 2 3.736830 -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 -ma x 0.881098 0.561966 0.881098 x 2.797460 0.561966 2.797460 x -ej 0.028985 3 2 -en 0.028985 2 0.287184 -ema 0.028985 3 x 7.293140 x 7.293140 x x x x x -ej 0.197963 2 1 -en 0.303501 1 1 |bgzip -c >${MSMS}

MSMS=fst_folded/input/msout.gz
gunzip -c ${MSMS} >fst_folded/input/msout
MSMS=fst_folded/input/msout
ODIR=fst_folded/output/


echo "Generating genotype likelihood based on the haplotypes" >>${LOG} 2>&1
${WDIR}/misc/msToGlf -in ${MSMS} -out ${ODIR}/glout -err 0.005 -depth 8 -singleOut 1 -regLen 0 -simpleRand 0  >>${LOG} 2>&1
echo "Splitting gl file into different populations" >>${LOG} 2>&1
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 45 1 13 >${ODIR}/pop1.glf.gz 2>>${LOG}
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 45 14 28 >${ODIR}/pop2.glf.gz 2>>${LOG}
${WDIR}/misc/splitgl ${ODIR}/glout.glf.gz 45 29 45 >${ODIR}/pop3.glf.gz 2>>${LOG}
echo 1 250000000 >${ODIR}/fai.fai
echo "Calculating perpop saf" >>${LOG} 2>&1
${WDIR}/angsd -glf ${ODIR}/pop1.glf.gz -nind 13 -doSaf 1 -out ${ODIR}/pop1 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop2.glf.gz -nind 15 -doSaf 1 -out ${ODIR}/pop2 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}
${WDIR}/angsd -glf ${ODIR}/pop3.glf.gz -nind 17 -doSaf 1 -out ${ODIR}/pop3 -fai ${ODIR}/fai.fai -issim 1  2>>${LOG}

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


echo "Calculating parwise 2d sfs (folded)" >>${LOG} 2>&1
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx -seed -1 -fold 1 >${ODIR}/pop1.pop2.saf.idx.ml.fold 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx -seed -1 -fold 1 >${ODIR}/pop1.pop3.saf.idx.ml.fold 2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -seed -1 -fold 1 >${ODIR}/pop2.pop3.saf.idx.ml.fold 2>>${LOG}
##need to fix below still
echo "Calculating fst index for 3 pairwise, and multi fst" >>${LOG} 2>&1
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx -fstout ${ODIR}/pop1.pop2.fold -sfs ${ODIR}/pop1.pop2.saf.idx.ml.fold -fold 1 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop3.fold -sfs ${ODIR}/pop1.pop3.saf.idx.ml.fold -fold 1 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop2.pop3.fold -sfs ${ODIR}/pop2.pop3.saf.idx.ml.fold -fold 1 2>>${LOG}
${WDIR}/misc/realSFS fst index ${ODIR}/pop1.saf.idx ${ODIR}/pop2.saf.idx ${ODIR}/pop3.saf.idx -fstout ${ODIR}/pop1.pop2.pop3.fold -sfs ${ODIR}/pop1.pop2.saf.idx.ml.fold -sfs ${ODIR}/pop1.pop3.saf.idx.ml.fold -sfs ${ODIR}/pop2.pop3.saf.idx.ml.fold -fold 1 2>>${LOG}
echo "Calculating fst stats for 3 pairwise, and multi fst" >>${LOG} 2>&1
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.fold.fst.idx >${ODIR}/pop1.pop2.fst.idx.res.fold 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop3.fold.fst.idx >${ODIR}/pop1.pop3.fst.idx.res.fold 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop2.pop3.fold.fst.idx >${ODIR}/pop2.pop3.fst.idx.res.fold 2>>${LOG}
${WDIR}/misc/realSFS fst stats ${ODIR}/pop1.pop2.pop3.fold.fst.idx >${ODIR}/pop1.pop2.pop3.fst.idx.res.fold 2>>${LOG}




##when generated the results:
#md5sum fst_folded/output/*|grep -v -P "arg$" >fst_folded/md5/fst.md5sum
${MD5}  -c fst_folded/md5/fst.md5sum >>${LOG} 2>&1 || exit 1
