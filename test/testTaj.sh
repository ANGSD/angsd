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
echo "Cleaning old output dir taj/" &>${LOG}
rm -rf tajima/output ${LOG}
mkdir -p tajima/output

MSMS=tajima/input/msout.gz
gunzip -c ${MSMS} >tajima/input/msout
MSMS=tajima/input/msout
ODIR=tajima/output/
##below is for simulating haplotypes
#mkdir -p tajima/input
#msms -ms 40 100 -t 900 -r 400 -SAA 1000 -SaA 500 -N 10000 -SF 0 -Sp .5 -oTPi -oAFS >${MSMS}
#grep -i afs ${MSMS} |grep -i sum|cut -f2 -d":"|tr " " "\n" >${MSMS}.sfs

echo "Generating genotype likelihood based on the haplotypes" >>${LOG} 2>&1
${WDIR}/misc/msToGlf -in ${MSMS} -out ${ODIR}/glout -err 0.005 -depth 8 -singleOut 1 -simpleRand 0 >>${LOG} 2>&1
${WDIR}/angsd -isSim 1 -glf ${ODIR}/glout.glf.gz -out ${ODIR}/norm -doSaf 1 -nInd 20 -fai hg19.fa.fai  2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/norm.saf.idx -P 24 -nSites 1000000 -m 0 -seed -1 >${ODIR}/norm.saf.em.ml 2>>${LOG}
${WDIR}/angsd -isSim 1 -glf ${ODIR}/glout.glf.gz -out ${ODIR}/norm -nInd 20 -doThetas 1 -doSaf 1 -pest ${ODIR}/norm.saf.em.ml -fai hg19.fa.fai  2>>${LOG}
${WDIR}/misc/realSFS saf2theta ${ODIR}/norm.saf.idx -sfs ${ODIR}/norm.saf.em.ml -outname ${ODIR}/thetaFromSaf 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetas.idx 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/norm.thetas.idx -outnames ${ODIR}/norm.thetas.idx.win -step 1000 -win 5000 2>>${LOG}

echo "2) Will do folded analysis" >>${LOG} 2>&1

${WDIR}/angsd -isSim 1 -glf ${ODIR}/glout.glf.gz -out ${ODIR}/fold -doSaf 1 -nInd 20 -fold 1 -fai hg19.fa.fai  2>>${LOG}
${WDIR}/misc/realSFS ${ODIR}/fold.saf.idx -P 24 -nSites 1000000 -m 0 -seed -1 >${ODIR}/fold.saf.em.ml 2>>${LOG}
${WDIR}/angsd -isSim 1 -glf ${ODIR}/glout.glf.gz -out ${ODIR}/fold -nInd 20 -doThetas 1 -doSaf 1 -pest ${ODIR}/fold.saf.em.ml -fold 1  -fai hg19.fa.fai 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/fold.thetas.idx -nChr 40 2>>${LOG}
${WDIR}/misc/thetaStat do_stat ${ODIR}/fold.thetas.idx -nChr 40 -outnames ${ODIR}/fold.thetas.idx.win -step 1000 -win 5000 2>>${LOG}


echo -e "\tThe theta estimates from msms simulated haplotypes for 8x 0.5% errorrate"  >>${LOG} 2>&1
echo -e "\tWatterson\tpairwise\ttajimasD" >>${LOG} 2>&1
grep ThetaW ${MSMS}|grep Sum >>${LOG} 2>&1
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (unfolded)" >>${LOG} 2>&1
cat ${ODIR}/norm.thetas.idx.pestPG >>${LOG} 2>&1
echo -e "-------\nThe estimated thetas (sum of 100 reps) and the tajima value (unfolded)" >>${LOG} 2>&1
cat ${ODIR}/fold.thetas.idx.pestPG >>${LOG} 2>&1
echo "You should eyeball the above and see if they are comparable column (1-5,9) (not all defined for unfold)" >>${LOG} 2>&1

##when generated the results:
##md5sum tajima/output/*pestPG tajima/output/*.ml tajima/output/*.saf tajima/output/*.gz  >tajima/md5/pestPG.md5sum
${MD5}  -c tajima/md5/pestPG.md5sum >>${LOG} 2>&1 || exit 1
